nextflow.enable.dsl=2

// Parameters
params.samplesheet = "${workflow.projectDir}/data/samplesheet.csv"
params.bed         = "${workflow.projectDir}/data/merged_output.bed"
params.outdir      = "results"

params.vcfdir      = "${params.outdir}/vcf"
params.qcdir       = "${params.outdir}/qc"
params.reportdir   = "${params.outdir}/reports"
params.scriptdir   = "${workflow.projectDir}/scripts"

// VEP parameters
def HOME = System.getenv('HOME') ?: '.'
params.run_vep   = params.run_vep ?: true
params.vep_cache = params.vep_cache ?: "${HOME}/vep_data"
params.vep_fasta = params.vep_fasta ?: "${HOME}/vep_data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

// Input channel
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.sample, file(row.vcf), file(row.bam)) }
    .set { sample_ch }

// Step 1: Filter VCF by BED regions
process BedFilterVCF {
    tag "$sample"
    publishDir "${params.vcfdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(vcf), path(bam)
    path bed

    output:
    tuple val(sample), path("${sample}.bed_filtered.vcf.gz"), path(bam)

    script:
    """
    # Ensure index is created for the staged copy
    tabix -p vcf $vcf || bcftools index -t $vcf
    
    bcftools view -R $bed $vcf -Oz -o ${sample}.bed_filtered.vcf.gz
    tabix -p vcf ${sample}.bed_filtered.vcf.gz
    """
}

// Step 2: Normalize VCF (split multi-allelics)
process NormalizeVCF {
    tag "$sample"
    publishDir "${params.vcfdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(vcf), path(bam)

    output:
    tuple val(sample), path("${sample}.normalized.vcf.gz"), path(bam)

    script:
    """
    bcftools norm -m -any $vcf -Oz -o ${sample}.normalized.vcf.gz
    tabix -p vcf ${sample}.normalized.vcf.gz
    """
}

// Step 3: Filter by Depth and Quality
process FilterVCF {
    tag "$sample"
    publishDir "${params.vcfdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(vcf), path(bam)

    output:
    tuple val(sample), path("${sample}.filtered.vcf.gz"), path(bam)

    script:
    """
    bcftools view -i 'INFO/DP >= 20 && QUAL >= 30' $vcf -Oz -o ${sample}.filtered.vcf.gz
    tabix -p vcf ${sample}.filtered.vcf.gz
    """
}

// Step 4: Add VAF tags (no filtering)
process AddVAF {
    tag "$sample"
    publishDir "${params.vcfdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(vcf), path(bam)

    output:
    tuple val(sample), path("${sample}.vaf_added.vcf.gz")

    script:
    """
    bcftools +fill-tags $vcf -Oz -o ${sample}.vaf_added.vcf.gz -- -t FORMAT/VAF
    tabix -p vcf ${sample}.vaf_added.vcf.gz
    """
}

// Step 5: Filter BAM by BED regions
process BedFilterBAM {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(bam)
    path bed

    output:
    tuple val(sample), path(vcf), path("${sample}.bed_filtered.bam"), path("${sample}.bed_filtered.bam.bai")

    script:
    """
    samtools view -L $bed -b -@ 4 $bam -o tmp.bam
    samtools sort -o ${sample}.bed_filtered.bam tmp.bam
    samtools index ${sample}.bed_filtered.bam
    """
}

// Step 6: Basic coverage summary
process CoverageSummary {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)
    path bed

    output:
    tuple val(sample), path("${sample}_coverage_summary.sorted.txt")

    script:
    """
    bedtools coverage -a $bed -b $bam -d > ${sample}_coverage_per_base.txt
    
    awk '{
        key=\$1":"\$2"-"\$3
        total[key]++
        if(\$5>=20) c20[key]++
        if(\$5>=30) c30[key]++
        if(\$5>=50) c50[key]++
        if(\$5>=100) c100[key]++
    }
    END {
        for (k in total) {
            printf "%s\\t>=20x:%.2f%%\\t>=30x:%.2f%%\\t>=50x:%.2f%%\\t>=100x:%.2f%%\\n", \
                   k, (c20[k]/total[k])*100, (c30[k]/total[k])*100, (c50[k]/total[k])*100, (c100[k]/total[k])*100
        }
    }' ${sample}_coverage_per_base.txt > ${sample}_coverage_summary.txt

    sort -k1,1V -k1.2,1n -t: -k2,2n ${sample}_coverage_summary.txt > ${sample}_coverage_summary.sorted.txt
    """
}

// Step 7: R1/R2 read pair ratio
process R1R2Ratio {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai)
    path bed

    output:
    tuple val(sample), path("${sample}_r1r2_per_exon.tsv")

    script:
    """
    echo "BED preview:" >&2
    head -n 5 $bed >&2

    echo "Running on BAM: $bam" >&2
    samtools view -H $bam | head -n 3 >&2
    while read chrom start end; do
        region="\${chrom}:\${start}-\${end}"
        counts=\$(samtools view -F 0x4 $bam "\$region" | \
          awk '{flag=\$2; if(and(flag,64)) r1++; if(and(flag,128)) r2++} END {if(r1+r2>0) printf("%d\\t%d\\t%.3f\\n", r1, r2, r1/(r1+r2)); else print "0\\t0\\tNA"}')
        echo -e "\${chrom}\\t\${start}\\t\${end}\\t\${counts}"
    done < $bed > ${sample}_r1r2_per_exon.tsv
    """
}

// Step 8: Forward/Reverse strand balance
process ForwardReverseRatio {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai)
    path bed

    output:
    tuple val(sample), path("${sample}_frstrand_per_exon.tsv")

    script:
    """
    while read chrom start end; do
        region="\${chrom}:\${start}-\${end}"
        echo "Processing region: \$region" >&2
        counts=\$(samtools view -F 0x4 $bam "\$region" | \
          awk '{
              flag=\$2
              if(and(flag,16)) rev++
              else fwd++
          } END {
              if(fwd+rev>0) {
                  frac = rev/(fwd+rev)
                  balance = (fwd/(fwd+rev) < rev/(fwd+rev) ? fwd/(fwd+rev) : rev/(fwd+rev))
                  printf("%d\\t%d\\t%.3f\\t%.3f\\n", fwd, rev, frac, balance)
              } else {
                  print "0\\t0\\tNA\\tNA"
              }
          }')
        echo -e "\${chrom}\\t\${start}\\t\${end}\\t\${counts}"
    done < $bed > ${sample}_frstrand_per_exon.tsv
    """
}

// Step 10: VEP Annotation with plugins
process VEP_Annotate {
  tag "$sample"
  publishDir "${params.vcfdir}", mode: 'copy'

  input:
  tuple val(sample), path(vcf)

  output:
  tuple val(sample), path("${sample}.vep.vcf.gz")

  script:
  """
  set -euo pipefail

  VCF_PATH=\$(realpath "$vcf")
  VCF_DIR=\$(dirname "\$VCF_PATH")
  VCF_BASE=\$(basename "\$VCF_PATH")

  # Mount the process work dir so VEP can read/write there
  docker run --rm -t \\
    -u \$(id -u):\$(id -g) \\
    -v "${params.vep_cache}:/cache:ro" \\
    -v "\$VCF_DIR:/work" \\
    ensemblorg/ensembl-vep \\
      vep \\
        -i /work/\$VCF_BASE \\
        -o /work/${sample}.vep.vcf.gz \\
        --vcf --compress_output bgzip \\
        --offline --cache --dir_cache /cache \\
        --fasta /cache/\$(basename "${params.vep_fasta}") \\
        --assembly GRCh38 --species homo_sapiens \\
        --force_overwrite --no_stats

  # index the output for downstream tools
  # tabix -p vcf ${sample}.vep.vcf.gz
  """
}

process LeanReport {
    tag "$sample"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(exon_cov), path(r1r2), path(frstrand)

    output:
    tuple val(sample), path("${sample}_variants_lean.xlsx")

    script:
    """
    python ${params.scriptdir}/generate_lean_report.py \
        $vcf \
        $exon_cov \
        $r1r2 \
        $frstrand \
        ${sample}_variants_lean.xlsx
    """
}

// Workflow
workflow {

    bed_ch = Channel.of(file(params.bed))

    // Step 1: VCF processing
    BedFilterVCF(sample_ch, bed_ch)
    NormalizeVCF(BedFilterVCF.out)
    FilterVCF(NormalizeVCF.out)
    AddVAF(FilterVCF.out)

    // Step 2: Optional VEP annotation
    vep_ch = params.run_vep ? VEP_Annotate(AddVAF.out) : AddVAF.out

    // Step 3: BAM processing
    BedFilterBAM(sample_ch, bed_ch)
    bam_sample_ch = BedFilterBAM.out.map { sample, vcf, bam, bai -> tuple(sample, bam, bai) }

    CoverageSummary(bam_sample_ch.map { sample, bam, bai -> tuple(sample, bam) }, bed_ch)
    R1R2Ratio(bam_sample_ch, bed_ch)
    ForwardReverseRatio(bam_sample_ch, bed_ch)

    // Step 4: Join all outputs by sample name
    lean_input_ch = vep_ch
        .join(CoverageSummary.out)
        .join(R1R2Ratio.out)
        .join(ForwardReverseRatio.out)

    // Step 5: Generate lean report
    LeanReport(lean_input_ch)
}
