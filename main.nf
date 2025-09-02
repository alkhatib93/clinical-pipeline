nextflow.enable.dsl=2

// Parameters
params.samplesheet = "${workflow.projectDir}/data/samplesheet.csv"
params.bed         = "${workflow.projectDir}/data/merged_output.bed"
params.outdir      = "results"

params.vcfdir      = "${params.outdir}/vcf"
params.qcdir       = "${params.outdir}/qc"
params.reportdir   = "${params.outdir}/reports"
params.scriptdir   = "${workflow.projectDir}/scripts"

// Sample Summary parameters
//params.verifybamid2_svd_prefix = "/path/to/verifybamid2/wgs"  // Update this path
//params.run_verifybamid2 = params.run_verifybamid2 ?: false
//params.run_somalier = params.run_somalier ?: false

// VEP parameters
def HOME = System.getenv('HOME') ?: '.'
params.run_vep   = params.run_vep ?: true
//params.vep_cache = params.vep_cache ?: "${HOME}/vep_data"
//params.vep_fasta = params.vep_fasta ?: "${HOME}/vep_data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
//params.vep_cache = params.vep_cache ?: "${System.getenv('HOME')}/vep_data"
//params.vep_fasta = params.vep_fasta ?: "${System.getenv('HOME')}/vep_data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

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
    bcftools view -i 'FORMAT/DP >= 20 && QUAL >= 30' $vcf -Oz -o ${sample}.filtered.vcf.gz
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
    tuple val(sample), path("${sample}_coverage_summary.sorted.txt"), path("${sample}_coverage_per_base.txt")

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

// Step 9: Mark duplicates - samtools markdup
process MarkDuplicates {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_dedup.bam"), path("${sample}_dedup.bam.bai"), path("${sample}_markdup.txt")

    script:
    """
    set -euo pipefail

    samtools sort -n -@ ${task.cpus ?: 4} -o ${sample}.namesort.bam $bam
    samtools fixmate -m -@ ${task.cpus ?: 4} ${sample}.namesort.bam ${sample}.fixmate.bam
    samtools sort -@ ${task.cpus ?: 4} -o ${sample}.positionsort.bam ${sample}.fixmate.bam

    samtools markdup -@ ${task.cpus ?: 4} \\
      -f ${sample}_markdup.stats.txt \\
      ${sample}.positionsort.bam ${sample}.dedup.bam

    samtools index -@ ${task.cpus ?: 4} ${sample}.dedup.bam

    rm -f ${sample}.namesort.bam ${sample}.fixmate.bam ${sample}.positionsort.bam
    """
}

// Step 9: Sample Summary - samtools flagstat
process SamtoolsFlagstat {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_flagstat.txt")

    script:
    """
    samtools flagstat $bam > ${sample}_flagstat.txt
    """
}

// Step 10: Sample Summary - samtools stats
process SamtoolsStats {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_stats.txt")

    script:
    """
    samtools stats $bam > ${sample}_stats.txt
    """
}

process SamtoolsDepth {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)
    path bed

    output:
    tuple val(sample),
        path("${sample}.acmg.depth.txt"),
        path("${sample}.acmg_gaps_lt20.bed"),
        path("${sample}.acmg_gaps_lt30.bed")

    script:
    """
    set -euo pipefail

    # Per-base depths within ACMG exons; -a emits zeros inside BED intervals
    # Adjust filters if you want callable bases only: e.g., add -q 20 (MAPQ) and/or -Q 20 (BQ)
    samtools depth -a -b $bed $bam > ${sample}.acmg.depth.txt

    # Make 1-bp BEDs for bases below thresholds and merge into contiguous gaps
    awk '(\$3 < 20){print \$1"\t"(\$2-1)"\t"\$2}' ${sample}.acmg.depth.txt \
      | bedtools sort -i - | bedtools merge -i - > ${sample}.acmg_gaps_lt20.bed

    awk '(\$3 < 30){print \$1"\t"(\$2-1)"\t"\$2}' ${sample}.acmg.depth.txt \
      | bedtools sort -i - | bedtools merge -i - > ${sample}.acmg_gaps_lt30.bed
    """
}

// Step 11: Sample Summary - Picard CollectAlignmentSummaryMetrics
process PicardMetrics {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_alignment_metrics.txt"), path("${sample}_insert_size_metrics.txt"), path("${sample}_insert_size_histogram.pdf")

    script:
    """
    set -euo pipefail

    picard CollectAlignmentSummaryMetrics \
        I=$bam \
        O=${sample}_alignment_metrics.txt \
        R=${params.vep_fasta}

    picard CollectInsertSizeMetrics \
        I=$bam \
        O=${sample}_insert_size_metrics.txt \
        R=${params.vep_fasta}
    """
}

// Step 1: Sample Summary - verifyBamID2 for contamination
process VerifyBamID2 {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_verifybamid2.txt")

    script:
    """
    verifyBamID2 \\
        --SVDPrefix ${params.verifybamid2_svd_prefix} \\
        --BamFile $bam \\
        --OutFile ${sample}_verifybamid2.txt
    """
}

// Step 14: Sample Summary - mosdepth for coverage metrics
process MosdepthCoverage {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai)
    path bed

    output:
    tuple val(sample),
        path("${sample}.mosdepth.summary.txt"),
        path("${sample}.regions.bed.gz"),
        path("${sample}.thresholds.bed.gz"),
        path("${sample}_coverage_summary.overall.txt")

    script:
    """
    set -euo pipefail

    # 1) Run mosdepth on WES targets with thresholds
    mosdepth --by $bed --thresholds 10,20,30,50,100 --fast-mode $sample $bam

    python ${params.scriptdir}/summarize_mosdepth.py \
        --prefix $sample \
        --summary ${sample}.mosdepth.summary.txt \
        --thresholds ${sample}.thresholds.bed.gz \
        --out ${sample}_coverage_summary.overall.txt
    """
}

// Step 15: Sample Summary - Sex check using read depth ratio
process SexCheck {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_sex_check.txt")

    script:
    """
    # Calculate read depth for X and Y chromosomes
    x_depth=\$(samtools idxstats $bam | grep -E "^X\\s" | awk '{print \$3}')
    y_depth=\$(samtools idxstats $bam | grep -E "^Y\\s" | awk '{print \$3}')
    
    # Calculate ratio
    if [ "\$y_depth" -gt 0 ]; then
        ratio=\$(echo "scale=3; \$x_depth / \$y_depth" | bc)
    else
        ratio="NA"
    fi
    
    # Determine sex based on ratio (X/Y ratio > 4 typically indicates female)
    if [ "\$ratio" != "NA" ]; then
        if (( \$(echo "\$ratio > 4" | bc -l) )); then
            sex="Female"
        else
            sex="Male"
        fi
    else
        sex="Unknown"
    fi
    
    echo "Sample: $sample" > ${sample}_sex_check.txt
    echo "X chromosome depth: \$x_depth" >> ${sample}_sex_check.txt
    echo "Y chromosome depth: \$y_depth" >> ${sample}_sex_check.txt
    echo "X/Y ratio: \$ratio" >> ${sample}_sex_check.txt
    echo "Predicted sex: \$sex" >> ${sample}_sex_check.txt
    """
}

// Step 16: Sample Summary - bcftools stats for Ti/Tv and Het/Hom ratios
process BcftoolsStats {
    tag "$sample"
    publishDir "${params.qcdir}", mode: 'copy'
    
    input:
    tuple val(sample), path(vcf)

    output:
    tuple val(sample), path("${sample}_bcftools_stats.txt")

    script:
    """
    bcftools stats $vcf > ${sample}_bcftools_stats.txt
    """
}

process VEP_Annotate {
  tag "$sample"
  publishDir "${params.vcfdir}", mode: 'copy'

  input:
  tuple val(sample), path(vcf)

  output:
  tuple val(sample), path("${sample}.vep.vcf")

  script:
  """
  set -euo pipefail

  # make a plain VCF for VEP (avoid STDIN quirks)
  if [[ "$vcf" == *.vcf.gz ]]; then
    gunzip -c "$vcf" > INPUT_FOR_VEP.vcf
  else
    cp "$vcf" INPUT_FOR_VEP.vcf
    #bgzip INPUT_FOR_VEP.vcf
  fi

  vep \
    -i INPUT_FOR_VEP.vcf \
    -o ${sample}.vep.vcf \
    --offline --cache --dir_cache /cache \
    --fasta /cache/\$(basename "${params.vep_fasta}") \
    --assembly GRCh38 --species homo_sapiens \
    --hgvs --symbol  \
    --vcf \
    --everything \
    --canonical \
    --plugin REVEL,/cache/\$(basename "${params.revel_vcf}") \
    --plugin ClinVar,/cache/\$(basename "${params.clinvar_vcf}") \
    #--custom /cache/\$(basename "${params.gnomad_vcf}"),gnomAD,vcf,exact,0,AF,AC
    #--tab \
    #--fields "Uploaded_variation,Location,Allele,IMPACT,SYMBOL,Gene,Feature,Feature_type,Consequence,HGVSc,HGVSp,Protein_position,Amino_acids,Existing_variation,CLIN_SIG,SYMBOL_SOURCE,CANONICAL" \
    #--force_overwrite --no_stats --no_progress

  #rm -f INPUT_FOR_VEP.vcf

  # ensure output exists
  #test -s "${sample}.vep.tsv"
  """
}

process LeanReport {
  tag "$sample"
  publishDir "${params.outdir}/reports", mode: 'copy'

  input:
  tuple val(sample),
        path(vcf), path(exon_cov), path(r1r2), path(frstrand),
        path(flagstat), path(stats),
        path(mosdepth_summary),
        path(sex_check),
        path(gaps20), path(gaps30),
        path(thresholds)

  output:
  tuple val(sample), path("${sample}_variants_lean.xlsx")

  script:
  """
  python ${params.scriptdir}/generate_lean_report_org.py \
    $vcf \
    $exon_cov \
    $r1r2 \
    $frstrand \
    ${sample}_variants_lean.xlsx \
    --sample-id ${sample} \
    --assay "WES" \
    --build "GRCh38" \
    --flagstat ${flagstat} \
    --stats ${stats} \
    --picard-align ${sample}.picard_alignment_metrics.txt \
    --picard-insert ${sample}.picard_insert_metrics.txt \
    --mosdepth-summary ${mosdepth_summary} \
    --acmg-thresholds ${thresholds} \
    --verifybamid ${sample}.verifybamid.selfSM \
    --sexcheck ${sex_check} \
    --sf-genes ${params.data}/acmg_sf_gene_list.txt \
    --gaps20 ${gaps20} \
    --gaps30 ${gaps30}
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

    // Step 4: Sample Summary QC processes
    //MarkDuplicates(sample_ch.map { sample, vcf, bam -> tuple(sample, bam) })
    SamtoolsFlagstat(bam_sample_ch.map { sample, bam, bai -> tuple(sample, bam) })
    SamtoolsStats(bam_sample_ch.map { sample, bam, bai -> tuple(sample, bam) })
    SamtoolsDepth(bam_sample_ch.map { sample, bam, bai -> tuple(sample, bam) }, bed_ch)
    //PicardMetrics(sample_ch.map { sample, vcf, bam -> tuple(sample, bam) })
    MosdepthCoverage(bam_sample_ch, bed_ch)
    SexCheck(bam_sample_ch.map { sample, bam, bai -> tuple(sample, bam) })
    //BcftoolsStats(vep_ch)
    BcftoolsStats(vep_ch.map { sample, vcf -> tuple(sample, vcf) })
    
    // Optional: verifyBamID2 for contamination (if enabled)
    //verifybamid2_ch = params.run_verifybamid2 ? VerifyBamID2(sample_ch.map { sample, vcf, bam -> tuple(sample, bam) }) : Channel.empty()
    exon_cov_ch = CoverageSummary.out.map { sample, summary, per_base ->tuple(sample, summary)}
    // SamtoolsDepth/ACMGGaps: split to gaps20 / gaps30 (use *annotated* gaps)
    //SamtoolsDepth.out.into { gapsA; gapsB }
    gaps20_ch = SamtoolsDepth.out.map { sample, depth, gaps20, gaps30 -> tuple(sample, gaps20)}
    gaps30_ch = SamtoolsDepth.out.map { sample, depth, gaps20, gaps30 -> tuple(sample, gaps30)}

    // MosdepthCoverage.out is: (sample, summary.txt, regions.bed.gz, thresholds.bed.gz, overall.txt)
    mosdepth_summary_ch = MosdepthCoverage.out.map { sample, summary, regions, thresholds, overall -> tuple(sample, overall)}
    thresholds_ch = MosdepthCoverage.out.map { sample, summary, regions, thresholds, overall -> tuple(sample, thresholds)}

    // Step 6: Join all outputs by sample name for lean report
    lean_input_ch = vep_ch
        .join(exon_cov_ch)
        .join(R1R2Ratio.out)
        .join(ForwardReverseRatio.out)
        .join(SamtoolsFlagstat.out)
        .join(SamtoolsStats.out)
        .join(mosdepth_summary_ch)
        .join(SexCheck.out)
        .join(gaps20_ch)
        .join(gaps30_ch)
        .join(thresholds_ch)
        //.join(BcftoolsStats.out)

    // Step 7: Generate lean report
    LeanReport(lean_input_ch)

}
