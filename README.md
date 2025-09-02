# Clinical Pipeline

A **Nextflow pipeline** for clinical variant analysis, combining variant filtering, quality control, and comprehensive annotation into a streamlined workflow.

## 🧬 Overview

This pipeline performs comprehensive clinical variant analysis:

1. **BED Filtering**: Keep only variants in predefined regions of interest.
2. **VCF Normalization**: Split multiallelic variants and left-align indels.
3. **Quality Filtering**: Keep variants with `DP ≥ 20` and `QUAL ≥ 30`.
4. **VAF Annotation**: Add VAF (Variant Allele Frequency) tags.
5. **VEP Annotation**: Comprehensive variant annotation with pathogenicity scores and clinical databases.
6. **Enhanced Coverage Analysis**: Calculate per-base and per-region coverage at multiple thresholds (10x, 20x, 30x, 50x, 100x) with gap analysis.
7. **Comprehensive QC Metrics**: BAM statistics, read balance analysis, sex determination, and coverage gaps.
8. **Read Balance QC**: R1/R2 read pair ratios and Forward/Reverse strand balance per region.
9. **Lean Reporting**: Produce a comprehensive Excel report with variant stratification by confidence and clinical significance.

---

## ⚙️ Requirements

- Nextflow (v20+)
- Python ≥3.7
- `bcftools`, `samtools`, `bedtools`, `tabix`, `vep`, `mosdepth`
- Python packages: `pandas`, `openpyxl`, `cyvcf2`
- Docker (for VEP annotation)

---

## 📁 Project Structure

```
clinical-pipeline/
├── main.nf                 # Main Nextflow script
├── nextflow.config         # Pipeline configuration
├── nextflow_schema.json    # Parameter schema
├── README.md              # You are here
├── CHANGELOG.md           # Version history
├── data/                  # Input data
│   ├── samplesheet.csv    # Sample definitions
│   ├── merged_output.bed  # Regions of interest
│   ├── vep_cache/         # VEP cache directory
│   └── vep_plugins/       # VEP plugins directory
├── bin/                   # Helper scripts
│   └── test_setup.sh      # Setup validation script
├── scripts/               # Python scripts
│   ├── generate_lean_report.py      # Report generation script
│   ├── generate_lean_report_org.py  # Enhanced report generation
│   └── summarize_mosdepth.py        # Coverage summary script
├── results/               # Output folder (auto-created)
│   ├── vcf/               # Filtered/normalized/annotated VCFs
│   ├── qc/                # Coverage and balance summaries
│   └── reports/           # Final lean variant reports
└── docs/                  # Documentation
    └── USAGE.md           # Detailed usage guide
```

---

## 🧾 Input Files

### `samplesheet.csv`
CSV file with sample name, VCF path, and BAM path:
```csv
sample,vcf,bam
HG003,data/HG003.vcf.gz,data/HG003.bam
```

### `regions.bed`
BED file with three columns: `chrom`, `start`, `end`. Example:
```
1   17018701   17018978
1   17022587   17022750
```

### VEP Data Requirements
The pipeline requires VEP cache and reference data:

#### VEP Cache Setup
The pipeline supports both Docker and local VEP installations:

**Option 1: Docker (Recommended)**
```bash
# Docker is enabled by default
nextflow run main.nf --use_docker true
```

**Option 2: Local VEP Installation**
```bash
# Use the provided setup script
./scripts/setup_vep_cache.sh --cache-dir ~/vep_data --version 109 --assembly GRCh38

# Run pipeline with local VEP
nextflow run main.nf --use_docker false --vep_cache ~/vep_data --vep_fasta ~/vep_data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

#### VEP Cache Requirements
- Assembly: GRCh38
- Cache directory: Contains VEP cache files for human genome
- Reference fasta: `Homo_sapiens.GRCh38.dna.toplevel.fa.gz` with index

#### VEP Plugins (Optional)
For enhanced annotation, you can add plugin files to your VEP cache:
- `revel_all_chromosomes.tsv.gz` - REVEL scores
- `spliceai_scores.raw.snv.hg38.vcf.gz` - SpliceAI SNV scores
- `spliceai_scores.raw.indel.hg38.vcf.gz` - SpliceAI indel scores
- `clinvar.vcf.gz` - ClinVar database
- `whole_genome_SNVs.tsv.gz` - CADD scores
- `dbNSFP4.3a_grch38.gz` - dbNSFP database

---

## 🚀 Running the Pipeline

### Minimal Run
```bash
nextflow run main.nf
```

### Custom Parameters
```bash
nextflow run main.nf \
  --samplesheet data/samplesheet.csv \
  --bed data/merged_output.bed \
  --outdir results \
  --vep_cache ~/vep_data \
  --vep_fasta ~/vep_data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

### Parameter Summary
| Parameter       | Default                 | Description                          |
|----------------|-------------------------|--------------------------------------|
| `--samplesheet`| `data/samplesheet.csv`  | CSV with sample, vcf, and bam paths  |
| `--bed`        | `data/merged_output.bed`| BED file of regions to filter        |
| `--outdir`     | `results`               | Output directory                     |
| `--run_vep`    | `true`                  | Enable VEP annotation                |
| `--use_docker` | `true`                  | Use Docker for VEP (vs local install)|
| `--vep_cache`  | `~/vep_data`            | VEP cache directory                  |
| `--vep_fasta`  | `~/vep_data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz` | Reference fasta |

---

## 📤 Output Overview

### 🔍 VCF Outputs (`results/vcf/`)
- `{sample}.bed_filtered.vcf.gz`: Only variants in BED regions
- `{sample}.normalized.vcf.gz`: Left-aligned, multiallelics split
- `{sample}.filtered.vcf.gz`: Variants with `DP ≥ 20`, `QUAL ≥ 30`
- `{sample}.vaf_added.vcf.gz`: VCF with VAF tags added
- `{sample}.vep_annotated.vcf.gz`: VCF with comprehensive VEP annotations

### 📊 Enhanced QC Summaries (`results/qc/`)
- `{sample}_coverage_summary.sorted.txt`: % of base positions covered at 20x, 30x, 50x, 100x
- `{sample}_coverage_summary.overall.txt`: Overall coverage statistics with thresholds
- `{sample}_r1r2_per_exon.tsv`: Read 1 / Read 2 balance per exon
- `{sample}_frstrand_per_exon.tsv`: Forward/Reverse strand balance
- `{sample}_flagstat.txt`: BAM alignment statistics
- `{sample}_stats.txt`: Detailed BAM statistics
- `{sample}_bcftools_stats.txt`: VCF statistics (Ti/Tv ratios, etc.)
- `{sample}_sex_check.txt`: Sex determination from read depth ratios
- `{sample}.acmg.depth.txt`: Per-base depth within ACMG regions
- `{sample}.acmg_gaps_lt20.bed`: Coverage gaps below 20x
- `{sample}.acmg_gaps_lt30.bed`: Coverage gaps below 30x
- `{sample}.mosdepth.summary.txt`: Mosdepth coverage summary
- `{sample}.regions.bed.gz`: Per-region coverage data
- `{sample}.thresholds.bed.gz`: Coverage at multiple thresholds

### 📓 Enhanced Lean Reports (`results/reports/`)
Comprehensive Excel report with multiple analysis sheets:

**Main Variant Sheets:**
- **All Variants**: Complete variant list with all annotations and QC metrics
- **Unique HighConf Rare**: High-confidence rare variants (VAF 35-65% for het, ≥90% for hom, gnomAD AF <1%)
- **Population HighConf**: High-confidence common variants (VAF 35-65% for het, ≥90% for hom, gnomAD AF ≥1%)
- **Medium Confidence**: Variants with VAF 20-34% or 65-89%
- **Low Confidence**: Variants with VAF <20%
- **ClinVar Pathogenic**: Variants with ClinVar pathogenic annotations

**Summary Sheets:**
- **Coverage Summary**: Statistical summary of coverage metrics (10x, 20x, 30x, 50x, 100x)
- **Gene Summary**: Per-gene variant counts, VAF ranges, and consequence summaries

**Enhanced Variant Information Included:**
- **Basic Info**: Chromosome, position, reference/alternate alleles, gene, transcript
- **Annotations**: VEP consequences, impact, HGVS notation (cDNA/protein), ClinVar status
- **Quality Metrics**: Depth, VAF, zygosity, coverage at multiple thresholds
- **Read Balance**: R1/R2 read pair ratios, forward/reverse strand balance
- **Population Data**: gnomAD allele frequencies
- **Pathogenicity Scores**: REVEL, SpliceAI scores (when available)
- **Coverage Gaps**: ACMG region coverage gaps below 20x and 30x thresholds

### 📈 Pipeline Reports
- `pipeline_report.html`: DAG and summary
- `timeline_report.html`: Step timing
- `trace.txt`: Detailed process trace

---

## 🧬 VEP Annotation Features

The pipeline includes comprehensive VEP annotation with:

### Core Annotations
- **Gene & Transcript**: Gene symbols, transcript IDs, and feature information
- **Consequences**: Variant effects (missense, nonsense, splice site, etc.)
- **Impact**: Severity assessment (HIGH, MODERATE, LOW, MODIFIER)
- **HGVS Notation**: cDNA (HGVSc) and protein (HGVSp) change descriptions
- **Population Data**: gnomAD allele frequencies (genome and exome)
- **Pathogenicity Scores**: REVEL and SpliceAI scores (when available)
- **Clinical Data**: ClinVar clinical significance and star ratings

### Annotation Output
Each variant includes:
- **Functional Impact**: Detailed consequence and impact predictions
- **Molecular Changes**: Precise cDNA and protein change descriptions
- **Population Frequency**: gnomAD allele frequencies for variant assessment
- **Pathogenicity Prediction**: REVEL scores for missense variants
- **Splicing Impact**: SpliceAI scores for splice site predictions
- **Clinical Significance**: ClinVar annotations and confidence levels

---

## 🎯 Enhanced Variant Confidence Stratification

The pipeline automatically stratifies variants by confidence level based on VAF and population frequency:

### **High Confidence Variants**
- **Heterozygous**: VAF 35-65% (expected range for true heterozygotes)
- **Homozygous**: VAF ≥90% (expected for true homozygotes)
- **Rare**: gnomAD AF <1% (likely disease-relevant)
- **Common**: gnomAD AF ≥1% (population variants)

### **Medium Confidence Variants**
- VAF 20-34% or 65-89% (borderline ranges)
- May represent true variants with technical artifacts
- Require additional validation

### **Low Confidence Variants**
- VAF <20% (likely sequencing artifacts)
- Often filtered out in clinical analysis
- May represent somatic variants or contamination

### **Clinical Prioritization**
- **ClinVar Pathogenic**: Variants with established clinical significance
- **Rare High-Confidence**: Primary candidates for clinical interpretation
- **Population Variants**: Common variants for reference

---

## 🔬 Enhanced Quality Control Features

### **Multi-Threshold Coverage Analysis**
- Coverage analysis at 10x, 20x, 30x, 50x, and 100x thresholds
- Gap identification below critical coverage levels (20x, 30x)
- Per-region coverage statistics for clinical interpretation

### **Comprehensive BAM QC**
- Alignment statistics and quality metrics
- Read pair balance analysis
- Strand balance assessment
- Sex determination from chromosome coverage ratios

### **VCF Quality Assessment**
- Transition/transversion ratios
- Heterozygosity/homozygosity ratios
- Quality score distributions
- Depth and VAF distributions

---

## 🛠 Example Setup
```bash
mkdir -p data/vep_cache data/vep_plugins

# Create example samplesheet
cat <<EOF > data/samplesheet.csv
sample,vcf,bam
HG003,data/HG003.vcf.gz,data/HG003.bam
EOF

# Add your BED file
cp my_exons.bed data/merged_output.bed

# Set up VEP data (see VEP documentation for details)
# Download and install VEP cache and plugins

# Run the pipeline
nextflow run main.nf
```

---

## 🔧 Performance Optimization

### Resource Configuration
Modify `nextflow.config` to adjust resource allocation:
```groovy
process {
    cpus = 4
    memory = '8 GB'
    time = '2h'
}
```

### Parallel Execution
The pipeline automatically parallelizes across samples. Use `-profile` for different execution environments:
```bash
nextflow run main.nf -profile docker
nextflow run main.nf -profile singularity
nextflow run main.nf -profile conda
```

---

## 📚 Documentation

For detailed usage instructions, troubleshooting, and advanced configuration, see:
- [USAGE.md](docs/USAGE.md) - Comprehensive usage guide
- [CHANGELOG.md](CHANGELOG.md) - Version history

---

## 🐛 Troubleshooting

Common issues and solutions:

1. **VEP Cache Not Found**: 
   - Ensure VEP cache is properly installed using `./scripts/setup_vep_cache.sh`
   - Check that `--vep_cache` and `--vep_fasta` paths are correct
   - For Docker: ensure Docker is running and the image is available

2. **Docker Issues**:
   - If Docker fails, try `--use_docker false` with local VEP installation
   - Ensure Docker daemon is running: `docker ps`
   - Pull the VEP image: `docker pull ensemblorg/ensembl-vep:latest`

3. **BAM Index Missing**: Index BAM files with `samtools index`

4. **Memory Issues**: Increase memory allocation in `nextflow.config`

5. **VEP Process Fails**:
   - Check VEP cache directory exists and contains required files
   - Verify reference fasta file exists and is indexed
   - Ensure input VCF files are properly formatted and indexed

6. **Coverage Analysis Issues**:
   - Ensure mosdepth is installed: `conda install -c bioconda mosdepth`
   - Check BED file format and coordinates
   - Verify BAM files have proper indexing

For detailed troubleshooting, see [USAGE.md](docs/USAGE.md).

---

## 📖 Citation

If you use this pipeline in your research, please cite:
- **Nextflow**: Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **VEP**: McLaren, W. et al. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122.
- **REVEL**: Ioannidis, N. M. et al. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense variants. The American Journal of Human Genetics, 99(4), 877-885.
- **SpliceAI**: Jaganathan, K. et al. (2019). Predicting splicing from primary sequence with deep learning. Cell, 176(3), 535-548.
- **Mosdepth**: Pedersen, B. S. & Quinlan, A. R. (2018). Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867-868.

---

Happy variant hunting! 🧬✨
