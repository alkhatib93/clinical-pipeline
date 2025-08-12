# Clinical Pipeline

A **Nextflow pipeline** for clinical variant analysis, combining variant filtering, quality control, and comprehensive annotation into a streamlined workflow.

## 🧬 Overview

This pipeline performs:

1. **BED Filtering**: Keep only variants in predefined regions of interest.
2. **VCF Normalization**: Split multiallelic variants and left-align indels.
3. **Quality Filtering**: Keep variants with `DP ≥ 20` and `QUAL ≥ 30`.
4. **VAF Annotation**: Add VAF (Variant Allele Frequency) tags.
5. **VEP Annotation**: Comprehensive variant annotation with REVEL, SpliceAI, ClinVar, CADD, and dbNSFP.
6. **Coverage Analysis**: Calculate per-base and per-region coverage.
7. **Read Balance QC**: R1/R2 and Forward/Reverse strand ratios.
8. **Lean Reporting**: Produce an Excel report with comprehensive variant summaries.

---

## ⚙️ Requirements

- Nextflow (v20+)
- Python ≥3.7
- `bcftools`, `samtools`, `bedtools`, `tabix`, `vep`
- Python packages: `pandas`, `openpyxl`

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
│   └── generate_lean_report.py  # Report generation script
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
The pipeline requires VEP cache and plugin data:

#### VEP Cache
- Assembly: GRCh38
- Cache directory: `data/vep_cache/`

#### VEP Plugins
Required plugin files in `data/vep_plugins/`:
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
  --vep_cache data/vep_cache \
  --vep_plugins data/vep_plugins
```

### Parameter Summary
| Parameter       | Default                 | Description                          |
|----------------|-------------------------|--------------------------------------|
| `--samplesheet`| `data/samplesheet.csv`  | CSV with sample, vcf, and bam paths  |
| `--bed`        | `data/merged_output.bed`| BED file of regions to filter        |
| `--outdir`     | `results`               | Output directory                     |
| `--vep_cache`  | `data/vep_cache`        | VEP cache directory                  |
| `--vep_plugins`| `data/vep_plugins`      | VEP plugins directory                |
| `--min_depth`  | `20`                    | Minimum read depth for filtering     |
| `--min_qual`   | `30`                    | Minimum quality score for filtering  |

---

## 📤 Output Overview

### 🔍 VCF Outputs (`results/vcf/`)
- `{sample}.bed_filtered.vcf.gz`: Only variants in BED regions
- `{sample}.normalized.vcf.gz`: Left-aligned, multiallelics split
- `{sample}.filtered.vcf.gz`: Variants with `DP ≥ 20`, `QUAL ≥ 30`
- `{sample}.vaf_added.vcf.gz`: VCF with VAF tags added
- `{sample}.vep_annotated.vcf.gz`: VCF with comprehensive VEP annotations

### 📊 QC Summaries (`results/qc/`)
- `{sample}_coverage_summary.sorted.txt`: % of base positions covered at 20x, 30x, 50x, 100x
- `{sample}_r1r2_per_exon.tsv`: Read 1 / Read 2 balance per exon
- `{sample}_frstrand_per_exon.tsv`: Forward/Reverse strand balance

### 📓 Lean Reports (`results/reports/`)
Excel report with multiple sheets:
- **Variants**: All variants with annotations and QC metrics
- **Summary**: Statistical summary of variants and coverage
- **Quality_Metrics**: Per-region quality metrics

### 📈 Pipeline Reports
- `pipeline_report.html`: DAG and summary
- `timeline_report.html`: Step timing
- `trace.txt`: Detailed process trace

---

## 🧬 VEP Annotation Features

The pipeline includes comprehensive VEP annotation with:

### Core Annotations
- **REVEL**: Rare Exome Variant Ensemble Learner scores for pathogenicity prediction
- **SpliceAI**: Splice site prediction scores for splicing impact
- **ClinVar**: Clinical variant database for clinical significance
- **CADD**: Combined Annotation Dependent Depletion scores
- **dbNSFP**: Database for nonsynonymous SNPs' functional predictions

### Annotation Output
Each variant includes:
- Consequence predictions
- Impact assessments
- Gene and transcript information
- Protein change predictions
- Pathogenicity scores
- Clinical significance

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

1. **VEP Cache Not Found**: Ensure VEP cache is properly installed and path is correct
2. **BAM Index Missing**: Index BAM files with `samtools index`
3. **Memory Issues**: Increase memory allocation in `nextflow.config`

For detailed troubleshooting, see [USAGE.md](docs/USAGE.md).

---

## 📖 Citation

If you use this pipeline in your research, please cite:
- **Nextflow**: Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **VEP**: McLaren, W. et al. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122.
- **REVEL**: Ioannidis, N. M. et al. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense variants. The American Journal of Human Genetics, 99(4), 877-885.
- **SpliceAI**: Jaganathan, K. et al. (2019). Predicting splicing from primary sequence with deep learning. Cell, 176(3), 535-548.

---

Happy variant hunting! 🧬✨
