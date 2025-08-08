# Clinical Pipeline

A **Nextflow pipeline** for clinical variant analysis, combining variant filtering, quality control, and annotation into a streamlined workflow.

## 🧬 Overview

This pipeline performs:

1. **BED Filtering**: Keep only variants in predefined regions of interest.
2. **VCF Normalization**: Split multiallelic variants and left-align indels.
3. **Quality Filtering**: Keep variants with `DP ≥ 20` and `QUAL ≥ 30`.
4. **VAF Annotation**: Annotate VAF and filter based on thresholds.
5. **Coverage Analysis**: Calculate per-base and per-region coverage.
6. **Read Balance QC**: R1/R2 and Forward/Reverse strand ratios.
7. **Lean Reporting**: Produce an Excel report with tab-separated QC summaries.

---

## ⚙️ Requirements

- Nextflow (v20+)
- Python ≥3.7
- `bcftools`, `samtools`, `bedtools`, `tabix`
- Python packages: `pandas`, `cyvcf2`, `openpyxl`

---

## 📁 Project Structure

```
clinical-pipeline/
├── main.nf                 # Main Nextflow script
├── nextflow.config         # Pipeline configuration
├── nextflow_schema.json    # Parameter schema (optional)
├── README.md              # You are here
├── CHANGELOG.md           # Version history
├── data/                  # Input data
│   ├── samplesheet.csv    # Sample definitions
│   └── regions.bed        # Regions of interest
├── scripts/               # Helper Python scripts
├── results/               # Output folder (auto-created)
│   ├── vcf/               # Filtered/normalized VCFs
│   ├── qc/                # Coverage and balance summaries
│   └── reports/           # Final lean variant reports
└── docs/                  # Documentation (usage, etc.)
    └── USAGE.md
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
  --bed data/regions.bed \
  --outdir results \
  --scriptdir scripts
```

### Parameter Summary
| Parameter       | Default                 | Description                          |
|----------------|-------------------------|--------------------------------------|
| `--samplesheet`| `data/samplesheet.csv`  | CSV with sample, vcf, and bam paths  |
| `--bed`        | `data/regions.bed`      | BED file of regions to filter        |
| `--outdir`     | `results`               | Output directory                     |
| `--scriptdir`  | `scripts`               | Directory for helper Python scripts  |

---

## 📤 Output Overview

### 🔍 VCF Outputs (`results/vcf/`)
- `{sample}.bed_filtered.vcf.gz`: Only variants in BED regions
- `{sample}.normalized.vcf.gz`: Left-aligned, multiallelics split
- `{sample}.filtered.vcf.gz`: Variants with `DP ≥ 20`, `QUAL ≥ 30`
- `{sample}.vaf_added.vcf.gz`: VCF with `VAF` added to FORMAT field

### 📊 QC Summaries (`results/qc/`)
- `{sample}_coverage_summary.txt`: % of base positions covered at 20x, 30x, 50x, 100x
- `{sample}_r1r2_per_exon.tsv`: Read 1 / Read 2 balance per exon
- `{sample}_frstrand_per_exon.tsv`: Forward/Reverse strand balance

### 📓 Lean Reports (`results/reports/`)
Excel report with multiple tabs:
- **High Confidence Variants**: `(GT=0/1 & 0.35≤VAF≤0.65)` or `(GT=1/1 & VAF≥0.90)`
- **Intermediate Confidence**: VAF between 0.1–0.34 or 0.65–0.89
- **Low Confidence**: VAF < 0.1
- **Coverage Tab** *(optional)*: Summary of regions with low exon coverage

### 📈 Pipeline Reports
- `pipeline_report.html`: DAG and summary
- `timeline_report.html`: Step timing
- `trace.txt`: Detailed process trace

---

## 🛠 Example Setup
```bash
mkdir -p data

# Create example samplesheet
cat <<EOF > data/samplesheet.csv
sample,vcf,bam
HG003,data/HG003.vcf.gz,data/HG003.bam
EOF

# Add your BED file
cp my_exons.bed data/regions.bed

# Run the pipeline
nextflow run main.nf --scriptdir scripts
```

---

## 🔧 Configuration

The `nextflow.config` file defines:
- Executor settings (e.g., local, SLURM, Docker)
- Default resources per process
- Parameters and defaults
- PublishDir settings for each output

---

## 🧪 Troubleshooting

- ❌ **Missing Index**: Make sure BAMs are indexed (`.bai` present)
- ❌ **Permissions**: Ensure read/write access to all files
- ❌ **VAF is NA**: Check that VCF has AD or VAF FORMAT fields
- 💡 **Not enough variants?** Try adjusting filters or BED content

---

## ➕ Extending the Pipeline

You can easily add:
- Additional annotations (e.g., dbNSFP, ClinVar)
- New summary tabs in the report script
- Plots or interactive HTML reports

---

## 📄 License

This pipeline is provided under MIT License for research and clinical prototyping.

---

## 👨‍🔬 Contact
For questions or contributions, open an issue or contact [your.name@domain.org](mailto:your.name@domain.org).

---

Happy variant hunting! 🧬✨
