# Clinical Pipeline

A Nextflow pipeline for clinical variant calling with basic filtering and quality control.

## Overview

This simplified pipeline performs the following steps:
1. **BED Filtering**: Filter VCF files to regions of interest defined in a BED file
2. **VCF Normalization**: Split multi-allelic variants
3. **Quality Filtering**: Filter variants by depth (≥20x) and quality (≥30)
4. **Coverage Analysis**: Generate coverage summaries for regions of interest

## Requirements

- Nextflow (version 20.0 or later)
- bcftools
- bedtools
- samtools
- tabix

## Project Structure

```
clinical-pipeline/
├── main.nf                 # Main pipeline script
├── nextflow.config         # Pipeline configuration
├── README.md              # This file
├── data/                  # Input data directory
│   ├── samplesheet.csv    # Sample information
│   └── regions.bed        # BED file with regions of interest
├── results/               # Output directory (created by pipeline)
│   ├── vcf/              # VCF files
│   └── qc/               # Quality control files
└── bin/                   # Helper scripts (optional)
```

## Input Files

### Samplesheet
Create a CSV file in `data/samplesheet.csv` with the following columns:
```csv
sample,vcf,bam
sample1,/path/to/sample1.vcf.gz,/path/to/sample1.bam
sample2,/path/to/sample2.vcf.gz,/path/to/sample2.bam
```

### BED File
Place your BED file as `data/regions.bed` defining regions of interest (chromosome, start, end).

## Usage

### Basic Run
```bash
nextflow run main.nf
```

### With Custom Parameters
```bash
nextflow run main.nf \
    --samplesheet data/my_samples.csv \
    --bed data/my_regions.bed \
    --outdir my_results
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `data/samplesheet.csv` | Path to samplesheet CSV file |
| `--bed` | `data/regions.bed` | Path to BED file with regions of interest |
| `--outdir` | `results` | Output directory |

## Output

The pipeline generates the following outputs:

### VCF Files (`results/vcf/` directory)
- `{sample}.bed_filtered.vcf.gz`: VCF filtered to BED regions
- `{sample}.normalized.vcf.gz`: Normalized VCF with split multi-allelics
- `{sample}.filtered.vcf.gz`: Final filtered VCF (depth ≥20x, quality ≥30)

### QC Files (`results/qc/` directory)
- `{sample}_coverage_summary.txt`: Coverage summary per region

### Reports
- `pipeline_report.html`: Pipeline execution report
- `timeline_report.html`: Timeline of pipeline execution
- `trace.txt`: Detailed execution trace

## Example Setup

```bash
# Create data directory
mkdir -p data

# Create samplesheet
echo "sample,vcf,bam" > data/samplesheet.csv
echo "NA12878,/data/NA12878.vcf.gz,/data/NA12878.bam" >> data/samplesheet.csv

# Copy your BED file
cp /path/to/your/regions.bed data/regions.bed

# Run pipeline
nextflow run main.nf
```

## Configuration

The pipeline uses `nextflow.config` for configuration. You can modify:
- Resource requirements (CPU, memory, time)
- Error handling strategy
- Output directory structure
- Executor settings

## Troubleshooting

1. **Missing tools**: Ensure all required tools (bcftools, bedtools, etc.) are installed and in PATH
2. **File permissions**: Ensure read access to input files and write access to output directory
3. **Memory issues**: Adjust memory requirements in `nextflow.config` if needed

## Extending the Pipeline

To add more steps, you can:
1. Add new processes to `main.nf`
2. Modify the workflow section to include new steps
3. Update the configuration as needed

## License

This pipeline is provided as-is for educational and research purposes. 