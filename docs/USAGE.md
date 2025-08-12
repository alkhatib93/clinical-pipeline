# Clinical Pipeline Usage Guide

## Overview

The Clinical Pipeline is a Nextflow-based workflow for comprehensive variant analysis and quality control of genomic data. It performs VCF filtering, normalization, annotation, and generates quality control reports.

## Quick Start

### Basic Usage

```bash
# Run with default parameters
nextflow run main.nf

# Run with custom parameters
nextflow run main.nf \
    --samplesheet data/samplesheet.csv \
    --bed data/merged_output.bed \
    --outdir results
```

### Required Input Files

#### 1. Samplesheet (CSV format)
Create a CSV file with the following columns:
- `sample`: Sample identifier
- `vcf`: Path to VCF file
- `bam`: Path to BAM file

Example:
```csv
sample,vcf,bam
HG003,/path/to/HG003.vcf.gz,/path/to/HG003.bam
HG004,/path/to/HG004.vcf.gz,/path/to/HG004.bam
```

#### 2. BED File
A BED file defining regions of interest for analysis:
```bed
chr1    1000000 1001000 exon1
chr1    2000000 2001000 exon2
chr2    500000  501000  exon3
```

## Pipeline Steps

### VCF Processing Pipeline
1. **BedFilterVCF**: Filter VCF by BED regions
2. **NormalizeVCF**: Split multi-allelic variants
3. **FilterVCF**: Apply depth and quality filters
4. **AddVAF**: Add VAF (Variant Allele Frequency) tags
5. **VEP_Annotate**: Annotate variants with VEP and plugins

### BAM Processing Pipeline
6. **BedFilterBAM**: Filter BAM by BED regions
7. **CoverageSummary**: Calculate coverage statistics
8. **R1R2Ratio**: Calculate R1/R2 read pair ratios
9. **ForwardReverseRatio**: Calculate strand balance

### Reporting
10. **LeanReport**: Generate comprehensive variant report

## Parameters

### Input/Output Parameters
- `--samplesheet`: Path to samplesheet CSV file (default: `data/samplesheet.csv`)
- `--bed`: Path to BED file (default: `data/merged_output.bed`)
- `--outdir`: Output directory (default: `results`)

### VEP Parameters
- `--vep_cache`: Path to VEP cache directory (default: `data/vep_cache`)
- `--vep_plugins`: Path to VEP plugins directory (default: `data/vep_plugins`)

### Filtering Parameters
- `--min_depth`: Minimum read depth (default: 20)
- `--min_qual`: Minimum quality score (default: 30)

## Output Structure

```
results/
├── vcf/                    # VCF files
│   ├── *.bed_filtered.vcf.gz
│   ├── *.normalized.vcf.gz
│   ├── *.filtered.vcf.gz
│   ├── *.vaf_added.vcf.gz
│   └── *.vep_annotated.vcf.gz
├── qc/                     # Quality control files
│   ├── *.bed_filtered.bam
│   ├── *.bed_filtered.bam.bai
│   ├── *_coverage_summary.sorted.txt
│   ├── *_r1r2_per_exon.tsv
│   └── *_frstrand_per_exon.tsv
└── reports/                # Final reports
    └── *_variants_lean.xlsx
```

## VEP Annotation

The pipeline includes comprehensive VEP annotation with the following plugins:

### Core Annotations
- **REVEL**: Rare Exome Variant Ensemble Learner scores
- **SpliceAI**: Splice site prediction scores
- **ClinVar**: Clinical variant database
- **CADD**: Combined Annotation Dependent Depletion scores
- **dbNSFP**: Database for nonsynonymous SNPs' functional predictions

### Required VEP Data
Ensure the following files are available in your VEP plugins directory:
- `revel_all_chromosomes.tsv.gz`
- `spliceai_scores.raw.snv.hg38.vcf.gz`
- `spliceai_scores.raw.indel.hg38.vcf.gz`
- `clinvar.vcf.gz`
- `whole_genome_SNVs.tsv.gz`
- `dbNSFP4.3a_grch38.gz`

## Quality Control Metrics

### Coverage Analysis
- Coverage at 20x, 30x, 50x, and 100x thresholds
- Per-region coverage statistics

### Read Pair Analysis
- R1/R2 read pair ratios per region
- Forward/reverse strand balance

### Variant Quality
- Depth and quality filtering
- VAF calculations
- Comprehensive annotation scores

## Troubleshooting

### Common Issues

1. **VEP Cache Not Found**
   ```
   Error: VEP cache directory not found
   ```
   Solution: Ensure VEP cache is properly installed and path is correct

2. **BAM Index Missing**
   ```
   Error: Could not retrieve index file
   ```
   Solution: Ensure BAM files are indexed with `samtools index`

3. **Memory Issues**
   ```
   Error: Process killed due to memory limit
   ```
   Solution: Increase memory allocation in `nextflow.config`

### Debug Mode
Run with debug information:
```bash
nextflow run main.nf -debug
```

## Performance Optimization

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

## Citation

If you use this pipeline in your research, please cite:
- Nextflow: Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- VEP: McLaren, W. et al. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122.

## Support

For issues and questions:
1. Check the troubleshooting section above
2. Review the Nextflow logs in `.nextflow.log`
3. Open an issue on the GitHub repository
