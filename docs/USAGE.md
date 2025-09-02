# Clinical Pipeline Usage Guide

## Overview

The Clinical Pipeline is a Nextflow-based workflow for comprehensive variant analysis and quality control of genomic data. It performs VCF filtering, normalization, annotation, and generates comprehensive quality control reports with enhanced sample summary metrics.

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

## Enhanced Pipeline Steps

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

### Enhanced Sample Summary QC Pipeline
10. **SamtoolsFlagstat**: Generate BAM alignment statistics
11. **SamtoolsStats**: Generate detailed BAM statistics
12. **SamtoolsDepth**: Calculate per-base depth and identify coverage gaps
13. **MosdepthCoverage**: Multi-threshold coverage analysis (10x, 20x, 30x, 50x, 100x)
14. **SexCheck**: Determine sex from chromosome coverage ratios
15. **BcftoolsStats**: Generate VCF quality statistics

### Reporting
16. **LeanReport**: Generate comprehensive variant report with enhanced QC metrics

## Enhanced Parameters

### Input/Output Parameters
- `--samplesheet`: Path to samplesheet CSV file (default: `data/samplesheet.csv`)
- `--bed`: Path to BED file (default: `data/merged_output.bed`)
- `--outdir`: Output directory (default: `results`)

### VEP Parameters
- `--vep_cache`: Path to VEP cache directory (default: `data/vep_cache`)
- `--vep_plugins`: Path to VEP plugins directory (default: `data/vep_plugins`)
- `--revel_vcf`: Path to REVEL scores file
- `--clinvar_vcf`: Path to ClinVar database file

### Filtering Parameters
- `--min_depth`: Minimum read depth (default: 20)
- `--min_qual`: Minimum quality score (default: 30)

### Docker Configuration
- `--use_docker`: Enable Docker for VEP annotation (default: true)

## Enhanced Output Structure

```
results/
├── vcf/                    # VCF files
│   ├── *.bed_filtered.vcf.gz
│   ├── *.normalized.vcf.gz
│   ├── *.filtered.vcf.gz
│   ├── *.vaf_added.vcf.gz
│   └── *.vep_annotated.vcf.gz
├── qc/                     # Enhanced quality control files
│   ├── *.bed_filtered.bam
│   ├── *.bed_filtered.bam.bai
│   ├── *_coverage_summary.sorted.txt
│   ├── *_coverage_summary.overall.txt
│   ├── *_r1r2_per_exon.tsv
│   ├── *_frstrand_per_exon.tsv
│   ├── *_flagstat.txt
│   ├── *_stats.txt
│   ├── *_bcftools_stats.txt
│   ├── *_sex_check.txt
│   ├── *.acmg.depth.txt
│   ├── *.acmg_gaps_lt20.bed
│   ├── *.acmg_gaps_lt30.bed
│   ├── *.mosdepth.summary.txt
│   ├── *.regions.bed.gz
│   └── *.thresholds.bed.gz
└── reports/                # Final reports
    └── *_variants_lean.xlsx
```

## Enhanced VEP Annotation

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

## Enhanced Quality Control Metrics

### Multi-Threshold Coverage Analysis
- Coverage analysis at 10x, 20x, 30x, 50x, and 100x thresholds
- Per-region coverage statistics with gap identification
- Coverage gaps below critical thresholds (20x, 30x) for clinical interpretation

### Comprehensive BAM QC
- **Alignment Statistics**: Total reads, mapped reads, duplicate rates
- **Quality Metrics**: Insert size distributions, alignment quality scores
- **Coverage Analysis**: Per-base depth within target regions
- **Gap Identification**: Contiguous regions below coverage thresholds

### Read Pair and Strand Analysis
- R1/R2 read pair ratios per region
- Forward/reverse strand balance assessment
- Read orientation patterns for quality assessment

### Sex Determination
- X/Y chromosome coverage ratios
- Automated sex prediction from read depth patterns
- Quality metrics for sex determination accuracy

### VCF Quality Assessment
- Transition/transversion ratios
- Heterozygosity/homozygosity ratios
- Quality score distributions
- Depth and VAF distributions

## Enhanced Lean Report Features

The pipeline generates comprehensive Excel reports with multiple analysis sheets:

### Variant Analysis Sheets
- **All Variants**: Complete variant list with all annotations and QC metrics
- **High Confidence Variants**: Stratified by VAF and population frequency
- **Medium Confidence Variants**: Borderline variants requiring validation
- **Low Confidence Variants**: Likely artifacts or somatic variants
- **ClinVar Pathogenic**: Variants with established clinical significance

### Quality Control Sheets
- **Coverage Summary**: Statistical summary of coverage metrics across thresholds
- **Gene Summary**: Per-gene variant counts and quality metrics
- **QC Metrics**: Comprehensive quality control statistics

### Enhanced Variant Information
- **Basic Variant Data**: Chromosome, position, alleles, gene information
- **Functional Annotations**: VEP consequences, impact, HGVS notation
- **Quality Metrics**: Depth, VAF, coverage at multiple thresholds
- **Coverage Gaps**: ACMG region coverage gaps below critical thresholds
- **Population Data**: gnomAD allele frequencies and population statistics

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

### Docker Integration
The pipeline supports Docker for VEP annotation:
```bash
# Enable Docker (default)
nextflow run main.nf --use_docker true

# Use local VEP installation
nextflow run main.nf --use_docker false
```

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

4. **Coverage Analysis Issues**
   ```
   Error: mosdepth command not found
   ```
   Solution: Install mosdepth: `conda install -c bioconda mosdepth`

5. **Docker Issues**
   ```
   Error: Docker daemon not running
   ```
   Solution: Start Docker daemon or use local VEP installation

### Debug Mode
Run with debug information:
```bash
nextflow run main.nf -debug
```

### Log Analysis
Check Nextflow logs for detailed error information:
```bash
# View recent logs
tail -f .nextflow.log

# Check specific process logs
find work/ -name "*.command.log" -exec grep -l "ERROR" {} \;
```

## Advanced Configuration

### Custom VEP Plugins
Add custom VEP plugins by modifying the VEP_Annotate process:
```groovy
process VEP_Annotate {
    // ... existing configuration ...
    
    script:
    """
    vep \
        -i INPUT_FOR_VEP.vcf \
        -o ${sample}.vep.vcf \
        --plugin CustomPlugin,/path/to/plugin \
        # ... other options ...
    """
}
```

### Coverage Threshold Customization
Modify coverage thresholds in the MosdepthCoverage process:
```groovy
process MosdepthCoverage {
    // ... existing configuration ...
    
    script:
    """
    mosdepth --by $bed --thresholds 5,10,15,20,25,30,50,100 --fast-mode $sample $bam
    """
}
```

## Citation

If you use this pipeline in your research, please cite:
- **Nextflow**: Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **VEP**: McLaren, W. et al. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122.
- **REVEL**: Ioannidis, N. M. et al. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense variants. The American Journal of Human Genetics, 99(4), 877-885.
- **SpliceAI**: Jaganathan, K. et al. (2019). Predicting splicing from primary sequence with deep learning. Cell, 176(3), 535-548.
- **Mosdepth**: Pedersen, B. S. & Quinlan, A. R. (2018). Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867-868.

## Support

For issues and questions:
1. Check the troubleshooting section above
2. Review the Nextflow logs in `.nextflow.log`
3. Check the pipeline reports in `results/`
4. Open an issue on the GitHub repository

## Version History

See [CHANGELOG.md](../CHANGELOG.md) for detailed version information and updates.
