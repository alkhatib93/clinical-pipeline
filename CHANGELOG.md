# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-01-01

### Added
- Initial release of clinical pipeline
- BED filtering of VCF files
- VCF normalization (split multi-allelics)
- Quality filtering (depth ≥20x, quality ≥30)
- Coverage analysis
- Basic reporting and timeline generation

### Features
- 4 main processes: BedFilterVCF, NormalizeVCF, FilterVCF, CoverageSummary
- Configurable parameters for input files and output directories
- Comprehensive documentation and examples
- Nextflow configuration with resource management 