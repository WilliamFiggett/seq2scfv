# Fork-Specific Changes

This document describes the changes made in this forked version of [`ngs-ai-org/seq2scfv`](https://github.com/ngs-ai-org/seq2scfv), originally designed as a collection of tools for analysis of full-length sequenced scFvs from in vitro display experiments.

## Purpose of the Fork

This fork was created to:
- Improve compatibility and usability through updates to the Dockerfile.
- Apply small corrections and refinements to the Python scripts.
- Address minor issues affecting reproducibility and ease of installation.

## Summary of Changes

- **Dockerfile updates**:
  - Improved compatibility for diverse environments.
  - Updated package installation steps and base image.
  - Resolved conflicting dependencies.

- **Python script corrections**:
  - Minor changes in analysis scripts.
  - Improved error handling and output formatting.
  - Enhanced code readability and maintainability.

- **Other modifications**:
  - Updated documentation as needed to reflect these changes.
  - Minor refinements to ensure reproducibility.

## Usage Differences

- The updated Dockerfile may change the installation or build steps. Please refer to the [Dockerfile](./Dockerfile) and documentation for new instructions, including how R scripts are handled.
- Python scripts may differ slightly in usage or output formatting due to corrections.

## Attribution & License

This fork is based on the original [`seq2scfv`](https://github.com/ngs-ai-org/seq2scfv) project. Please refer to the original repository and license for additional information.
