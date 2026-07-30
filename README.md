# Reconstructing-East-Asian-Data-Assimilation 

Code to reproduce figures from the paper "Reconstructing East Asian Temperatures from 1368 to 1911 Using Historical Documents, Climate Models, and Data Assimilation" by Eric Sun, Kuan-hui Elaine Lin, Wan-Ling Tseng, Pao K. Wang and Hsin-Cheng Huang ([arxiv link](http://arxiv.org/abs/2410.21790)).

This repository contains the code and prepared data used to reproduce the computational results in the paper:

> **Reconstructing East Asian Temperatures from 1368 to 1911 Using Historical Documents, Climate Models, and Data Assimilation**  
> Eric Sun, Kuan-hui Elaine Lin, Wan-Ling Tseng, Pao K. Wang, and Hsin-Cheng Huang  
> [arXiv preprint](http://arxiv.org/abs/2410.21790)

The repository includes implementations of:

- interval-censored kriging for the REACHES documentary temperature index;
- uncertainty-aware quantile mapping;
- ridge-penalized nonstationary AR(1) estimation;
- Kalman filtering and Rauch–Tung–Striebel smoothing;
- functional clustering and spatial regionalization;
- documentary-coverage and historical-population analysis; and
- held-out interval-compatibility validation.

The workflow generates Figures 2–10 in the main manuscript, Figures S1–S5 in the Supplementary Material, and the associated computational tables and diagnostic outputs.

## Repository organization

```text
.
├── Code/                   Analysis and figure-generation scripts
├── Data/                   Input data and documented prepared datasets
├── Output/                 Generated figures, tables, and intermediate files
├── renv/                   Project environment infrastructure
├── run_all.R               Master reproducibility workflow
├── renv.lock               Locked R package environment
├── output_manifest.csv     Mapping from manuscript items to code and outputs
└── README.md


## Organization

### Code 

Detailed descriptions of the code and usage instructions are available in `Code/README.md`.

### Data

Detailed descriptions of the data and usage instructions are available in `Data/README.md`.

### Figure  

The provided code generates Figures 2–10 as presented in the paper and Figures S1–S5 as presented in the supplementary.


