# Exploratory GWAS for Circadian-Calcium Signaling in Canine Cognitive Dysfunction

Analysis code for:

> *Exploratory Genome-Wide Association Study for Circadian-Calcium Signaling to Canine Cognitive Dysfunction in Companion Dogs*
> Younghan Song

This repository contains the full analysis pipeline used in the manuscript, from genotype QC through GWAS, functional annotation, linkage disequilibrium characterization, and polygenic score (PRS) construction. Scripts are numbered in pipeline order and grouped below by manuscript section.

## Data availability

This analysis uses genotype and Health and Life Experience Survey (HLES) data from the Dog Aging Project (DAP) 2022 Curated Data Open Access Release. **Raw genotype and survey data are not included in this repository** and require a separate data access request through the Dog Aging Project (https://dogagingproject.org/data-access). See the manuscript's Data availability statement for full details.

## Software

- PLINK v1.9 and v2.0
- R (version — *fill in*), packages: *fill in key packages (e.g., lme4, performance, DAVID query tools, etc.)*
- Python (version — *fill in*), for `.ipynb` notebooks
- DAVID Bioinformatics (https://davidbioinformatics.nih.gov/)

## Pipeline overview

Several scripts appear twice (a base version and a `_local` version). Both were used: the base scripts were run on a lab server, and the `_local` versions re-run the same analysis on a local machine with restricted memory/core allocation (used as an execution-environment check, not a separate analysis). Where a script name includes a MAF-related suffix (e.g., `001`, `005`, `01`), this refers to the MAF threshold tested (0.001 / 0.005 / 0.01) as part of the sensitivity analysis described in Methods B.2 (Table S5); **MAF ≥ 0.01 is the primary model reported in the main text**.

*Note: file-to-section mapping below is inferred from filenames and script order — please correct any misattributions before publishing.*

### B.1 – Study Population / data preparation
| File | Purpose |
|---|---|
| `00_datacheck.r` | Initial data checks on raw genotype/phenotype files |
| `00_Make_Bed.sh` | Build PLINK binary fileset (.bed/.bim/.fam) from raw data |

### B.2 – QC and Association Testing
| File | Purpose |
|---|---|
| `01_MAF.sh` | MAF filtering (PLINK `--maf`) |
| `02-1_Missing_Rate.sh` | Compute sample/variant missingness |
| `02_missing_threshold.r` | Derive the prespecified missingness threshold, (Q3+MAX)/2, at each candidate MAF cutoff |
| `02-2_Missing_Rate_Filtered.sh` | Apply missingness filter (`--mind`/`--geno`) |
| `02-3_MR_Filtered_Check.sh` | Verify post-filter sample/variant counts |
| `03-1_HWE.sh` | Hardy-Weinberg equilibrium filtering (PLINK `--hwe`) |
| `04-1_DuplicateID_Filter.sh` | Remove duplicate/ineligible IDs prior to regression |
| `04-2_Regression.sh` | Base association regression script (template) |
| `04-2_Regression_01_1e6.sh` / `04-2_Regression_01_1e6_local.sh` | **Primary model**: MAF ≥ 0.01 (server / local re-run) |
| `04-2_Regression_005_1e6.sh` / `04-2_Regression_005_1e6_local.sh` | Sensitivity: MAF ≥ 0.005 (server / local re-run) |
| `04-2_Regression_001_1e6.sh` / `04-2_Regression_001_1e6_local.sh` | Sensitivity: MAF ≥ 0.001 (server / local re-run) |
| `04-2_Regression_1e6_parm_local.sh` | Parameter-sweep version of the regression (local) |
| `04-2_Regression_check.sh` | Post-regression output verification |
| `05_Regression_Result.r`, `05_Regression_Result2.r`, `05_Regression_Result3.R` | Compile/summarize regression outputs across MAF×PC configurations (Table 1, Table S5) |

### B.3 – Population Structure (PCA)
| File | Purpose |
|---|---|
| `03-2_PCA.sh` | Run PCA (PLINK `--pca`) |
| `03_PCA.r` | Evaluate PCs 5–30, variance explained, PC selection (PC1–8) |

### B.4 – Covariates
| File | Purpose |
|---|---|
| `03-2_Add_Phenotype.r` | Merge CSLB/phenotype data with genotype IDs |
| `03-3_Select_Phenotype.r` | Construct covariates: age, sex (X:autosome ratio), breed dummy variables, neuter status, activity PC |

### B.6 – Exploratory Functional Annotation
| File | Purpose |
|---|---|
| `06_match_gene.ipynb` | Map exploratory-threshold SNPs to genes by intragenic overlap (Ensembl GTF) |
| `08_David_Tables.R` | Process DAVID GO/KEGG output tables |
| `08_GO_KEGG.ipynb` | GO (BP/CC/MF) and KEGG enrichment analysis (Figs. 3–6) |
| `10_barplot.r` | Generate enrichment bar/dot plots |
| `11_Biomart.r` | Gene annotation lookup (biomaRt) |

### B.7 – Linkage Disequilibrium
| File | Purpose |
|---|---|
| `09_Linkage_Disequilibrium.sh` | Regional LD calculation (PLINK `--r2`), pooled sample |
| `09_Linkage_Disequilibrium_local.sh` | Same, local re-run |
| `09_Breed_Keep_Files.R` | Filter breeds with ≥10 individuals for breed-stratified LD |
| `09_Linkage_Disequilibrium_local_BreedStratified.sh` | Within-breed LD calculation |
| `09_LD_Comparison.R` | Compare pooled vs. within-breed LD estimates (Table S3, Fig. S4) |
| `08_LD_SNP.R` | Per-SNP LD tagging variant extraction |
| `09_Regional_Plot.r`, `09_Regional_Plot2.R` | Regional association plots (Figs. 7–11) |
| `10-1_Extract_SNP.sh` | Extract tagging/candidate SNPs for downstream analysis |

### B.8 – Cross-Species Functional Context
| File | Purpose |
|---|---|
| `14_eQTL.R` | Lift-over to GRCh38 and GTEx brain-tissue eQTL query |

### B.9 – Polygenic Score (PRS) Construction and Evaluation
| File | Purpose |
|---|---|
| `12_CalculateGeno_Beta.sh` / `12_CalculateGeno_Beta_local.sh` | Extract genotype dosages and effect sizes for PRS input (server / local) |
| `13-1_Geno.sh` | Genotype extraction for PRS scoring |
| `13_PRS-1.R` | Method 1: standard additive PRS |
| `13_PRS-2.R` | Method 2: LD-block PRS |
| `13_PRS-3.R` | Method 3: correlation-adjusted (eigen-decomposition) PRS |
| `13_PRS-4.R` | Method 4: PC-based PRS |

