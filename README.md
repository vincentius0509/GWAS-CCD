# GWAS of canine cognitive dysfunction scores in 7,115 Dog Aging Project dogs

Analysis code for: Song Y. and Lee Y. *Calcium-response genes and shared genetic effects underlie cognitive dysfunction scores in 7,115 companion dogs*. Animal Genetics ([year], [volume], [DOI]).

- Data (summary statistics, gene and gene-set results, analysed dog and variant IDs): Zenodo [data DOI]
- This code, archived: Zenodo [code DOI]

## Data source

All analyses use the Dog Aging Project 2023 Curated Data Open Access Release. Individual-level genotype and survey data are **not** included in this repository; they are available to researchers with an approved Data Use Agreement at https://dogagingproject.org/data-access.

> Dog Aging Project. (2024). Dog Aging Project – 2023 Curated Data Open Access Release, version 1.0 [Data file and codebook]. https://app.terra.bio/ : Terra at the Broad Institute of MIT and Harvard.

This research is based on publicly available data collected by the Dog Aging Project, under U19 grant AG057377 (PI: Daniel Promislow) from the National Institute on Aging, a part of the National Institutes of Health, and by additional grants and private donations, including generous support from the Glenn Foundation for Medical Research, the Tiny Foundation at Myriad Canada, and the WoodNext Foundation. These data are housed on the Terra platform at the Broad Institute of MIT and Harvard. Use of these data does not imply collaboration with, or endorsement by, the Dog Aging Project.

## Repository contents

`scripts/` holds the 75 R and shell scripts of the analysis. Script numbers follow the order in which the scripts were run (`07-2` runs after `07-1`; R and shell scripts share one sequence). Each script reads the outputs of earlier steps from the results directory and writes its own outputs there with its step number as a prefix (for example, `04-5_gwas/`, `06-5_gsa_all.csv.gz`).

Not included:

- scripts that only draw the paper figures (Fig. 1, Fig. 2, Figs S1–S3); the results they plot are produced by the scripts listed below;
- small helper scripts that only monitored long-running jobs.

Comments were written in Korean during the analysis and have been removed; each script starts with a one-line English description. Some console messages and a few output labels remain in Korean, because changing them would change the files that later scripts read.

## Software

| Software | Version | Used for |
|---|---|---|
| PLINK | 1.90b7.1 (18 Oct 2023) | QC steps, LD, regional analyses |
| PLINK 2 | v2.00a6LM AVX2 (23 Nov 2023) | QC, PCA, KING kinship, association checks |
| GCTA | 1.95.4 | GRM, MLMA-LOCO, GREML, COJO, dominance GRM |
| REGENIE | v4.1.3 | SNP × age and SNP × body weight interaction, binary models |
| MAGMA | v1.10 | gene and gene-set analysis |
| minimap2 / samtools | 2.28 / 1.19.2 | read alignment and depth counts (10-3) |
| SRA Toolkit | 3.4.1 | download of raw reads (06-11, 10-3) |
| STAR | 2.7.2a (bundled with Cell Ranger 7.2.0) | RNA-seq quantification (06-11) |
| R | 4.4.2 | all R scripts |
| R packages | data.table 1.18.6.1, ggplot2 3.4.4, Matrix 1.6.5, mgcv 1.9.1, gridExtra 2.3, susieR 0.14.2, parallel (base) |  |
| pandoc | 2.9.2.1 | Word output of supplementary tables (14-11) |

## Before running

- **Paths.** Paths are written for the original server: the project directory `/data/02_ByProject/DogAgingProject` and the tool directory `/data/00_tools`. Replace them with your own locations, for example:

  ```bash
  sed -i 's#/data/02_ByProject/DogAgingProject#/path/to/project#g; s#/data/00_tools#/path/to/tools#g' scripts/*
  ```

- **Inputs from the Dog Aging Project:** genotype files (PLINK format) with metadata, genetic ancestry proportions, the 2023 codebook, and the survey tables used in `00_datacheck.R` (CSLB, Health and Life Experience Survey, Annual Follow-Up Survey, Brain Health Study enrolment).
- **Public inputs** (downloaded by the scripts or described in their headers): NCBI CanFam3.1 gene annotation, KEGG and Gene Ontology gene sets, GWAS Catalog gene lists, human–dog orthologues, Human Protein Atlas consensus tissue RNA data, published body-size markers (Plassais et al. 2019), SRA PRJNA800779 (DAP raw reads) and PRJNA815057.
- Long steps (GWAS of all models, GREML, read alignment) take hours to days; most shell scripts accept thread and job settings as environment variables at the top of the file.

## Scripts

### 00 Data

| Script | What it does | Paper |
|---|---|---|
| `00_datacheck.R` | Checks the DAP 2023 release files and defines the analysis sample, covariates and phenotype (most recent CSLB score) | Methods; Fig. S1 |

### 01 Variant and sample missingness

| Script | What it does | Paper |
|---|---|---|
| `01-1_Missing_Rate.sh` | Extracts the 7,166 dogs with phenotype and covariates, corrects sex, and computes unfiltered missingness and MAF | Table S2 |
| `01-1_missing_threshold.R` | Examines variant missingness distributions to choose the --geno threshold | Table S2 |
| `01-2_Geno_Filter.sh` | Applies variant missingness < 5% and recomputes missingness | Table S2 |
| `01-3_Group_Missing.sh` | Counts missing genotypes per SNP within each ancestry group | Table S2 |
| `01-3_group_missing_threshold.R` | Flags SNPs with missingness concentrated in one ancestry group (> 10% and binomial p < 1e-6) | Table S2 |
| `01-4_Group_Exclude.sh` | Removes the group-specific missingness SNPs and recomputes missingness | Table S2 |
| `01-4_mind_threshold.R` | Chooses the sample missingness (--mind) threshold | Table S2 |
| `01-4_Mind_Filter.sh` | Applies sample missingness < 5% (51 dogs removed) | Table S2 |

### 02 Minor allele frequency

| Script | What it does | Paper |
|---|---|---|
| `02-1_MAF_Filter.sh` | Re-applies variant missingness and builds MAF 0.01 / 0.05 / 0.1 sets with allele counts per ancestry group | Table S2 |
| `02-2_maf_group_threshold.R` | Compares the MAF sets on QC metrics and removes SNPs whose minor alleles are confined to one ancestry group | Table S2 |

### 03 HWE, PCA

| Script | What it does | Paper |
|---|---|---|
| `03-1_HWE.sh` | Hardy-Weinberg tests in the pooled sample and within ancestry groups | Table S2 |
| `03-2_hwe_threshold.R` | Compares pooled vs. within-group HWE and defines the exclusion rule (p < 1e-6 in >= 2 groups) | Table S2 |
| `03-3_HWE_Exclude.sh` | Removes HWE failures; final QC genotype set (4,842,736 autosomal SNPs, 7,115 dogs) | Table S2 |
| `03-4_PCA.sh` | Genotype principal components: PC-AiR (primary) and two alternatives (all unrelated dogs; ancestry-proportion PCA) | Methods; Appendix S1.4 |
| `03-5_pca_check.R` | Compares the PC methods and writes the GCTA covariate files | Appendix S1.4 |

### 04 Primary GWAS

| Script | What it does | Paper |
|---|---|---|
| `04-1_GRM.sh` | Builds the genetic relationship matrix (whole genome and per chromosome for LOCO) | Methods |
| `04-2_pheno.R` | Writes GCTA phenotype files (CSLB rank-based inverse-normal, raw score, simulated phenotypes) | Methods |
| `04-3_Sim_MLMA.sh` | Runs GCTA MLMA on 20 simulated breed-structure phenotypes with 0-30 PCs | Appendix S1.5 |
| `04-4_pc_choice.R` | Summarises lambda of the simulated phenotypes to guide the number of PCs | Appendix S1.5 |
| `04-5_MLMA_LOCO.sh` | Primary GWAS: GCTA MLMA-LOCO of CSLB scores with 0-30 PCs | Fig. 1a, b |
| `04-6_sim_check.R` | Checks the simulated-phenotype experiment | Appendix S1.5 |
| `04-7_model_compare.R` | Compares the 0-30 PC models on the real phenotype (lambda, loci, QQ tail) to choose five PCs | Appendix S1.5 |

### 05 Strata and alternative models

| Script | What it does | Paper |
|---|---|---|
| `05-1_sens_inputs.R` | Defines the stratified and alternative covariate models and writes their input files | Methods; Table S5 |
| `05-2_Sens_MLMA.sh` | Runs the GCTA MLMA-LOCO models for strata and alternative covariates | Fig. 2a; Table S5 |
| `05-3_Sens_REGENIE.sh` | REGENIE models: SNP x age interaction and CSLB >= 50 case-control (all dogs and dogs >= 8 years) | Fig. 2a; Table S6 |
| `05-4_sens_compare.R` | Compares every model with the primary model (lambda, loci, effect sizes) | Table S5; Table S6a |
| `05-5_model_catalog.R` | Catalogue of models (what differs from the primary model) and quality checks | Appendix S1.5 |

### 06 Post-GWAS

| Script | What it does | Paper |
|---|---|---|
| `06-1_threshold_overview.R` | Counts SNPs below each significance threshold and draws Manhattan/QQ plots for every model | supporting check |
| `06-2_grm_het.R` | (A) GRM over-correction diagnostics; (B) heterogeneity of SNP effects between strata | Table S5a |
| `06-3_magma_gene.R` | MAGMA gene-level analysis for each model, permutation null and gene-level heterogeneity | Methods |
| `06-4_gene_sets.R` | Builds the 14 prespecified gene sets and the KEGG/GO sets (before any gene-set result is read) | Table S3 |
| `06-5_magma_geneset.R` | MAGMA competitive gene-set tests (prespecified and genome-wide) | Table S3; Fig. 2a |
| `06-6_h2_rg_power.R` | SNP heritability (GREML), genetic correlations between strata, and power | Fig. 1c; Table S4a |
| `06-7_cslb_measurement.R` | CSLB measurement properties in all DAP surveys: test-retest reliability, dementia-diagnosis AUC | Fig. S2; Table S10 |
| `06-8_age8_h2.R` | Heritability in dogs >= 8 years with mean scores, and genotype-by-age GREML | Table S4 |
| `06-9_network.R` | Network propagation of gene-level signals on a protein interaction network | not reported |
| `06-11_sandor_quant.sh` | Downloads and quantifies a public dog prefrontal-cortex ageing RNA-seq data set (SRA PRJNA815057) | not reported |
| `06-12_sandor_aging.R` | Ageing differential expression from 06-11 and enrichment of CSLB signals | not reported |
| `06-13_gblup_prs.R` | Cross-validated genome-wide polygenic prediction (GBLUP) of CSLB scores | not reported |

### 07-12 Alternative explanations

| Script | What it does | Paper |
|---|---|---|
| `07-1_roh.R` | Runs of homozygosity and F_ROH per dog | Table S9 |
| `07-2_roh_assoc.R` | Tests F_ROH against CSLB in a GREML mixed model | Table S9 |
| `08-1_ancestry.R` | Builds a breed-ancestry similarity matrix from DAP ancestry proportions | Table S9 |
| `08-2_ancestry_test.R` | Variance in CSLB explained by breed ancestry (REML, LRT) | Table S9 |
| `09-1_size_pgs.R` | Body-size polygenic score from published body-size markers (Plassais et al. 2019) | Table S9 |
| `09-2_size_age.R` | Tests body-size polygenic score x age in a GREML mixed model | Table S9 |
| `10-1_chrx_qc.R` | X-chromosome genotype QC and X-specific GRM | Table S9 |
| `10-2_chrx_test.R` | X-chromosome variance test (not testable: almost no heterozygous females after imputation) | Table S9 |
| `10-3_depth_counts.sh` | Aligns DAP raw reads (SRA PRJNA800779) and counts reads in 10 kb bins | Table S9 |
| `10-4_depth_test.R` | mtDNA copy number and Y-chromosome dosage from read depth, and their association with CSLB | Table S9 |
| `11-1_pheno_purify.R` | Flags dogs with a sensory or neurological diagnosis before the CSLB survey | Table S9 |
| `11-2_purified_h2.R` | SNP heritability after excluding the flagged dogs | Table S9 |
| `11-3_premise_check.R` | Checks whether the flagged conditions raise CSLB scores | Appendix S1.12 |
| `12-1_dominance_grm.R` | Additive and dominance GRMs | Table S9 |
| `12-2_dominance_greml.R` | GREML with additive and dominance components | Table S9 |

### 13 Loci and pathways

| Script | What it does | Paper |
|---|---|---|
| `13-1_suggestive_loci.R` | Lead loci (p < 1e-5 and 1e-4) per model, cross-stratum consistency, robustness across covariate models | Table S5 |
| `13-2_regional_ld.R` | Regional plots across models and LD comparison between purebred-like and mixed-like dogs | not reported |
| `13-3_finemap.R` | Conditional analysis (COJO) and fine-mapping (SuSiE) at lead loci | not reported |
| `13-4_breed_ld.R` | LD around lead SNPs within breeds | not reported |
| `13-5_covariate_spec.R` | Explains stratum-specific effects by age, body weight, sex and ancestry interactions | Appendix S1.11 |
| `13-6_int_wt_REGENIE.sh` | Genome-wide SNP x body weight interaction GWAS (REGENIE) | Results (CFA19); Table S7 |
| `13-7_magma.R` | MAGMA gene and gene-set analysis for age-floor, binary and SNP x body weight models | Table S6 |
| `13-8_pathway.R` | Pathway trends across age and body-weight models | not reported |
| `13-9_locus_pathway_map.R` | Maps lead loci to genes and pathways for each model (input to 13-10) | Table S6b |
| `13-10_pathway_test.R` | Tests locus-nominated pathways after removing nominating loci, in complementary strata | Table S6b; Fig. 2 |
| `13-11_pathway_by_model.R` | Tabulates the 13-10 results by model (no new computation) | Table S6b |
| `13-12_calcium_check.R` | Calcium-response set: removal of nominating loci, leave-one-gene-out, cumulative removal, age comparison | Fig. 2b; Fig. S3; Table S8a |
| `13-13_prev_lookup.R` | Looks up the loci of an earlier analysis in the present data | cover letter only |

### 14 Paper tables

| Script | What it does | Paper |
|---|---|---|
| `14-3_power_table.R` | Sample sizes for 80% power under alternative survey designs | Table S10 |
| `14-4_brain_expression.R` | Brain-elevated expression (Human Protein Atlas) of the calcium-response set and of genes near loci | Table S8b |
| `14-5_age_pathway_scan.R` | Genome-wide comparison of gene-set effects between dogs < 8 and >= 8 years | Table S6c |
| `14-6_cfa19_breed_check.R` | CFA19 SNP x body weight locus: breed adjustment, breed exclusion, neighbouring SNPs | Table S7 |
| `14-9_cfa38_breed_check.R` | Same checks for the CFA38 locus | Table S7 |
| `14-11_supp_tables.R` | Assembles Supplementary Tables S1-S10 into a Word file (pandoc) | Tables S1-S10 |
| `14-12_zenodo_data.R` | Builds the Zenodo data folder (summary statistics, gene and gene-set results, dog and variant IDs) | Data availability |

## Paper items and the scripts that produce them

| Paper item | Scripts |
|---|---|
| Fig. 1a, b (Manhattan, QQ) | 04-5 |
| Fig. 1c (SNP heritability) | 06-6, 06-8 |
| Fig. 2a (gene sets × models) | 05-2, 05-3, 06-3, 06-4, 06-5, 13-7 |
| Fig. 2b (calcium-response set after locus removal) | 13-12 |
| Fig. S1 (workflow) | 01-1 to 06-7 (numbers in the diagram) |
| Fig. S2 (CSLB test–retest, AUC) | 06-7 |
| Fig. S3 (leave-one-gene-out, cumulative removal) | 13-12 |
| Table S1 (breed composition) | 00_datacheck, 14-11 |
| Table S2 (genotype QC) | 01-1 to 03-3 |
| Table S3 (prespecified gene sets) | 06-4, 06-5, 13-7 |
| Table S4 (SNP heritability, G × age) | 06-6, 06-8 |
| Table S5 (lead loci, cross-stratum consistency) | 05-4, 06-2, 13-1 |
| Table S6 (older-dog analyses) | 05-4, 13-7, 13-9, 13-10, 14-5 |
| Table S7 (CFA19 and CFA38) | 13-6, 14-6, 14-9 |
| Table S8 (calcium-response genes, brain expression) | 13-12, 14-4 |
| Table S9 (alternative explanations) | 07-1 to 12-2 |
| Table S10 (power and sample size) | 06-7, 14-3 |
| All supplementary tables (Word) | 14-11 |
| Zenodo data folder | 14-12 |

## Licence

[licence]
