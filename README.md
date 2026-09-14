# Age at asthma onset and shared genetic liability with depression and anxiety

Analysis code for a genetic-epidemiology study of the shared genetic architecture and
Mendelian randomization (MR) relationships between **childhood-onset (COA)** and
**adult-onset (AOA)** asthma and **depression (MDD)** and **anxiety (ANX)**, using
European-ancestry GWAS summary statistics, with replication in FinnGen.

> **Summary of findings.** Genome-wide genetic correlations with depression and anxiety are
> substantially stronger for adult-onset than childhood-onset asthma, and this onset difference
> exceeds that for bipolar disorder and schizophrenia. Genetic liability to depression is
> associated with adult-onset asthma (robust to adjustment for BMI and smoking initiation) and,
> in FinnGen, with asthma overall and with childhood asthma, so the association is not confined
> to adult-onset disease. Candidate immune genes (IL33, IL1RL1, IL4R, TSLP, TYK2) showed no
> colocalizing regulatory signal.

This repository contains **code only**. It does not contain GWAS summary statistics or any
individual-level data; all inputs are publicly available from the sources listed below.

## Analysis pipeline

Scripts are numbered in execution order. Each reads harmonized inputs and writes to
sub-directories of the project folder (`harmonized/`, `ldsc/`, `lava/`, `placo/`, `mr/`,
`coloc/`, `Result/`).

| Step | Script | Purpose |
|---|---|---|
| 0 | `R/Step0_harmonize.R` | Harmonize asthma, psychiatric and lung-function GWAS to a common schema |
| 0b | `R/Step0b_harmonize_lung.R` | Harmonize FEV1 / FVC |
| 0c | `R/Step0c_harmonize_immune.R` | Harmonize eosinophil / CRP (GWAS-VCF) |
| 1 | `R/Step1_ldsc.R` | Genome-wide genetic correlation (LDSC via GenomicSEM); liability-scale h²; cross-trait intercepts; subtype contrasts (delta method) |
| 2 | `R/Step2_lava.R` | Local genetic correlation (LAVA) |
| 3 | `R/Step3_placo.R` | Cross-trait pleiotropy screen (PLACO / PLACO+) |
| 4 | `R/Step4_mr.R` | Bidirectional two-sample MR between asthma subtypes and depression / anxiety |
| 5 | `R/Step5_mr_lung.R` | MR involving lung function (FEV1, FVC, FEV1/FVC) |
| 6 | `R/Step6_figures.R` | Exploratory composite figures (superseded by Step 9) |
| 7 | `R/Step7_immune_cis.R` | Whole-blood cis-eQTL MR (eQTLGen) and polygenic eosinophil / CRP MR |
| 8a | `R/Step8a_fetch_tissue_eqtl.R` | Fetch GTEx v8 lung / brain-cortex cis-eQTL instruments |
| 8b | `R/Step8b_tissue_cis_mr.R` | Tissue-specific cis-MR (Wald ratio / IVW) |
| 8c | `R/Step8c_coloc.R` | Colocalization (coloc.abf), full window and ±100 kb |
| 8d | `R/Step8d_coloc_susie.R` | SuSiE-based colocalization allowing multiple causal variants |
| 11 | `R/Step11_mvmr.R` | Multivariable MR: depression + BMI + smoking initiation → asthma subtypes |
| 11b | `R/Step11b_mvmr_qhet.R` | MVMR Q-minimization sensitivity (weak instruments) |
| 12 | `R/Step12_finngen_replication.R` | Replication in FinnGen R13 (Howard 2019 depression exposure); provisional Table S11 |
| 13 | `R/Step13_final_sensitivity.R` | Final sensitivity analyses: MR-PRESSO (seed 20260911, 10,000 simulations), liability-scale Steiger across prevalence grids, final Tables S8 / S11 / S12 |
| 9 | `R/Step9_figures.R` | Manuscript Figures 1–4 |
| 10 | `R/Step10_tables.R` | Manuscript tables (main Tables 1–5, supplementary S1–S12) and assembled workbook |

Run order for a full reproduction: 0 → 0b → 0c → 1 → 2 → 3 → 4 → 5 → 7 → 8a → 8b → 8c → 8d →
11 → 11b → 12 → 13 → 9 → 10. Step 13 must be run after Step 12; running Step 12 alone
overwrites Table S11 with provisional values.

## Reproducing the analysis

1. **Set paths.** Each script defines the project root near the top as a placeholder
   (`PATH/TO/PROJECT`). Replace it with your local project directory. Also set
   `PATH/TO/ld_ref` (1000 Genomes EUR PLINK files), `PATH/TO/eur_w_ld_chr` (HapMap3 LD scores),
   and `PATH/TO/ldref_ascii` (a copy of the EUR PLINK files on a path without spaces or
   non-ASCII characters, which PLINK requires).
2. **Obtain inputs.** Download the GWAS summary statistics listed below into `raw/`.
3. **Run in order** with `source()` from R. Steps 9–10 read the result CSVs from earlier steps
   and regenerate all figures and tables under `Result/`.

### Software

| Component | Version |
|---|---|
| R | 4.5.1 (Steps 0–8), 4.6.1 (Steps 11–13) |
| TwoSampleMR | 0.7.9 |
| MRPRESSO | 1.0 |
| GenomicSEM | 0.0.5 |
| LAVA | 0.1.5 |
| MendelianRandomization | 0.10.0 |
| MVMR | 0.4.8 |
| ieugwasr | 1.1.0.9000 (local clumping) |
| genetics.binaRies | 0.1.2 (PLINK 1.9 binary) |
| data.table | 1.17.8 / 1.18.6.1 |
| ggplot2 | 4.0.0 / 4.0.3 |

Other packages: coloc, susieR, Rsamtools, patchwork, openxlsx, psych. Exact versions for each
step are recorded in the `_stepN_sessionInfo_*.txt` files written by the scripts.

> `Step3_placo.R` requires the `PLACO` function source (`PLACO_v0.2.0.R`) obtained from the
> PLACO authors (Ray & Chatterjee, *PLoS Genet* 2020); it is not redistributed here.

## Data availability (input GWAS)

| Trait | Source |
|---|---|
| Childhood-/adult-onset asthma | Ferreira et al., *Am J Hum Genet* 2019 (PMID 30929738) |
| Depression | PGC MDD (Adams et al., *Cell* 2025; EUR, no UKB, no 23andMe) |
| Anxiety | PGC-ANX (Strom et al., *Nat Genet* 2026; without Utah cohort) |
| Bipolar disorder | PGC3 (Mullins et al., *Nat Genet* 2021; no UKB) |
| Schizophrenia | PGC3 (Trubetskoy et al., *Nature* 2022; EUR) |
| FEV1 / FVC / FEV1:FVC | SpiroMeta (Shrine et al., *Nat Genet* 2019) |
| Eosinophil count | Astle et al., *Cell* 2016 |
| C-reactive protein | Ligthart et al., *Am J Hum Genet* 2018 |
| Body mass index | GIANT (Locke et al., *Nature* 2015; no UKB) |
| Smoking initiation | GSCAN (Liu et al., *Nat Genet* 2019; without UKB and 23andMe) |
| Replication: asthma, childhood asthma | FinnGen R13 (J10_ASTHMA_EXMORE, ASTHMA_CHILD_EXMORE) |
| Replication: depression | Howard et al., *Nat Neurosci* 2019 (PGC + UKB, no 23andMe) |
| Tissue cis-eQTL | GTEx v8 / eQTL Catalogue; whole blood: eQTLGen |

## Citation

If you use this code, please cite: **[manuscript citation — to be added on acceptance]**.

## License

Code released under the MIT License (see `LICENSE`). Third-party GWAS and eQTL datasets remain
under their original licenses.

## Contact

Hye Jin Lee, MD, PhD — Division of Allergy and Pulmonology, Department of Pediatrics,
Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea.
[e-mail — to be added]
