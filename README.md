
***

# Single-Cell Gene Expression and Dental Hard Tissue Disease: A Mendelian Randomization Study

![Status: Under Review](https://img.shields.io/badge/Status-Under_Review_at_JDR-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green)
![R >= 4.3.3](https://img.shields.io/badge/R-%3E%3D%204.3.3-orange)

> **📝 Manuscript Status:**  
> The manuscript detailing this single-cell MR analytical framework and its extensive findings is currently **under review** at the *Journal of Dental Research (JDR)*. The full datasets and detailed extended results described in the Appendix are generated using the following pipeline.

## 👥 Authors and Affiliations
**Tianyi Wang<sup>1,2†</sup>, Qi Xing<sup>1,3†</sup>, Liu Liu<sup>1,3,#</sup>, Yongwen Guo<sup>1,2</sup>, Yang Song<sup>1</sup>, Dingming Huang<sup>1,3</sup>, Lan Zhang<sup>1,3,#</sup>**

<sup>†</sup> These authors contributed equally to this article.  
<sup>#</sup> **Correspondence:** Prof. Lan Zhang (zlnancy914@sina.com) and Dr. Liu Liu (liul@scu.edu.cn, Twitter: @LiuLiu_DDS).  

<sup>1</sup> *State Key Laboratory of Oral Diseases & National Center for Stomatology & National Clinical Research Center for Oral Diseases, West China Hospital of Stomatology, Sichuan University, Chengdu 610041, China*
<sup>2</sup> *Department of Orthodontics, West China Hospital of Stomatology, Sichuan University.*  
<sup>3</sup> *Department of Conservative Dentistry and Endodontics, West China Hospital of Stomatology, Sichuan University.*

---

## 📖 Introduction
Welcome to the official repository for our large-scale integration of single-cell expression quantitative trait loci (sc-eQTLs) and GWAS data to elucidate the genetic architecture of dental hard tissue traits. 

To overcome the analytical challenges of multi-cohort single-cell dimensions, we structured our reproducible MR pipeline into **Four Core Phases**, progressing from initial discovery to rigorous validation across 4 major sc-eQTL databases (1MscBloodNL, OneK1K, DICE, CD4dynamic), culminating in a meticulously designed **Scenario-Adaptive Scoring System** and bulk-tissue benchmarking.

## ⚙️ Analytical Pipeline & Scripts Overview

Our pipeline is controlled by master bash scripts that sequence the R-based data processing steps.

* **Preparation:** `installpackages.r` sets up the environment. `00datacleaning_outcome.R` and `01datacleaning_exposure.R` format raw GWAS and eQTL datasets.
* **Phase 1 (LD Clumping & Exploratory MR):** Scripts like `01.1dynamicTclump0.2.r` to `01.6onek1kclump0.01.r` perform rigorous LD pruning, utilizing adaptive R² thresholds (e.g., 0.2 for phase 1, 0.01 for phase 2) to maximize instrumental variant capture without compromising independence. The Discovery pipeline, utilizing R² = 0.2 clumping data and primary FDR corrections, Executed via `01_run_phase1_discovery.sh`, utilizing `02analysispipline.R`. 
* **Phase 2 (Replication Prep & Confirmatory MR):** Controlled by `02_process_and_prep_phase2.sh` （utilizing `02.1merge_chunks_calculate_fdr.R`） and `03_run_phase2_replication.sh`（utilizing `02analysispipline.R` to run TwoSampleMR algorithms under R² = 0.01 clumping data and Bonferonni corrections）  .
* **Phase 3 (Multi-Cohort Harmonization):** Executed via `00.1_master_validation.sh` and `03-05 R and bash scripts`. Aligns and heavily cross-validates genetic instruments across the 4 core databases (1MscBloodNL, OneK1K, Validation using DICE & CD4dynamic).
* **Phase 4 (Bulk Assessment & Scoring System):** Executed via `00.2_master_phase4_bulk.sh` and `06_run_phase4_bulk.sh`. Leverages eQTLGen data (`06validation_analysis_lite_bulkeqtl.R`) for systemic validation. Finally, `07postphase4cleaning_updated.R` strictly calculates the quantitative **Scenario-Adaptive Scores**, tiering targets into robust Confidence classifications.
* **Downstream Analysis:** `08colocalization.R` computes Bayesian posterior probabilities to confirm shared causal variants.
* Environmental parameters are stored in `aa_discovery_config.sh`
* Utility script `Rterminator.sh` is provided for dynamic memory and parallel process management.

## 📂 Repository Structure
```text
scEqtl-DentalHT-MR/
├── aa_discovery_config.sh              # Global configuration variables
├── 00_master_runner.sh                 # Master execution shell
├── 00.1_master_validation.sh           # Phase 3 Multi-database harmonizer
├── 00.2_master_phase4_bulk.sh          # Phase 4 Bulk & Scoring runner
├── *.sh                                # Sub-step runners (01_run_**, 02_process_**)
├── Rterminator.sh                      # Memory management utility
│
├── installpackages.r                   # Dependencies installation
├── 00datacleaning_outcome.R
├── 01datacleaning_exposure.R
│
├── 01.*clump*.r                        # Phase 1: LD Clumping scripts (0.01 & 0.2 thresholds)
├── 02analysispipline.R                 # Core TwoSampleMR pipeline
├── 02.1merge_chunks_calculate_fdr.R
├── 03filter_and_mapping_*.R            # Phase 3: DICE/CD4dynamic mappings
├── 04validation_analysis_lite.R        
├── 05postphase3cleaning.R              
├── 06validation_analysis_*.R           # Phase 4: eQTLGen Bulk referencing
├── 07postphase4cleaning_updated.R      # Phase 4: SCENARIO-ADAPTIVE SCORING SYSTEM
└── 08colocalization.R                  # Bayesian Colocalization
│
├── exposure/                           # Processed eGenes instruments
├── outcome/                            # Formatted Dental Traits GWAS
└── rawinstruments/                     # Pre-clumped raw variants
```

## 📊 Data Availability
* **Dental Phenotypes:** [FinnGen Consortium Data Freeze 12](https://www.finngen.fi/en)
* **Single-Cell eQTLs:** OneK1K, 1MscBloodNL, DICE, and CD4dynamic.
* **Bulk eQTLs:** eQTLGen Consortium.

*(Please refer to `Appendix Tables 1 & 2` in our manuscript for full data accession protocols).*

## 📚 Please Cite
If you utilize our pipeline architecture, scoring system, codes, or conceptually build upon this workflow, please cite our pending manuscript alongside the foundational frameworks:

**Primary Reference:**
1. [Tianyi Wang, Qi Xing], Liu Liu, *et al.* (2024). Single-Cell Gene Expression and Dental Hard Tissue Disease: A Mendelian Randomization Study. *Journal of Dental Research* [Under Review].
2. Liu, L., *et al.* (2024). Genetically Supported Drug Targets and Dental Traits: A Mendelian Randomization Study. *Journal of Dental Research*. [https://doi.org/10.1177/00220345241272045](https://doi.org/10.1177/00220345241272045)

**sc-eQTL Data Sources & Methodological Foundations:**
* **1mscblood:** Oelen, R. K., *et al.* (2022). *Nature Communications*, 13(1):3267.
* **OneK1K:** Yazar, S., *et al.* (2022). *Science*, 376(6589), eabf3041.
* **CD4dynamic:** Soskic, B., *et al.* (2022). *Nature Genetics*, 54(6), 817-826.
* **DICE:** Schmiedel, B. J., *et al.* (2018). *Cell*, 175(6), 1701-1715.

## 🤝 Acknowledgments
We profoundly acknowledge the original contributors of the QTL working pipeline readapted in our study. Kindly reference the pivotal methodological articles mentioned above. Furthermore, we extend our greatest gratitude to the thousands of participants and the dedicated researchers of the participating consortia for making their invaluable summary statistics publicly available.

## 📜 License
This methodology and code pipeline is provided under the [MIT License](LICENSE).
