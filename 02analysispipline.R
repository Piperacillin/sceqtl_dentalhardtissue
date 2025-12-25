####敏感性分析版本 ####
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(stringr)
  library(qs)
  library(data.table)
  library(MendelianRandomization) # For cML, Robust, dIVW
  library(mr.raps)                # For MR-RAPS
  library(MRPRESSO)               # For MR-PRESSO
  library(plyr)
  library(plinkbinr)
  library(ieugwasr)
})

# ==============================================================================
# 1. Setup & Data Loading (Same as before)
# ==============================================================================

EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA")
OUTCOME       <- Sys.getenv("OUTCOME")
START         <- as.numeric(Sys.getenv("START"))
END           <- as.numeric(Sys.getenv("END"))
R2_THRESH     <- Sys.getenv("R2_THRESH")
# 从环境变量读取是否只跑主分析 (Phase 2用)
ONLY_PRIMARY <- Sys.getenv("ONLY_PRIMARY")

# Paths
EXP_FILE <- paste0("exposure/", EXPOSURE_DATA,"/",EXPOSURE_DATA, "_r2_", R2_THRESH, ".qs")
OUT_FILE <- paste0("outcome/", OUTCOME, ".qs")
BASE_DIR <- paste0(EXPOSURE_DATA, "_", OUTCOME)
RES_DIR  <- file.path(BASE_DIR, "results")
SENS_DIR <- file.path(BASE_DIR, "sensitivity") # New separate folder for extra analysis
DATA_DIR <- file.path(BASE_DIR, "data")
# ---  QC 结果文件夹 ---
QC_DIR   <- file.path(BASE_DIR, "eggercochranesq") 
# Modify directory creation loop to include QC_DIR
for(d in c(RES_DIR, SENS_DIR, DATA_DIR, QC_DIR)) {
  if(!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Load Data & Slice
out <- qs::qread(OUT_FILE)
exp0 <- qs::qread(EXP_FILE)
all_exposures <- unique(exp0$exposure)
real_end <- min(END, length(all_exposures))

if (START > length(all_exposures)) {
  message("idx out of range")
  quit(save="no")
}

exp_to_keep <- all_exposures[START:real_end]
exp <- subset(exp0, exp0$exposure %in% exp_to_keep)

# Optional remove list logic here (omitted for brevity, same as previous)

# Harmonise
dat0 <- harmonise_data(exposure_dat = exp, outcome_dat = out)
dat <- subset(dat0, mr_keep == TRUE)
if(nrow(dat) == 0) quit(save="no")

# Steiger
if(!"r.exposure" %in% names(dat)) dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
if(!"r.outcome" %in% names(dat))  dat$r.outcome  <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
steiger_pass <- subset(directionality_test(dat), correct_causal_direction == TRUE)
dat_steiger <- subset(dat, dat$exposure %in% steiger_pass$exposure)

if(nrow(dat_steiger) == 0) quit(save="no")

# ==============================================================================
# 2. Primary MR Analysis (Generating the Nominal P-values)
# ==============================================================================

snp_counts <- as.data.frame(table(dat_steiger$exposure))
exposures_le2 <- as.character(snp_counts$Var1[snp_counts$Freq <= 2])
exposures_gt2 <- as.character(snp_counts$Var1[snp_counts$Freq > 2])

dat_le2 <- subset(dat_steiger, exposure %in% exposures_le2)
dat_gt2 <- subset(dat_steiger, exposure %in% exposures_gt2)

# --- A. Uncorrelated / Small SNP Count (<=2) ---
res_le2 <- data.frame()
if(nrow(dat_le2) > 0) {
  try({
    temp <- mr(dat_le2, method_list = c("mr_ivw", "mr_wald_ratio"))
    if(!plyr::empty(temp)){
      res_le2 <- temp[, c("exposure", "outcome", "nsnp", "method", "b", "se", "pval")]
      names(res_le2) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")
      res_le2$clump_thresh <- R2_THRESH
    }
  })
}

# --- B. Correlated / Large SNP Count (>2) ---
res_gt2 <- data.frame()
qc_res_container <- data.frame()
if(nrow(dat_gt2) > 0){
  unique_gt2 <- unique(dat_gt2$exposure)
  source('./rawinstruments/dat_to_MRInput_local_mod.R') # Ensure this path is correct
  
  for(exp_id in unique_gt2){
    dat_sub <- dplyr::filter(dat_gt2, exposure == exp_id)
    
    # Error handling for individual primary analysis
    tryCatch({
      dat_mrinput <- dat_to_MRInput_local(dat_sub, get_correlation=TRUE)
      
      # Try IVW first (Primary Method)
      ivw_mod <- MendelianRandomization::mr_ivw(dat_mrinput[[1]], correl=TRUE)
      
      # 计算 Egger Intercept & Heterogeneity
      egger_mod <- MendelianRandomization::mr_egger(dat_mrinput[[1]], correl=TRUE)
      
      # 提取指标
      temp_qc <- data.frame(
        exposure = egger_mod@Exposure,
        outcome  = egger_mod@Outcome,
        nsnp     = egger_mod@SNPs,
        # 多效性相关
        egger_intercept        = egger_mod@Intercept,
        egger_intercept_95_ci  = paste(round(egger_mod@CILower.Int, 4), round(egger_mod@CIUpper.Int, 4), sep=", "),
        egger_intercept_pvalue = egger_mod@Pvalue.Int,
        # 异质性相关 (MendelianRandomization包中 Heter.Stat[1]是Q, [2]是P)
        cochrans_q      = egger_mod@Heter.Stat[1],
        cochrans_q_pval = egger_mod@Heter.Stat[2],
        stringsAsFactors = FALSE
      )
      
      # 计算 I2 (异质性程度)
      # 逻辑: I2 = (Q - df) / Q; df = nsnp - 2 (对于Egger) 或 nsnp - 1 (对于IVW)
      # 为了保守起见，通常用 IVW 的 Q 或 Egger 的 Q 均可，这里沿用你之前提供的逻辑计算
      q_val <- egger_mod@Heter.Stat[1]
      df    <- egger_mod@SNPs - 1 # 近似自由度
      
      if (q_val >= df) {
        temp_qc$i2 <- (q_val - df) / q_val
      } else {
        temp_qc$i2 <- 0
      }
      
      # 存入容器
      qc_res_container <- rbind(qc_res_container, temp_qc)

      temp_row <- data.frame(
        exposure = ivw_mod@Exposure, outcome = ivw_mod@Outcome, nsnp = ivw_mod@SNPs,
        method = ivw_mod@class[1], beta = ivw_mod@Estimate, se = ivw_mod@StdError,
        p = ivw_mod@Pvalue, clump_thresh = R2_THRESH, stringsAsFactors = F
      )
      res_gt2 <- rbind(res_gt2, temp_row)
      
      # Capture Egger Intercept here if needed (omitted for brevity to focus on sensitivity logic)
      
    }, error = function(e) {
      message(paste("[Error] Primary MR failed for:", exp_id, "-", e$message))
    })
  }
}

# --- [新增] 保存 QC (Egger Intercept & Cochran's Q) 结果 ---
if(nrow(qc_res_container) > 0) {
  qc_file_name <- paste0("results_egger_cochransq_", EXPOSURE_DATA, "_", OUTCOME, "_",R2_THRESH , "_chunk_", START, ".txt")
  write.table(qc_res_container, file.path(QC_DIR, qc_file_name), sep = "\t", row.names = F, quote = F)
}

# Combine Primary Results
all_res <- rbind(res_le2, res_gt2)

# Save Primary Results
write.table(all_res, file.path(RES_DIR, paste0("results_primary_r2_", R2_THRESH, "_chunk_", START, ".txt")), sep = "\t", row.names = F, quote = F)

# ==============================================================================
# 3. Conditional Sensitivity Analysis (Only if P < 0.05)
# ==============================================================================
if(ONLY_PRIMARY == "TRUE") {
  message("ONLY_PRIMARY flag is set. Skipping sensitivity analysis.")
  message("Chunk Processing Complete.")
  quit(save = "no")
}

# Identify significant exposures
sig_exposures <- unique(all_res$exposure[all_res$p < 0.05])

if(length(sig_exposures) > 0) {
  
  message(paste("Running sensitivity analysis for", length(sig_exposures), "significant exposures..."))
  
  # Containers for sensitivity results
  sens_mr_res    <- data.frame() # Weighted Median/Mode
  sens_raps_res  <- data.frame() # RAPS
  sens_presso_res<- data.frame() # PRESSO
  sens_mrb_res   <- data.frame() # cML, Robust, dIVW
  
  # Loop through each significant exposure
  for(exp_id in sig_exposures) {
    print(paste(exp_id, "in", length(sig_exposures)))
    # Subset data once for this exposure
    dat_sub <- subset(dat_steiger, exposure == exp_id)
    nsnp <- nrow(dat_sub)
    if(nsnp == 0) next
    
    message(paste("  > Processing:", exp_id, "(nSNP =", nsnp, ")"))
    
    # 定义一个通用的 sample size (用于 cML)
    # 如果有多个SNP，取平均或中位数样本量
    n_sample <- median(dat_sub$samplesize.exposure, na.rm=TRUE)
    if(is.na(n_sample)) n_sample <- 1000 # Fallback if missing
    
    # ------------------------------------------------------
    # Method 1: Weighted Median & Weighted Mode (via TwoSampleMR)
    # ------------------------------------------------------
    message(paste("Weighted Median & Weighted Mode."))
    tryCatch({
      if(nsnp >= 3) {
        wm_res <- mr(dat_sub, method_list = c("mr_weighted_median", "mr_weighted_mode"))
        if(nrow(wm_res) > 0) {
          sens_mr_res <- rbind(sens_mr_res, wm_res)
        }
      }
    }, error = function(e) message(paste("    [Fail] Weighted methods:", exp_id)))
    
    # ------------------------------------------------------
    # Method 2: MR-RAPS
    # ------------------------------------------------------
    message(paste("MR-RAPS."))
    tryCatch({
      if(nsnp >= 2) {
      raps_fit <- mr.raps.all(b_exp = dat_sub$beta.exposure, b_out = dat_sub$beta.outcome,
                              se_exp = dat_sub$se.exposure, se_out = dat_sub$se.outcome)
      diag <- mr.raps(b_exp = dat_sub$beta.exposure, b_out = dat_sub$beta.outcome,
                      se_exp = dat_sub$se.exposure, se_out = dat_sub$se.outcome, 
                      diagnosis = FALSE)
      # Extracting RAPS results (usually index 1 is robust)
      raps_row <- data.frame(
        exposure = exp_id, outcome = OUTCOME, nsnp = nsnp, method = "MR-RAPS",
        beta = raps_fit[6,3], se = raps_fit[6,4], p = diag$beta.p.value
      )
      sens_raps_res <- rbind(sens_raps_res, raps_row)
      }
    }, error = function(e) message(paste("    [Fail] RAPS:", exp_id)))
    
    # ------------------------------------------------------
    # Method 3: MR-PRESSO (Require >= 4 SNPs)
    # ------------------------------------------------------
    message(paste("MR-PRESSO."))
    tryCatch({
      if(nsnp >= 4) {
        # MR-PRESSO is strict about data frames, ensuring independent SNPs usually
        # Note: PRESSO expects global variables or direct columns. 
        # We map 'dat_sub' columns to function args.
        presso_out <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                                SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                data = dat_sub, NbDistribution = 1000, SignifThreshold = 0.05)
        
        # Extract Main Result (Raw) and Outlier Corrected Result
        main_res <- presso_out$`Main MR results`
        
        # Raw
        if(nrow(main_res) >= 1){
          sens_presso_res <- rbind(sens_presso_res, data.frame(
            exposure = exp_id, outcome = OUTCOME, nsnp = nsnp, method = "MR-PRESSO (Raw)",
            beta = main_res[1, "Causal Estimate"], se = main_res[1, "Sd"], 
            p = main_res[1, "P-value"], global_pval = presso_out$`MR-PRESSO results`$`Global Test`$Pvalue
          ))
        }
        # Outlier Corrected (if exists)
        if(nrow(main_res) >= 2){
          sens_presso_res <- rbind(sens_presso_res, data.frame(
            exposure = exp_id, outcome = OUTCOME, nsnp = nsnp, method = "MR-PRESSO (Outlier-corrected)",
            beta = main_res[2, "Causal Estimate"], se = main_res[2, "Sd"], 
            p = main_res[2, "P-value"], global_pval = NA
          ))
        }
      }
    }, error = function(e) message(paste("    [Fail] PRESSO:", exp_id, e$message)))
    
    # ------------------------------------------------------
    # Method 4: MendelianRandomization Pkg (cML, Robust IVW, dIVW)
    # ------------------------------------------------------
    # Need to convert to MRInput object first
    message(paste("cML, Robust IVW, dIVW."))
    tryCatch({
      mr_input_obj <- dat_to_MRInput(dat_sub) # Use standard unc-correlated input for generalized sensitivity
      # Note: If you specifically want Correlated cML, use dat_to_MRInput_local logic, 
      # but most sensitivity checks (like cML-DP) handle validity internally.
      
      # 4a. cML
      try({
        cml_fit <- mr_cML(mr_input_obj[[1]], n = n_sample) # Need sample size!
        sens_mrb_res <- rbind(sens_mrb_res, data.frame(
          exposure = exp_id, outcome = OUTCOME, nsnp = cml_fit@SNPs, method = "cML",
          beta = cml_fit@Estimate, se = cml_fit@StdError, p = cml_fit@Pvalue
        ))
      }, silent = TRUE)
      
      # 4b. Robust IVW
      try({
        rob_fit <- mr_ivw(mr_input_obj[[1]], robust = TRUE)
        sens_mrb_res <- rbind(sens_mrb_res, data.frame(
          exposure = exp_id, outcome = OUTCOME, nsnp = rob_fit@SNPs, method = "Robust IVW",
          beta = rob_fit@Estimate, se = rob_fit@StdError, p = rob_fit@Pvalue
        ))
      }, silent = TRUE)
      
      # 4c. dIVW
      try({
        # mr_divw usually implies debiased IVW
        divw_fit <- mr_divw(mr_input_obj[[1]])
        sens_mrb_res <- rbind(sens_mrb_res, data.frame(
          exposure = exp_id, outcome = OUTCOME, nsnp = divw_fit@SNPs, method = "dIVW",
          beta = divw_fit@Estimate, se = divw_fit@StdError, p = divw_fit@Pvalue
        ))
      }, silent = TRUE)
      
    }, error = function(e) message(paste("    [Fail] MR Pkg Methods:", exp_id)))
    
  } # End Loop
  
  # ==========================================================================
  # 4. Save Sensitivity Results Files
  # ==========================================================================
  
  suffix <- paste0("_r2_", R2_THRESH, "_chunk_", START, ".txt")
  
  # ==========================================
  # 41. 提取元数据 (ID信息)
  # ==========================================
  # 默认设为 NA，以防 sens_mr_res 为空
  fill_id_exp <- NA
  fill_id_out <- NA
  # 如果 sens_mr_res 有数据，从中提取 ID
  if (exists("sens_mr_res") && nrow(sens_mr_res) > 0) {
    fill_id_exp <- unique(sens_mr_res$id.exposure)[1]
    fill_id_out <- unique(sens_mr_res$id.outcome)[1]
  }
  # ==========================================
  # 42. 定义标准化的列名列表
  # ==========================================
  # 我们希望最后保留这些列
  target_cols <- c("id.exposure", "id.outcome", "outcome", "exposure", 
                   "method", "nsnp", "b", "se", "pval", "global_pval")
  # ==========================================
  # 43. 分别处理每个数据框
  # ==========================================
  # --- 处理 sens_mr_res ---
  # 特点: 已经是 b, pval; 有 id; 缺 global_pval
  df1 <- if (exists("sens_mr_res") && nrow(sens_mr_res) > 0) {
    sens_mr_res %>%
      mutate(global_pval = NA) # 补上缺少的列
  } else { NULL }
  # --- 处理 sens_raps_res ---
  # 特点: beta -> b; p -> pval; 缺 ids; 缺 global_pval
  df2 <- if (exists("sens_raps_res") && nrow(sens_raps_res) > 0) {
    sens_raps_res %>%
      dplyr::rename(b = beta, pval = p) %>%
      dplyr::mutate(id.exposure = fill_id_exp,
             id.outcome = fill_id_out,
             global_pval = NA)
  } else { NULL }
  # --- 处理 sens_presso_res ---
  # 特点: beta -> b; p -> pval; 缺 ids; 有 global_pval
  df3 <- if (exists("sens_presso_res") && nrow(sens_presso_res) > 0) {
    sens_presso_res %>%
      dplyr::rename(b = beta, pval = p) %>%
      dplyr::mutate(id.exposure = fill_id_exp,
             id.outcome = fill_id_out) # global_pval 已经有了，不用补
  } else { NULL }
  # --- 处理 sens_mrb_res ---
  # 特点: beta -> b; p -> pval; 缺 ids; 缺 global_pval
  df4 <- if (exists("sens_mrb_res") && nrow(sens_mrb_res) > 0) {
    sens_mrb_res %>%
      dplyr::rename(b = beta, pval = p) %>%
      dplyr::mutate(id.exposure = fill_id_exp,
             id.outcome = fill_id_out,
             global_pval = NA)
  } else { NULL }
  # ==========================================
  # 44. 合并并重排
  # ==========================================
  final_sensitivity_res <- bind_rows(df1, df2, df3, df4) %>%
    # 确保只保留需要的列，并按指定顺序排列
    dplyr::select(any_of(target_cols)) 
  
  if(nrow(final_sensitivity_res) > 0)    write.table(final_sensitivity_res, file.path(SENS_DIR, paste0("sensitivity_res", suffix)), sep="\t", row.names=F, quote=F)

  
} else {
  message("No significant results (p<0.05) found in this chunk. Skipping sensitivity analysis.")
}

message("Chunk Processing Complete.")


