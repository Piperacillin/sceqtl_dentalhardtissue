# 04validation_analysis_lite.R
# 专为验证阶段设计的轻量级 MR 分析脚本
# 功能：读取子集 QS -> 协调数据 -> 仅运行 IVW 和 Wald Ratio -> 贴回基因/细胞信息

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

# ================= 接收环境变量 =================
EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA") # 例如: CD4_Dynamic
OUTCOME       <- Sys.getenv("OUTCOME")       # 例如: K11_ABRASION
START         <- as.numeric(Sys.getenv("START"))
END           <- as.numeric(Sys.getenv("END"))
R2_THRESH     <- Sys.getenv("R2_THRESH")
INPUT_QS      <- Sys.getenv("INPUT_QS")      # 直接传入文件路径

# 路径定义
BASE_DIR <- paste0(EXPOSURE_DATA, "_", OUTCOME)
QC_DIR   <- file.path(BASE_DIR, "eggercochranesq") 
RES_DIR  <- file.path(BASE_DIR, "results")
if(!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
if(!dir.exists(QC_DIR)) dir.create(QC_DIR, recursive = TRUE)
# ================= 1. 数据读取与切片 =================
tryCatch({
  # 读取 Outcome
  out_file <- paste0("outcome/", OUTCOME, ".qs")
  out <- qread(out_file)
  
  # 读取 Exposure (用户预处理好的子集 QS)
  exp_full <- qread(INPUT_QS)
  
  # 获取当前 Chunk 的基因列表
  all_genes <- unique(exp_full$exposure)
  
  # 防止索引越界
  if (START > length(all_genes)) quit(save="no")
  real_end <- min(END, length(all_genes))
  
  target_genes <- all_genes[START:real_end]
  exp_sub <- subset(exp_full, exposure %in% target_genes)
  
  # --- [关键步骤] 提取元数据 (Metadata) ---
  # 使用大括号 {} 在管道中进行条件判断：
  # 如果 act_time 列不存在，就 mutate 创建为 "none"，否则保持原样
  meta_info <- exp_sub %>%
    { if (!"act_time" %in% names(.)) mutate(., act_time = "none") else . } %>%
    select(exposure, gene_symbol, cell_type, act_time) %>%
    distinct()
  
}, error = function(e) {
  message("[Error] Data loading failed: ", e$message)
  quit(save="no")
})

# ================= 2. 协调数据 (Harmonise) =================
dat <- harmonise_data(exposure_dat = exp_sub, outcome_dat = out)
dat <- subset(dat, mr_keep == TRUE)

if(nrow(dat) == 0) {
  message("No valid SNPs after harmonisation.")
  quit(save="no")
}

# ================= 3. 运行 MR (仅 IVW & Wald) =================
# Steiger
if(!"r.exposure" %in% names(dat)) dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
if(!"r.outcome" %in% names(dat))  dat$r.outcome  <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
steiger_pass <- subset(directionality_test(dat), correct_causal_direction == TRUE)
dat_steiger <- subset(dat, dat$exposure %in% steiger_pass$exposure)

if(nrow(dat_steiger) == 0) quit(save="no")

# Primary MR Analysis (Generating the Nominal P-values)

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


# ================= 4. 贴回拆分信息并保存 =================
if(nrow(all_res) > 0) {
  
  # 使用 left_join 将 gene_symbol, cell_type 等信息贴回来
  # 这样就不需要再用正则去猜了，直接使用输入文件里已经做好的拆分
  res_annotated <- all_res %>%
    left_join(meta_info, by = "exposure") %>%
    # 调整列顺序，把重要的放前面
    select(gene_symbol, cell_type, any_of("act_time"), everything())
  
  # 保存结果
  out_name <- paste0("results_val_chunk_", START, ".txt")
  
  write.table(res_annotated, file.path(RES_DIR, out_name), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message(paste("Success:", nrow(res_annotated), "associations saved."))
  
} else {
  message("No results generated.")
}