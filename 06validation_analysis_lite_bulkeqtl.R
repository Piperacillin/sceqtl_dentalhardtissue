# 06validation_analysis_bulk.R
# Phase 4 专用：Bulk eQTL 验证分析脚本
# 功能：适配 eQTLGen/GTEx 数据格式 -> 运行 IVW/Wald -> 格式化输出
head(qs::qread('./exposure/eqtlgen/eqtlgen_r2_0.2_validation.qs'))
head(qs::qread('./exposure/GTEx/GTEx_r2_0.2_validation.qs'))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(qs)
  library(data.table)
})

# ================= 接收环境变量 =================
EXPOSURE_DATA <- Sys.getenv("EXPOSURE_DATA") # e.g., eqtlgen
OUTCOME       <- Sys.getenv("OUTCOME")       # e.g., K11_ABRASION
START         <- as.numeric(Sys.getenv("START"))
END           <- as.numeric(Sys.getenv("END"))
INPUT_QS      <- Sys.getenv("INPUT_QS")      

# 路径定义
BASE_DIR <- paste0(EXPOSURE_DATA, "_", OUTCOME)
RES_DIR  <- file.path(BASE_DIR, "results")
if(!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)

# ================= 1. 数据读取与预处理 =================
tryCatch({
  # 读取 Outcome
  out_file <- paste0("outcome/", OUTCOME, ".qs")
  out <- qread(out_file)
  
  # 读取 Exposure (Script 12 生成的 validation.qs)
  exp_full <- qread(INPUT_QS)
  
  # --- [Bulk 适配] 标准化列名 ---
  # Bulk 数据通常没有 'exposure' 列作为唯一ID，通常是 'gene.exposure' 或 'symbol.exposure'
  # 我们需要构造一个标准的 'exposure' 列用于分块
  
  if(!"exposure" %in% names(exp_full)) {
    if("phenotype.exposure" %in% names(exp_full)) {
      exp_full$exposure <- exp_full$phenotype.exposure
    } else if("gene.exposure" %in% names(exp_full)) {
      exp_full$exposure <- exp_full$gene.exposure
    } else {
      stop("Cannot determine Gene ID column in Bulk data")
    }
  }
  
  # --- [Bulk 适配] 补全缺失的 Metadata ---
  # 确保后续贴标签时不会报错
  if(!"gene_symbol" %in% names(exp_full)) exp_full$gene_symbol <- exp_full$gene.exposure
  if(!"cell_type" %in% names(exp_full))   exp_full$cell_type <- paste0("Bulk_", EXPOSURE_DATA)
  if(!"act_time" %in% names(exp_full))    exp_full$act_time <- "Static"
  
  # 切片逻辑
  all_genes <- unique(exp_full$exposure)
  
  if (START > length(all_genes)) quit(save="no")
  real_end <- min(END, length(all_genes))
  
  target_genes <- all_genes[START:real_end]
  exp_sub <- subset(exp_full, exposure %in% target_genes)
  
  # 提取元数据用于最后合并
  meta_info <- exp_sub %>%
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
# 这里的逻辑与 Phase 3 完全一致
res_final <- data.frame()
unique_exps <- unique(dat$exposure)

for(e in unique_exps) {
  dat_single <- subset(dat, exposure == e)
  try({
    mr_res <- mr(dat_single, method_list = c("mr_ivw", "mr_wald_ratio"))
    if(nrow(mr_res) > 0) {
      clean_res <- mr_res %>%
        select(exposure, outcome, method, nsnp, b, se, pval) %>%
        rename(beta = b, p = pval)
      res_final <- rbind(res_final, clean_res)
    }
  }, silent = TRUE)
}

# ================= 4. 保存结果 =================
if(nrow(res_final) > 0) {
  
  # 贴回元数据
  res_annotated <- res_final %>%
    left_join(meta_info, by = "exposure") %>%
    select(gene_symbol, cell_type, act_time, everything())
  
  out_name <- paste0("results_val_chunk_", START, ".txt")
  write.table(res_annotated, file.path(RES_DIR, out_name), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste("Success:", nrow(res_annotated), "associations saved."))
  
} else {
  message("No results generated.")
}
