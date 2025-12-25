#!/usr/bin/env Rscript

# Rscript 01.10eqtlgenclump0.01.r 2>&1 | tee 01.10eqtlgenclump0.01.r.log

# ==============================================================================
# 配置区 
# ==============================================================================
WORK_DIR       <- "/home/publicdata/32445978_bookliu11/eqtlgen"
INPUT_FILE     <- "eqtlgen_cis_1000kb.qs"
OUTPUT_FILE    <- "eqtlgen_r2_0.01.qs"
#CHAIN_PATH     <- "/home/piperacillin/publicdata/rawinstruments/hg38ToHg19.over.chain"
BIM_FILE       <- "/opt/biosoft/rawinstruments/EUR.bim"
PLINK_BFILE    <- "/opt/biosoft/rawinstruments/EUR"
N_CORES        <- 30
CLUMPR <- 0.01

# !!! 【重要】请在这里填入你刚才用 `which plink` 查到的绝对路径 !!!
# 如果你的服务器上没有 plink，需要先去官网下载 plink 1.9 并解压
#plinkbinr::get_plink_exe()
PLINK_BIN      <- "/usr/local/lib/R/library/plinkbinr/bin/plink_Linux" 
# 或者可能是: "/home/piperacillin/bin/plink" 

# ==============================================================================
# 初始化与日志函数
# ==============================================================================
log_msg <- function(msg) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"))
}

log_msg("正在检查环境...")

# 0. 检查 PLINK 是否存在 (关键步骤)
if (!file.exists(PLINK_BIN)) {
  stop(paste("严重错误：PLINK 可执行文件未找到！路径:", PLINK_BIN, "\n请修改脚本中的 PLINK_BIN 变量。"))
} else {
  log_msg(paste("已定位 PLINK:", PLINK_BIN))
}

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
  library(ieugwasr)
  library(data.table)
  library(parallel)
})

setwd(WORK_DIR)

# ==============================================================================
# 1. 加载与 LiftOver (与之前相同，略去详细注释)
# ==============================================================================
# # log_msg("正在加载数据并进行 LiftOver...")
# dt <- qs::qread(INPUT_FILE)
# chain <- import.chain(CHAIN_PATH)
# dt_unique <- dt %>% dplyr::distinct(SNP, chr.exposure, pos.exposure)
# 
# gr_hg38 <- GRanges(
#   seqnames = paste0("chr", dt_unique$chr.exposure),
#   ranges = IRanges(start = dt_unique$pos.exposure, end = dt_unique$pos.exposure),
#   snp_id = dt_unique$SNP
# )
# results_gr <- unlist(liftOver(gr_hg38, chain))
# df_hg19_map <- data.frame(SNP = mcols(results_gr)$snp_id, pos.exposure.hg19 = start(results_gr), stringsAsFactors = FALSE)
# dt_lifted <- dt %>% inner_join(df_hg19_map, by = "SNP") %>%
#   mutate(pos.exposure = pos.exposure.hg19) %>% select(-pos.exposure.hg19)
# 
# # 过滤 BIM
# valid_snps <- fread(BIM_FILE, select = 2, header = FALSE)$V2
# dt_valid <- dt_lifted %>% filter(SNP %in% valid_snps)
# 
# log_msg(paste("预处理完成，进入 Clump 阶段。待处理 SNP:", nrow(dt_valid)))
dt_valid <- qs::qread(INPUT_FILE)
# ==============================================================================
# 2. 定义 Clump 函数 (修正版)
# ==============================================================================
safe_clump <- function(sub_df) {
  if(nrow(sub_df) == 0) return(NULL)
  
  tryCatch({
    res <- ld_clump(
      dplyr::tibble(rsid=sub_df$SNP, pval=sub_df$pval.exposure, id=sub_df$id.exposure),
      clump_r2 = CLUMPR,
      clump_kb = 10000,
      clump_p = 0.99, 
      plink_bin = PLINK_BIN,  # <--- 这里显式调用了绝对路径
      bfile = PLINK_BFILE
    )
    return(res)
  }, error = function(e) {
    # 打印错误信息到标准错误输出
    message(paste0("ERROR in ID ", unique(sub_df$id.exposure)[1], ": ", e$message))
    return(NULL)
  })
}

# ==============================================================================
# 3. 冒烟测试 (Smoke Test) - 防止跑完几万个才发现全是错的
# ==============================================================================
log_msg("正在进行 Clump 冒烟测试 (测试前 20 个 ID)...")

test_ids <- unique(dt_valid$id.exposure)[1:20]
test_data_list <- split(dt_valid %>% filter(id.exposure %in% test_ids), 
                        dt_valid %>% filter(id.exposure %in% test_ids) %>% pull(id.exposure))

# 单核运行测试，确保能捕获错误
test_res <- lapply(test_data_list, safe_clump)
test_df <- bind_rows(test_res)

if (nrow(test_df) == 0) {
  stop("严重错误：测试运行结果为空！\n这通常意味着 PLINK 运行失败。\n请检查日志中是否有 'Command not found' 或权限错误。\n请检查 EUR 参考面板路径是否正确。")
} else {
  log_msg(paste("测试通过！测试样本 Clump 后保留 SNP:", nrow(test_df)))
  log_msg("开始进行全量并行计算...")
}

# ==============================================================================
# 4. 全量并行计算
# ==============================================================================
dt_list <- split(dt_valid, dt_valid$id.exposure)

# 并行
clumped_results_list <- mclapply(dt_list, safe_clump, mc.cores = N_CORES)

clumped_df <- bind_rows(clumped_results_list)

log_msg(paste("Clump 完成。保留独立 SNP 总数:", nrow(clumped_df)))

log_msg(paste("Clump 完成。原有独立exposure 总数:", length(unique(dt_valid$id.exposure))))
# 保存
dt_final <- dt_valid %>%
  filter(SNP %in% clumped_df$rsid & id.exposure %in% clumped_df$id)
log_msg(paste("Clump 完成。保留独立exposure 总数:", length(unique(dt_final$id.exposure))))

qs::qsave(dt_final, paste0(WORK_DIR ,"/",OUTPUT_FILE))
log_msg("任务完成。")
