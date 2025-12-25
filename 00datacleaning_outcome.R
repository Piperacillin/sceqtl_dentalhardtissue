####meta PAP####
rm(list = ls())
library(TwoSampleMR)
library(data.table)

setwd('/home/piperacillin/publicdata/datasets/DENTALTRAITS/PAP')
# # 读取命令行参数
# args <- commandArgs(trailingOnly = TRUE)
file <- 'meta_analysis_mvp_ukbb_summary_stats_K11_PULP_PERIAPICAL_meta_out.tsv.gz'
output_file <- 'PULP_PERIAPICAL.qs'

finn_info <- fread('./meta_analysis_mvp_ukbb_FinnGen_R12_MVP_UKBB_mapping.tsv',data.table = F) %>%
  dplyr::filter(fg_phenotype =='K11_PULP_PERIAPICAL')

#paste0(tools::file_path_sans_ext(basename(gz_files)), ".gz") basename截取最后一段

# 循环读取并保存成rdata文件

finn <- fread(file,data.table = F)
head(finn)
finn2 <- finn
names(finn2)[names(finn2) == "#CHR"] <- "CHR"


cols_needed <- c("rsid",
                 "leave_MVP_AFR_inv_var_meta_beta",
                 "leave_MVP_AFR_inv_var_meta_sebeta",
                 "leave_MVP_AFR_inv_var_meta_p",
                 "ALT","REF","meta_eaf","CHR","POS","leave_MVP_AFR_N")

finn_sub <- finn2[, intersect(cols_needed, names(finn2))]

# 去除明显缺失或无效行
finn_sub <- subset(finn_sub,
                   !is.na(rsid) &
                     !is.na(leave_MVP_AFR_inv_var_meta_beta) &
                     !is.na(leave_MVP_AFR_inv_var_meta_sebeta) &
                     !is.na(leave_MVP_AFR_inv_var_meta_p))

# 去重：对重复 rsid 先根据 p 值选最显著一条，防止 format_data 内部递归去重
finn_sub <- finn_sub[order(finn_sub$rsid, finn_sub$leave_MVP_AFR_inv_var_meta_p), ]
finn_sub <- finn_sub[!duplicated(finn_sub$rsid), ]

outcome_format <- format_data(
  dat = finn_sub,
  type = "outcome",
  snp_col = "rsid",
  beta_col = "leave_MVP_AFR_inv_var_meta_beta",
  se_col = "leave_MVP_AFR_inv_var_meta_sebeta",
  pval_col = "leave_MVP_AFR_inv_var_meta_p",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "meta_eaf",
  chr_col = "CHR",
  pos_col = "POS",
  samplesize_col = "leave_MVP_AFR_N"
)

outcome_format$outcome <- 'PULP_PERIAPICAL'
outcome_format$samplesize.outcome <- finn_info$fg_n_cases + finn_info$fg_n_controls + finn_info$ukbb_n_cases + finn_info$ukbb_n_controls + finn_info$MVP_EUR_n_cases  + finn_info$MVP_EUR_n_controls   + finn_info$MVP_AMR_n_cases   + finn_info$MVP_AMR_n_controls  

head(outcome_format)
unique(outcome_format$samplesize.outcome)

data <- outcome_format
rm(finn,finn_info,finn_sub,finn2,outcome_format)
####lifetover to hg19####
library(rtracklayer) ## BiocManager::install("rtracklayer")
library(GenomicRanges)


# 0. 准备路径和数据
# ===============================
# 假设你的数据框名为 data (已经读入)
CHAIN_PATH <- "/home/piperacillin/publicdata/rawinstruments/hg38ToHg19.over.chain" # 请修改为你的实际路径
# BIM_FILE <- "path/to/eur.bim" # 如果你有参考面板的BIM文件
head(data)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

# 1) 构建 GRanges（确保染色体名带有 'chr' 前缀）
dtu <- data %>% 
  filter(!is.na(chr.outcome), !is.na(pos.outcome), !is.na(SNP)) %>% 
  distinct(SNP, chr.outcome, pos.outcome)

seqs <- ifelse(grepl("^chr", as.character(dtu$chr.outcome)), 
               as.character(dtu$chr.outcome), 
               paste0("chr", dtu$chr.outcome))

gr_hg38 <- GRanges(
  seqnames = seqs,
  ranges   = IRanges(start = dtu$pos.outcome, end = dtu$pos.outcome),
  snp_id   = dtu$SNP
)

# 2) 执行 liftover
chain <- import.chain(CHAIN_PATH)
res_gr <- unlist(liftOver(gr_hg38, chain))

# 3) 合并回 data，更新 pos.outcome 为 hg19，并备份原位置
map_hg19 <- data.frame(
  SNP = mcols(res_gr)$snp_id,
  pos.outcome.hg19 = start(res_gr),
  stringsAsFactors = FALSE
)

data_hg19 <- data %>%
  inner_join(map_hg19, by = "SNP") %>%
  mutate(pos.outcome.hg38_backup = pos.outcome,
         pos.outcome = pos.outcome.hg19) %>%
  dplyr::select(-pos.outcome.hg19)

head(data_hg19)

# data_hg19 即为 liftover 后的数据框

log_msg(paste("LiftOver 完成。原始数量:", nrow(data), "-> 转换后数量:", nrow(data_hg19)))

# 保存
qs::qsave(data_hg19, file = output_file)


####Periapical####
rm(list = ls())
library(TwoSampleMR)
library(data.table)

setwd('/home/piperacillin/publicdata/datasets/DENTALTRAITS/PAP')
# # 读取命令行参数
# args <- commandArgs(trailingOnly = TRUE)
file <- 'necrosis_or_apical_periodontitis.gz'
output_file <- 'NECROSIS_PERIAPICAL.qs'

finn_info <- fread(file,data.table = T) 
head(finn_info)

library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38) 

# 你的列名带有 #，fread 可能会将其读为 "#chrom" 或者自动去掉 #
raw_dt <- finn_info
# 标准化列名：将 "#chrom" 改为 "CHR"，"pos" 改为 "POS"
# 这里处理以此防止列名带 '#' 导致的引用错误
colnames(raw_dt)[1] <- "CHR" 
if("pos" %in% colnames(raw_dt)) setnames(raw_dt, "pos", "POS")
# 确保 CHR 列没有 "chr" 前缀 (dbSNP 包通常需要无前缀或标准化前缀，我们手动统一处理)
head(raw_dt)
raw_dt[, CHR := as.character(CHR)]
raw_dt[, CHR := sub("^chr", "", CHR, ignore.case = TRUE)]

# ==========================================
# 1. 修正染色体编码 (把 23 改为 X)
# ==========================================

# 手动替换非标准编码
# 注意：R data.table 的引用更新写法
raw_dt[CHR == "23", CHR := "X"]
raw_dt[CHR == "24", CHR := "Y"]  # 如果有 Y 染色体
raw_dt[CHR == "25" | CHR == "M" | CHR == "MT", CHR := "MT"] # 线粒体

# ==========================================
# 2. 重新构建 GRanges 对象
# ==========================================
gr_query <- GRanges(
  seqnames = raw_dt$CHR,
  ranges = IRanges(start = raw_dt$POS, end = raw_dt$POS)
)

# ==========================================
# 3. 设定 seqlevelsStyle 并查询
# ==========================================
# 这一步告诉 GRanges 把染色体名统一成 NCBI 标准 (1, 2, ..., X, Y)
seqlevelsStyle(gr_query) <- "NCBI"

# 再次运行查询
# 使用 suppressWarnings 是因为数据库包含很多 contig (如 alt scafolds)，
# 而你的数据只有主染色体，这会导致"sequence levels not in the other"警告，这是无害的。
message("正在重新查询...")
found_snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr_query)

message(paste("查询成功！找到 SNP 数量:", length(found_snps)))
head(found_snps)
# 将查询结果转换为 data.table
# found_snps 包含 seqnames, pos, RefSNP_id (rsID), alleles_as_ambig
dt_snps <- as.data.table(found_snps)
dt_snps <- dt_snps[, .(CHR = as.character(seqnames), POS = pos, rsid = RefSNP_id)]
# 将查到的 rsID 合并回原始数据
# 注意：一个位置可能对应多个 rsID (极其罕见的多等位位点重叠)，这里为了简化取合并
# 如果你非常在意 Ref/Alt 匹配，这里需要更复杂的校验，但在 MR 中通常直接基于位置匹配即可
head(dt_snps)

dt_annotated <- merge(raw_dt, dt_snps, by = c("CHR", "POS"), all.x = TRUE)
# 检查有多少没匹配到
missing_count <- sum(is.na(dt_annotated$rsid))
message(paste("已补全 rsID。未找到 rsID 的行数:", missing_count))
# 如果必须要有 rsID 才能做 MR，可以过滤掉 NA
dt_annotated <- dt_annotated[!is.na(rsid)]
# ===========================
# 3. LiftOver (hg38 -> hg19)
# ===========================
message("正在进行 LiftOver (hg38 -> hg19)...")
# 再次构建 GRanges，这次是为了 LiftOver
# LiftOver 的 Chain 文件通常要求 seqnames 为 "chr1", "chr2" 格式
gr_hg38 <- GRanges(
  seqnames = paste0("chr", dt_annotated$CHR),
  ranges = IRanges(start = dt_annotated$POS, end = dt_annotated$POS),
  snp_id = dt_annotated$rsid # 携带刚补全的 rsID
)
# 导入 Chain 文件 (请确保路径正确)
# 下载地址: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
CHAIN_PATH <- "/home/publicdata/rawinstruments/hg38ToHg19.over.chain" 
chain <- import.chain(CHAIN_PATH)
results_gr <- unlist(liftOver(gr_hg38, chain))
# 提取转换后的结果
df_hg19_map <- data.frame(
  rsid = mcols(results_gr)$snp_id,
  POS_hg19 = start(results_gr),
  stringsAsFactors = FALSE
)
# ===========================
# 4. 合并最终结果
# ===========================
final_data <- dt_annotated %>%
  inner_join(df_hg19_map, by = "rsid") %>%
  mutate(
    POS_hg38 = POS,    # 备份旧坐标
    POS = POS_hg19     # 更新为 hg19 坐标，作为主要 POS 列
  ) %>%
  select(-POS_hg19)    # 移除临时列
# 此时 final_data 拥有了 rsID，且 POS 已经是 hg19
head(final_data)

# 去重：对重复 rsid 先根据 p 值选最显著一条，防止 format_data 内部递归去重
final_data <- final_data[order(final_data$rsid, final_data$pval), ]
final_data <- final_data[!duplicated(final_data$rsid), ]
final_data$ncase <- 103832
final_data$nctrl <- 353106
head(final_data)

outcome_format <- format_data(
  dat = final_data,
  type = "outcome",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "sebeta",
  pval_col = "pval",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  chr_col = "CHR",
  pos_col = "POS",
  ncase_col = "ncase",  #48,120
  ncontrol_col = "nctrl"
)

outcome_format$outcome <- 'NECROSIS_PERIAPICALL'
#outcome_format$samplesize.outcome <- finn_info$fg_n_cases + finn_info$fg_n_controls + finn_info$ukbb_n_cases + finn_info$ukbb_n_controls + finn_info$MVP_EUR_n_cases  + finn_info$MVP_EUR_n_controls   + finn_info$MVP_AMR_n_cases   + finn_info$MVP_AMR_n_controls  

head(outcome_format)
unique(outcome_format$samplesize.outcome)
output_file
# 保存
qs::qsave(outcome_format, file = output_file)

####PULPITIS####
rm(list = ls())
library(TwoSampleMR)
library(data.table)

setwd('/home/piperacillin/publicdata/datasets/DENTALTRAITS/PAP')
# # 读取命令行参数
# args <- commandArgs(trailingOnly = TRUE)
file <- 'pulpitis.gz'
output_file <- 'PULPITIS.qs'

finn_info <- fread(file,data.table = T) 
head(finn_info)

library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38) 

# 你的列名带有 #，fread 可能会将其读为 "#chrom" 或者自动去掉 #
raw_dt <- finn_info
# 标准化列名：将 "#chrom" 改为 "CHR"，"pos" 改为 "POS"
# 这里处理以此防止列名带 '#' 导致的引用错误
colnames(raw_dt)[1] <- "CHR" 
if("pos" %in% colnames(raw_dt)) setnames(raw_dt, "pos", "POS")
# 确保 CHR 列没有 "chr" 前缀 (dbSNP 包通常需要无前缀或标准化前缀，我们手动统一处理)
head(raw_dt)
raw_dt[, CHR := as.character(CHR)]
raw_dt[, CHR := sub("^chr", "", CHR, ignore.case = TRUE)]

# ==========================================
# 1. 修正染色体编码 (把 23 改为 X)
# ==========================================

# 手动替换非标准编码
# 注意：R data.table 的引用更新写法
raw_dt[CHR == "23", CHR := "X"]
raw_dt[CHR == "24", CHR := "Y"]  # 如果有 Y 染色体
raw_dt[CHR == "25" | CHR == "M" | CHR == "MT", CHR := "MT"] # 线粒体

# ==========================================
# 2. 重新构建 GRanges 对象
# ==========================================
gr_query <- GRanges(
  seqnames = raw_dt$CHR,
  ranges = IRanges(start = raw_dt$POS, end = raw_dt$POS)
)

# ==========================================
# 3. 设定 seqlevelsStyle 并查询
# ==========================================
# 这一步告诉 GRanges 把染色体名统一成 NCBI 标准 (1, 2, ..., X, Y)
seqlevelsStyle(gr_query) <- "NCBI"

# 再次运行查询
# 使用 suppressWarnings 是因为数据库包含很多 contig (如 alt scafolds)，
# 而你的数据只有主染色体，这会导致"sequence levels not in the other"警告，这是无害的。
message("正在重新查询...")
found_snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr_query)

message(paste("查询成功！找到 SNP 数量:", length(found_snps)))
head(found_snps)
# 将查询结果转换为 data.table
# found_snps 包含 seqnames, pos, RefSNP_id (rsID), alleles_as_ambig
dt_snps <- as.data.table(found_snps)
dt_snps <- dt_snps[, .(CHR = as.character(seqnames), POS = pos, rsid = RefSNP_id)]
# 将查到的 rsID 合并回原始数据
# 注意：一个位置可能对应多个 rsID (极其罕见的多等位位点重叠)，这里为了简化取合并
# 如果你非常在意 Ref/Alt 匹配，这里需要更复杂的校验，但在 MR 中通常直接基于位置匹配即可
head(dt_snps)

dt_annotated <- merge(raw_dt, dt_snps, by = c("CHR", "POS"), all.x = TRUE)
# 检查有多少没匹配到
missing_count <- sum(is.na(dt_annotated$rsid))
message(paste("已补全 rsID。未找到 rsID 的行数:", missing_count))
# 如果必须要有 rsID 才能做 MR，可以过滤掉 NA
dt_annotated <- dt_annotated[!is.na(rsid)]
# ===========================
# 3. LiftOver (hg38 -> hg19)
# ===========================
message("正在进行 LiftOver (hg38 -> hg19)...")
# 再次构建 GRanges，这次是为了 LiftOver
# LiftOver 的 Chain 文件通常要求 seqnames 为 "chr1", "chr2" 格式
gr_hg38 <- GRanges(
  seqnames = paste0("chr", dt_annotated$CHR),
  ranges = IRanges(start = dt_annotated$POS, end = dt_annotated$POS),
  snp_id = dt_annotated$rsid # 携带刚补全的 rsID
)
# 导入 Chain 文件 (请确保路径正确)
# 下载地址: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
CHAIN_PATH <- "/home/publicdata/rawinstruments/hg38ToHg19.over.chain" 
chain <- import.chain(CHAIN_PATH)
results_gr <- unlist(liftOver(gr_hg38, chain))
# 提取转换后的结果
df_hg19_map <- data.frame(
  rsid = mcols(results_gr)$snp_id,
  POS_hg19 = start(results_gr),
  stringsAsFactors = FALSE
)
# ===========================
# 4. 合并最终结果
# ===========================
final_data <- dt_annotated %>%
  inner_join(df_hg19_map, by = "rsid") %>%
  mutate(
    POS_hg38 = POS,    # 备份旧坐标
    POS = POS_hg19     # 更新为 hg19 坐标，作为主要 POS 列
  ) %>%
  select(-POS_hg19)    # 移除临时列
# 此时 final_data 拥有了 rsID，且 POS 已经是 hg19
head(final_data)

# 去重：对重复 rsid 先根据 p 值选最显著一条，防止 format_data 内部递归去重
final_data <- final_data[order(final_data$rsid, final_data$pval), ]
final_data <- final_data[!duplicated(final_data$rsid), ]
final_data$ncase <- 48120
final_data$nctrl <- 353106
head(final_data)

outcome_format <- format_data(
  dat = final_data,
  type = "outcome",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "sebeta",
  pval_col = "pval",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  chr_col = "CHR",
  pos_col = "POS",
  ncase_col = "ncase",  #48,120
  ncontrol_col = "nctrl"
)

outcome_format$outcome <- 'PULPITIS'
#outcome_format$samplesize.outcome <- finn_info$fg_n_cases + finn_info$fg_n_controls + finn_info$ukbb_n_cases + finn_info$ukbb_n_controls + finn_info$MVP_EUR_n_cases  + finn_info$MVP_EUR_n_controls   + finn_info$MVP_AMR_n_cases   + finn_info$MVP_AMR_n_controls  

head(outcome_format)
unique(outcome_format$samplesize.outcome)
output_file
# 保存
qs::qsave(outcome_format, file = output_file)