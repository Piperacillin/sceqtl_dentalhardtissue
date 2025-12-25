getwd()
library(data.table)
library(readr)
library(stringr)
library(dplyr)
library(TwoSampleMR)

####1mscblood####
rm(list = ls())
setwd('/home/publicdata/datasets/eqtl/1mscblood')

# 根目录
base_dir <- "/home/publicdata/datasets/eqtl/1mscblood"
in_dir   <- file.path(base_dir, "filtered_p1e5")
# 列出所有条件文件夹（24hCA, 24hMTB, ...）
cell_dirs <- list.dirs(in_dir, full.names = TRUE, recursive = FALSE)
# 从ProbeName提取第一个'_'之前的部分
probe_prefix <- function(x) {
  # 对没有'_'的情况保持原值
  sub("^([^_]+)_.*$", "\\1", x)
}
# 从文件名提取细胞类型（B/DC/NK/CD8T/monocyte/CD4T）
cell_type_from_filename <- function(fn) {
  sub("^([^_]+)_.*$", "\\1", basename(fn))
}
# 读取一个文件并添加exposure列
read_one <- function(fpath, folder_name) {
  # fread自动识别.gz
  dt <- fread(fpath, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  
  # 确保存在ProbeName列（你已确认每个文件有该列）
  if (!"HGNCName" %in% names(dt)) {
    stop(sprintf("File %s is missing 'HGNCName' column", fpath))
  }
  
  # 生成exposure
  # 组成：HGNCName的前缀 + 文件名第一个'_'前的细胞类型 + 文件夹名（如24hCA）
  file_cell  <- cell_type_from_filename(fpath)
  probe_cell <- probe_prefix(dt[["HGNCName"]])
  folder_lab <- folder_name
  
  dt[, exposure := paste0(probe_cell, "_", file_cell, "_", folder_lab)]
  
  # 可选的追踪信息
  dt[, source_file := basename(fpath)]
  dt[, source_dir  := folder_lab]
  dt[, source_cellfile := file_cell]
  
  dt
}
# 遍历并合并
all_list <- list()
for (dir_path in cell_dirs) {
  folder_name <- basename(dir_path)
  files <- list.files(dir_path, pattern = "_filtered\\.txt\\.gz$", full.names = TRUE)
  if (length(files) == 0) next
  for (f in files) {
    all_list[[length(all_list) + 1]] <- read_one(f, folder_name)
  }
}
# 合并为data.table（若没有文件会报错，先判断）
if (length(all_list) == 0) {
  stop("No filtered files found under: ", in_dir)
}
dt_all <- rbindlist(all_list, use.names = TRUE, fill = TRUE)

head(dt_all)

dt <- copy(dt_all)
# 1) 解析“Beta (SE)”列计算两数据集固定效应meta beta/se
parse_beta_se_pair <- function(x) {
  split_pairs <- strsplit(x, ";", fixed = TRUE)
  res <- lapply(split_pairs, function(parts) {
    if (length(parts) < 2) {
      return(list(beta1 = NA_real_, se1 = NA_real_, beta2 = NA_real_, se2 = NA_real_,
                  meta_beta = NA_real_, meta_se = NA_real_))
    }
    parse_one <- function(s) {
      s <- trimws(s)
      m <- stringr::str_match(s, "^\\s*([+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?)\\s*\\(\\s*([0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?)\\s*\\)\\s*$")
      if (is.na(m[1,1])) c(NA_real_, NA_real_) else c(as.numeric(m[1,2]), as.numeric(m[1,4]))
    }
    b1s1 <- parse_one(parts[1]); b2s2 <- parse_one(parts[2])
    beta1 <- b1s1[1]; se1 <- b1s1[2]
    beta2 <- b2s2[1]; se2 <- b2s2[2]
    if (is.na(beta1) || is.na(se1) || is.na(beta2) || is.na(se2) || se1 <= 0 || se2 <= 0) {
      meta_beta <- NA_real_; meta_se <- NA_real_
    } else {
      w1 <- 1 / (se1^2); w2 <- 1 / (se2^2)
      meta_beta <- (beta1*w1 + beta2*w2) / (w1 + w2)
      meta_se   <- sqrt(1 / (w1 + w2))
    }
    list(beta1 = beta1, se1 = se1, beta2 = beta2, se2 = se2, meta_beta = meta_beta, meta_se = meta_se)
  })
  data.table::rbindlist(res)
}
bt <- parse_beta_se_pair(dt[["Beta (SE)"]])
dt[, `:=`(beta1 = bt$beta1, se1 = bt$se1, beta2 = bt$beta2, se2 = bt$se2,
          beta = bt$meta_beta, se = bt$meta_se)]
# 2) 整理基本列
if ("PValue" %in% names(dt)) data.table::setnames(dt, "PValue", "pval")
if ("SNPName" %in% names(dt)) data.table::setnames(dt, "SNPName", "rsid")
# 解析other_allele
parse_other_allele <- function(snp_type, effect_allele) {
  mapply(function(v, ea) {
    alleles <- toupper(unlist(strsplit(toupper(v), "/", fixed = TRUE)))
    ea <- toupper(ea)
    other <- alleles[alleles != ea]
    if (length(other) == 0) NA_character_ else other[1]
  }, snp_type, effect_allele, USE.NAMES = FALSE)
}
dt[, alt := toupper(AlleleAssessed)]
dt[, ref := toupper(parse_other_allele(SNPType, AlleleAssessed))]
# 样本量求和
dt[, ma_samples := sapply(strsplit(DatasetsNrSamples, ";", fixed = TRUE), function(v) {
  v <- suppressWarnings(as.numeric(v))
  if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE)
})]
# 3) 用 eQTLGen 的频率文件补齐 eaf
# 读取 eQTLGen 频率文件：包含SNP, AlleleB, AlleleB_all（AlleleB的频率）
alleles_eqtlgen <- read_tsv(
  gzfile("/home/piperacillin/publicdata/eQTL/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"),
  col_types = cols(.default = "c")
) %>% as.data.table()
alleles_eqtlgen <- alleles_eqtlgen[, .(SNP, AlleleB, AlleleB_all)]
alleles_eqtlgen[, `:=`(SNP = as.character(SNP),
                       AlleleB = toupper(AlleleB),
                       AlleleB_all = suppressWarnings(as.numeric(AlleleB_all)))]
# 对齐key（你的dt使用rsid列）
dt[, rsid := as.character(rsid)]
# 合并
dt_m <- merge(dt, alleles_eqtlgen, by.x = "rsid", by.y = "SNP", all.x = TRUE)
# 仅保留标准碱基并构造eaf（以alt=效应等位基因为准）
is_base <- function(x) !is.na(x) & x %in% c("A","C","G","T")
dt_m[, `:=`(alt = toupper(alt), ref = toupper(ref))]
dt_m[!(is_base(alt) & is_base(ref) & is_base(AlleleB)), AlleleB_all := NA_real_]
# 方向对齐：如果 alt == AlleleB -> eaf = AlleleB_all；若 alt == 对应另一等位基因 -> eaf = 1 - AlleleB_all
dt_m[, eaf_from_eqtlgen := NA_real_]
dt_m[!is.na(AlleleB_all) & alt == AlleleB, eaf_from_eqtlgen := AlleleB_all]
dt_m[!is.na(AlleleB_all) & alt != AlleleB & !is.na(ref) & ref == AlleleB, eaf_from_eqtlgen := 1 - AlleleB_all]
# 范围检查
dt_m[!(eaf_from_eqtlgen > 0 & eaf_from_eqtlgen < 1), eaf_from_eqtlgen := NA_real_]
# 若你的dt已有eaf列（部分来源于别处），优先已有；否则用eQTLGen补齐
if (!"eaf" %in% names(dt_m)) dt_m[, eaf := NA_real_]
dt_m[is.na(eaf) & !is.na(eaf_from_eqtlgen), eaf := eaf_from_eqtlgen]
# 4) 基本QC：可选剔除回文位点，或仅在eaf远离0.5时保留
is_palindromic <- function(a1, a2) {
  sapply(seq_along(a1), function(i) {
    pair <- paste0(sort(c(a1[i], a2[i])), collapse = "/")
    pair %in% c("A/T", "C/G")
  })
}
dt_m[, pal := is_palindromic(alt, ref)]

# 你可以选择保留全部，或剔除回文
# dt_m <- dt_m[pal == FALSE | (eaf <= 0.42 | eaf >= 0.58)]
# 5) 继续你的筛选（例如UT来源）
dt_ut <- dt_m[source_dir %in% "UT"]
# 6) 使用TwoSampleMR格式化（现在提供eaf_col）
exposure_dat <- format_data(
  dt_m,
  type = "exposure",
  header = TRUE,
  phenotype_col = "exposure",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  samplesize_col = "ma_samples"
)
head(exposure_dat)

qs::qsave(exposure_dat, '1mscblood_filtered_exposure.qs')
length(unique(exposure_dat$exposure))

####ONEK1K####
rm(list = ls())
setwd('~/publicdata/datasets/eqtl/onek1k')
onek1k <- fread('./eqtl_p1e5.tsv.gz')
head(onek1k)

N_onek1k <- 982  # 请用真实 N 替换

#构造格式化前的数据：effect_allele = A2；other_allele = A1；eaf 用 A2 的频率
dt_m <- onek1k[
  !is.na(RSID),
  .(
    exposure = paste(GENE, CELL_TYPE, sep = "_"),
    rsid     = RSID,
    beta     = SPEARMANS_RHO,
    se       = (1 - SPEARMANS_RHO^2) / sqrt(N_onek1k - 3),
    eaf      = fcoalesce(A2_FREQ_ONEK1K, A2_FREQ_HRC),
    alt      = A2,
    ref      = A1,
    pval     = P_VALUE,
    chr = CHR,
    pos = POS,
    ma_samples = N_onek1k,
    gene_query_id = GENE
  )
]
# 统一大写（正确语法）
dt_m[, c("alt", "ref") := .(toupper(alt), toupper(ref))]
# 或：dt_m[, `:=`(alt = toupper(alt), ref = toupper(ref))]
# 基础过滤
dt_m <- dt_m[!is.na(beta) & !is.na(se) & !is.na(eaf) & !is.na(pval)]
dt_m <- dt_m[eaf > 0 & eaf < 1]
head(dt_m)
# TwoSampleMR 格式
exposure_dat <- format_data(
  dt_m,
  type = "exposure",
  header = TRUE,
  phenotype_col = "exposure",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  samplesize_col = "ma_samples",
  chr_col = "chr",
  pos_col = "pos",
  gene_col = "gene_query_id"
)


head(exposure_dat)
qs::qsave(exposure_dat, 'onek1k_filtered_exposure.qs')

####dynamicCD4T####
rm(list=ls())
library(arrow)
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(TwoSampleMR)
library(GenomicRanges)
library(biomaRt) # 用于获取基因坐标

# 确保安装了必要的注释包
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

# 1. 批量读取并合并数据

setwd("~/publicdata/datasets/eqtl/CD4Tdynamic")
# 获取文件列表
files <- list.files(pattern = "_500kb_combined\\.parquet$")
cat("检测到", length(files), "个文件，开始批量读取...\n")
# 使用 lapply 批量读取并处理
data_list <- lapply(files, function(f) {
  # 1. 提取条件后缀 (去除 _500kb_combined.parquet)
  condition_suffix <- stringr::str_remove(f, "_500kb_combined\\.parquet")
  
  # 2. 读取数据 (仅读取必要列以节省内存)
  # 假设 parquet 中包含: phenotype_id, variant_id, slope, slope_se, pval_nominal, maf (或类似列)
  df <- read_parquet(f) %>%
    select(phenotype_id, variant_id, slope, slope_se, pval_nominal, ma_samples, maf) %>% # 确保这里包含了 maf 列
    mutate(
      # 修改 phenotype_id
      phenotype_id = paste0(phenotype_id, "_", condition_suffix)
    )
  
  return(df)
})
# 合并为一个巨大的 Data Table
cat("正在合并所有数据框...\n")
full_dat <- rbindlist(data_list)
cat("合并完成。总行数:", nrow(full_dat), "\n")

# ==============================================================================
# 2. 预处理坐标 (拆分 variant_id)
# ==============================================================================
# 假设 variant_id 格式为 chr_pos_ref_alt (例如 1_12345_A_G)
full_dat <- full_dat %>%
  separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = "_", remove = FALSE) %>%
  mutate(
    chr = as.character(chr),
    pos = as.numeric(pos)
  )

# ==============================================================================
# 3. 高效查询 RSID (仅对 Unique 位点查询)
# ==============================================================================
cat("提取唯一位点进行 RSID 注释...\n")
# 提取不重复的 chr:pos 组合
unique_sites <- full_dat %>%
  select(chr, pos) %>%
  distinct()
cat("需注释的唯一位点数:", nrow(unique_sites), "\n")
# 创建 GRanges 对象
# 注意：SNPlocs.Hsapiens.dbSNP155.GRCh38 使用 "ch" 或 "1" 格式取决于版本，
# 通常 SNPlocs 输入需要与 seqlevels(snplocs) 匹配。
# 这里先尝试直接转换，如果报错可能需要加 "chr" 前缀。
gr <- GRanges(seqnames = unique_sites$chr, 
              ranges = IRanges(start = unique_sites$pos, end = unique_sites$pos))
# 修正 seqlevels (如果 SNPlocs 需要 "chr1" 而你的是 "1")
#seqlevelsStyle(gr) <- "UCSC" # 强制转为 chr1, chr2... 格式以匹配 SNPlocs
cat("正在查询 SNPlocs (这可能需要几分钟)...\n")
snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)
head(snps)
# 提取映射表
rsid_map <- as.data.frame(snps) %>%
  dplyr::select(seqnames, pos, RefSNP_id) %>%
  dplyr::rename(chr_std = seqnames, rsid = RefSNP_id) %>%
  mutate(
    # 将 chr_std (如 chr1) 转回原始格式 (如 1) 以便合并
    chr = gsub("chr", "", chr_std),
    pos = as.numeric(pos)
  ) %>%
  select(chr, pos, rsid)
head(rsid_map)
# ==============================================================================
# 4. 将 RSID 合并回总数据
# ==============================================================================
cat("将 RSID 合并回原始大数据...\n")
# 使用 data.table 的合并更高效
full_dat_dt <- as.data.table(full_dat)
rsid_map_dt <- as.data.table(rsid_map)
# Left join
final_dat <- merge(full_dat_dt, rsid_map_dt, by = c("chr", "pos"), all.x = TRUE)

# 计算以 "rs" 开头的数量
n_match <- sum(grepl("^rs", final_dat$rsid))
n_total <- nrow(final_dat)
# 输出数量和百分比
cat("RSIDs starting with 'rs':", n_match, "/", n_total, 
    sprintf("(%.2f%%)\n", (n_match / n_total) * 100))


# 填充缺失的 RSID
final_dat[is.na(rsid), rsid := paste0("chr", chr, ":", pos)]
cat("RSID 注释完成。\n")
head(final_dat)

# ==============================================================================
# 1. 初步过滤 (P值 和 F统计量)
# ==============================================================================
# 这一步放在最前面，可以极大减少后续需要注释基因坐标的行数
cat("正在进行 P 值和 F 统计量过滤...\n")
# 确保 final_dat 是 data.table 格式
setDT(final_dat)
# 计算 F 统计量并过滤
# 注意：OneK1K 数据通常已经是 cis-eQTL (文件名含 500kb)，但我们需要确保满足你的 1000kb 要求
final_dat[, F_stat := (slope / slope_se)^2]
# 执行过滤
filtered_dat <- final_dat[pval_nominal < 1e-5 & F_stat >= 10]
cat("过滤前行数:", nrow(final_dat), "\n")
cat("过滤后行数:", nrow(filtered_dat), "\n")
# ==============================================================================
# 2. 获取基因坐标 (用于 Cis 定义)
# ==============================================================================
cat("正在获取基因坐标以验证 Cis 距离...\n")
# 2.1 提取纯净的 Gene ID
# 你的 phenotype_id 现在长这样: "ENSG00000123456_TN_NFKB"
# 我们需要提取前面的 "ENSG00000123456" 去查询坐标
# 假设原 phenotype_id 不包含下划线，或者它是 Ensembl ID。
# 如果原 ID 是 Gene Symbol (如 STAT4)，逻辑也类似。
# 提取第一个 "_" 之前的部分作为查询 ID (根据你的实际 ID 格式调整)
# 这里假设原 phenotype_id 是 Ensembl ID (OneK1K 常用格式)
filtered_dat[, gene_query_id := sub("_.*", "", phenotype_id)] 
unique_genes <- unique(filtered_dat$gene_query_id)
cat("待查询基因数:", length(unique_genes), "\n")
# 2.2 使用 biomaRt 获取 hg38 坐标
# 如果网络不好，这一步可能会慢。也可以使用 AnnotationDbi
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # 默认是 hg38

gene_coords <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique_genes,
  mart = ensembl
)
setDT(gene_coords)
# 2.3 计算 TSS (转录起始位点)
# 正链(+) TSS = start_position
# 负链(-) TSS = end_position
gene_coords[, tss := ifelse(strand == 1, start_position, end_position)]
gene_coords[, chromosome_name := as.character(chromosome_name)]
saveRDS(gene_coords, "gene_coords_hg38.rds") ####保存联网部分####
# ==============================================================================
# 3. 合并坐标并计算距离
# ==============================================================================
cat("合并坐标并计算 SNP-Gene 距离...\n")
# 合并
# 注意：filtered_dat 中的 chr 可能是 "1", "2"，biomaRt 也是 "1", "2" (无chr前缀)
# 需确保格式一致
cis_check <- merge(
  filtered_dat, 
  gene_coords[, .(ensembl_gene_id, chr_gene = chromosome_name, tss)], 
  by.x = "gene_query_id", 
  by.y = "ensembl_gene_id",
  all.x = FALSE # 丢弃查不到坐标的基因
)
# 计算距离并判定 Cis
# 条件1: 染色体相同
# 条件2: 距离 <= 1000kb (1,000,000 bp)
cis_check[, dist := abs(pos - tss)]
cis_check[, is_cis := (chr == chr_gene) & (dist <= 1000000)]
# 保留 Cis-eQTLs
final_cis_dat <- cis_check[is_cis == TRUE]
cat("Cis 过滤完成。剩余 SNP 数量:", nrow(final_cis_dat), "\n")
# ==============================================================================
# 4. 最终格式化 (输出给 TwoSampleMR)
# ==============================================================================
# 此时 final_cis_dat 已经满足：
# 1. P < 1e-5
# 2. F >= 10
# 3. 距离 TSS 1Mb 以内 (Cis)
# 4. 拥有 RSID (来自上一步)
exposure_dat <- format_data(
  final_cis_dat,
  type = "exposure",
  header = TRUE,
  phenotype_col = "phenotype_id",
  snp_col = "rsid",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",             # 再次提醒：这里用的是 MAF，后续需注意方向
  effect_allele_col = "alt",
  other_allele_col  = "ref",
  pval_col = "pval_nominal",
  samplesize_col = "ma_samples",
  chr_col = "chr",
  pos_col = "pos",
  gene_col = "gene_query_id"   # 可选：保留基因名列
)
# 添加 F 统计量到输出对象中 (TwoSampleMR 默认不带 F 列，手动加上方便查看)
exposure_dat$F_stat <- final_cis_dat$F_stat
head(exposure_dat)

# 1. 准备映射表 (去除 Symbol 为空的行)
# 假设 gene_coords 是上一步 biomaRt 下载的结果
gene_map <- gene_coords %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol) %>%
  dplyr::distinct() %>%
  dplyr::filter(hgnc_symbol != "") # 过滤掉没有 Symbol 的记录
# 2. 执行转换
# 逻辑：将 exposure_dat 与映射表合并。
# 如果能匹配到 Symbol，则替换；如果匹配不到（NA或空），则保留原来的 Ensembl ID。
exposure_dat <- exposure_dat %>%
  left_join(gene_map, by = c("gene.exposure" = "ensembl_gene_id")) %>%
  mutate(
    # 优先使用 hgnc_symbol，如果为空则保持原样(Ensembl ID)
    gene.exposure = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol, gene.exposure)
  ) %>%
  dplyr::select(-hgnc_symbol) # 删除辅助列

head(exposure_dat)
cat("unique eQTL:", length(unique(exposure_dat$exposure)), "\n")
cat("unique gene:", length(unique(exposure_dat$gene.exposure)), "\n")
# 保存最终的大文件 (建议保存为 qs 或 rds 格式，比 csv 快且小)
qs::qsave(exposure_dat, "CD4_Dynamic_Exposure.qs")

setwd("~/publicdata/datasets/eqtl/CD4Tdynamic")
dt <- qs::qread("CD4_Dynamic_Exposure.qs")
library(rtracklayer)
library(GenomicRanges)
head(dt)

# 1. 加载 Chain 文件 (R会自动处理 .gz 解压)
chain_path <- "/home/piperacillin/publicdata/rawinstruments/hg38ToHg19.over.chain"
chain <- import.chain(chain_path)

# 2. 构建 GRanges 对象 (注意：Chain文件通常要求染色体带 "chr" 前缀)
# 确保你的 chr.exposure 没有 "chr" 前缀，如果有的话需要调整代码
gr_hg38 <- GRanges(
  seqnames = paste0("chr", dt$chr.exposure), 
  ranges = IRanges(start = dt$pos.exposure, end = dt$pos.exposure),
  snp_id = dt$SNP  # 携带 SNP ID 以便后续合并
)

# 3. 执行 LiftOver
# liftOver 返回的是一个列表，因为某些位置可能映射失败或映射到多个位置
results_list <- liftOver(gr_hg38, chain)
results_gr <- unlist(results_list)

# 4. 提取转换后的数据并整理
# 提取 SNP ID 和 新的 hg19 位置
df_hg19 <- data.frame(
  SNP = mcols(results_gr)$snp_id,
  pos.exposure.hg19 = start(results_gr),
  stringsAsFactors = FALSE
)
head(df_hg19)
# 5. 将新坐标合并回原始数据框 dt
# inner_join 会自动去除那些 LiftOver 失败（即在 hg19 中找不到对应位置）的 SNP
dt_lifted <- dt %>%
  # 使用 inner_join 自动过滤掉无法转换坐标的 SNP
  dplyr::inner_join(df_hg19, by = "SNP") %>%
  dplyr::mutate(
    pos.exposure.hg38 = pos.exposure,      # 备份原始 hg38 坐标
    pos.exposure = pos.exposure.hg19       # 更新主坐标列为 hg19 (给 Clumping 用)
  ) %>%
  dplyr::select(-pos.exposure.hg19)  

# 6. 查看结果
head(dt_lifted)
# 此时 dt_lifted 中的 pos.exposure 是 hg19 的，pos.exposure.hg38 是 hg38 的

# clump liberal analysis
library(ieugwasr)
library(dplyr)
library(data.table)
library(parallel) # 用于并行计算

# ==============================================================================
# 1. 预处理：读取参考面板的 SNP 列表，进行预过滤
# ==============================================================================
# 这一步极其关键！它能防止因 SNP ID 不匹配（如 chr:pos 格式）导致的 PLINK 报错
# 假设你的 .bim 文件路径是 /home/piperacillin/publicdata/rawinstruments/EUR.bim
bim_file <- "/home/piperacillin/publicdata/rawinstruments/EUR.bim"

message("正在读取参考面板 BIM 文件以校对 SNP ID...")
# 只读取第二列（SNP ID），速度很快
valid_snps <- fread(bim_file, select = 2, header = FALSE)$V2

# 仅保留存在于参考面板中的 SNP
dt_valid <- dt_lifted %>%
  filter(SNP %in% valid_snps)

message(paste0("过滤完成：原始 SNP ", nrow(dt_lifted), " 个，保留有效 SNP ", nrow(dt_valid), " 个"))
# 如果这里保留的太少，说明你的 rsID 和参考面板对不上，需要检查数据

# ==============================================================================
# 2. 定义一个安全的 Clump 函数
# ==============================================================================
safe_clump <- function(sub_df) {
  # 如果子数据框为空或行数太少，直接返回 NULL
  if(nrow(sub_df) == 0) return(NULL)
  
  tryCatch({
    res <- ld_clump(
      dplyr::tibble(rsid=sub_df$SNP, pval=sub_df$pval.exposure, id=sub_df$id.exposure),
      clump_r2 = 0.2,
      clump_kb = 10000,
      clump_p = 0.99, # 设为 0.99 保证不过滤 P 值，完全依赖输入的 sub_df
      plink_bin = get_plink_exe(),
      bfile = "/home/piperacillin/publicdata/rawinstruments/EUR"
    )
    return(res)
  }, error = function(e) {
    # 如果单个 ID 出错，只打印警告，不中断程序
    # message(paste("ID:", unique(sub_df$id.exposure)[1], "Clump 失败")) 
    return(NULL)
  })
}

# ==============================================================================
# 3. 并行计算执行 Clump
# ==============================================================================

# 将大表按 id.exposure 拆分成列表
# 这一步可能会消耗一些内存，如果几万个 ID 导致内存溢出，需要分批次处理
dt_list <- split(dt_valid, dt_valid$id.exposure)

message(paste("开始并行处理", length(dt_list), "个独立的 exposure ID..."))

# 设置并行核心数（根据你的服务器配置调整，建议 10-20 核）
n_cores <- 20 

# 使用 mclapply 进行并行运算 (Linux/Mac 专用，Windows需用 parLapply)
clumped_results_list <- mclapply(dt_list, safe_clump, mc.cores = n_cores)
# ==============================================================================
# 4. 合并结果并提取最终数据
# ==============================================================================
# 剔除 NULL 结果并合并
clumped_df <- bind_rows(clumped_results_list)

message(paste("Clump 完成。保留 SNP 数量:", nrow(clumped_df)))

# 将 Clump 后的结果回贴到原始数据（保留所有详细列）
dt_final <- dt_valid %>%
  dplyr::filter(SNP %in% clumped_df$rsid & id.exposure %in% clumped_df$id)

head(dt_final)

qs::qsave(dt_final, "CD4_Dynamic_Exposure_clump0.2.qs")


####dice####
rm(list = ls())
library(VariantAnnotation)
library(TwoSampleMR)
library(dplyr)
library(tibble)
library(tools) # 用于处理文件名

# 设置工作目录
setwd('/home/publicdata/datasets/eqtl/dice')

# =======================================================
# 1. 定义核心处理函数 (读取单个文件并处理)
# =======================================================
process_single_vcf <- function(vcf_filename) {
  
  message(paste0("\n>>> 正在处理文件: ", vcf_filename))
  
  # --- A. 读取 VCF ---
  # 使用 tryCatch 防止某个文件损坏导致整个循环中断
  vcf_obj <- tryCatch({
    VariantAnnotation::readVcf(vcf_filename)
  }, error = function(e) {
    message(paste("读取失败:", vcf_filename, "-", e$message))
    return(NULL)
  })
  
  if(is.null(vcf_obj)) return(NULL)
  
  # --- B. 筛选 P < 1e-5 ---
  raw_pvalues <- info(vcf_obj)$Pvalue
  pvalues_numeric <- as.numeric(as.character(raw_pvalues))
  
  # 处理 NA
  if(any(is.na(pvalues_numeric))) {
    pvalues_numeric[is.na(pvalues_numeric)] <- 1
  }
  
  threshold <- 1e-5
  keep_idx <- which(pvalues_numeric < threshold)
  
  n_total <- length(pvalues_numeric)
  n_keep <- length(keep_idx)
  message(paste0("    原始位点: ", n_total, " -> 筛选后: ", n_keep))
  
  if(n_keep == 0) {
    message("    警告: 该文件没有符合阈值的位点，跳过。")
    return(NULL)
  }
  
  # 创建子集
  vcf_subset <- vcf_obj[keep_idx, ]
  
  # --- C. 解析数据 (使用你之前的逻辑) ---
  rr <- rowRanges(vcf_subset)
  chrom <- as.character(seqnames(rr))
  pos <- start(rr)
  
  rsid <- rownames(vcf_subset)
  if(is.null(rsid) || length(rsid) == 0) { rsid <- paste0("chr", chrom, ":", pos) }
  
  ref_allele <- as.character(rr$REF)
  alt_allele <- as.character(unlist(lapply(rr$ALT, function(x) as.character(x[1]))))
  
  info_dat <- info(vcf_subset)
  beta <- as.numeric(info_dat$Beta)
  pval <- as.numeric(info_dat$Pvalue)
  
  # 基因符号
  gene_symbol <- if("GeneSymbol" %in% names(info_dat)) as.character(info_dat$GeneSymbol) else rep("Unknown", length(beta))
  
  # SE 反向计算
  safe_pval <- ifelse(pval < 1e-300, 1e-300, pval)
  z_score <- qnorm(safe_pval / 2, lower.tail = FALSE)
  se <- abs(beta) / z_score
  se[is.infinite(se) | is.nan(se)] <- 0
  
  # --- D. 构建特定的 Phenotype 名称 ---
  # 需求：'gene_symbol'_'文件名去除.vcf'
  
  # 1. 获取纯文件名 (去除路径和扩展名)
  # 注意：file_path_sans_ext 只去除最后一个后缀，如果是 .vcf.gz 需处理两次或用 gsub
  file_tag <- basename(vcf_filename)
  file_tag <- gsub("\\.vcf$", "", file_tag, ignore.case = TRUE) # 去除 .vcf
  file_tag <- gsub("\\.gz$", "", file_tag, ignore.case = TRUE)  # 如果有 .gz 也去除
  file_tag <- gsub("\\.vcf$", "", file_tag, ignore.case = TRUE) # 再次去除 (针对 .vcf.gz)
  
  # 2. 组合名称
  # 处理 GeneSymbol 为空的情况
  gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <- "UnknownGene"
  new_phenotype <- paste(gene_symbol, file_tag, sep = "_")
  
  # --- E. 返回 Tibble ---
  res_df <- tibble(
    SNP = rsid,
    chr = chrom,
    pos = pos,
    beta = beta,
    se = se,
    pval = pval,
    eaf = NA, 
    effect_allele = alt_allele,
    other_allele = ref_allele,
    Phenotype = new_phenotype # 这里是我们构造的新列
  )
  
  # 显式进行垃圾回收，这在处理循环大文件时非常重要
  rm(vcf_obj, vcf_subset)
  gc() 
  
  return(res_df)
}

# =======================================================
# 2. 批量执行
# =======================================================

# 获取目录下所有的 VCF 文件
vcf_files <- list.files(pattern = "\\.vcf", full.names = TRUE)

message(paste("共检测到", length(vcf_files), "个 VCF 文件。开始处理..."))

# 使用 lapply 循环处理所有文件，结果保存为一个 List
# 如果文件非常多，这里可能会跑一段时间
df_list <- lapply(vcf_files, process_single_vcf)

# 移除 List 中的 NULL 元素 (处理失败或无显著位点的文件)
df_list <- df_list[!sapply(df_list, is.null)]

message("所有文件处理完毕，正在合并数据...")

# =======================================================
# 3. 合并并格式化为 TwoSampleMR 对象
# =======================================================

if(length(df_list) > 0) {
  
  # 合并将所有小表变成一个总表
  master_df <- bind_rows(df_list)
  
  message(paste("合并完成！总行数:", nrow(master_df)))
  print(head(master_df))
  
  message("正在转换为 TwoSampleMR 格式...")
  
  final_exposure_dat <- format_data(
    master_df,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    chr_col = "chr",
    pos_col = "pos",
    phenotype_col = "Phenotype" # 关键：使用自定义合成的列
  )
  
  message("转换成功！final_exposure_dat 已准备好。")
  print(head(final_exposure_dat))
  
  # 查看一下有多少个独特的表型 (Gene_CellType)
  message("包含的独特 Phenotype 数量: ", length(unique(final_exposure_dat$phenotype)))
  
} else {
  stop("这就尴尬了：没有任何数据符合筛选条件 (P < 1e-5)。请检查文件或放宽阈值。")
}


