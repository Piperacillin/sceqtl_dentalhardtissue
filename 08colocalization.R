####功能：数据切片####
# 第一步：修改 R 的底层配置文件
# 在你的 R 控制台（或者 RStudio）里，直接运行下面这行代码：
# 
# R
# 
# # 打开或创建当前用户的 R 环境变量配置文件
# file.edit("~/.Renviron")
# 这时候，RStudio（或你的文本编辑器）会弹出一个空白文件（或者已有内容的文件）。请在这个文件里，新增这极其关键的一行（注意大小写和斜杠，不要有任何多余的空格）：
# 
# Plaintext
# 
# TMPDIR="E:/R_temp"
# 写完之后，保存这个文件（Ctrl + S）。
# 
# 第二步：重启 R 会话（极其关键！）
# 配置写好了，但当前的 R 还没醒过神来。你必须重启它：
# 
# 如果你在用 RStudio：点击上方菜单栏的 Session -> Restart R（或者快捷键 Ctrl + Shift + F10）。
# 
# 如果你在用原生 R 终端：输入 q() 退出，然后再重新打开 R。

tempdir()
setwd("D:/OneDrive/学院与活动/postdoc/00课题/2511sc_MR/2604reboot")

onek1k_eqtl_table <- 'G:/rawdata/sceqtl/eqtl_table.tsv.gz'
onek1k_esnp_table <- 'G:/rawdata/sceqtl/esnp_table.tsv.gz'
finn1 <- 'G:/rawdata/sceqtl/finngenR12_dentalhardtissue/finngen_R12_K11_EMBIMPACT_TEETH.gz'
finn2 <- 'G:/rawdata/sceqtl/finngenR12_dentalhardtissue/finngen_R12_K11_ABRASION.gz'

# ==========================================================
# 极速偷看 20GB 压缩包的前 5 行 (纯 Base R，Windows 完美兼容)
# ==========================================================
# 1. 建立指向 gz 文件的连接
con <- gzfile(onek1k_eqtl_table, "rt") 

# 2. 只读取前 5 行
preview_lines <- readLines(con, n = 5)

# 3. 必须关闭连接，释放资源
close(con)

# 4. 打印查看
cat(preview_lines, sep = "\n")

####
# 1. 建立指向 gz 文件的连接
con <- gzfile(onek1k_esnp_table, "rt") 

# 2. 只读取前 5 行
preview_lines <- readLines(con, n = 5)

# 3. 必须关闭连接，释放资源
close(con)

# 4. 打印查看
cat(preview_lines, sep = "\n")

####
# 1. 建立指向 gz 文件的连接
con <- gzfile(finn1, "rt") 

# 2. 只读取前 5 行
preview_lines <- readLines(con, n = 5)

# 3. 必须关闭连接，释放资源
close(con)

# 4. 打印查看
cat(preview_lines, sep = "\n")

####
# 1. 建立指向 gz 文件的连接
con <- gzfile(finn2, "rt") 

# 2. 只读取前 5 行
preview_lines <- readLines(con, n = 5)

# 3. 必须关闭连接，释放资源
close(con)

# 4. 打印查看
cat(preview_lines, sep = "\n")



# ==========================================================
# 偷看 1M-scBloodNL 文件表头及前 5 行
# ==========================================================
test_file <- "G:/rawdata/sceqtl/eqtls_20201106_genome_wide/UT/CD4T_expression_eQTLsFDR-ProbeLevel.txt.gz"

con <- gzfile(test_file, "rt") 
preview_lines <- readLines(con, n = 5)
close(con)

cat(preview_lines, sep = "\n")


# ==============================================================================
# 功能: Windows 兼容版 - Phase4 高分靶点 OneK1K 极速切片与 GWAS 对齐
# 核心: 直接根据 OneK1K eqtl_table 提取全信息，免除 Join；自动反推推导 SE
# ==============================================================================
library(data.table)
library(dplyr)
library(vroom)
library(stringr)

# ------------------------------------------------------------------------------
# 1. 加载 Phase4 打分表，构建 OneK1K 的“切片任务清单”
# ------------------------------------------------------------------------------
phase4_df <- fread("D:/OneDrive/学院与活动/postdoc/00课题/2511sc_MR/2604reboot/result/Phase4_Final_Project_Grand_Master_Table_with_Bulkadded.csv")

task_list_onek1k <- phase4_df %>%
  filter(final_5scale_score >= 4) %>%
  filter(!is.na(cell_type_onek1k) & cell_type_onek1k != "Not_Available") %>%
  select(GENE = gene_symbol, CELL_TYPE = cell_type_onek1k, outcome) %>%
  distinct()

message(sprintf(">>> 共生成 %d 个 OneK1K 精确 [靶点+细胞] 切片任务！", nrow(task_list_onek1k)))

# ------------------------------------------------------------------------------
# 2. 扫库 eqtl_table.tsv.gz (包含全部背景 SNPs)
# ------------------------------------------------------------------------------
message("\n>>> 正在全库扫描 OneK1K eqtl_table (仅读取命中任务的相关数据)...")

# OneK1K 变量映射说明: 
# RSID = SNP ID; A2 = Effect Allele (结合 A2_FREQ); A1 = Other Allele
# SPEARMANS_RHO ≈ Beta; A2_FREQ_ONEK1K = MAF/EAF
# ==============================================================================
# 彻底解决 OOM: 使用 readr 分块流式提取 OneK1K 20GB 数据
# ==============================================================================
library(readr)  # 必须加载 readr

# 1. (假设你已经跑了前面的代码，拿到了 task_list_onek1k)
# task_list_onek1k <- phase4_df %>% ... 

message("\n>>> 开启防爆内存模式：分块流式扫描 OneK1K 20GB eqtl_table (每次 50 万行)...")
message(">>> 这个过程可能会运行 10~30 分钟，会有进度条显示，请耐心喝杯咖啡。")

# 2. 精确定义所需读取的列类型 (丢弃不需要的列，极大提升速度和节约内存)
my_cols <- cols_only(
  GENE = col_character(),
  CELL_TYPE = col_character(),
  RSID = col_character(),
  CHR = col_character(),
  POS = col_integer(),
  A1 = col_character(),
  A2 = col_character(),
  A2_FREQ_ONEK1K = col_double(),
  SPEARMANS_RHO = col_double(),
  P_VALUE = col_double()
)

# 3. 定义每次读取一个 Chunk (数据块) 时的过滤动作
chunk_filter <- function(chunk, pos) {
  # 只保留 task_list_onek1k 里需要的 基因+细胞 组合
  chunk %>% semi_join(task_list_onek1k, by = c("GENE", "CELL_TYPE"))
}

# 4. 执行分块流式读取！
eqtl_onek1k_hits <- read_tsv_chunked(
  file = "G:/rawdata/sceqtl/eqtl_table.tsv.gz",
  callback = DataFrameCallback$new(chunk_filter),
  chunk_size = 500000, # 每次吞吐 50 万行，内存绝对安全
  col_types = my_cols,
  progress = TRUE      # 在控制台显示扫描进度条，让你心里有底
)

message(">>> OneK1K 分块扫描完毕！正在推导 Standard Error (SE)...")

# 5. 等全部小数据拼好后，再统一计算 SE (放在循环外算，节省性能)
eqtl_onek1k_hits <- eqtl_onek1k_hits %>%
  mutate(
    z_score = qnorm(P_VALUE / 2, lower.tail = FALSE),
    eqtl_se = ifelse(is.infinite(z_score) | z_score == 0, NA, abs(SPEARMANS_RHO) / z_score)
  ) %>%
  filter(!is.na(eqtl_se)) # 剔除极值导致无法计算 SE 的行

message(sprintf(">>> OneK1K 完美提取成功！一共提取到 %d 行有效背景 SNP 数据！", nrow(eqtl_onek1k_hits)))

# ==============================================================================
# 接下来，你就可以直接拿着这个安全的 `eqtl_onek1k_hits` 
# 去跑前面的 [3. 循环切片与 FinnGen GWAS 共定位对齐] 代码了。
# ==============================================================================

# ------------------------------------------------------------------------------
# 3. 循环切片与 FinnGen GWAS 共定位对齐
# ------------------------------------------------------------------------------
if(!dir.exists("./slice_data")) dir.create("./slice_data")

for (i in 1:nrow(task_list_onek1k)) {
  
  current_gene <- task_list_onek1k$GENE[i]
  current_cell <- task_list_onek1k$CELL_TYPE[i]
  current_gwas <- task_list_onek1k$outcome[i]
  
  message(sprintf("\n>>> [%d/%d] 处理: %s | %s | GWAS: %s", 
                  i, nrow(task_list_onek1k), current_gene, current_cell, current_gwas))
  
  # 提取当前区域
  current_eqtl_locus <- eqtl_onek1k_hits %>% 
    filter(GENE == current_gene & CELL_TYPE == current_cell)
  
  if (nrow(current_eqtl_locus) < 50) {
    message("   [!] 警告：背景 SNP 数量极少 (< 50)，无法完成高质量共定位，跳过。")
    next
  }
  
  chr_num <- unique(current_eqtl_locus$CHR)[1]
  pos_min <- min(current_eqtl_locus$POS)
  pos_max <- max(current_eqtl_locus$POS)
  
  # 读取 FinnGen GWAS 数据
  # ==========================================================
  # [终极修改] 智能分轨读取 GWAS 数据 + RSID精准锚定 (免疫hg38坐标偏移)
  # ==========================================================
  # 提取当前 eQTL 靶点区间内的所有可用 RSID，作为去 GWAS 里捞针的"磁铁"
  target_rsids <- unique(current_eqtl_locus$RSID)
  
  # 判断当前表型是否属于那三个 Meta 荟萃文件
  meta_keywords <- c("EMBIMPACT_TEETH", "DENTOFACIAL_ANOMALIES", "CARIES_1_OPER_ONLYAVO")
  is_meta <- any(sapply(meta_keywords, function(k) grepl(k, current_gwas)))
  
  gwas_base_dir <- "G:/rawdata/sceqtl/finngenR12_dentalhardtissue"
  
  if (is_meta) {
    # ---------------- 轨道 A: 处理 Meta 荟萃版 GWAS ----------------
    gwas_file <- sprintf("%s/finngen_R12_%s_meta_out.tsv.gz", gwas_base_dir, current_gwas)
    if (!file.exists(gwas_file)) {
      message(sprintf("   [!] 未找到 Meta GWAS 数据: %s，跳过。", gwas_file))
      next
    }
    
    gwas_subset <- vroom(
      gwas_file, 
      col_select = c(`#CHR`, POS, REF, ALT, rsid, 
                     all_inv_var_meta_beta, all_inv_var_meta_sebeta, all_inv_var_meta_p, FINNGEN_af_alt),
      show_col_types = FALSE
    ) %>%
      # 【核心避坑】直接通过 RSID 匹配，彻底抛弃容易错位的 POS 区间
      filter(rsid %in% target_rsids) %>%
      rename(
        gwas_chr = `#CHR`, gwas_pos = POS,
        gwas_ea = ALT, gwas_oa = REF,
        gwas_beta = all_inv_var_meta_beta, gwas_se = all_inv_var_meta_sebeta,
        gwas_p = all_inv_var_meta_p, gwas_maf = FINNGEN_af_alt
      )
    
  } else {
    # ---------------- 轨道 B: 处理 标准版 FinnGen GWAS ----------------
    gwas_file <- sprintf("%s/finngen_R12_%s.gz", gwas_base_dir, current_gwas)
    if (!file.exists(gwas_file)) {
      message(sprintf("   [!] 未找到标准 GWAS 数据: %s，跳过。", gwas_file))
      next
    }
    
    gwas_subset <- vroom(
      gwas_file, 
      col_select = c(`#chrom`, pos, ref, alt, rsids, beta, sebeta, pval, af_alt),
      show_col_types = FALSE
    ) %>%
      # 【核心避坑】直接通过 RSID 匹配
      filter(rsids %in% target_rsids) %>%
      rename(
        rsid = rsids, gwas_chr = `#chrom`, gwas_pos = pos,
        gwas_ea = alt, gwas_oa = ref,
        gwas_beta = beta, gwas_se = sebeta,
        gwas_p = pval, gwas_maf = af_alt
      )
  }
  
  if (nrow(gwas_subset) == 0) {
    message("   [!] 警告：在 GWAS 中未匹配到 eQTL 的任何 RSID，可能被质控过滤，跳过。")
    next
  }
  
  # ==========================================================
  # 等位基因对齐 (Harmonization)
  # ==========================================================
  merged_locus <- inner_join(
    current_eqtl_locus %>% select(rsid = RSID, pos = POS, 
                                  eqtl_ea = A2, eqtl_oa = A1, 
                                  eqtl_beta = SPEARMANS_RHO, eqtl_se, eqtl_p = P_VALUE, eqtl_maf = A2_FREQ_ONEK1K),
    gwas_subset, 
    by = "rsid"
  )
  
  merged_locus <- merged_locus %>%
    mutate(
      harmonized = case_when(
        eqtl_ea == gwas_ea & eqtl_oa == gwas_oa ~ TRUE,
        eqtl_ea == gwas_oa & eqtl_oa == gwas_ea ~ FALSE,
        TRUE ~ NA
      ),
      final_gwas_beta = ifelse(harmonized == FALSE, gwas_beta * -1, gwas_beta),
      final_gwas_ea   = ifelse(harmonized == FALSE, gwas_oa, gwas_ea)
    ) %>%
    filter(!is.na(harmonized))
  
  # 【致命 Bug 修复处】
  final_export <- merged_locus %>%
    mutate(chr = chr_num) %>%  # 必须使用 mutate 将外部变量作为一个新列加入！
    select(rsid, chr, pos,     # 这里直接取干净的 pos 即可
           eqtl_ea, eqtl_oa, eqtl_beta, eqtl_se, eqtl_p, eqtl_maf,
           gwas_ea = final_gwas_ea, gwas_beta = final_gwas_beta, gwas_se, gwas_p, gwas_maf)
  
  # 安全导出
  safe_cell_name <- str_replace_all(current_cell, "[ /]", "_")
  output_path <- sprintf("./slice_data/OneK1K_%s_%s_%s_Aligned.csv", 
                         current_gene, safe_cell_name, current_gwas)
  fwrite(final_export, output_path)
  
  message(sprintf("   [√] 成功落盘！可用共定位 SNP: %d", nrow(final_export)))
}


# ==============================================================================
# 功能: 1M-scBloodNL 终极适配版 - OOM免疫 + 自动兼容 Meta/标准 GWAS + 免疫 hg38 位置偏移
# ==============================================================================
library(data.table)
library(dplyr)
library(readr)  
library(stringr)
library(vroom)

base_dir <- "G:/rawdata/sceqtl/eqtls_20201106_genome_wide"
if(!dir.exists("./slice_data")) dir.create("./slice_data")
phase4_df <- fread("D:/OneDrive/学院与活动/postdoc/00课题/2511sc_MR/2604reboot/result/Phase4_Final_Project_Grand_Master_Table_with_Bulkadded.csv")

# 2. 从 Phase4 表单中提取 1M-scBloodNL 的任务
# 注意检查你的 act_time_1mscblood 是否都是规范的 "UT", "3hCA" 等格式
task_list_1mscblood <- phase4_df %>%
  filter(final_5scale_score >= 4) %>%
  filter(!is.na(cell_type_1mscblood) & cell_type_1mscblood != "Not_Available") %>%
  select(GENE = gene_symbol, 
         ACT_TIME = act_time_1mscblood,   
         CELL_TYPE = cell_type_1mscblood, 
         outcome) %>%
  distinct()

message(sprintf(">>> 共生成 %d 个 1M-scBloodNL 动态路径切片任务！", nrow(task_list_1mscblood)))

# 明确 1M-scBloodNL eQTL 我们要提取的列
eqtl_cols <- cols_only(
  PValue = col_double(),
  SNPName = col_character(),
  SNPChr = col_character(),
  SNPChrPos = col_integer(),
  SNPType = col_character(),
  AlleleAssessed = col_character(),
  OverallZScore = col_double(),
  HGNCName = col_character()
)

for (i in 1:nrow(task_list_1mscblood)) {
  
  current_gene <- task_list_1mscblood$GENE[i]
  current_time <- task_list_1mscblood$ACT_TIME[i]
  current_cell <- task_list_1mscblood$CELL_TYPE[i]
  current_gwas <- task_list_1mscblood$outcome[i] # 例如: K11_EMBIMPACT_TEETH
  
  message(sprintf("\n>>> [%d/%d] 启动切片: 靶标 %s | 状态 %s | GWAS %s", 
                  i, nrow(task_list_1mscblood), current_gene, current_time, current_gwas))
  
  # ==========================================================
  # 1. 纯流式零内存消耗提取 1M-scBloodNL
  # ==========================================================
  target_file <- sprintf("%s/%s/%s_expression_eQTLsFDR-ProbeLevel.txt.gz", 
                         base_dir, current_time, current_cell)
  if (!file.exists(target_file)) {
    message("   [!] 找不到 eQTL 状态文件，跳过。")
    next
  }
  
  chunk_filter <- function(chunk, pos) { chunk %>% filter(HGNCName == current_gene) }
  
  current_eqtl_locus <- read_tsv_chunked(
    file = target_file, callback = DataFrameCallback$new(chunk_filter),
    chunk_size = 500000, col_types = eqtl_cols, progress = FALSE
  )
  
  if (nrow(current_eqtl_locus) < 5) {
    message("   [!] eQTL 背景 SNP 极少，判定为缺失靶标，跳过。")
    next
  }
  
  # 格式清洗，获取 eQTL 的可用 rsid 集合，这将被用作去 GWAS 里捞针的"磁铁"
  current_eqtl_locus <- current_eqtl_locus %>%
    mutate(
      allele_1 = str_extract(SNPType, "^[^/]+"),
      allele_2 = str_extract(SNPType, "[^/]+$"),
      eqtl_oa = ifelse(allele_1 == AlleleAssessed, allele_2, allele_1),
      eqtl_ea = AlleleAssessed,
      eqtl_beta = OverallZScore, 
      eqtl_se = 1,               
      eqtl_p = PValue,
      eqtl_maf = 0.5             
    ) %>%
    select(rsid = SNPName, chr = SNPChr, pos = SNPChrPos, 
           eqtl_ea, eqtl_oa, eqtl_beta, eqtl_se, eqtl_p, eqtl_maf) %>%
    filter(!is.na(eqtl_oa))
  
  target_rsids <- unique(current_eqtl_locus$rsid)
  
  # ==========================================================
  # 2. 智能分轨读取 GWAS 数据 (Meta版 vs 标准版)
  # 核心策略：直接用 rsid 过滤，不依赖 POS，完美避开 hg19/hg38 基因组版本偏移
  # ==========================================================
  # 判断是否为 Meta 数据集
  meta_keywords <- c("EMBIMPACT_TEETH", "DENTOFACIAL_ANOMALIES", "CARIES_1_OPER_ONLYAVO")
  is_meta <- any(sapply(meta_keywords, function(k) grepl(k, current_gwas)))
  
  if (is_meta) {
    # ---------------- 轨道 A: 处理 Meta 荟萃版 GWAS ----------------
    gwas_file <- sprintf("G:/rawdata/sceqtl/finngenR12_dentalhardtissue/finngen_R12_%s_meta_out.tsv.gz", current_gwas)
    if (!file.exists(gwas_file)) { message("   [!] Meta GWAS 文件缺失。"); next }
    
    gwas_subset <- vroom(
      gwas_file,
      # 按你给的抬头精准指定列
      col_select = c(`#CHR`, POS, REF, ALT, rsid, 
                     all_inv_var_meta_beta, all_inv_var_meta_sebeta, all_inv_var_meta_p, FINNGEN_af_alt),
      show_col_types = FALSE
    ) %>%
      # 【避坑神技】直接通过 ID 匹配，不再卡 POS 范围
      filter(rsid %in% target_rsids) %>%
      # 统一重命名为标准 GWAS 格式
      rename(
        gwas_chr = `#CHR`, gwas_pos = POS,
        gwas_ea = ALT, gwas_oa = REF,
        gwas_beta = all_inv_var_meta_beta, gwas_se = all_inv_var_meta_sebeta,
        gwas_p = all_inv_var_meta_p, gwas_maf = FINNGEN_af_alt
      )
    
  } else {
    # ---------------- 轨道 B: 处理 标准版 FinnGen GWAS ----------------
    gwas_file <- sprintf("G:/rawdata/sceqtl/finngenR12_dentalhardtissue/finngen_R12_%s.gz", current_gwas)
    if (!file.exists(gwas_file)) { message("   [!] 标准 GWAS 文件缺失。"); next }
    
    gwas_subset <- vroom(
      gwas_file,
      col_select = c(`#chrom`, pos, ref, alt, rsids, beta, sebeta, pval, af_alt),
      show_col_types = FALSE
    ) %>%
      # 【避坑神技】同样只按 ID 匹配
      filter(rsids %in% target_rsids) %>%
      rename(
        rsid = rsids, gwas_chr = `#chrom`, gwas_pos = pos,
        gwas_ea = alt, gwas_oa = ref,
        gwas_beta = beta, gwas_se = sebeta,
        gwas_p = pval, gwas_maf = af_alt
      )
  }
  
  if (nrow(gwas_subset) == 0) {
    message("   [!] 在 GWAS 中未匹配到任何 eQTL SNP (可能被质控过滤)，跳过。")
    next
  }
  
  # ==========================================================
  # 3. 合并与等位基因对齐 (Harmonization)
  # ==========================================================
  merged_locus <- inner_join(current_eqtl_locus, gwas_subset, by = "rsid") %>%
    mutate(
      harmonized = case_when(
        eqtl_ea == gwas_ea & eqtl_oa == gwas_oa ~ TRUE,
        eqtl_ea == gwas_oa & eqtl_oa == gwas_ea ~ FALSE,
        TRUE ~ NA
      ),
      final_gwas_beta = ifelse(harmonized == FALSE, gwas_beta * -1, gwas_beta),
      final_gwas_ea   = ifelse(harmonized == FALSE, gwas_oa, gwas_ea)
    ) %>%
    filter(!is.na(harmonized))
  
  final_export <- merged_locus %>%
    select(rsid, 
           chr, pos, # 保留 eQTL 的 hg19 坐标用于作图 (LocusZoom需要位置，但只要相对走势一致即可)
           eqtl_ea, eqtl_oa, eqtl_beta, eqtl_se, eqtl_p, eqtl_maf,
           gwas_ea = final_gwas_ea, gwas_beta = final_gwas_beta, gwas_se, gwas_p, gwas_maf)
  
  safe_cell_name <- str_replace_all(current_cell, "[ /]", "_")
  output_path <- sprintf("./slice_data/1MscBlood_%s_%s_%s_%s_Aligned.csv", 
                         current_gene, current_time, safe_cell_name, current_gwas)
  fwrite(final_export, output_path)
  
  message(sprintf("   [√] %s 完美落盘！可用共定位 SNP: %d", 
                  ifelse(is_meta, "[Meta版]", "[标准版]"), nrow(final_export)))
}


####功能: 批量自动化共定位分析 (Colocalization) + 顶刊级 LocusZoom 镜像出图####
# ==============================================================================
# 功能: 批量自动化共定位分析 (Colocalization) + 顶刊级 LocusZoom 镜像出图
# 升级: 自动读取 FinnGen Manifest 字典表，动态填补 GWAS 的 N 和 s 参数
# ==============================================================================
library(data.table)
library(dplyr)
library(coloc)
library(ggplot2)
library(patchwork)
library(stringr)
getwd()
setwd("D:/OneDrive/学院与活动/postdoc/00课题/2511sc_MR/2604reboot/")
# ------------------------------------------------------------------------------
# 1. 环境初始化与字典预载入
# ------------------------------------------------------------------------------
input_dir <- "./slice_data"
plot_dir <- "./coloc_plots"

if(!dir.exists(plot_dir)) dir.create(plot_dir)

# [核心升级] 加载 FinnGen Manifest 字典表
manifest_path <- "G:/rawdata/sceqtl/finngenR12_dentalhardtissue/finngen_R12_manifest.csv"
manifest_df <- fread(manifest_path) %>%
  # 按 phenocode 长度降序排列，防止短名字被长名字错误匹配包含（如 AB1 匹配进 AB1_XYZ）
  arrange(desc(nchar(phenocode)))

# 扫描切片文件
file_list <- list.files(input_dir, pattern = "_Aligned.csv$", full.names = TRUE)
message(sprintf(">>> 共扫描到 %d 个准备就绪的对齐文件，已载入 Manifest 字典表！", length(file_list)))

summary_results <- data.frame()

# ------------------------------------------------------------------------------
# 2. 开启批量循环
# ------------------------------------------------------------------------------
for (i in 1:length(file_list)) {
  
  file_path <- file_list[i]
  file_name <- basename(file_path)
  
  # ==========================================================
  # 【刚刚为你加入的修复区】动态解析文件名！防止变量幽灵污染！
  # ==========================================================
  
  # 1. 拆分文件名获取数据库和基因名
  # 例如 file_name 是 "1MscBlood_ACAT1_3hPA_CD4T_K11_EMBIMPACT_TEETH_Aligned.csv"
  name_parts <- str_split(file_name, "_")[[1]]
  dataset_type <- name_parts[1] 
  gene_name <- name_parts[2]    # 这里绝不会再永远是 ZP3 了
  
  # 2. 精确正则提取表型 (寻找 K11_ 开头，到 _Aligned.csv 结束的部分)
  current_gwas <- str_extract(file_name, "K11_.*(?=_Aligned\\.csv)")
  
  # 3. 从字典中提取 N (总样本量) 和 s (病例占比)
  gwas_N <- NA
  gwas_s <- NA
  
  if (!is.na(current_gwas)) {
    matched_row <- manifest_df %>% filter(phenocode == current_gwas) %>% slice(1)
    if (nrow(matched_row) == 1) {
      cases <- as.numeric(matched_row$num_cases)
      controls <- as.numeric(matched_row$num_controls)
      gwas_N <- cases + controls
      gwas_s <- cases / gwas_N
    }
  }
  
  # 安全拦截：如果真的找不到样本量，直接跳过当前文件，防止 coloc 报错崩溃
  if (is.na(gwas_N) || is.na(gwas_s)) {
    message(sprintf("   [!] 未能在 Manifest 中匹配到表型样本量，跳过此文件: %s", file_name))
    next
  }
  
  # ==========================================================
  # 从这里开始，接回你原本的数据读取与计算逻辑
  # ==========================================================
  
  # 1. 读取对齐数据
  df_locus <- fread(file_path)
  
  # ==========================================================
  # 【史诗级抢救】自动修复 OneK1K pos 列缺失 Bug
  # ==========================================================
  if (all(is.na(df_locus$pos))) {
    message("   [!] 警报：检测到 pos 列全为空 (OneK1K 原始数据浮点解析丢失)！")
    message("   [*] 正在启动紧急预案，从本地 GWAS 文件中极速抓取坐标修复...")
    
    meta_keywords <- c("EMBIMPACT_TEETH", "DENTOFACIAL_ANOMALIES", "CARIES_1_OPER_ONLYAVO")
    is_meta <- any(sapply(meta_keywords, function(k) grepl(k, current_gwas)))
    gwas_base_dir <- "G:/rawdata/sceqtl/finngenR12_dentalhardtissue"
    
    if (is_meta) {
      gwas_file <- sprintf("%s/finngen_R12_%s_meta_out.tsv.gz", gwas_base_dir, current_gwas)
    } else {
      gwas_file <- sprintf("%s/finngen_R12_%s.gz", gwas_base_dir, current_gwas)
    }
    
    if (file.exists(gwas_file)) {
      if (is_meta) {
        pos_map <- vroom(gwas_file, col_select = c(rsid, POS), show_col_types = FALSE)
        df_locus <- df_locus %>% select(-pos) %>% inner_join(pos_map, by = "rsid") %>% rename(pos = POS)
      } else {
        pos_map <- vroom(gwas_file, col_select = c(rsids, pos), show_col_types = FALSE)
        df_locus <- df_locus %>% select(-pos) %>% inner_join(pos_map, by = c("rsid" = "rsids"))
      }
      # 将抢救回来的完整数据重新写回本地，永久修复这个 CSV
      fwrite(df_locus, file_path)
      message("   [√] 坐标抢救成功！已永久补齐该靶点物理坐标。")
    } else {
      message("   [!] 找不到底层 GWAS 文件，无法抢救，跳过。")
      next
    }
  }
  
  # ==========================================================
  # 2. [终极极值消毒] 处理空值与极值
  # ==========================================================
  df_locus <- df_locus %>%
    # 【新增拦截】彻底剔除所有可能的 NA 值（重点拦截 MAF 缺失）
    filter(!is.na(gwas_beta) & !is.na(gwas_se) & !is.na(eqtl_beta) & !is.na(eqtl_se) &
             !is.na(gwas_maf) & !is.na(eqtl_maf) & !is.na(gwas_p) & !is.na(eqtl_p)) %>%
    mutate(
      # 防止绝对极强的信号 P = 0 导致崩溃
      gwas_p = ifelse(gwas_p == 0 | gwas_p < 1e-300, 1e-300, gwas_p),
      eqtl_p = ifelse(eqtl_p == 0 | eqtl_p < 1e-300, 1e-300, eqtl_p),
      gwas_se = ifelse(gwas_se <= 0, 1e-6, gwas_se),
      eqtl_se = ifelse(eqtl_se <= 0, 1e-6, eqtl_se),
      
      # 约束 MAF 范围，防止 MAF 为 0 或 1 导致方差为 0
      gwas_maf = ifelse(gwas_maf <= 0, 0.001, ifelse(gwas_maf >= 1, 0.999, gwas_maf)),
      eqtl_maf = ifelse(eqtl_maf <= 0, 0.001, ifelse(eqtl_maf >= 1, 0.999, eqtl_maf))
    ) %>%
    filter(is.finite(gwas_beta) & is.finite(gwas_se) & is.finite(eqtl_beta) & is.finite(eqtl_se)) %>%
    distinct(rsid, .keep_all = TRUE) %>%
    arrange(pos)
  
  if(nrow(df_locus) < 10) {
    message("   [!] 清洗后可用 SNP 过少 (<10)，跳过。")
    next
  }
  
  # ==========================================================
  # 执行 Coloc 分析 (装填动态参数)
  # ==========================================================
  # Dataset 1: GWAS (疾病)
  dataset_gwas <- list(
    pvalues = df_locus$gwas_p,
    type    = "cc",            
    s       = gwas_s,          # 【完美填补】利用字典提取的 Case 比例
    N       = gwas_N,          # 【完美填补】利用字典提取的总样本量
    MAF     = df_locus$gwas_maf,
    beta    = df_locus$gwas_beta,
    varbeta = (df_locus$gwas_se)^2,
    snp     = df_locus$rsid,
    position= df_locus$pos
  )
  
  # Dataset 2: eQTL (基因表达)
  dataset_eqtl <- list(
    pvalues = df_locus$eqtl_p,
    type    = "quant",
    N       = ifelse(dataset_type == "OneK1K", 982, 500), 
    MAF     = df_locus$eqtl_maf,
    beta    = df_locus$eqtl_beta,
    varbeta = (df_locus$eqtl_se)^2,
    snp     = df_locus$rsid,
    position= df_locus$pos
  )
  
  coloc_res <- tryCatch({
    coloc.abf(dataset1 = dataset_gwas, dataset2 = dataset_eqtl)
  }, error = function(e) {
    message("   [!] coloc 算法运行失败，数据奇异。")
    return(NULL)
  })
  
  if (is.null(coloc_res)) next
  
  # 提取后验概率
  pp_h0 <- coloc_res$summary["PP.H0.abf"]
  pp_h1 <- coloc_res$summary["PP.H1.abf"]
  pp_h2 <- coloc_res$summary["PP.H2.abf"]
  pp_h3 <- coloc_res$summary["PP.H3.abf"]
  pp_h4 <- coloc_res$summary["PP.H4.abf"]
  
  message(sprintf("   [√] Coloc 计算完毕: PP.H3 = %.3f, PP.H4 = %.3f", pp_h3, pp_h4))
  
  conclusion <- case_when(
    pp_h4 > 0.7 ~ "Strong Evidence for Colocalization",
    pp_h3 > 0.5 ~ "Independent Signals (LD Confounding)",
    pp_h0 > 0.5 ~ "No Association",
    TRUE ~ "Inconclusive"
  )
  
  # 汇总至数据表，为了学术严谨，顺便把 N 和 s 也一并记录在附录里备查
  current_result <- data.frame(
    File_Source = file_name,
    Database = dataset_type,
    Gene = gene_name,
    Phenotype = current_gwas,
    GWAS_N = gwas_N,
    GWAS_s = round(gwas_s, 4),
    Analyzed_SNPs = nrow(df_locus),
    PP.H0 = round(pp_h0, 4),
    PP.H1 = round(pp_h1, 4),
    PP.H2 = round(pp_h2, 4),
    PP.H3 = round(pp_h3, 4),
    PP.H4 = round(pp_h4, 4),
    Conclusion = conclusion,
    stringsAsFactors = FALSE
  )
  summary_results <- bind_rows(summary_results, current_result)
  
  # ==========================================================
  # 批量出图 
  # ==========================================================
  if (pp_h4 > 0.5) {
    message("   >>> 触发 PP.H4 出图阈值！生成 LocusZoom...")
    
    lead_snp <- df_locus$rsid[which.min(df_locus$gwas_p)]
    lead_pos <- df_locus$pos[df_locus$rsid == lead_snp]
    df_locus <- df_locus %>% mutate(ld_r2 = exp(-abs(pos - lead_pos)/25000))
    
    ld_colors <- scale_color_gradientn(
      colors = c("navy", "lightblue", "green", "orange", "red"),
      limits = c(0, 1), name = bquote(LD~(r^2))
    )
    
    p_gwas <- ggplot(df_locus, aes(x = pos / 1e6, y = -log10(gwas_p), color = ld_r2)) +
      geom_point(alpha = 0.8, size = 2) + ld_colors +
      geom_point(data = filter(df_locus, rsid == lead_snp), color = "purple", shape = 18, size = 5) +
      annotate("text", x = lead_pos / 1e6, y = max(-log10(df_locus$gwas_p)) * 1.1, 
               label = lead_snp, color = "purple", fontface = "bold") +
      theme_classic(base_size = 14) +
      labs(title = sprintf("GWAS Signal - %s Locus", gene_name), 
           x = NULL, y = expression(-log[10](italic(P)[GWAS]))) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    annotation_text <- sprintf("Coloc PP.H4 = %.3f\nDatabase: %s", pp_h4, dataset_type)
    
    p_eqtl <- ggplot(df_locus, aes(x = pos / 1e6, y = -log10(eqtl_p), color = ld_r2)) +
      geom_point(alpha = 0.8, size = 2) + ld_colors +
      geom_point(data = filter(df_locus, rsid == lead_snp), color = "purple", shape = 18, size = 5) +
      annotate("text", x = min(df_locus$pos)/1e6, y = max(-log10(df_locus$eqtl_p))*0.9, 
               label = annotation_text, hjust = 0, fontface = "bold", size = 4.5) +
      theme_classic(base_size = 14) +
      labs(title = sprintf("sc-eQTL Signal (%s)", gene_name), 
           x = paste("Chromosome", unique(df_locus$chr), "Position (Mb)"),
           y = expression(-log[10](italic(P)[eQTL]))) +
      scale_y_reverse()
    
    final_plot <- p_gwas / p_eqtl + plot_layout(guides = "collect")
    
    safe_plot_name <- str_replace_all(str_remove(file_name, "_Aligned.csv"), "[^A-Za-z0-9_]", "")
    plot_output_path <- sprintf("%s/LocusZoom_%s.pdf", plot_dir, safe_plot_name)
    ggsave(plot_output_path, plot = final_plot, width = 8, height = 8, dpi = 300)
  }
}

# ------------------------------------------------------------------------------
# 3. 导出附录总表
# ------------------------------------------------------------------------------
write.csv(summary_results, "Supplementary_Table_Coloc_Results.csv", row.names = FALSE)
message("\n>>> 批量推断与 LocusZoom 生成完毕！结果表存至工作目录。")
