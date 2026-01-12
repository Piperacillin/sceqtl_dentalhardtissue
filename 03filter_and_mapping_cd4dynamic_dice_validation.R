# ==============================================================================
# 脚本名称: 03_filter_and_intersect.R
# 功能: 
# 1. 读取 result 中的 1mscblood 结果，筛选 FDR<0.05 & CD4T 的基因
# 2. 读取 exposure/CD4_Dynamic 下的 qs 文件
# 3. 对 qs 数据进行特殊的列拆分 (第一个和最后一个下划线)
# 4. 取交集并输出
# ==============================================================================

# 1. 加载必要的包
# 如果没有安装，请取消注释下一行进行安装
# install.packages(c("data.table", "tidyverse", "qs"))
rm(list=ls())
library(data.table)
library(tidyverse)
library(qs)

# ================= 配置路径 =================
# 输入文件路径
msc_result_file <- "result/1mscblood_Phase1_All_Results.txt"
cd4_dynamic_file <- "exposure/CD4_Dynamic/CD4_Dynamic_r2_0.2.qs"

# 输出文件路径 (可自定义)
output_file <- "exposure/CD4_Dynamic/CD4_Dynamic_r2_0.2_validation.qs"

# ================= Step 1: 读取并筛选 1mscblood 结果 =================
message(">>> Step 1: Reading and filtering 1mscblood results...")

if (!file.exists(msc_result_file)) {
  stop(paste("Error: File not found:", msc_result_file))
}

# 读取 txt 文件 (假设是 Tab 分隔，根据之前的 bash 脚本逻辑)
msc_data <- fread(msc_result_file)
head(msc_data)
unique(msc_data$method)

# 1. & 2. 分组并计算 FDR
msc_data <- msc_data %>%
  group_by(outcome) %>%
  mutate(FDR = p.adjust(p, method = "fdr")) %>%
  ungroup() # 记得取消分组

# 3. 检查并创建输出目录
dir.create(dirname(msc_result_file), showWarnings = FALSE, recursive = TRUE)
# 4. 保存文件
write_delim(msc_data, msc_result_file, delim = "\t")

# 检查列名是否存在
required_cols <- c("FDR", "cell_type", "gene_symbol")
if (!all(required_cols %in% colnames(msc_data))) {
  stop(paste("Error: Missing columns in 1mscblood file. Found:", paste(colnames(msc_data), collapse=", ")))
}

# 筛选: FDR < 0.05 且 cell_type 为 CD4T
# 注意: 之前的 bash 脚本生成的 cell_type 列内容可能带有引号或空格，这里做一下清理比较稳妥
target_genes_df <- msc_data %>%
  filter(as.numeric(FDR) < 0.05) %>% ####筛选标准是FDR<0.05####
  filter(str_detect(cell_type, "CD4T")) %>% # 使用 str_detect 容错性更好 (e.g. "CD4T" or "CD4T_activated")
  select(gene_symbol) %>%
  distinct()

target_gene_list <- target_genes_df$gene_symbol

message(paste("Found", length(target_gene_list), "target genes (FDR<0.05 & CD4T)."))

if (length(target_gene_list) == 0) {
  warning("No genes passed the filter! The output will be empty.")
}

# ================= Step 2: 读取并处理 CD4 Dynamic 数据 =================
message(">>> Step 2: Reading CD4 Dynamic data and splitting columns...")

if (!file.exists(cd4_dynamic_file)) {
  stop(paste("Error: File not found:", cd4_dynamic_file))
}

# 读取 QS 文件
cd4_data <- qread(cd4_dynamic_file)

# 确保 exposure 列存在
if (!"exposure" %in% colnames(cd4_data)) {
  stop("Error: 'exposure' column not found in qs file.")
}

# --- 核心逻辑: 按照第一个和最后一个下划线拆分 ---
# 正则表达式解释:
# ^([^_]+)   : 捕获组1 (gene_symbol) -> 从头开始，非下划线的字符，直到遇到第一个下划线
# _          : 匹配第一个下划线
# (.+)       : 捕获组2 (cell_type)   -> 贪婪匹配中间所有内容，直到遇到最后一个下划线
# _          : 匹配最后一个下划线 (因为前面的 .+ 是贪婪的，它会吃到只剩最后一个 _)
# ([^_]+)$   : 捕获组3 (act_time)    -> 最后一个下划线后的非下划线字符，直到结尾

pattern <- "^([^_]+)_(.+)_(.+)$"

cd4_processed <- cd4_data %>%
  mutate(
    # 使用 str_match 提取三个捕获组
    matches = str_match(exposure, pattern),
    gene_symbol = matches[, 2],
    cell_type   = matches[, 3],
    act_time    = matches[, 4]
  ) %>%
  select(-matches) # 移除临时列
head(cd4_processed)
# 检查拆分是否有 NA (即不符合规则的行)
invalid_rows <- sum(is.na(cd4_processed$gene_symbol))
if (invalid_rows > 0) {
  warning(paste("Warning:", invalid_rows, "rows in exposure column did not match the split pattern."))
}

# ================= Step 3: 取交集 =================
message(">>> Step 3: Intersecting data...")

# 筛选: 保留 cd4_processed 中 gene_symbol 在 target_gene_list 里的行
final_result <- cd4_processed %>%
  filter(gene_symbol %in% target_gene_list)

message(paste("    Original rows in CD4 Dynamic:", nrow(cd4_processed)))
message(paste("    Rows after intersection:", nrow(final_result)))
message(paste("    Unique genes in intersection:", length(unique(final_result$gene_symbol))))

# ================= Step 4: 输出结果 =================
message(paste(">>> Step 4: Saving result to", output_file))

# 创建目录(如果不存在)
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

qs::qsave(final_result, output_file)

message("Done.")


# 脚本名称: 06_onek1k_dice_validation_prep.R
# 功能: 
# 1. 读取 OneK1K 结果，筛选 FDR<0.05 的阳性结果
# 2. 基于"谱系宽泛映射表" (Lineage-based Broad Mapping) 将 OneK1K 细胞映射到 DICE 细胞
# 3. 读取 DICE 数据 (qs格式)
# 4. 拆分 DICE exposure 列 (静态格式：Gene_CellType)
# 5. 取交集并生成验证子集


# 加载必要的包
library(data.table)
library(tidyverse)
library(qs)

# ================= 1. 路径配置 =================
# 输入文件
onek1k_result_file <- "result/onek1k_Phase1_All_Results.txt"
dice_data_file     <- "exposure/dice/dice_r2_0.2.qs"

# 输出文件
output_file        <- "exposure/dice/dice_r2_0.2_validation.qs"

# ================= 2. 定义免疫细胞谱系映射表 =================
# 逻辑：OneK1K 的输入(Input) -> 对应 DICE 的所有可能的验证靶点(Targets)
# 只要属于同一谱系，我们就保留用于验证

mapping_rules <- list(
  # --- CD4+ T Cell Family ---
  "CD4" = list(
    inputs = c("CD4 Naive/Central memory T cell", 
               "CD4 Effector memory/TEMRA", 
               "CD4 SOX4 T cell"),
    targets = c("CD4_NAIVE", "CD4_STIM", "TREG_NAIVE", "TREG_MEM", 
                "TFH", "TH1", "TH2", "TH17", "THSTAR")
  ),
  
  # --- CD8+ T Cell Family ---
  "CD8" = list(
    inputs = c("CD8 Naive/Central memory T cell", 
               "CD8 Effector memory", 
               "CD8 S100B T cell"),
    targets = c("CD8_NAIVE", "CD8_STIM")
  ),
  
  # --- Myeloid Family (Monocytes & DC) ---
  "Myeloid" = list(
    inputs = c("Classic Monocyte", 
               "Non-classic Monocyte", 
               "Dendritic Cell"),
    targets = c("MONOCYTES", "M2")
  ),
  
  # --- B Cell Family ---
  "B_Cell" = list(
    inputs = c("Naïve/Immature B Cell", 
               "Memory B Cell", 
               "Plasma Cell"),
    targets = c("B_CELL_NAIVE")
  ),
  
  # --- NK Cell Family ---
  "NK" = list(
    inputs = c("Natural Killer Cell", 
               "Natural Killer Recruiting Cell"),
    targets = c("NK")
  )
)

# ================= 3. 处理 OneK1K 发现集 =================
message(">>> Step 1: Processing OneK1K results...")

if (!file.exists(onek1k_result_file)) stop("OneK1K result file not found!")

# 读取并筛选 FDR < 0.05
onek1k_data <- fread(onek1k_result_file)

# 读取 txt 文件 (假设是 Tab 分隔，根据之前的 bash 脚本逻辑)
msc_data <- fread(msc_result_file)
head(onek1k_data)
unique(onek1k_data$method)

# 1. & 2. 分组并计算 FDR
onek1k_data <- onek1k_data %>%
  group_by(outcome) %>%
  mutate(FDR = p.adjust(p, method = "fdr")) %>%
  ungroup() # 记得取消分组

# 3. 检查并创建输出目录
dir.create(dirname(onek1k_result_file), showWarnings = FALSE, recursive = TRUE)
# 4. 保存文件
write_delim(onek1k_data, onek1k_result_file, delim = "\t")

# 确保列名一致 (脚本04生成的列名通常是小写下划线)
# 如果之前的脚本生成的列名不同，请在此处 rename
sig_onek1k <- onek1k_data %>%
  filter(as.numeric(FDR) < 0.05) %>%
  select(gene_symbol, cell_type) %>% # 保留 SNP 以便后续可能的同向性检查
  distinct()

message(paste("    Found", nrow(sig_onek1k), "significant gene-cell pairs in OneK1K."))

# --- 构建查询清单 (Query List) ---
# 将 OneK1K 的 (Gene, Cell_OneK1K) 转换为 (Gene, Cell_DICE_List)

generate_queries <- function(row_idx, df) {
  gene <- df$gene_symbol[row_idx]
  cell <- df$cell_type[row_idx]
  
  # 查找该 cell 属于哪个谱系
  target_dice_cells <- character(0)
  
  for (lineage in names(mapping_rules)) {
    # 模糊匹配：因为 OneK1K 名字可能有细微变体，或者完全匹配
    # 使用 %in% 进行精确匹配 (基于我们定义的 inputs 列表)
    if (cell %in% mapping_rules[[lineage]]$inputs) {
      target_dice_cells <- mapping_rules[[lineage]]$targets
      break
    }
  }
  
  # 如果没找到映射（比如该细胞不在我们的清单里），返回 NULL
  if (length(target_dice_cells) == 0) return(NULL)
  
  # 返回扩展的数据框: Gene x All_Target_Cells
  return(data.frame(gene_symbol = gene, 
                    cell_type = target_dice_cells,
                    onek1k_source_cell = cell,
                    stringsAsFactors = FALSE))
}

# 批量构建查询表
message("    Expanding OneK1K hits to DICE targets (Lineage Expansion)...")
query_list <- lapply(1:nrow(sig_onek1k), generate_queries, df = sig_onek1k)
query_df <- do.call(rbind, query_list) %>% distinct()

message(paste("    Generated", nrow(query_df), "validation queries for DICE."))

if (nrow(query_df) == 0) {
  stop("No valid queries generated. Check mapping names or OneK1K cell type names.")
}

# ================= 4. 处理 DICE 验证集 =================
message(">>> Step 2: Processing DICE data...")

if (!file.exists(dice_data_file)) stop("DICE qs file not found!")

dice_data <- qread(dice_data_file)
head(dice_data)
# 确保 exposure 列存在
if (!"exposure" %in% colnames(dice_data)) stop("Column 'exposure' missing in DICE data")

# --- 拆分 exposure 列 (静态格式：Gene_Cell) ---
# 规则：第一个下划线前是 Gene，第一个下划线后是 Cell
# 示例: GAPDH_CD4_NAIVE -> Gene: GAPDH, Cell: CD4_NAIVE

message("    Splitting DICE exposure column...")

# 使用 separate 或 extract (tidyr)
# 这里使用 extract 正则更加稳健
# ^([^_]+)  : 开头直到第一个下划线 (Gene)
# _         : 下划线分隔符
# (.*)$     : 剩余所有内容 (Cell Type)
dice_processed <- dice_data %>%
  extract(exposure, into = c("gene_symbol", "cell_type"), 
          regex = "^([^_]+)_(.*)$", remove = FALSE)
unique(dice_processed$cell_type)
# ================= 5. 取交集 (Extraction) =================
message(">>> Step 3: Extracting validation subset (Fixed Many-to-Many)...")

# 1. 预处理 query_df：解决 Many-to-Many 问题
# 目标：保证每个 (gene_symbol, dice_cell_type) 组合在 query_df 中只出现一次
# 策略：如果同一个 DICE 靶点对应多个 OneK1K 来源，将来源细胞合并为字符串

query_unique <- query_df %>%
  group_by(gene_symbol, cell_type) %>%
  summarise(
    # 将来源细胞用分号连接，例如: "CD4 Naive/Central memory; CD4 Effector memory"
    onek1k_source_cell = paste(unique(onek1k_source_cell), collapse = "; "),
    .groups = "drop"
  )

message(paste("    Collapsed query rows from", nrow(query_df), "to", nrow(query_unique)))

# 2. 执行 Join
# 注意：确保你的 DICE 数据列名也是 "dice_cell_type"。
# 如果之前代码生成的列名是 "cell_type"，请将下方的 "dice_cell_type" 改为 "cell_type"

validation_subset <- dice_processed %>%
  # 此时变成了 Many-to-One 关系 (DICE 有多行 SNP, Query 只有一行映射)
  # 这样是安全的，不会产生数据膨胀
  inner_join(query_unique, by = c("gene_symbol", "cell_type"))

# 3. 最终检查
message(paste("    DICE total rows:", nrow(dice_processed)))
message(paste("    Validation subset rows:", nrow(validation_subset)))
message(paste("    Unique genes covered:", length(unique(validation_subset$gene_symbol))))

# 检查是否有重复 SNP (以防万一)
# 假设 SNP 列名为 "SNP" 或 "rsid"，如果没有该列可跳过此检查
if ("SNP" %in% colnames(validation_subset)) {
  dup_count <- sum(duplicated(validation_subset[, c("gene_symbol", "cell_type", "SNP")]))
  if (dup_count > 0) {
    warning(paste("Warning: Still found", dup_count, "duplicated SNPs! Check input data unique keys."))
  } else {
    message("    Success: No duplicated SNP records found.")
  }
}

# ================= 6. 保存结果 =================
message(paste(">>> Step 4: Saving to", output_file))

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
qs::qsave(validation_subset, output_file)

message("Done.")
