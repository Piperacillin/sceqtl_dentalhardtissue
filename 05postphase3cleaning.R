rm(list=ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

onek1k_p1_ar <- fread('./result/onek1k_Phase1_All_Results.txt')
onemscblood_p1_ar <- fread('./result/1mscblood_Phase1_All_Results.txt')
head(onemscblood_p1_ar)
# ==============================================================================
# 脚本名称: clean_1mscblood_dynamic.R
# 功能: 
#### 1. 读取 1mscblood 主结果 ####
# 2. 依据 UT 组计算动态特异性 (Z-test)
# 3. 生成符合要求的汇总表格
# ==============================================================================
rm(list=ls())

# 加载包
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# 1. 读取数据
input_file <- "./result/1mscblood_Phase1_All_Results.txt"
if (!file.exists(input_file)) stop("Input file not found!")

message(">>> Reading input file: ", input_file)
dt <- fread(input_file)

# 确保数值列类型正确
cols_to_numeric <- c("beta", "se", "FDR", "p")
dt[, (cols_to_numeric) := lapply(.SD, as.numeric), .SDcols = cols_to_numeric]

# ==============================================================================
# 2. 准备基线数据 (UT Reference)
# ==============================================================================
message(">>> Preparing Baseline (UT) references...")

# 提取 UT 组的数据，保留用于匹配的关键列 (Gene, Cell, Outcome) 和 统计量 (Beta, SE)
# 注意：同一个 Gene-Cell 可能对应多个 Outcome，必须按 Outcome 匹配
dt_ut <- dt[act_time == "UT", .(gene_symbol, cell_type, outcome, beta_ut = beta, se_ut = se)]

# 检查是否有重复 (即同一个 Gene-Cell-Outcome 在 UT 里有多条记录，理论上不应发生，除非有多个 SNP)
# 如果是基于 Gene-Level 的汇总，这里假设每个 Gene-Cell-Outcome 只有一条最佳 SNP 记录
# 如果原始数据包含多个 SNP，我们取 P 值最小的或者去重
dt_ut <- unique(dt_ut, by = c("gene_symbol", "cell_type", "outcome"))

# ==============================================================================
# 3. 计算动态差异 (Z-test)
# ==============================================================================
message(">>> Calculating Dynamic Specificity (Z-test vs UT)...")

# 将原始数据与 UT 数据合并 (Left Join)
dt_merged <- merge(dt, dt_ut, by = c("gene_symbol", "cell_type", "outcome"), all.x = TRUE)

# 计算 Z-score 和 P_diff
# Z = (b_stim - b_ut) / sqrt(se_stim^2 + se_ut^2)
dt_merged[, z_score := (beta - beta_ut) / sqrt(se^2 + se_ut^2)]
dt_merged[, p_diff := 2 * pnorm(-abs(z_score))]

# ==============================================================================
# 4. 生成最终状态列
# ==============================================================================
message(">>> Generating final status columns...")

# 定义逻辑:
# 1. is_sig_fdr: FDR < 0.05
# 2. dynamic_status:
#    - 如果 FDR >= 0.05 -> "Not_Sig" (不显著，无需讨论动态)
#    - 如果 act_time == "UT" -> "Baseline" (这是基线信号)
#    - 如果 p_diff < 0.05 -> "PASS" (是动态 eQTL)
#    - 如果 p_diff >= 0.05 -> "FAIL" (虽然显著，但跟基线没区别，是通用 eQTL)
#    - 如果 p_diff is NA -> "No_UT_Data" (没有对应的 UT 数据可比)

dt_final <- dt_merged %>%
  mutate(
    # 第三列：是否 FDR < 0.05
    is_fdr_sig = if_else(FDR < 0.05, TRUE, FALSE),
    
    # 第四列：是否是动态 eQTL
    dynamic_eqtl_status = case_when(
      FDR >= 0.05 ~ "Not_Sig",      # 根本不显著
      act_time == "UT" ~ "Baseline", # 是基线数据
      is.na(p_diff) ~ "No_UT_Data",  # 缺少基线对照
      p_diff < 0.05 ~ "PASS",        # 显著且与基线不同 -> 动态
      TRUE ~ "FAIL"                  # 显著但与基线相同 -> 静态/通用
    )
  ) %>%
  # 按要求选择和排序列
  # 第一列 Gene, 第二列 Cell, 第三列 Sig, 第四列 Dynamic
  # 同时也保留 Outcome, Time, SNP, Beta 等信息供后续分析
  select(gene_symbol, cell_type, act_time,  is_fdr_sig, dynamic_eqtl_status, 
         FDR, p_diff, outcome, beta, se, everything()) %>%
  arrange(gene_symbol, cell_type, act_time)


input_file_e <- "./result/1mscblood_Phase1_Egger_Cochransq_Results.txt"
if (!file.exists(input_file_e)) stop("Input file not found!")

message(">>> Reading input file: ", input_file)
eggercochranesq <- fread(input_file_e)

head(eggercochranesq)


# ==============================================================================
# 5. 输出与检查
# ==============================================================================
output_file <- "./result/1mscblood_Cleaned_Dynamic_Master.csv"
message(">>> Saving cleaned master table to: ", output_file)
write.csv(dt_final, output_file, row.names = FALSE)

# 打印一些统计信息
message("\n--- Summary Statistics ---")
message("Total rows: ", nrow(dt_final))
message("Significant Hits (FDR<0.05): ", sum(dt_final$is_fdr_sig))
message("\nDynamic Status counts (for Significant Hits):")
print(table(dt_final[dt_final$is_fdr_sig == TRUE, ]$dynamic_eqtl_status))


# ==============================================================================
# 脚本名称: clean_1mscblood_qc.R
# 功能: 
#### 1. 读取 Egger & Cochran's Q 结果####
# 2. 添加 QC 状态标签 (PASS/FAIL)
# 3. 输出清洗后的 QC 表格
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# 1. 读取数据
input_file_e <- "./result/1mscblood_Phase1_Egger_Cochransq_Results.txt"
if (!file.exists(input_file_e)) stop("Input file not found!")

message(">>> Reading QC input file: ", input_file_e)
egger_dt <- fread(input_file_e)

# 2. 数据清洗与打标签
# 逻辑：P值 >= 0.05 为 PASS (无显著多效性/异质性)，否则为 FAIL
# 注意处理 NA 值的情况 (有些行可能因为SNP少无法计算 Egger)

qc_final <- egger_dt %>%
  mutate(
    # --- Egger 多效性检查 ---
    # 如果 P >= 0.05 -> PASS (Good)
    # 如果 P < 0.05  -> FAIL (Bad, 多效性存在)
    # 如果是 NA      -> NA
    egger_status = case_when(
      is.na(egger_intercept_pvalue) ~ "NA", 
      egger_intercept_pvalue >= 0.05 ~ "PASS",
      TRUE ~ "FAIL"
    ),
    
    # --- Cochran's Q 异质性检查 ---
    # 如果 P >= 0.05 -> PASS (Homogeneous)
    # 如果 P < 0.05  -> FAIL (Heterogeneous)
    heterogeneity_status = case_when(
      is.na(cochrans_q_pval) ~ "NA",
      cochrans_q_pval >= 0.05 ~ "PASS",
      TRUE ~ "FAIL"
    )
  ) %>%
  # 调整列顺序，把重要的状态放在前面
  select(gene_symbol, cell_type, act_time, outcome, nsnp, 
         egger_status, heterogeneity_status,
         egger_intercept_pvalue, cochrans_q_pval, i2,
         everything()) %>%
  arrange(gene_symbol, cell_type, act_time)

# 3. 保存
output_file <- "./result/1mscblood_Cleaned_QC_Master.csv"
message(">>> Saving cleaned QC table to: ", output_file)
write.csv(qc_final, output_file, row.names = FALSE)

# 4. 简单统计
message("\n--- QC Summary ---")
message("Total QC records: ", nrow(qc_final))
message("\nEgger Status Counts:")
print(table(qc_final$egger_status))
message("\nHeterogeneity Status Counts:")
print(table(qc_final$heterogeneity_status))


# ==============================================================================
# 脚本名称: 04_clean_sensitivity_to_wide.R
# 功能: 
#### 1. 读取灵敏度分析结果 (cML, Robust IVW, MR-RAPS, Weighted Median等)####
# 2. 定义 PASS/FAIL (P < 0.05 为 PASS，表示验证成功)
# 3. 长表转宽表: 每种算法一列
# ==============================================================================

# 1. 读取数据
input_file_s <- "./result/1mscblood_Phase1_sensitivity_test_Results.txt"
if (!file.exists(input_file_s)) stop("Input file not found!")

message(">>> Reading Sensitivity input file: ", input_file_s)
sens_dt <- fread(input_file_s)

# 2. 数据清洗与状态定义
# 注意：MR-PRESSO Global Test 通常 P>0.05 是好(无多效性)，但这里 sens_dt 里的 pval 列
# 对于 "MR-PRESSO (Raw)" 或 "Outlier-corrected" 通常是指因果效应的 P 值。
# 既然是验证因果性是否稳健，我们统一标准：P < 0.05 即认为该方法支持因果结论 (PASS)

sens_clean <- sens_dt %>%
  # 过滤掉 method 为 NA 的行
  filter(!is.na(method)) %>%
  mutate(
    # 定义状态
    validation_status = if_else(as.numeric(pval) < 0.05, "PASS", "FAIL"),
    
    # 规范化方法名称 (变成适合做列名的格式)
    # 例如: "Robust IVW" -> "Robust_IVW", "MR-PRESSO (Raw)" -> "MR_PRESSO_Raw"
    method_clean = str_replace_all(method, "[ \\(\\)-]+", "_") %>% 
      str_remove("_$") # 去除末尾可能的下划线
  ) %>%
  # 只保留需要的列
  select(gene_symbol, cell_type, act_time, outcome, method_clean, validation_status)

# 3. 长表转宽表 (Pivoting)
# 每一行变成: Gene, Cell, Time, Outcome, cML_status, Robust_IVW_status ...
sens_wide <- sens_clean %>%
  pivot_wider(
    names_from = method_clean,
    values_from = validation_status,
    names_prefix = "Sens_" # 加上前缀区分，例如 Sens_cML
  )

sens_wide_processed <- sens_wide %>%
  # 第一步：根据逻辑生成新列 Sens_MR_PRESSO
  # 使用 case_when 按顺序判断，一旦满足条件即停止
  mutate(Sens_MR_PRESSO = case_when(
    # 1. 任意一列是 PASS -> 记作 PASS
    Sens_MR_PRESSO_Raw == "PASS" | Sens_MR_PRESSO_Outlier_corrected == "PASS" ~ "PASS",
    
    # 2. 两列都是 NA -> 记作 NA (这一步先保持为R的NA格式，最后一步统一替换)
    is.na(Sens_MR_PRESSO_Raw) & is.na(Sens_MR_PRESSO_Outlier_corrected) ~ NA_character_,
    
    # 3. 剩余的所有情况 -> 记作 FAIL
    # 这里包含了：都是FAIL、以及一个FAIL一个NA的情况
    TRUE ~ "FAIL"
  )) %>%
  # 2. 删除原始的两列
  select(-Sens_MR_PRESSO_Raw, -Sens_MR_PRESSO_Outlier_corrected) %>%
  # 第二步：将整个数据框中所有的 NA (包括原有和新生成的) 填充为 "Not Available"
  # 使用 across(everything()) 遍历所有列
  mutate(across(everything(), ~ replace_na(as.character(.), "Not_Available")))

# 查看结果
head(sens_wide_processed)

# 4. 保存结果
output_file <- "./result/1mscblood_Cleaned_Sensitivity_Wide.csv"
message(">>> Saving cleaned Sensitivity table to: ", output_file)
write.csv(sens_wide_processed, output_file, row.names = FALSE)

# 5. 打印统计信息
message("\n--- Sensitivity Methods Summary ---")
message("Unique methods found: ", paste(unique(sens_clean$method_clean), collapse = ", "))
message("Total rows in wide table: ", nrow(sens_wide))

# ==============================================================================
# 脚本名称: 05_clean_phase2_bonferroni.R
# 功能: 
#### 1. 读取 Phase 2 (R2=0.01 + Bonferroni) 结果 ####
# 2. 将 bonferroni_pass 列转换为 PASS/FAIL 状态
# 3. 生成简化列名 R2_01_Bonf
# ==============================================================================


# 1. 读取数据
input_file <- "./result/1mscblood_Phase2_Final_Report.csv"
if (!file.exists(input_file)) stop("Input file not found!")

message(">>> Reading Phase 2 input file: ", input_file)
dt_p2 <- fread(input_file)

# 2. 数据清洗与状态定义
# 简化列名建议: R2_01_Bonf (代表 R2=0.01 with Bonferroni correction)

dt_clean <- dt_p2 %>%
  mutate(
    # 将逻辑值转换为字符状态
    # TRUE -> PASS
    # FALSE -> FAIL
    R2_01_Bonf = if_else(as.logical(bonferroni_pass), "PASS", "FAIL")
  ) %>%
  # 只保留用于合并的索引列和新的状态列
  # 注意：这里假设 gene_symbol, cell_type, act_time 已经存在于文件中
  select(gene_symbol, cell_type, act_time, outcome, R2_01_Bonf)

# 3. 保存清洗后的宽表
output_file <- "./result/1mscblood_Cleaned_Phase2_Bonf.csv"
message(">>> Saving cleaned Phase 2 table to: ", output_file)
write.csv(dt_clean, output_file, row.names = FALSE)

# 4. 打印统计
message("\n--- Phase 2 Bonferroni Summary ---")
print(table(dt_clean$R2_01_Bonf))

# ==============================================================================
# 脚本名称: 06_merge_1mscblood_final.R
# 功能: 
# 1. 以 1mscblood_Cleaned_Dynamic_Master 为核心
# 2. 合并 QC (Egger/Q), Sensitivity (Wide), Phase 2 (Bonf)
# 3. 输出 1mscblood_Final_Integrated_Master.csv
# ==============================================================================

rm(list = ls()) # 清空环境
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# ==============================================================================
# 1. 定义所有输入文件路径
# ==============================================================================
# A. 主表 (核心: 含 Beta, P, FDR, Dynamic Status)
f_main <- "./result/1mscblood_Cleaned_Dynamic_Master.csv"

# B. QC 表 (含 Egger/Heterogeneity Status)
f_qc   <- "./result/1mscblood_Cleaned_QC_Master.csv"

# C. 灵敏度表 (含 Sens_cML, Sens_Robust_IVW 等)
f_sens <- "./result/1mscblood_Cleaned_Sensitivity_Wide.csv"

# D. Phase 2 表 (含 R2_01_Bonf)
f_p2   <- "./result/1mscblood_Cleaned_Phase2_Bonf.csv"

# ==============================================================================
# 2. 读取并检查核心主表
# ==============================================================================
message(">>> [Step 1] Loading Main Master Table...")
if (!file.exists(f_main)) stop("Error: Main table not found at ", f_main)

dt_final <- fread(f_main)
message(paste0("    Rows: ", nrow(dt_final)))

# 定义合并的“锚点” (Keys)
# 对于 1mscblood，必须包含 act_time
join_keys <- c("gene_symbol", "cell_type", "act_time", "outcome")

# ==============================================================================
# 3. 合并 QC 数据
# ==============================================================================
message(">>> [Step 2] Merging QC Data...")

if (file.exists(f_qc)) {
  dt_qc <- fread(f_qc)
  
  # 挑选需要的列 (避免重复 beta, se 等列)
  # 我们只需要 QC 的状态和具体的 P 值
  cols_qc_keep <- c(join_keys, "egger_status", "heterogeneity_status", 
                    "egger_intercept_pvalue", "cochrans_q_pval", "i2")
  
  # 为了防止 QC 表缺某些列导致报错，取交集
  cols_to_use <- intersect(colnames(dt_qc), cols_qc_keep)
  
  dt_final <- dt_final %>%
    left_join(dt_qc %>% select(all_of(cols_to_use)), by = join_keys)
  
  message("    + QC data merged.")
} else {
  warning("    ! QC file not found. Skipping.")
}

# ==============================================================================
# 4. 合并 Sensitivity 数据
# ==============================================================================
message(">>> [Step 3] Merging Sensitivity Data...")

if (file.exists(f_sens)) {
  dt_sens <- fread(f_sens)
  
  # 灵敏度表已经是宽表了，直接合并即可
  # 同样只保留 key 和 以 Sens_ 开头的列，以防万一
  sens_cols <- colnames(dt_sens)[grep("^Sens_", colnames(dt_sens))]
  cols_to_use <- c(join_keys, sens_cols)
  
  dt_final <- dt_final %>%
    left_join(dt_sens %>% select(all_of(cols_to_use)), by = join_keys)
  
  message("    + Sensitivity data merged.")
} else {
  warning("    ! Sensitivity file not found. Skipping.")
}

# ==============================================================================
# 5. 合并 Phase 2 Bonferroni 数据
# ==============================================================================
message(">>> [Step 4] Merging Phase 2 Bonferroni Data...")

if (file.exists(f_p2)) {
  dt_p2 <- fread(f_p2)
  
  # 只需要 R2_01_Bonf 这一列状态
  cols_p2_keep <- c(join_keys, "R2_01_Bonf")
  
  dt_final <- dt_final %>%
    left_join(dt_p2 %>% select(all_of(cols_p2_keep)), by = join_keys)
  
  message("    + Phase 2 data merged.")
} else {
  warning("    ! Phase 2 file not found. Skipping.")
}

# ==============================================================================
# 6. 最终清洗 (填充 NA)
# ==============================================================================
message(">>> [Step 5] Finalizing and Filling NAs...")

# 定义需要填充 "Not_Tested" 或 "Not_Calculated" 的状态列
status_cols <- c("egger_status", "heterogeneity_status", "R2_01_Bonf")
# 动态获取所有 Sens_ 开头的列
sens_cols_final <- colnames(dt_final)[grep("^Sens_", colnames(dt_final))]
all_cols_to_fill <- c(status_cols, sens_cols_final)

# 仅对表中实际存在的列进行填充
cols_exist <- intersect(all_cols_to_fill, colnames(dt_final))

if(length(cols_exist) > 0) {
  dt_final <- dt_final %>%
    mutate(across(all_of(cols_exist), ~replace_na(., "Not_Avaliable")))
}
head(dt_final)
# ==============================================================================
# 7. 保存结果
# ==============================================================================
out_file <- "./result/1mscblood_Final_Integrated_Master.csv"
message(">>> [Output] Saving Final Table to: ", out_file)

write.csv(dt_final, out_file, row.names = FALSE)

message(">>> Done! 1mscblood integration complete.")


# ==============================================================================
# 脚本名称: 07_merge_cd4dynamic_validation.R
# 功能: 
# 1. 读取 1mscblood 总表
# 2. 读取 CD4 Dynamic Phase 3 验证结果
# 3. 聚合 CD4 Dynamic 结果 (按 Gene+Outcome 维度，取最小 P 值判断 PASS/FAIL)
# 4. 将结果合并回总表 (仅针对 CD4T 细胞行)
# ==============================================================================

rm(list = ls()) # 清空环境
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ==============================================================================
# 1. 读取数据
# ==============================================================================
# A. 1mscblood 总表
file_master <- "./result/1mscblood_Final_Integrated_Master.csv"
if (!file.exists(file_master)) stop("Error: Master table not found at ", file_master)
message(">>> Loading Master Table...")
dt_master <- fread(file_master)

# B. CD4 Dynamic Phase 3 结果
file_cd4dy <- "./result/CD4_Dynamic_Phase3_Validation_All_Results.txt"
if (!file.exists(file_cd4dy)) stop("Error: CD4 Dynamic result not found at ", file_cd4dy)
message(">>> Loading CD4 Dynamic Phase 3 Results...")
dt_cd4dy <- fread(file_cd4dy)

# ==============================================================================
# 2. 聚合 CD4 Dynamic 数据 (核心逻辑)
# ==============================================================================
message(">>> Aggregating CD4 Dynamic Data...")

# 逻辑：
# 1. 按 gene_symbol 和 outcome 分组 (忽略 act_time 和 cell_type 的差异)
# 2. 计算最小 P 值 (min_p)
# 3. 判定状态: 如果 min_p < 0.05 -> PASS, 否则 FAIL

dt_cd4dy_agg <- dt_cd4dy %>%
  group_by(gene_symbol, outcome) %>%
  summarise(
    min_p = min(p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    temp_cd4dy_status = if_else(min_p < 0.05, "PASS", "FAIL")
  ) %>%
  select(gene_symbol, outcome, temp_cd4dy_status)

message(paste0("    Unique Gene-Outcome pairs in validation set: ", nrow(dt_cd4dy_agg)))

# ==============================================================================
# 3. 合并到总表
# ==============================================================================
message(">>> Merging into Master Table...")

# 1. 先把聚合后的状态 Left Join 上去
#    这里会匹配所有 gene_symbol 和 outcome，不管 cell_type 是什么
dt_merged <- dt_master %>%
  left_join(dt_cd4dy_agg, by = c("gene_symbol", "outcome"))

# 2. 生成最终列 cd4dynamic (应用业务逻辑)
#    逻辑：
#    - 如果 cell_type 不包含 "CD4T" -> "Not_Applicable" (非CD4细胞不适用此验证)
#    - 如果 cell_type 包含 "CD4T":
#         - 如果 temp_cd4dy_status 有值 (匹配上了) -> 使用该值 (PASS/FAIL)
#         - 如果 temp_cd4dy_status 为 NA (没匹配上) -> "Not_Available"

dt_final <- dt_merged %>%
  mutate(
    cd4dynamic = case_when(
      # 条件 1: 不是 CD4T 细胞 (使用 str_detect 增强容错)
      !str_detect(cell_type, "CD4T") ~ "Not_Applicable",
      
      # 条件 2: 是 CD4T 细胞，且匹配到了验证结果
      !is.na(temp_cd4dy_status) ~ temp_cd4dy_status,
      
      # 条件 3: 是 CD4T 细胞，但没匹配到结果
      TRUE ~ "Not_Available"
    )
  ) %>%
  # 移除临时列
  select(-temp_cd4dy_status)

# ==============================================================================
# 4. 保存与检查
# ==============================================================================
output_file <- "./result/1mscblood_Final_Integrated_Master_v2.csv"
message(">>> [Output] Saving Updated Master Table to: ", output_file)
write.csv(dt_final, output_file, row.names = FALSE)

# 打印统计信息
message("\n--- Summary of 'cd4dynamic' column ---")
# 只统计 CD4T 细胞行的分布情况
cd4_rows <- dt_final %>% filter(str_detect(cell_type, "CD4T"))
print(table(cd4_rows$cd4dynamic))

message("\nDone.")

####onek1k####
# 脚本名称: 10_merge_onek1k_final_static.R
# 功能: OneK1K 静态数据清洗、验证与汇总 (集成 DICE 谱系验证)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# Part 1: 读取并清洗 OneK1K 主结果 (Static Master) =============================
message("\n>>> [Part 1] Cleaning OneK1K Main Results...")

input_file <- "./result/onek1k_Phase1_All_Results.txt"
if (!file.exists(input_file)) stop("OneK1K Main Input file not found!")

dt <- fread(input_file)

# 1.1 确保数值列类型正确 -------------------------------------------------------
cols_to_numeric <- c("beta", "se", "FDR", "p")
cols_exist <- intersect(cols_to_numeric, colnames(dt))
dt[, (cols_exist) := lapply(.SD, as.numeric), .SDcols = cols_exist]

# 1.2 生成基础主表 -------------------------------------------------------------
# 注意：静态数据不计算动态差异 (Z-test)，直接根据 FDR 判断
dt_main <- dt %>%
  mutate(
    is_fdr_sig = if_else(FDR < 0.05, TRUE, FALSE),
    study_type = "Static_OneK1K"
  ) %>%
  # 调整列顺序
  select(gene_symbol, cell_type, is_fdr_sig, study_type, 
         FDR, p, outcome, beta, se, everything())

message("   Main table prepared. Rows: ", nrow(dt_main))


# Part 2: 清洗 QC 结果 (Egger & Heterogeneity) =================================
message("\n>>> [Part 2] Cleaning QC Results...")
f_qc <- "./result/onek1k_Phase1_Egger_Cochransq_Results.txt"
dt_qc_clean <- NULL

if (file.exists(f_qc)) {
  egger_dt <- fread(f_qc)
  
  dt_qc_clean <- egger_dt %>%
    mutate(
      # Egger Status: P >= 0.05 is PASS
      egger_status = case_when(
        is.na(egger_intercept_pvalue) ~ "NA", 
        egger_intercept_pvalue >= 0.05 ~ "PASS",
        TRUE ~ "FAIL"
      ),
      # Heterogeneity Status: P >= 0.05 is PASS
      heterogeneity_status = case_when(
        is.na(cochrans_q_pval) ~ "NA",
        cochrans_q_pval >= 0.05 ~ "PASS",
        TRUE ~ "FAIL"
      )
    ) %>%
    # 只保留关键列，OneK1K 无 act_time
    select(gene_symbol, cell_type, outcome, 
           egger_status, heterogeneity_status,
           egger_intercept_pvalue, cochrans_q_pval, i2)
}


# Part 3: 清洗 Sensitivity 结果 (转宽表) =======================================
message("\n>>> [Part 3] Cleaning Sensitivity Results...")
f_sens <- "./result/onek1k_Phase1_sensitivity_test_Results.txt"
dt_sens_clean <- NULL

if (file.exists(f_sens)) {
  sens_dt <- fread(f_sens)
  
  # 3.1 清洗并转宽表 -----------------------------------------------------------
  sens_wide <- sens_dt %>%
    filter(!is.na(method)) %>%
    mutate(
      validation_status = if_else(as.numeric(pval) < 0.05, "PASS", "FAIL"),
      # 规范化列名
      method_clean = str_replace_all(method, "[ \\(\\)-]+", "_") %>% str_remove("_$")
    ) %>%
    select(gene_symbol, cell_type, outcome, method_clean, validation_status) %>%
    pivot_wider(
      names_from = method_clean, values_from = validation_status, names_prefix = "Sens_"
    )
  
  # 3.2 整合 MR-PRESSO 逻辑 ----------------------------------------------------
  if("Sens_MR_PRESSO_Raw" %in% colnames(sens_wide) | "Sens_MR_PRESSO_Outlier_corrected" %in% colnames(sens_wide)) {
    # 补全缺失列以防报错
    if(!"Sens_MR_PRESSO_Raw" %in% colnames(sens_wide)) sens_wide$Sens_MR_PRESSO_Raw <- NA
    if(!"Sens_MR_PRESSO_Outlier_corrected" %in% colnames(sens_wide)) sens_wide$Sens_MR_PRESSO_Outlier_corrected <- NA
    
    sens_wide <- sens_wide %>%
      mutate(Sens_MR_PRESSO = case_when(
        Sens_MR_PRESSO_Raw == "PASS" | Sens_MR_PRESSO_Outlier_corrected == "PASS" ~ "PASS",
        is.na(Sens_MR_PRESSO_Raw) & is.na(Sens_MR_PRESSO_Outlier_corrected) ~ NA_character_,
        TRUE ~ "FAIL"
      )) %>%
      select(-Sens_MR_PRESSO_Raw, -Sens_MR_PRESSO_Outlier_corrected)
  }
  
  dt_sens_clean <- sens_wide
}


# Part 4: 清洗 Phase 2 (Bonferroni) ============================================
message("\n>>> [Part 4] Cleaning Phase 2 Results...")
f_p2 <- "./result/onek1k_Phase2_Final_Report.csv"
dt_p2_clean <- NULL

if (file.exists(f_p2)) {
  dt_p2 <- fread(f_p2)
  dt_p2_clean <- dt_p2 %>%
    mutate(R2_01_Bonf = if_else(as.logical(bonferroni_pass), "PASS", "FAIL")) %>%
    select(gene_symbol, cell_type, outcome, R2_01_Bonf)
}


# Part 5: DICE 验证逻辑 (Lineage-based Broad Validation) =======================
message("\n>>> [Part 5] Processing DICE Validation (Lineage Mapping)...")
f_dice <- "./result/dice_Phase3_Validation_All_Results.txt" # 确保这是合并后的DICE结果
dt_dice_agg <- NULL

# 5.1 定义谱系映射函数 ---------------------------------------------------------
get_dice_lineage <- function(cell_type) {
  case_when(
    cell_type %in% c("CD4_NAIVE", "CD4_STIM", "TREG_NAIVE", "TREG_MEM", "TFH", "TH1", "TH2", "TH17", "THSTAR") ~ "CD4",
    cell_type %in% c("CD8_NAIVE", "CD8_STIM") ~ "CD8",
    cell_type %in% c("MONOCYTES", "M2") ~ "Myeloid",
    cell_type %in% c("B_CELL_NAIVE") ~ "B_Cell",
    cell_type %in% c("NK") ~ "NK",
    TRUE ~ "Unknown"
  )
}

get_onek1k_lineage <- function(cell_type) {
  case_when(
    str_detect(cell_type, "CD4") | str_detect(cell_type, "SOX4") ~ "CD4",
    str_detect(cell_type, "CD8") | str_detect(cell_type, "S100B") ~ "CD8",
    str_detect(cell_type, "Monocyte") | str_detect(cell_type, "Dendritic") ~ "Myeloid",
    str_detect(cell_type, "B Cell") | str_detect(cell_type, "Plasma") ~ "B_Cell",
    str_detect(cell_type, "Natural Killer") ~ "NK",
    TRUE ~ "Unknown"
  )
}

# 5.2 聚合 DICE 数据 -----------------------------------------------------------
if (file.exists(f_dice)) {
  dt_dice <- fread(f_dice)
  
  # 按 Gene + Lineage + Outcome 取最小 P 值 (任意子集命中即命中)
  dt_dice_agg <- dt_dice %>%
    mutate(lineage = get_dice_lineage(cell_type)) %>%
    filter(lineage != "Unknown") %>%
    group_by(gene_symbol, outcome, lineage) %>%
    summarise(
      min_p = min(p, na.rm = TRUE),
      # 记录具体是哪个子集命中的，方便写文章
      hit_subtypes = paste(unique(cell_type[p < 0.05]), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(temp_dice_status = if_else(min_p < 0.05, "PASS", "FAIL"))
  
  message("   DICE aggregated rows: ", nrow(dt_dice_agg))
} else {
  warning("   DICE validation file not found. DICE columns will be Not_Available.")
}


# Part 6: 终极合并 (Integration) ===============================================
message("\n>>> [Part 6] Merging All into dt_final_stable...")

# 6.1 初始化主表 (打上 Lineage 标签以便合并 DICE) ------------------------------
dt_final_stable <- dt_main %>%
  mutate(lineage = get_onek1k_lineage(cell_type))

# 定义静态连接键 (无 act_time)
join_keys <- c("gene_symbol", "cell_type", "outcome")

# 6.2 依次合并各模块 -----------------------------------------------------------
# 合并 QC
if (!is.null(dt_qc_clean)) {
  dt_final_stable <- dt_final_stable %>% left_join(dt_qc_clean, by = join_keys)
}

# 合并 Sensitivity
if (!is.null(dt_sens_clean)) {
  dt_final_stable <- dt_final_stable %>% left_join(dt_sens_clean, by = join_keys)
}

# 合并 Phase 2
if (!is.null(dt_p2_clean)) {
  dt_final_stable <- dt_final_stable %>% left_join(dt_p2_clean, by = join_keys)
}

# 合并 DICE (使用 Gene + Outcome + Lineage)
if (!is.null(dt_dice_agg)) {
  dt_final_stable <- dt_final_stable %>%
    left_join(dt_dice_agg %>% select(gene_symbol, outcome, lineage, temp_dice_status, hit_subtypes), 
              by = c("gene_symbol", "outcome", "lineage"))
}


# Part 7: 最终清洗与保存 =======================================================
message("\n>>> [Part 7] Finalizing Columns...")

dt_final_output <- dt_final_stable %>%
  mutate(
    # 生成最终 dice_validation 列
    dice_validation = case_when(
      lineage == "Unknown" ~ "Not_Applicable", # 非5大谱系的细胞
      !is.na(temp_dice_status) ~ temp_dice_status, # 匹配上了
      TRUE ~ "Not_Available" # 没匹配上
    )
  ) %>%
  select(-temp_dice_status) %>% # 移除临时列
  # 填充所有 NA 为 Not_Available
  mutate(across(
    any_of(c("egger_status", "heterogeneity_status", "R2_01_Bonf", "dice_validation")), 
    ~replace_na(., "Not_Available")
  )) %>%
  # 填充所有 Sensitivity 列的 NA
  mutate(across(starts_with("Sens_"), ~replace_na(., "Not_Available")))

# 保存结果
out_file <- "./result/onek1k_Final_Integrated_Master.csv"
write.csv(dt_final_output, out_file, row.names = FALSE)

message(">>> Success! Saved to: ", out_file)
message("--- Summary ---")
message("Total Rows: ", nrow(dt_final_output))
message("DICE Validated (PASS): ", sum(dt_final_output$dice_validation == "PASS"))
head(dt_final_output)
head(dt_final)
message(">>> Script Complete.")



####Final_FUSION(Kaaaaaaaaa)####
# 脚本名称: 11_merge_two_studies_final.R ---------------------------------------
# 功能: 合并 1mscblood (动态) 和 OneK1K (静态) 生成终极宽表
# RStudio Note: 使用 trailing dashes (----) 以确保大纲正确显示

rm(list = ls()) # 清空环境
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# 1. 读取数据 -------------------------------------------------------------------
# A. OneK1K (Static + DICE Validation)
file_k1k <- "./result/onek1k_Final_Integrated_Master.csv"
if (!file.exists(file_k1k)) stop("Error: OneK1K Final table not found!")
message(">>> Loading OneK1K Final Table...")
dt_k1k <- fread(file_k1k)

# B. 1mscblood (Dynamic + CD4 Validation)
# 优先读取 v2 版本
file_msc <- "./result/1mscblood_Final_Integrated_Master_v2.csv"
if (!file.exists(file_msc)) {
  file_msc <- "./result/1mscblood_Final_Integrated_Master.csv"
}
if (!file.exists(file_msc)) stop("Error: 1mscblood Final table not found!")
message(">>> Loading 1mscblood Final Table...")
dt_msc <- fread(file_msc)


# 2. 统一 Lineage 定义 (合并的关键锚点) -----------------------------------------
message(">>> Standardizing Lineages for Merge...")

# OneK1K 已经有 lineage 列了，给 1mscblood 也加上
get_msc_lineage <- function(cell_type) {
  case_when(
    str_detect(cell_type, "CD4T") ~ "CD4",
    str_detect(cell_type, "CD8T") ~ "CD8",
    str_detect(cell_type, "monocyte") | str_detect(cell_type, "DC") ~ "Myeloid",
    str_detect(cell_type, "NK") ~ "NK",
    str_detect(cell_type, "B") ~ "B_Cell", 
    TRUE ~ "Other"
  )
}

dt_msc <- dt_msc %>%
  mutate(lineage = get_msc_lineage(cell_type))


# 3. 准备合并 (列名重命名) ------------------------------------------------------
message(">>> Renaming columns to avoid conflicts...")

# 定义连接键 (保持不变)
join_keys <- c("gene_symbol", "outcome", "lineage")

# 处理 1mscblood: 除了 join_keys 外，所有列加后缀 _1mscblood
dt_msc_renamed <- dt_msc %>%
  rename_with(~ paste0(., "_1mscblood"), .cols = -all_of(join_keys))

# 处理 OneK1K: 除了 join_keys 外，所有列加后缀 _onek1k
dt_k1k_renamed <- dt_k1k %>%
  rename_with(~ paste0(., "_onek1k"), .cols = -all_of(join_keys))


# 4. 执行全连接 (Full Join) -----------------------------------------------------
message(">>> Performing Full Join (One-to-Many)...")

# OneK1K (静态) 自动匹配 1mscblood (动态) 的多个时间点
dt_merged <- full_join(dt_msc_renamed, dt_k1k_renamed, by = join_keys)

# 5. [核心新增] 清洗与去重 (Scoring & Deduplication) ----------------------------
message(">>> Cleaning Gene Symbols and Deduplicating...")

dt_cleaned <- dt_merged %>%
  mutate(
    # 5.1 备份原始基因名 (方便追溯)
    raw_gene_symbol = gene_symbol,
    
    # 5.2 清洗基因名：去除 .数字 后缀
    gene_symbol = str_remove(gene_symbol, "\\.\\d+$"),
    
    # 5.3 计算验证强度得分 (Score)
    # 逻辑：优先保留"阳性结果最多"的行
    # 注意：列名现在带有后缀，需对应修改
    validation_score = 
      (if_else(is_fdr_sig_1mscblood == TRUE, 2, 0, missing = 0)) +
      (if_else(dynamic_eqtl_status_1mscblood == "PASS", 3, 0, missing = 0)) +
      (if_else(is_fdr_sig_onek1k == TRUE, 2, 0, missing = 0)) +
      (if_else(dice_validation_onek1k == "PASS", 1, 0, missing = 0)) +
      (if_else(cd4dynamic_1mscblood == "PASS", 1, 0, missing = 0)),
    
    # 5.4 辅助排序：得分相同时，优先保留 FDR 更小的 (用 1-FDR 模拟)
    fdr_tie_breaker = coalesce(1 - FDR_1mscblood, 0) + coalesce(1 - FDR_onek1k, 0)
  ) %>%
  
  # 5.5 执行去重
  # 按照 清洗后的基因名 + Outcome + Lineage 分组
  # (同一个基因在同一个谱系同一个结局下，只保留一行)
  group_by(gene_symbol, outcome, lineage) %>%
  arrange(desc(validation_score), desc(fdr_tie_breaker)) %>%
  slice(1) %>%
  ungroup()

message(paste0("    Rows after cleaning & deduplication: ", nrow(dt_cleaned)))
message(paste0("    Removed duplicates: ", nrow(dt_merged) - nrow(dt_cleaned)))


# 6. 生成分析辅助列 (回答验证问题) ----------------------------------------------
message(">>> Generating Summary Columns...")

dt_final_project <- dt_cleaned %>%
  mutate(
    # Q1: 哪些基因两个表都有？---------------------------------------------------
    presence_status = case_when(
      !is.na(FDR_1mscblood) & !is.na(FDR_onek1k) ~ "Both_Detected",
      !is.na(FDR_1mscblood) & is.na(FDR_onek1k)  ~ "Unique_to_1mscblood",
      is.na(FDR_1mscblood) & !is.na(FDR_onek1k)  ~ "Unique_to_OneK1K",
      TRUE ~ "Unknown"
    ),
    
    # Q2: 动态基因验证情况 ------------------------------------------------------
    # 寻找那些在 1mscblood 动态显著，但在 OneK1K 静态不显著的黄金靶点
    dynamic_validation_in_static = case_when(
      # --- 动态 (1mscblood) 为主视角的分类 ---
      dynamic_eqtl_status_1mscblood == "PASS" & is_fdr_sig_onek1k == TRUE ~ "Dynamic_&_Static_Sig",
      dynamic_eqtl_status_1mscblood == "PASS" & is_fdr_sig_onek1k == FALSE ~ "Dynamic_Specific (Unique)",
      dynamic_eqtl_status_1mscblood == "PASS" & is.na(is_fdr_sig_onek1k) ~ "Dynamic_Only (Not in OneK1K)",
      
      # --- 静态 (OneK1K) 为主视角的分类 (新增) ---
      # 1. Static_Specific: OneK1K 显著，但 1mscblood 有数据且未通过动态验证
      is_fdr_sig_onek1k == TRUE & (dynamic_eqtl_status_1mscblood %in% c("FAIL", "Not_Sig", "Baseline")) ~ "Static_Specific (Unique)",
      
      # 2. Static_Only: OneK1K 显著，但 1mscblood 完全没数据
      is_fdr_sig_onek1k == TRUE & (is.na(dynamic_eqtl_status_1mscblood) | dynamic_eqtl_status_1mscblood == "No_UT_Data") ~ "Static_Only (Not in 1mscblood)",
      
      TRUE ~ "Other"
    ),
    
    # Q3: 双重外部验证结果 (Dice & CD4) -----------------------------------------
    double_validation_status = case_when(
      dice_validation_onek1k == "PASS" & cd4dynamic_1mscblood == "PASS" ~ "Dice_CD4_Double_Pass",
      dice_validation_onek1k == "PASS" ~ "Dice_Pass_Only",
      cd4dynamic_1mscblood == "PASS" ~ "CD4_Pass_Only",
      TRUE ~ "No_External_Validation"
    )
  ) %>%
  # 调整列顺序 ------------------------------------------------------------------
select(
  gene_symbol, outcome, lineage, presence_status, 
  # 1mscblood 核心
  act_time_1mscblood, dynamic_eqtl_status_1mscblood, cd4dynamic_1mscblood, 
  is_fdr_sig_1mscblood, FDR_1mscblood, beta_1mscblood, cell_type_1mscblood,
  # OneK1K 核心
  is_fdr_sig_onek1k, dice_validation_onek1k, FDR_onek1k, beta_onek1k, 
  hit_subtypes_onek1k, cell_type_onek1k,
  # 分析列
  dynamic_validation_in_static, double_validation_status,
  # 其他列
  everything()
)


# 7. 保存与统计 -----------------------------------------------------------------
out_file <- "./result/Phase3_Final_Project_Grand_Master_Table.csv"
message(">>> [Output] Saving Grand Master Table to: ", out_file)
write.csv(dt_final_project, out_file, row.names = FALSE)

message("\n--- Final Integration Summary ---")
message("Total Rows: ", nrow(dt_final_project))

message("\nPresence Status:")
print(table(dt_final_project$presence_status))

message("\nDynamic vs Static Validation:")
print(table(dt_final_project$dynamic_validation_in_static))

message("\nReplication Validation Status:")
print(table(dt_final_project$double_validation_status))

message(">>> Script Complete.")

# 脚本名称: 12_filter_and_subset_external_data.R ------------------------------
# 功能: 读取已去重的总表 -> 筛选目标基因 -> 截取 eQTLGen/GTEx -> 保存
# RStudio Note: Simplified version

rm(list = ls()) # 清空环境
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(qs)
})

# 1. 读取总表 (已去重) ----------------------------------------------------------
file_grand_master <- "./result/Phase3_Final_Project_Grand_Master_Table.csv"
if (!file.exists(file_grand_master)) stop("Grand Master Table not found!")

message(">>> Loading Grand Master Table...")
dt_master <- fread(file_grand_master)
unique(dt_master$dynamic_validation_in_static)
# 2. 筛选目标分类基因 -----------------------------------------------------------
target_categories <- c(
  "Dynamic_&_Static_Sig",
  "Dynamic_Specific (Unique)",
  "Dynamic_Only (Not in OneK1K)",
  "Static_Specific (Unique)",
  "Static_Only (Not in 1mscblood)"
)

message(">>> Filtering genes based on categories: ")
print(target_categories)

dt_filtered <- dt_master %>%
  filter(dynamic_validation_in_static %in% target_categories)

# 直接使用 gene_symbol，因为它已经在 Script 11 中清洗干净了 (无 .1 后缀)
target_gene_list <- unique(dt_filtered$gene_symbol)

message(paste0("    Total unique target genes found: ", length(target_gene_list)))

if(length(target_gene_list) == 0) stop("No genes found matching criteria!")

# 3. 处理 eQTLGen 数据集 --------------------------------------------------------
message("\n>>> Processing eQTLGen Data...")
eqtlgen_path_in  <- "./exposure/eqtlgen/eqtlgen_r2_0.2.qs"
eqtlgen_path_out <- "./exposure/eqtlgen/eqtlgen_r2_0.2_validation.qs"

if (file.exists(eqtlgen_path_in)) {
  eqtlgen_data <- qread(eqtlgen_path_in)
  
  if ("symbol.exposure" %in% colnames(eqtlgen_data)) {
    eqtlgen_sub <- eqtlgen_data %>% 
      filter(symbol.exposure %in% target_gene_list)
  } else if ("gene.exposure" %in% colnames(eqtlgen_data)) {
    warning("    'symbol.exposure' not found, trying 'gene.exposure'.")
    eqtlgen_sub <- eqtlgen_data %>% 
      filter(gene.exposure %in% target_gene_list)
  } else {
    stop("    Cannot find gene symbol column in eQTLGen data!")
  }
  eqtlgen_sub$tissue <- "Blood"
  eqtlgen_sub$phenotype.exposure <- paste0(eqtlgen_sub$phenotype.exposure,'_',eqtlgen_sub$tissue)
  message(paste0("    Rows matched: ", nrow(eqtlgen_sub)))
  qsave(eqtlgen_sub, eqtlgen_path_out)
  message(paste0("    Saved to: ", eqtlgen_path_out))
  rm(eqtlgen_data, eqtlgen_sub); gc()
} else {
  warning(paste0("    eQTLGen file not found: ", eqtlgen_path_in))
}

# 4. 处理 GTEx 数据集 -----------------------------------------------------------
message("\n>>> Processing GTEx Data...")
gtex_path_in  <- "./exposure/GTEx/GTEx_r2_0.2.qs" 
gtex_path_out <- "./exposure/GTEx/GTEx_r2_0.2_validation.qs"

if (file.exists(gtex_path_in)) {
  gtex_data <- qread(gtex_path_in)
  
  if ("gene.exposure" %in% colnames(gtex_data)) {
    gtex_sub <- gtex_data %>% 
      filter(gene.exposure %in% target_gene_list)
  } else {
    stop("    Cannot find 'gene.exposure' column in GTEx data!")
  }
  
  message(paste0("    Rows matched: ", nrow(gtex_sub)))
  qsave(gtex_sub, gtex_path_out)
  message(paste0("    Saved to: ", gtex_path_out))
  rm(gtex_data, gtex_sub); gc()
} else {
  warning(paste0("    GTEx file not found: ", gtex_path_in))
}

message("\n>>> Script Complete.")

