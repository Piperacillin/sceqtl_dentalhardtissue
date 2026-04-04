# ==============================================================================
# # 功能: 
# 1. 清洗 eQTLGen (Blood) 验证结果
# 2. 清洗 GTEx (Multi-tissue) 验证结果 (聚合处理)
# 3. 将两者合并到 Phase 3 Grand Master Table 中
# 4. 生成最终的 Phase 4 Grand Master Table
# ==============================================================================

rm(list = ls()) # 清空环境
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ==============================================================================
# 1. 配置路径
# ==============================================================================
# 输入：Phase 3 生成的总表
file_phase3_master <- "./result/Phase3_Final_Project_Grand_Master_Table.csv"
if (!file.exists(file_phase3_master)) stop("Phase 3 Master Table not found! Run Script 11 first.")

# 输入：Phase 4 原始结果 (由 bash 脚本合并生成)
file_eqtlgen_raw <- "./result/eqtlgen_Phase4_Bulk_All_Results.txt"
file_gtex_raw    <- "./result/GTEx_Phase4_Bulk_All_Results.txt"

# 输出
file_phase4_master <- "./result/Phase4_Final_Project_Grand_Master_Table_with_Bulkadded.csv"

message(">>> Loading Phase 3 Master Table...")
dt_master <- fread(file_phase3_master)
message(paste0("    Base Rows: ", nrow(dt_master)))


# ==============================================================================
# 2. 处理 eQTLGen (Blood) 验证结果 (修改版)
# ==============================================================================
message("\n>>> [Part 1] Processing eQTLGen Validation...")

dt_eqtlgen_clean <- NULL

if (file.exists(file_eqtlgen_raw)) {
  raw_eqtlgen <- fread(file_eqtlgen_raw)
  
  dt_eqtlgen_clean <- raw_eqtlgen %>%
    mutate(
      p = as.numeric(p), 
      beta = as.numeric(beta),
      
      # --- [核心修改] 清洗 exposure 列 ---
      # 逻辑：去除结尾的 _Blood 或 _blood，提取前面的部分作为 gene_symbol
      # (?i) 开启忽略大小写模式
      gene_symbol = str_remove(exposure, "(?i)_(blood|Blood)$")
    ) %>%
    # 按清洗后的 gene_symbol 和 outcome 分组
    group_by(gene_symbol, outcome) %>%
    summarise(
      eqtlgen_min_p = min(p, na.rm = TRUE),
      eqtlgen_beta = beta[which.min(p)], # 保留最显著那个SNP的beta
      .groups = "drop"
    ) %>%
    mutate(
      # 生成状态列
      eqtlgen_validation = if_else(eqtlgen_min_p < 0.05, "PASS", "FAIL")
    ) %>%
    select(gene_symbol, outcome, eqtlgen_validation, eqtlgen_min_p, eqtlgen_beta)
  
  message(paste0("    eQTLGen validated genes: ", sum(dt_eqtlgen_clean$eqtlgen_validation == "PASS")))
  
} else {
  warning("    eQTLGen result file not found. Columns will be NA.")
}

# ==============================================================================
# 3. 处理 GTEx (Multi-tissue) 验证结果
# ==============================================================================
message("\n>>> [Part 2] Processing GTEx Validation...")

dt_gtex_clean <- NULL

if (file.exists(file_gtex_raw)) {
  raw_gtex <- fread(file_gtex_raw)
  
  # 清洗逻辑：
  # GTEx 包含多个组织。我们需要聚合信息：
  # 1. 是否在 *任意* 组织中显著? (Broad Validation)
  # 2. 具体在哪些组织显著? (Hit Tissues)
  
  # 检查是否有 tissue 信息 (假设 cell_type 列被用作 tissue，或者从 exposure ID 中提取)
  # 如果 validation_analysis_bulk.R 里 cell_type 是 "Bulk_GTEx"，我们可能无法区分组织
  # 但通常 ID (id.exposure) 会包含组织名 (e.g., Gene_Adipose_Subcutaneous)
  
  dt_gtex_prep <- raw_gtex %>%
    mutate(
      p = as.numeric(p),
      # 尝试提取组织名: 假设 exposure 列格式为 Gene_Tissue
      # 如果没有下划线，就用 cell_type
      extracted_tissue = str_extract(exposure, "(?<=_).*") 
    ) %>%
    mutate(tissue_label = coalesce(extracted_tissue, cell_type, "Unknown_Tissue"))
  
  # 聚合
  dt_gtex_clean <- dt_gtex_prep %>%
    group_by(gene_symbol, outcome) %>%
    summarise(
      gtex_min_p = min(p, na.rm = TRUE),
      # 记录所有通过验证的组织 (P < 0.05)
      gtex_hit_tissues = paste(unique(tissue_label[p < 0.05]), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      gtex_validation = if_else(gtex_min_p < 0.05, "PASS", "FAIL"),
      # 如果 hit_tissues 为空字符串 (即全FAIL)，填 NA
      gtex_hit_tissues = if_else(gtex_hit_tissues == "", NA_character_, gtex_hit_tissues)
    ) %>%
    select(gene_symbol, outcome, gtex_validation, gtex_min_p, gtex_hit_tissues)
  
  message(paste0("    GTEx validated genes: ", sum(dt_gtex_clean$gtex_validation == "PASS")))
  
} else {
  warning("    GTEx result file not found. Columns will be NA.")
}


# ==============================================================================
# 4. 合并至总表 (Grand Integration)
# ==============================================================================
message("\n>>> [Part 3] Merging Phase 4 into Grand Master Table...")

# 合并键：Gene + Outcome
# 注意：不使用 Lineage，因为 Bulk 验证是针对 Gene 的通用验证
join_keys <- c("gene_symbol", "outcome")

# 1. 合并 eQTLGen
if (!is.null(dt_eqtlgen_clean)) {
  dt_master <- dt_master %>%
    left_join(dt_eqtlgen_clean, by = join_keys)
} else {
  dt_master$eqtlgen_validation <- "Not_Tested"
  dt_master$eqtlgen_min_p <- NA
}

# 2. 合并 GTEx
if (!is.null(dt_gtex_clean)) {
  dt_master <- dt_master %>%
    left_join(dt_gtex_clean, by = join_keys)
} else {
  dt_master$gtex_validation <- "Not_Tested"
  dt_master$gtex_min_p <- NA
  dt_master$gtex_hit_tissues <- NA
}

# 3. 填充 NA
dt_final_phase4 <- dt_master %>%
  mutate(
    eqtlgen_validation = replace_na(eqtlgen_validation, "Not_Available"),
    gtex_validation = replace_na(gtex_validation, "Not_Available")
  )


# ==============================================================================
# 5. 生成高级分析列 (Total Evidence Score)
# ==============================================================================
message(">>> Generating Evidence Score...")
head(dt_final_phase4)
#unique(dt_final_phase4$eqtlgen_validation)

dt <- dt_final_phase4
  
# 定义敏感性分析的方法列名 (方便后续计算)
sens_methods <- c("Sens_cML", "Sens_Robust_IVW", "Sens_dIVW", 
                    "Sens_Weighted_median", "Sens_Weighted_mode", "Sens_MR_PRESSO")

# 拼接后缀
sens_cols_1msc <- paste0(sens_methods, "_1mscblood")
sens_cols_k1k  <- paste0(sens_methods, "_onek1k")

# 2. 核心计算逻辑 --------------------------------------------------------------
message(">>> Calculating Scores and Classifying Scenarios...")

dt_scored <- dt %>%
  rowwise() %>% # 开启逐行计算模式 (为了计算敏感性PASS的数量)
  mutate(
    # --- A. 敏感性分析 "3 Pass Rule" 计算 ---
    # 统计 1mscblood 中 PASS 的算法数量
    n_sens_pass_1msc = sum(c_across(all_of(sens_cols_1msc)) == "PASS", na.rm = TRUE),
    # 统计 OneK1K 中 PASS 的算法数量
    n_sens_pass_k1k  = sum(c_across(all_of(sens_cols_k1k)) == "PASS", na.rm = TRUE)
  ) %>%
  ungroup() %>% # 结束逐行模式，回到向量化模式提高速度
  mutate(
    
    # --- B. 敏感性得分 (0 或 2 分) ---
    score_sens_1msc = if_else(n_sens_pass_1msc >= 3, 2, 0),
    score_sens_k1k  = if_else(n_sens_pass_k1k >= 3, 2, 0),
    
    # --- C. 单细胞内部评分 (Internal Score) ---
    # [1mscblood 维度]
    score_1msc_total = 
      (if_else(coalesce(is_fdr_sig_1mscblood, FALSE), 2, 0)) +   # FDR显著 (2分)
      (if_else(coalesce(cd4dynamic_1mscblood, "NA") == "PASS", 0.5, 0)) + # 动态验证 (0.5分)
      (if_else(coalesce(R2_01_Bonf_1mscblood, "NA") == "PASS", 1, 0)) +   # 内部复现 (1分)
      (if_else(coalesce(egger_status_1mscblood, "NA") == "PASS", 1, 0)) + # Egger QC (1分)
      (if_else(coalesce(heterogeneity_status_1mscblood, "NA") == "PASS", 1, 0)) + # Heterogeneity QC (1分)
      score_sens_1msc, # 敏感性奖励 (2分 or 0分)
    
    # [OneK1K 维度]
    score_onek1k_total = 
      (if_else(coalesce(is_fdr_sig_onek1k, FALSE), 2, 0)) +    # FDR显著 (0.5分, 静态权重略低)
      (if_else(coalesce(dice_validation_onek1k, "NA") == "PASS", 0.5, 0)) + # Dice验证 (1分)
      (if_else(coalesce(R2_01_Bonf_onek1k, "NA") == "PASS", 1, 0)) +      # 内部复现 (1分)
      (if_else(coalesce(egger_status_onek1k, "NA") == "PASS", 1, 0)) +    # Egger QC (1分)
      (if_else(coalesce(heterogeneity_status_onek1k, "NA") == "PASS", 1, 0)) + # Heterogeneity QC (1分)
      score_sens_k1k, # 敏感性奖励 (2分 or 0分)
    
    # [总分]
    sc_internal_score = score_1msc_total + score_onek1k_total,
    
    # --- D. 核心生物学场景分类 (Scenario_Class) ---
    Scenario_Class = case_when(
      # 1. Robust_Core (双显著): 动态 PASS + 静态 Significant
      dynamic_eqtl_status_1mscblood == "PASS" & is_fdr_sig_onek1k == TRUE ~ "Robust_Core",
      
      # 2. Stimulus_Specific (动态特异): 动态 PASS + 静态 Not Significant (FALSE or NA)
      dynamic_eqtl_status_1mscblood == "PASS" & (is_fdr_sig_onek1k == FALSE | is.na(is_fdr_sig_onek1k)) ~ "Stimulus_Specific",
      
      # 3. Resting_Specific (静态特异): 动态 Not PASS + 静态 Significant
      # 注: 动态状态包括 FAIL, Not_Sig, Baseline (即非PASS状态)
      (dynamic_eqtl_status_1mscblood %in% c("FAIL", "Not_Sig", "Baseline")) & is_fdr_sig_onek1k == TRUE ~ "Resting_Specific",
      
      # 4. 其他情况 (如静态也只在OneK1K里有但1mscblood没测到等)
      TRUE ~ "Other_Unclassified"
    ),
    
    # --- E. Bulk 外部标签 (External Label) ---
    bulk_label = case_when(
      # 优先级最高：如果分类不明，直接打上 Other 标签，防止污染 "Sc_Discovery"
      Scenario_Class == "Other_Unclassified" ~ "Other_Unclassified", 
      
      # 优先级第二：如果通过了 Bulk 验证
      coalesce(eqtlgen_validation, "FAIL") == "PASS" | coalesce(gtex_validation, "FAIL") == "PASS" ~ "Bulk_Validated",
      
      # 优先级第三：既明确了场景，又没被 Bulk 验证 -> 真正的单细胞发现
      TRUE ~ "Sc_Discovery" 
    )
  ) %>%
  # 调整列顺序
  select(
    gene_symbol, outcome, lineage, 
    Scenario_Class, sc_internal_score, bulk_label,
    score_1msc_total, score_onek1k_total, n_sens_pass_1msc, n_sens_pass_k1k,
    everything()
  )  

# 2. 定义理论满分常数 ----------------------------------------------------------
# 基于 Script 14 的权重设定:
# 1mscblood total max = 2(FDR) + 0.5(Dyn) + 1(Rep) + 1(Egger) + 1(Het) + 2(Sens) = 7.5
# OneK1K total max    = 2(FDR) + 0.5(Dice) + 1(Rep) + 1(Egger) + 1(Het) + 2(Sens) = 7.5

MAX_SCORE_1MSC <- 7.5
MAX_SCORE_K1K  <- 7.5

# 3. 归一化与评级计算 ----------------------------------------------------------
message(">>> Performing Scenario-Adaptive Normalization...")

dt_final <- dt_scored %>%
  mutate(
    # --- Step 1: 确定该基因的理论满分 (Denominator) ---
    theoretical_max = case_when(
      Scenario_Class == "Robust_Core"       ~ (MAX_SCORE_1MSC + MAX_SCORE_K1K), # 14.0
      Scenario_Class == "Stimulus_Specific" ~ MAX_SCORE_1MSC,                   # 7.5
      Scenario_Class == "Resting_Specific"  ~ MAX_SCORE_K1K,                    # 6.5
      TRUE ~ NA_real_ 
    ),
    
    # --- Step 2: 确定该基因的有效得分 (Numerator) ---
    # 注意：对于 Specific 场景，我们只取对应侧的分数，忽略另一侧的"噪音"
    effective_raw_score = case_when(
      Scenario_Class == "Robust_Core"       ~ sc_internal_score,       # 两边之和
      Scenario_Class == "Stimulus_Specific" ~ score_1msc_total,        # 只看 1msc 分
      Scenario_Class == "Resting_Specific"  ~ score_onek1k_total,      # 只看 K1K 分
      TRUE ~ 0
    ),
    
    # --- Step 3: 计算得分率 (Ratio) ---
    score_ratio = if_else(is.na(theoretical_max), 0, effective_raw_score / theoretical_max),
    
    # --- Step 4: 映射到 5 分制 (Final Score) ---
    # 保留 1 位小数，方便排序
    final_5scale_score = round(score_ratio * 5, 1),
    
    # --- Step 5: 赋予等级标签 (Evidence Grade) ---
    evidence_grade = case_when(
      Scenario_Class == "Other_Unclassified" ~ "Unclassified",
      score_ratio >= 0.85 ~ "Elite (5-Star)",
      score_ratio >= 0.70 ~ "Strong (4-Star)",
      score_ratio >= 0.50 ~ "Moderate (3-Star)",
      score_ratio >= 0.30 ~ "Weak (2-Star)",
      TRUE                ~ "Low_Confidence (1-Star)"
    )
  ) %>%
  # 调整列顺序，将最终评级放在最显眼的位置
  select(
    gene_symbol, outcome, lineage, 
    Scenario_Class, evidence_grade, final_5scale_score, bulk_label,
    score_ratio, effective_raw_score, theoretical_max, # 保留中间计算过程备查
    everything()
  )

# 2. 核心过滤: 剔除 Other_Unclassified
dt_core <- dt_final %>%
  filter(Scenario_Class != "Other_Unclassified")

# 3. 因子化排序 (为了作图好看，规定顺序)
# 将等级按照 星级从高到低 排序
grade_levels <- c("Elite (5-Star)", "Strong (4-Star)", "Moderate (3-Star)", 
                  "Weak (2-Star)", "Low_Confidence (1-Star)")

# 将场景按照 稳健 -> 动态 -> 静态 排序
scenario_levels <- c("Robust_Core", "Stimulus_Specific", "Resting_Specific")

dt_core$evidence_grade <- factor(dt_core$evidence_grade, levels = grade_levels)
dt_core$Scenario_Class <- factor(dt_core$Scenario_Class, levels = scenario_levels)


message(paste0("Original Rows: ", nrow(dt_final)))
message(paste0("Core Rows:     ", nrow(dt_core)))
message(">>> Core Dataset Ready.")

# ==============================================================================
# 6. 保存与展示
# ==============================================================================
message(">>> [Output] Saving Core Phase 4 Grand Master Table to: ", file_phase4_master)
write.csv(dt_core, file_phase4_master, row.names = FALSE)

message("========================================================")
message(">>> 核心数据体检报告 (Core Data Health Check)")
message("========================================================")

# --- 统计 1: 总览 ---
message("\n[1] 总体规模 (Total Numbers):")
message(sprintf("    Total Core Genes: %d", nrow(dt_core)))
message(sprintf("    Unique Outcomes:  %d", length(unique(dt_core$outcome))))
message(sprintf("    Cell Lineages:    %s", paste(unique(dt_core$lineage), collapse=", ")))

# --- 统计 2: 三大场景分布 ---
message("\n[2] 生物学场景分布 (Scenario Distribution):")
dt_core %>% 
  count(Scenario_Class) %>% 
  mutate(Percent = round(n / sum(n) * 100, 1)) %>% 
  print()

# --- 统计 3: 证据等级分布 ---
message("\n[3] 证据等级分布 (Evidence Grade):")
dt_core %>% 
  count(evidence_grade) %>% 
  mutate(Percent = round(n / sum(n) * 100, 1)) %>% 
  print()

# --- 统计 4: 核心看点 - 场景 vs 等级 (交叉表) ---
message("\n[4] 场景 x 等级 (Scenario vs Grade Matrix):")
table(dt_core$evidence_grade, dt_core$Scenario_Class) %>% print()

# --- 统计 5: 创新性分析 - 场景 vs Bulk验证 ---
message("\n[5] 创新性分析 (Novelty Check):")
message("    * Sc_Discovery = Bulk漏掉的 (新发现)")
message("    * Bulk_Validated = 已知的 (高置信)")
table(dt_core$bulk_label, dt_core$Scenario_Class) %>% print()

# --- 统计 6: 顶级 Elite 基因预览 ---
message("\n[6] 顶级 Elite (5-Star) 基因预览 (Top 5 per Scenario):")
dt_core %>%
  filter(evidence_grade == "Elite (5-Star)") %>%
  group_by(Scenario_Class) %>%
  arrange(desc(final_5scale_score)) %>%
  slice_head(n = 3) %>% # 每个场景取前3个
  select(Scenario_Class, gene_symbol, lineage, final_5scale_score, bulk_label) %>%
  print()

message("\n========================================================")


message(">>> Script Complete.")
