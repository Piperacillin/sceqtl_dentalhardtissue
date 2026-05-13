# ==============================================================================
# # 功能: 
# 1. 清洗 eQTLGen (Blood) 与 GTEx (Multi-tissue) 验证结果
# 2. 将两者合并到 Phase 3 Grand Master Table 中
# 3. [升级版] 生成基于“梯度打分”的最终 Phase 4 核心总表
# 4. [防守模块] 执行打分权重鲁棒性测试 (Robustness Check)
# ==============================================================================

rm(list = ls()) # 清空环境
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})
setwd("D:/OneDrive/学院与活动/postdoc/00课题/2511sc_MR/2604reboot")
# ==============================================================================
# 1. 配置路径
# ==============================================================================
file_phase3_master <- "./result/Phase3_Final_Project_Grand_Master_Table.csv"
if (!file.exists(file_phase3_master)) stop("Phase 3 Master Table not found! Run Script 11 first.")

file_eqtlgen_raw <- "./result/eqtlgen_Phase4_Bulk_All_Results.txt"
file_gtex_raw    <- "./result/GTEx_Phase4_Bulk_All_Results.txt"
file_phase4_master <- "./result/Phase4_Final_Project_Grand_Master_Table_with_Bulkadded.csv"

message(">>> Loading Phase 3 Master Table...")
dt_master <- fread(file_phase3_master)
message(paste0("    Base Rows: ", nrow(dt_master)))

# ==============================================================================
# 2. 处理 eQTLGen (Blood) 验证结果
# ==============================================================================
message("\n>>> [Part 1] Processing eQTLGen Validation...")
dt_eqtlgen_clean <- NULL

if (file.exists(file_eqtlgen_raw)) {
  raw_eqtlgen <- fread(file_eqtlgen_raw)
  dt_eqtlgen_clean <- raw_eqtlgen %>%
    mutate(
      p = as.numeric(p), 
      beta = as.numeric(beta),
      gene_symbol = str_remove(exposure, "(?i)_(blood|Blood)$")
    ) %>%
    group_by(gene_symbol, outcome) %>%
    summarise(
      eqtlgen_min_p = min(p, na.rm = TRUE),
      eqtlgen_beta = beta[which.min(p)], 
      .groups = "drop"
    ) %>%
    mutate(eqtlgen_validation = if_else(eqtlgen_min_p < 0.05, "PASS", "FAIL")) %>%
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
  dt_gtex_prep <- raw_gtex %>%
    mutate(
      p = as.numeric(p),
      extracted_tissue = str_extract(exposure, "(?<=_).*") 
    ) %>%
    mutate(tissue_label = coalesce(extracted_tissue, cell_type, "Unknown_Tissue"))
  
  dt_gtex_clean <- dt_gtex_prep %>%
    group_by(gene_symbol, outcome) %>%
    summarise(
      gtex_min_p = min(p, na.rm = TRUE),
      gtex_hit_tissues = paste(unique(tissue_label[p < 0.05]), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      gtex_validation = if_else(gtex_min_p < 0.05, "PASS", "FAIL"),
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
join_keys <- c("gene_symbol", "outcome")

if (!is.null(dt_eqtlgen_clean)) {
  dt_master <- dt_master %>% left_join(dt_eqtlgen_clean, by = join_keys)
} else {
  dt_master$eqtlgen_validation <- "Not_Tested"
  dt_master$eqtlgen_min_p <- NA
}

if (!is.null(dt_gtex_clean)) {
  dt_master <- dt_master %>% left_join(dt_gtex_clean, by = join_keys)
} else {
  dt_master$gtex_validation <- "Not_Tested"
  dt_master$gtex_min_p <- NA
  dt_master$gtex_hit_tissues <- NA
}

dt_final_phase4 <- dt_master %>%
  mutate(
    eqtlgen_validation = replace_na(eqtlgen_validation, "Not_Available"),
    gtex_validation = replace_na(gtex_validation, "Not_Available")
  )

# ==============================================================================
# 5. 生成高级分析列 (Total Evidence Score - 【梯度打分升级版】)
# ==============================================================================
message(">>> Generating Evidence Score with Gradient Weighting...")
dt <- dt_final_phase4

sens_methods <- c("Sens_cML", "Sens_Robust_IVW", "Sens_dIVW", 
                  "Sens_Weighted_median", "Sens_Weighted_mode", "Sens_MR_PRESSO")
sens_cols_1msc <- paste0(sens_methods, "_1mscblood")
sens_cols_k1k  <- paste0(sens_methods, "_onek1k")

dt_scored <- dt %>%
  rowwise() %>% 
  mutate(
    n_sens_pass_1msc = sum(c_across(all_of(sens_cols_1msc)) == "PASS", na.rm = TRUE),
    n_sens_pass_k1k  = sum(c_across(all_of(sens_cols_k1k)) == "PASS", na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  mutate(
    # --- [核心升级] B. 敏感性得分 (连续梯度打分) ---
    # 最高依然是 2.0 分，避免改变总分架构，但区分了结果的灰度区间
    score_sens_1msc = case_when(
      n_sens_pass_1msc >= 5 ~ 2.0,   # 极度稳健 (5-6个通过)
      n_sens_pass_1msc >= 3 ~ 1.5,   # 多数稳健 (3-4个通过)
      n_sens_pass_1msc >= 1 ~ 0.5,   # 弱稳健 (1-2个通过)
      TRUE                  ~ 0.0    # 不稳健 (0个)
    ),
    score_sens_k1k = case_when(
      n_sens_pass_k1k >= 5 ~ 2.0,
      n_sens_pass_k1k >= 3 ~ 1.5,
      n_sens_pass_k1k >= 1 ~ 0.5,
      TRUE                 ~ 0.0
    ),
    
    # --- C. 单细胞内部评分 ---
    score_1msc_total = 
      (if_else(coalesce(is_fdr_sig_1mscblood, FALSE), 2, 0)) +   
      (if_else(coalesce(cd4dynamic_1mscblood, "NA") == "PASS", 0.5, 0)) + 
      (if_else(coalesce(R2_01_Bonf_1mscblood, "NA") == "PASS", 1, 0)) +   
      (if_else(coalesce(egger_status_1mscblood, "NA") == "PASS", 1, 0)) + 
      (if_else(coalesce(heterogeneity_status_1mscblood, "NA") == "PASS", 1, 0)) + 
      score_sens_1msc, 
    
    score_onek1k_total = 
      (if_else(coalesce(is_fdr_sig_onek1k, FALSE), 2, 0)) +    
      (if_else(coalesce(dice_validation_onek1k, "NA") == "PASS", 0.5, 0)) + 
      (if_else(coalesce(R2_01_Bonf_onek1k, "NA") == "PASS", 1, 0)) +      
      (if_else(coalesce(egger_status_onek1k, "NA") == "PASS", 1, 0)) +    
      (if_else(coalesce(heterogeneity_status_onek1k, "NA") == "PASS", 1, 0)) + 
      score_sens_k1k, 
    
    sc_internal_score = score_1msc_total + score_onek1k_total,
    
    # --- D. 核心生物学场景分类 ---
    Scenario_Class = case_when(
      dynamic_eqtl_status_1mscblood == "PASS" & is_fdr_sig_onek1k == TRUE ~ "Robust_Core",
      dynamic_eqtl_status_1mscblood == "PASS" & (is_fdr_sig_onek1k == FALSE | is.na(is_fdr_sig_onek1k)) ~ "Stimulus_Specific",
      (dynamic_eqtl_status_1mscblood %in% c("FAIL", "Not_Sig", "Baseline")) & is_fdr_sig_onek1k == TRUE ~ "Resting_Specific",
      TRUE ~ "Other_Unclassified"
    ),
    
    # --- E. Bulk 外部标签 ---
    bulk_label = case_when(
      Scenario_Class == "Other_Unclassified" ~ "Other_Unclassified", 
      coalesce(eqtlgen_validation, "FAIL") == "PASS" | coalesce(gtex_validation, "FAIL") == "PASS" ~ "Bulk_Validated",
      TRUE ~ "Sc_Discovery" 
    )
  )

MAX_SCORE_1MSC <- 7.5
MAX_SCORE_K1K  <- 7.5

message(">>> Performing Scenario-Adaptive Normalization...")
dt_final <- dt_scored %>%
  mutate(
    theoretical_max = case_when(
      Scenario_Class == "Robust_Core"       ~ (MAX_SCORE_1MSC + MAX_SCORE_K1K), 
      Scenario_Class == "Stimulus_Specific" ~ MAX_SCORE_1MSC,                   
      Scenario_Class == "Resting_Specific"  ~ MAX_SCORE_K1K,                    
      TRUE ~ NA_real_ 
    ),
    effective_raw_score = case_when(
      Scenario_Class == "Robust_Core"       ~ sc_internal_score,       
      Scenario_Class == "Stimulus_Specific" ~ score_1msc_total,        
      Scenario_Class == "Resting_Specific"  ~ score_onek1k_total,      
      TRUE ~ 0
    ),
    score_ratio = if_else(is.na(theoretical_max), 0, effective_raw_score / theoretical_max),
    final_5scale_score = round(score_ratio * 5, 1),
    evidence_grade = case_when(
      Scenario_Class == "Other_Unclassified" ~ "Unclassified",
      score_ratio >= 0.85 ~ "Elite (5-Star)",
      score_ratio >= 0.70 ~ "Strong (4-Star)",
      score_ratio >= 0.50 ~ "Moderate (3-Star)",
      score_ratio >= 0.30 ~ "Weak (2-Star)",
      TRUE                ~ "Low_Confidence (1-Star)"
    )
  )

# 剔除并排序
dt_core <- dt_final %>% filter(Scenario_Class != "Other_Unclassified")
grade_levels <- c("Elite (5-Star)", "Strong (4-Star)", "Moderate (3-Star)", 
                  "Weak (2-Star)", "Low_Confidence (1-Star)")
scenario_levels <- c("Robust_Core", "Stimulus_Specific", "Resting_Specific")

dt_core$evidence_grade <- factor(dt_core$evidence_grade, levels = grade_levels)
dt_core$Scenario_Class <- factor(dt_core$Scenario_Class, levels = scenario_levels)

message(">>> Core Dataset Ready.")

# ==============================================================================
# 6. 保存与基础数据报告
# ==============================================================================
write.csv(dt_core, file_phase4_master, row.names = FALSE)
message(">>> Core Data Saved to: ", file_phase4_master)

# ==============================================================================
# 7. [终极防守模块] 全局打分权重扰动与 Jaccard 相似性分析
# ==============================================================================
message("\n========================================================")
message(">>> [Part 7] 全局参数扰动与 Jaccard 相似性分析 (Global Robustness)")
message("========================================================")
library(ggplot2)

# 提取当前基准状态的 5-Star 基因作为 Ground Truth
base_elite_genes <- dt_core %>% 
  filter(evidence_grade == "Elite (5-Star)") %>% 
  pull(gene_symbol) %>% 
  unique()

message(sprintf("   [Baseline] 基准 Elite 基因数量: %d", length(base_elite_genes)))

# 如果基准基因太少，提早警告
if(length(base_elite_genes) == 0) {
  warning("没有找到 Elite 基因，无法进行鲁棒性测试！请检查前面的打分逻辑。")
}

# ------------------------------------------------------------------------------
# 7.1 定义权重扰动矩阵 (Weight Perturbation Matrix)
# 顺序: c(FDR, 二次验证/Dyn, 内部复现R2, Egger, Het, 敏感性Sens)
# ------------------------------------------------------------------------------
scenarios <- list(
  "Baseline"         = c(2.0, 0.5, 1.0, 1.0, 1.0, 2.0),
  "Equal_Weights"    = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0), # 所有指标一视同仁
  "Heavy_FDR"        = c(4.0, 0.5, 1.0, 0.5, 0.5, 1.0), # 极度看重主效应 P 值
  "Heavy_QC"         = c(1.0, 0.5, 1.0, 2.0, 2.0, 3.0), # 极度看重质控(防假阳性)
  "No_QC"            = c(2.0, 1.0, 2.0, 0.0, 0.0, 1.0), # 完全不看多效性和异质性
  "Heavy_Validation" = c(2, 1, 2, 0.5, 0.5, 1.0)  # 极度看重独立验证
)

# ------------------------------------------------------------------------------
# 7.2 定义测试函数: 根据给定权重重新计算得分并提取 Elite 基因
# ------------------------------------------------------------------------------
simulate_weights <- function(dt_input, w_vec) {
  # 解析权重
  w_fdr <- w_vec[1]; w_val <- w_vec[2]; w_rep <- w_vec[3]; 
  w_egg <- w_vec[4]; w_het <- w_vec[5]; w_sen <- w_vec[6]
  
  # 计算新的理论满分 (单侧)
  max_single <- w_fdr + w_val + w_rep + w_egg + w_het + w_sen
  
  dt_sim <- dt_input %>%
    mutate(
      # 重新计算单侧总分
      sim_1msc = 
        (if_else(coalesce(is_fdr_sig_1mscblood, FALSE), w_fdr, 0)) +   
        (if_else(coalesce(cd4dynamic_1mscblood, "NA") == "PASS", w_val, 0)) + 
        (if_else(coalesce(R2_01_Bonf_1mscblood, "NA") == "PASS", w_rep, 0)) +   
        (if_else(coalesce(egger_status_1mscblood, "NA") == "PASS", w_egg, 0)) + 
        (if_else(coalesce(heterogeneity_status_1mscblood, "NA") == "PASS", w_het, 0)) + 
        (score_sens_1msc / 2.0 * w_sen), # 敏感性按比例折算新权重
      
      sim_k1k = 
        (if_else(coalesce(is_fdr_sig_onek1k, FALSE), w_fdr, 0)) +    
        (if_else(coalesce(dice_validation_onek1k, "NA") == "PASS", w_val, 0)) + 
        (if_else(coalesce(R2_01_Bonf_onek1k, "NA") == "PASS", w_rep, 0)) +      
        (if_else(coalesce(egger_status_onek1k, "NA") == "PASS", w_egg, 0)) +    
        (if_else(coalesce(heterogeneity_status_onek1k, "NA") == "PASS", w_het, 0)) + 
        (score_sens_k1k / 2.0 * w_sen),
      
      # 根据场景自适应汇总有效得分与满分
      sim_effective = case_when(
        Scenario_Class == "Robust_Core"       ~ sim_1msc + sim_k1k,
        Scenario_Class == "Stimulus_Specific" ~ sim_1msc,
        Scenario_Class == "Resting_Specific"  ~ sim_k1k,
        TRUE ~ 0
      ),
      sim_max = case_when(
        Scenario_Class == "Robust_Core"       ~ max_single * 2,
        Scenario_Class %in% c("Stimulus_Specific", "Resting_Specific") ~ max_single,
        TRUE ~ NA_real_
      ),
      # 计算新比例
      sim_ratio = sim_effective / sim_max
    )
  
  # 依然以 85% (0.85) 作为 5-Star Elite 的门槛
  return(dt_sim %>% filter(sim_ratio >= 0.85) %>% pull(gene_symbol) %>% unique())
}

# ------------------------------------------------------------------------------
# 7.3 计算所有场景的 Jaccard 系数
# ------------------------------------------------------------------------------
# Jaccard = 交集大小 / 并集大小
calc_jaccard <- function(setA, setB) {
  if (length(setA) == 0 && length(setB) == 0) return(1.0)
  if (length(setA) == 0 || length(setB) == 0) return(0.0)
  return(length(intersect(setA, setB)) / length(union(setA, setB)))
}

results_list <- list()

for (scene_name in names(scenarios)) {
  test_genes <- simulate_weights(dt_core, scenarios[[scene_name]])
  jaccard_val <- calc_jaccard(test_genes, base_elite_genes)
  
  # 计算保留率(传统Overlap)作为辅助信息
  overlap_val <- ifelse(length(base_elite_genes)==0, 0, 
                        length(intersect(test_genes, base_elite_genes)) / length(base_elite_genes))
  
  results_list[[scene_name]] <- data.frame(
    Scenario = scene_name,
    Jaccard_Index = jaccard_val,
    Overlap_Rate = overlap_val,
    Gene_Count = length(test_genes)
  )
}

df_robustness <- bind_rows(results_list)

message("\n>>> 权重扰动结果汇总:")
print(df_robustness %>% mutate(Jaccard_Index = round(Jaccard_Index, 3), Overlap_Rate = round(Overlap_Rate, 3)))

# ------------------------------------------------------------------------------
# 7.4 生成供发表的可视化图表 (Publication-ready Plot)
# ------------------------------------------------------------------------------
# 排除 Baseline 自己跟自己比 (肯定是1)
df_plot <- df_robustness %>% filter(Scenario != "Baseline")

# 为了作图美观，按 Jaccard 排序
df_plot$Scenario <- factor(df_plot$Scenario, levels = df_plot$Scenario[order(df_plot$Jaccard_Index)])

p <- ggplot(df_plot, aes(x = Scenario, y = Jaccard_Index)) +
  # 添加柱状图
  geom_col(aes(fill = Jaccard_Index), width = 0.6, color = "black", alpha = 0.8) +
  # 添加数值标签
  geom_text(aes(label = sprintf("%.2f", Jaccard_Index)), hjust = -0.2, size = 5, fontface = "bold") +
  # 翻转坐标轴，横向更易读
  coord_flip(ylim = c(0, 1.1)) +
  # 渐变配色，越接近1越偏蓝色（代表高稳健性）
  scale_fill_gradient(low = "#ff9999", high = "#66b3ff", limits = c(0, 1)) +
  # 添加 0.7 的虚线作为“高一致性”参考线
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", size = 1) +
  # 完善主题与标题
  theme_minimal(base_size = 14) +
  labs(
    #title = "Global Robustness of the Evidence Scoring System",
    subtitle = "Jaccard Similarity Index across Extreme Weighting Perturbations",
    x = "Weighting Perturbation Scenario",
    y = "Jaccard Index (Similarity to Baseline Elite Core)",
    fill = "Jaccard\nIndex"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(face = "bold", size = 12),
    panel.grid.major.y = element_blank()
  )
p
# 保存图表
file_plot <- "./result/Robustness_Jaccard_Plot(Global Robustness of the Evidence Scoring System).pdf"
ggsave(file_plot, plot = p, width = 9, height = 6, dpi = 300)

message(paste0(">>> 鲁棒性可视化图表已保存至: ", file_plot))
message("========================================================")

message(">>> Script Complete.")
