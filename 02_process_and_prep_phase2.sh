#!/bin/bash
# 02_process_and_prep_phase2.sh
# ================= 接收参数区域 =================
EXPOSURE_NAME="${1:-CD4_Dynamic}"
OUTCOME_NAME="${2:-PULPITIS}"
# ==============================================
# 路径定义
WORK_DIR="${EXPOSURE_NAME}_${OUTCOME_NAME}"

PHASE1_RES_DIR="${WORK_DIR}/results"
PHASE1_SENS_DIR="${WORK_DIR}/sensitivity" # 新增灵敏度分析目录

# 输出文件定义
MERGED_FILE="${WORK_DIR}/Phase1_All_Results.txt"
#注意：文件名中有空格，Shell中引用时必须加引号
MERGED_SENS_FILE="${WORK_DIR}/Phase1_sensitivity test_Results.txt" 
SIG_GENE_LIST="${WORK_DIR}/Phase1_Significant_Genes.txt"
R_SCRIPT="02.1merge_chunks_calculate_fdr.R"
# ===========================================

echo "[Processing] 1. Merging Phase 1 Primary Results..."

# --- Check if Primary Results exist ---
if ls ${PHASE1_RES_DIR}/results_primary_r2_0.2_*.txt 1> /dev/null 2>&1; then
    # 获取表头 (取第一个文件)
    ls ${PHASE1_RES_DIR}/results_primary_r2_0.2_*.txt | head -n 1 | xargs head -n 1 > "$MERGED_FILE"
    # 合并内容 (跳过表头，过滤包含 'exposure' 字样的标题行)
    cat ${PHASE1_RES_DIR}/results_primary_r2_0.2_*.txt | grep -v "exposure" >> "$MERGED_FILE"
    echo "   -> Primary results merged into: $MERGED_FILE"
else
    echo "   -> Warning: No primary result files found in $PHASE1_RES_DIR"
fi

echo "[Processing] 2. Merging Phase 1 Sensitivity Results..."

# --- 1. 新增：合并 Sensitivity 文件夹 ---
# 这里假设 sensitivity 文件夹下的文件格式一致（如果包含异质性、多效性等不同表格，建议分开合并）
if [ -d "$PHASE1_SENS_DIR" ] && ls ${PHASE1_SENS_DIR}/*0.2*.txt 1> /dev/null 2>&1; then
    # 获取灵敏度文件的表头
    ls ${PHASE1_SENS_DIR}/*0.2*.txt | head -n 1 | xargs head -n 1 > "$MERGED_SENS_FILE"
    # 合并内容
    cat ${PHASE1_SENS_DIR}/*0.2*.txt | grep -v "exposure" >> "$MERGED_SENS_FILE"
    echo "   -> Sensitivity results merged into: $MERGED_SENS_FILE"
else
    echo "   -> Warning: No sensitivity files found in $PHASE1_SENS_DIR, skipping."
fi

echo "[Processing] 2.1. Merging Phase 1 Egger & Cochran's Q Results..."
# 定义路径 (复用之前的 WORK_DIR)
EGGER_DIR="${WORK_DIR}/eggercochranesq"
MERGED_EGGER_FILE="${WORK_DIR}/Phase1_Egger_Cochransq_Results.txt"
if [ -d "$EGGER_DIR" ] && ls ${EGGER_DIR}/*0.2*.txt 1> /dev/null 2>&1; then
    # 获取表头 (取第一个文件)
    ls ${EGGER_DIR}/*0.2*.txt | head -n 1 | xargs head -n 1 > "$MERGED_EGGER_FILE"
    
    # 合并内容 (跳过表头，假设表头包含 'exposure' 或 'id.exposure')
    cat ${EGGER_DIR}/*0.2*.txt | grep -v "exposure" >> "$MERGED_EGGER_FILE"
    
    echo "   -> Egger & Cochran's Q results merged into: $MERGED_EGGER_FILE"
else
    echo "   -> Warning: No Egger/Cochran's Q files found in $EGGER_DIR, skipping."
fi


echo "[Processing] 3. Calculating FDR and Filtering (R Script)..."

# --- Check if Primary Merged file is not empty/exists before running R ---
# Checks before calling R script retained...
if [ -s "$MERGED_FILE" ]; then
    Rscript $R_SCRIPT "$MERGED_FILE" "$SIG_GENE_LIST" "$EXPOSURE_NAME"
else
    echo "Error: Primary merged file is missing or empty."
fi

echo "[Processing] Preparation complete."
