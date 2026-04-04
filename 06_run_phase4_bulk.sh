#!/bin/bash
# 06_run_phase4_bulk.sh
# Phase 4 Bulk 验证专用 Runner

EXPOSURE_NAME="$1"
OUTCOME_NAME="$2"
MAX_JOBS="${GLOBAL_MAX_JOBS:-30}" 

# ================= 配置区域 =================
# Bulk 数据通常基因数较多但SNP也多，Chunk 可以保持 20-50
CHUNK_SIZE=20  
# [修改点] 使用 Bulk 专用的 R 脚本
R_SCRIPT="06validation_analysis_lite_bulkeqtl.R"

# [修改点] 输入文件路径必须匹配 Script 12 的输出
# Script 12 输出路径示例: exposure/eqtlgen/eqtlgen_r2_0.2_validation.qs
INPUT_QS="exposure/${EXPOSURE_NAME}/${EXPOSURE_NAME}_r2_0.2_validation.qs"

# 检查文件
if [ ! -f "$INPUT_QS" ]; then
    echo "[Error] Input file not found: $INPUT_QS"
    echo "       Please run Script 12 (filter_and_subset) first."
    exit 1
fi

# 临时工作目录
WORK_DIR="${EXPOSURE_NAME}_${OUTCOME_NAME}"
LOG_DIR="./logs/phase4_${EXPOSURE_NAME}_${OUTCOME_NAME}"
mkdir -p "$LOG_DIR"

echo ">>> [Phase 4] Start: $EXPOSURE_NAME vs $OUTCOME_NAME"

# 1. 获取基因总数 (利用 R 脚本中的逻辑来读取 exposure 列，增加兼容性)
# 注意：这里我们临时用 R 读取一下文件长度，确保 exposure 列被正确识别
TOTAL_EXP=$(Rscript -e "
library(qs); dat <- qread('${INPUT_QS}'); 
col <- if('exposure' %in% names(dat)) 'exposure' else if('symbol.exposure' %in% names(dat)) 'symbol.exposure' else 'gene.exposure';
cat(length(unique(dat[[col]])))
")

echo "    Total Genes: $TOTAL_EXP | Chunk Size: $CHUNK_SIZE"

# 2. 定义并行函数
run_bulk_task() {
    local start=$1
    local chunk=$2
    local total=$3
    local exp=$4
    local out=$5
    local qs_path=$6
    
    local end=$((start + chunk - 1))
    if [ "$end" -gt "$total" ]; then end=$total; fi
    
    env \
        EXPOSURE_DATA="$exp" \
        OUTCOME="$out" \
        START="$start" \
        END="$end" \
        INPUT_QS="$qs_path" \
        Rscript --vanilla "06validation_analysis_lite_bulkeqtl.R" > "${LOG_DIR}/job_${start}.log" 2>&1
}
export -f run_bulk_task
export LOG_DIR

# 3. 执行并行
seq 1 "$CHUNK_SIZE" "$TOTAL_EXP" | parallel \
    --jobs "$MAX_JOBS" \
    --line-buffer \
    --tagstring "[Batch {}]" \
    run_bulk_task {} "$CHUNK_SIZE" "$TOTAL_EXP" "$EXPOSURE_NAME" "$OUTCOME_NAME" "$INPUT_QS"

# ================= 结果合并 =================
echo ">>> [Phase 4] Merging results..."

TEMP_RES_DIR="${WORK_DIR}/results"
MERGED_FILENAME="${EXPOSURE_NAME}_${OUTCOME_NAME}_Phase4_Bulk_Result.txt"
MERGED_PATH="${WORK_DIR}/${MERGED_FILENAME}"

if ls ${TEMP_RES_DIR}/results_val_chunk_*.txt 1> /dev/null 2>&1; then
    # 取表头
    ls ${TEMP_RES_DIR}/results_val_chunk_*.txt | head -n 1 | xargs head -n 1 > "$MERGED_PATH"
    # 合并内容
    cat ${TEMP_RES_DIR}/results_val_chunk_*.txt | grep -v "gene_symbol" >> "$MERGED_PATH"
    echo "    Merged to: $MERGED_PATH"
else
    echo "    Warning: No results generated for $EXPOSURE_NAME vs $OUTCOME_NAME"
fi
