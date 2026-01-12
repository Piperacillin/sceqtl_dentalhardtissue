#!/bin/bash
# 05_run_validation.sh
# 验证阶段专用 Runner：并行计算 -> 内部合并 -> 输出到 result 文件夹

EXPOSURE_NAME="$1"
OUTCOME_NAME="$2"
MAX_JOBS="${GLOBAL_MAX_JOBS:-30}" # 默认并发 30 (200个基因很快就能跑完)

# ================= 配置区域 =================
R2_VAL="0.2"
CHUNK_SIZE=20  # 验证阶段基因少，Chunk切小一点，让并行更充分
R_SCRIPT="04validation_analysis_lite.R"

# 定义输入文件路径 (假设您生成的子集文件在 exposure/名字/名字_r2_0.2.qs)
# 如果文件名不同，请在此处修改
INPUT_QS="exposure/${EXPOSURE_NAME}/${EXPOSURE_NAME}_r2_0.2_validation.qs"

# 检查文件
if [ ! -f "$INPUT_QS" ]; then
    echo "[Error] Input file not found: $INPUT_QS"
    exit 1
fi

# 临时工作目录 (存放 chunk 结果)
WORK_DIR="${EXPOSURE_NAME}_${OUTCOME_NAME}"
LOG_DIR="./logs/validation_${EXPOSURE_NAME}_${OUTCOME_NAME}"
mkdir -p "$LOG_DIR"

echo ">>> [Validation] Start: $EXPOSURE_NAME vs $OUTCOME_NAME"

# 1. 获取基因总数
TOTAL_EXP=$(Rscript -e "cat(length(unique(qs::qread('${INPUT_QS}')[['exposure']])))")
echo "    Total Genes: $TOTAL_EXP | Chunk Size: $CHUNK_SIZE"

# 2. 定义并行函数
run_val_task() {
    local start=$1
    local chunk=$2
    local total=$3
    local exp=$4
    local out=$5
    local qs_path=$6
    local r2_val=$7
    
    local end=$((start + chunk - 1))
    if [ "$end" -gt "$total" ]; then end=$total; fi
    
    # 传递环境变量给 R
    env \
        EXPOSURE_DATA="$exp" \
        OUTCOME="$out" \
        START="$start" \
        END="$end" \
        INPUT_QS="$qs_path" \
        R2_THRESH="$r2_val" \
        Rscript --vanilla "04validation_analysis_lite.R" > "${LOG_DIR}/job_${start}.log" 2>&1
}
export -f run_val_task
export LOG_DIR

# 3. 执行并行 (GNU Parallel)
seq 1 "$CHUNK_SIZE" "$TOTAL_EXP" | parallel \
    --jobs "$MAX_JOBS" \
    --line-buffer \
    --tagstring "[Batch {}]" \
    run_val_task {} "$CHUNK_SIZE" "$TOTAL_EXP" "$EXPOSURE_NAME" "$OUTCOME_NAME" "$INPUT_QS" "$R2_VAL"

# ================= 结果合并 =================
echo ">>> [Validation] Merging results..."

# 目标文件名 (例如: result/CD4_Dynamic_Phase2_Validation.txt)
# 这里我们先生成每个 Outcome 的汇总，最后 Master 脚本会再处理 (或者直接在这里生成最终文件)
# 为了符合您 "生成各自的总表" 的需求，我们在这里把该 Exposure-Outcome 对的结果合并
TEMP_RES_DIR="${WORK_DIR}/results"
MERGED_FILENAME="${EXPOSURE_NAME}_${OUTCOME_NAME}_Phase3_Validation_Result.txt"
MERGED_PATH="${WORK_DIR}/${MERGED_FILENAME}"

if ls ${TEMP_RES_DIR}/results_val_chunk_*.txt 1> /dev/null 2>&1; then
    # 取表头 (第一个文件)
    ls ${TEMP_RES_DIR}/results_val_chunk_*.txt | head -n 1 | xargs head -n 1 > "$MERGED_PATH"
    # 合并内容 (跳过表头)
    cat ${TEMP_RES_DIR}/results_val_chunk_*.txt | grep -v "gene_symbol" >> "$MERGED_PATH"
    
    echo "    Merged to: $MERGED_PATH"
else
    echo "    Warning: No results generated for $EXPOSURE_NAME vs $OUTCOME_NAME"
fi

# 清理临时 Log (可选)
# rm -rf "$LOG_DIR"

