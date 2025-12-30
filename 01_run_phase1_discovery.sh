#!/bin/bash
# 01_run_phase1_discovery.sh

# ================= 接收参数区域 =================
EXPOSURE_NAME="${1:-CD4_Dynamic_demo}" 
OUTCOME_NAME="${2:-NECROSIS_PERIAPICAL}"
# 这是硬上限，即使内存足够，也不会超过这个并发数
MAX_JOBS="${GLOBAL_MAX_JOBS:-50}" 

# ================= 内部配置 =================
R2_VAL="0.2"
CHUNK_SIZE=100
R_SCRIPT="02analysispipline.R"
# 设置安全内存缓冲：如果系统空闲内存少于 25GB，不启动新任务
# 假设大文件(OneK1K)最坏情况单任务能跑 8GB (读取+PLINK)
# 为了安全，我们要求剩余内存必须 大于 (单任务峰值 + 15GB缓冲)
# 这样即使当前任务还没到峰值，预留的空间也足够它膨胀
MIN_MEM_FREE="25G"

# [关键修改 2] 启动延时
# 每启动一个新任务，强制等待 30秒。
# 作用：让先启动的任务先跑一会，尽早把内存占住，
# 给 parallel 反射出真实的内存压力，避免瞬间并发造成的误判。
START_DELAY="15s" 

# 目录准备
LOG_DIR="./logs/phase1_${EXPOSURE_NAME}_${OUTCOME_NAME}"
mkdir -p "$LOG_DIR"

echo "[Phase 1] Init: ${EXPOSURE_NAME} -> ${OUTCOME_NAME}"

# 1. 准备文件路径
R2_FILE="${R2_VAL//./}" 
QS_FILE="exposure/${EXPOSURE_NAME}/${EXPOSURE_NAME}_r2_${R2_FILE}.qs"

if [ ! -f "$QS_FILE" ]; then
    ALT_FILE="exposure/${EXPOSURE_NAME}/${EXPOSURE_NAME}_r2_${R2_VAL}.qs"
    if [ -f "$ALT_FILE" ]; then
        QS_FILE=$ALT_FILE
    else
        echo "Error: Exposure file not found for $EXPOSURE_NAME" && exit 1
    fi
fi

# 2. 获取暴露总数
TOTAL_EXP=$(Rscript -e "cat(length(unique(qs::qread('${QS_FILE}')[['exposure']])))")

if ! [[ "$TOTAL_EXP" =~ ^[0-9]+$ ]]; then
    echo "Error: Failed to get exposure count" && exit 1
fi
echo "Total exposures: $TOTAL_EXP | Chunk size: $CHUNK_SIZE"
echo "Strategy: Parallel execution with --memfree ${MIN_MEM_FREE} (Max jobs: ${MAX_JOBS})"

# ================= 定义单个任务的包装函数 =================
# 这个函数会被 parallel 调用。它计算结束点并运行 R
run_chunk_task() {
    local start_idx=$1
    # 接收外部传入的变量
    local chunk_size=$2
    local total_exp=$3
    local exp_name=$4
    local out_name=$5
    local r_script=$6
    local r2_val=$7
    local log_dir=$8

    # 计算 End Index
    local end_idx=$((start_idx + chunk_size - 1))
    if [ "$end_idx" -gt "$total_exp" ]; then end_idx=$total_exp; fi

    # 运行 R 脚本 (保持环境变量传递方式)
    # 注意：这里不需要 nohup 和 &，因为 parallel 会管理进程
    env \
        EXPOSURE_DATA="$exp_name" \
        OUTCOME="$out_name" \
        START="$start_idx" \
        END="$end_idx" \
        R2_THRESH="$r2_val" \
        ONLY_PRIMARY="FALSE" \
        Rscript --vanilla "$r_script" > "${log_dir}/job_${start_idx}_${end_idx}.log" 2>&1
    
    # 返回状态码 (parallel 会利用这个判断成功与否)
    return $?
}

# 必须导出函数和变量，以便 parallel 的子shell能看见
export -f run_chunk_task

# ================= 使用 GNU Parallel 执行 =================
# seq 生成从 1 到 TOTAL，步长为 CHUNK_SIZE 的序列 (1, 101, 201...)
# {} 会被 seq 的输出(start_idx)替换
echo "Strategy: Parallel with '--delay ${START_DELAY}' and '--memfree ${MIN_MEM_FREE}'"
echo "This prevents memory spikes (like PLINK) from overlapping dangerously."
seq 1 "$CHUNK_SIZE" "$TOTAL_EXP" | parallel \
    --jobs "$MAX_JOBS" \
    --memfree "$MIN_MEM_FREE" \
    --delay "$START_DELAY" \
    --line-buffer \
    --tagstring "[Batch {}]" \
    run_chunk_task {} "$CHUNK_SIZE" "$TOTAL_EXP" "$EXPOSURE_NAME" "$OUTCOME_NAME" "$R_SCRIPT" "$R2_VAL" "$LOG_DIR"

# ========================================================
echo "[Phase 1] All jobs completed at $(date)"

# 错误检查
ERR_COUNT=$(grep -l "Error\|ERROR" "${LOG_DIR}"/*.log 2>/dev/null | wc -l)
if [ "$ERR_COUNT" -gt 0 ]; then
    echo "WARNING: Found $ERR_COUNT log files with errors. Check $LOG_DIR"
    exit 1 # 返回非零状态码给Master
else
    echo "Success: No errors detected in logs."
    exit 0
fi
