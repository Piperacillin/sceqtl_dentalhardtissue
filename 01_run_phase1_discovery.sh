#!/bin/bash
# 01_run_phase1_discovery.sh
# ================= 接收参数区域 =================
# 使用方式: bash 01_run_phase1_discovery.sh "暴露名" "结局名"
EXPOSURE_NAME="${1:-CD4_Dynamic_demo}" # 默认为原值，防止空参报错
OUTCOME_NAME="${2:-NECROSIS_PERIAPICAL}"
# 允许从外部环境变量覆盖最大并发数，默认20
MAX_JOBS="${GLOBAL_MAX_JOBS:-20}" 
# ==============================================
# 原有内部配置
R2_VAL="0.2"
CHUNK_SIZE=100
R_SCRIPT="02analysispipline.R"

# 目录准备
LOG_DIR="./logs/phase1_${EXPOSURE_NAME}_${OUTCOME_NAME}"
#ERROR_DIR="./errors/phase1_${EXPOSURE_NAME}_${OUTCOME_NAME}"  # ← 修复1
#mkdir -p $LOG_DIR $ERROR_DIR
mkdir -p $LOG_DIR

echo "[Phase 1] Starting Discovery Analysis (R2=${R2_VAL})..."

# 1. 准备文件路径（处理小数点）
R2_FILE="${R2_VAL//./}"  # 0.2 → 02  ← 修复2
QS_FILE="exposure/${EXPOSURE_NAME}/${EXPOSURE_NAME}_r2_${R2_FILE}.qs"

# 检查文件存在性
if [ ! -f "$QS_FILE" ]; then
    echo "Error: $QS_FILE not found"
    echo "Looking for alternative naming..."
    # 尝试替代命名
    ALT_FILE="exposure/${EXPOSURE_NAME}/${EXPOSURE_NAME}_r2_${R2_VAL}.qs"
    if [ -f "$ALT_FILE" ]; then
        QS_FILE=$ALT_FILE
        echo "Using: $ALT_FILE"
    else
        exit 1
    fi
fi

# 2. 获取暴露总数（增强错误检测）  ← 修复4
TOTAL_EXP=$(Rscript -e "cat(length(unique(qs::qread('${QS_FILE}')[['exposure']])))")

if ! [[ "$TOTAL_EXP" =~ ^[0-9]+$ ]]; then
    echo "Error: Failed to get valid exposure count"
    exit 1
fi
echo "Total exposures detected: $TOTAL_EXP"
echo "Running Phase 1 for: ${EXPOSURE_NAME} -> ${OUTCOME_NAME} with Max Jobs: ${MAX_JOBS}"

# 3. 循环提交任务
for (( i=1; i<=TOTAL_EXP; i+=CHUNK_SIZE )); do
    START=$i
    END=$((i + CHUNK_SIZE - 1))
    if [ $END -gt $TOTAL_EXP ]; then END=$TOTAL_EXP; fi
    
    echo "[$(date '+%H:%M:%S')] Submitting batch: ${START}-${END}"
    
    # 使用 env 传递环境变量（更明确）  ← 修复3
    nohup env \
        EXPOSURE_DATA=$EXPOSURE_NAME \
        OUTCOME=$OUTCOME_NAME \
        START=$START \
        END=$END \
        R2_THRESH=$R2_VAL \
        ONLY_PRIMARY="FALSE" \
        nohup Rscript --vanilla $R_SCRIPT > "${LOG_DIR}/job_${START}_${END}.log" 2>&1 &
    
    # 进程控制：动态监控  ← 修复5
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 2
    done
done

wait
echo "[Phase 1] All jobs completed at $(date)"

# 可选：检查失败任务
echo "Checking for errors..."
grep -l "Error\|ERROR" ${LOG_DIR}/*.log 2>/dev/null | wc -l
