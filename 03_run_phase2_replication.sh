#!/bin/bash
# 03_run_phase2_replication.sh
# ================= 接收参数区域 =================
EXPOSURE_NAME="${1:-CD4_Dynamic_demo}"
OUTCOME_NAME="${2:-NECROSIS_PERIAPICAL}"
# 允许从外部环境变量覆盖最大并发数
MAX_JOBS="${GLOBAL_MAX_JOBS:-20}"
# ==============================================
# 原有配置
R2_VAL="0.01_phase2"
CHUNK_SIZE=50
R_SCRIPT="02analysispipline.R"
LOG_DIR="./logs/phase2_${EXPOSURE_NAME}_${OUTCOME_NAME}"
# ===========================================

INPUT_FILE="exposure/${EXPOSURE_NAME}/${EXPOSURE_NAME}_r2_${R2_VAL}.qs"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Phase 2 input file not found. Assuming no significant genes in Phase 1."
    exit 0
fi

mkdir -p $LOG_DIR

echo "[Phase 2] Starting Replication Analysis (R2=0.01 subset)..."

# 1. 获取 Phase 2 基因数（增强错误检测）
EXP_COUNT=$(Rscript -e "cat(length(unique(qs::qread('${INPUT_FILE}')[['exposure']])))")

if ! [[ "$EXP_COUNT" =~ ^[0-9]+$ ]]; then
    echo "Error: Failed to get valid exposure count"
    exit 1
fi
echo "Replicating $EXP_COUNT genes."

# 2. 简单循环
for (( i=1; i<=EXP_COUNT; i+=CHUNK_SIZE )); do
    START=$i
    END=$((i + CHUNK_SIZE - 1))
    if [ $END -gt $EXP_COUNT ]; then END=$EXP_COUNT; fi

    echo "[$(date '+%H:%M:%S')] Submitting batch: ${START}-${END}"

    # 使用 env 传递环境变量
    nohup env \
        EXPOSURE_DATA=$EXPOSURE_NAME \
        OUTCOME=$OUTCOME_NAME \
        START=$START \
        END=$END \
        R2_THRESH=$R2_VAL \
        ONLY_PRIMARY="TRUE" \
        Rscript --vanilla $R_SCRIPT > "${LOG_DIR}/job_${START}_${END}.log" 2>&1 &

    # 进程控制：动态监控
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 2
    done
done

wait
echo "[Phase 2] Analysis complete. Merging results..."


# 3. Phase 2 结果合并与 Bonferroni 矫正 (校正基数为 Phase 2 的总数)
# 我们需要 Phase 2 的总基因数来进行loose校正
ORIGINAL_COUNT=$(Rscript -e "cat(length(unique(qs::qread('${INPUT_FILE}')[['exposure']])))")

PHASE2_RES_DIR="${EXPOSURE_NAME}_${OUTCOME_NAME}/results"
FINAL_FILE="${EXPOSURE_NAME}_${OUTCOME_NAME}/Final_Replication_Result.txt"

# 合并 Phase 2 结果
ls ${PHASE2_RES_DIR}/results_primary_r2_0.01_phase2_*.txt | head -n 1 | xargs head -n 1 > $FINAL_FILE
cat ${PHASE2_RES_DIR}/results_primary_r2_0.01_phase2_*.txt | grep -v "exposure" >> $FINAL_FILE

# 计算 Bonferroni
Rscript -e "
res <- read.table('${FINAL_FILE}', header=T, sep='\t')
N_total <- ${ORIGINAL_COUNT}

# 筛选 IVW/Wald
res_main <- subset(res, method %in% c('Inverse variance weighted', 'Wald ratio', 'IVW'))

# Bonferroni 阈值
bonf_thresh <- 0.05 / N_total

res_main\$bonferroni_pass <- res_main\$p < bonf_thresh
res_main\$bonferroni_threshold <- bonf_thresh

# 保存最终报告
write.csv(res_main, './${EXPOSURE_NAME}_${OUTCOME_NAME}/${EXPOSURE_NAME}_${OUTCOME_NAME}_Phase2_Final_Report.csv', row.names=F)
cat(paste0('Bonferroni Threshold: ', bonf_thresh, '\n'))
cat(paste0('Genes passed: ', sum(res_main\$bonferroni_pass), '\n'))
"

echo "[Success] Pipeline Finished. Check ${EXPOSURE_NAME}_${OUTCOME_NAME}_Phase2_Final_Report.csv"
