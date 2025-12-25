#!/bin/bash
#bash 00_master_runner.sh 2>&1 | tee master_runner.log

# ================= 全局配置 =================
# 1. 定义所有暴露 (Exposures)
# 可以是数组，也可以从文本文件读取
EXPOSURE_LIST=(
    "CD4_Dynamic"
    "onek1k"
    "1mscblood"
    "dice"
)

# 2. 定义所有结局 (Outcomes)
OUTCOME_LIST=(
    "NECROSIS_PERIAPICAL"
    "PULP_PERIAPICAL"
    "PULPITIS"
)

# 3. 高效利用 CPU 的关键设置
# 根据你的服务器核心数设置。例如：如果你有 64 核，设为 50-60。
# 这会覆盖子脚本里的 MAX_JOBS=20。
export GLOBAL_MAX_JOBS=50 

# ===========================================

# 获取开始时间
START_TIME=$(date)
echo "========================================================"
echo "Master Pipeline Started at: $START_TIME"
echo "Exposures: ${#EXPOSURE_LIST[@]}"
echo "Outcomes: ${#OUTCOME_LIST[@]}"
echo "Max Parallel Jobs per Pair: $GLOBAL_MAX_JOBS"
echo "========================================================"

# 双层循环遍历所有组合
for exp in "${EXPOSURE_LIST[@]}"; do
    for out in "${OUTCOME_LIST[@]}"; do
        
        echo "" 
        echo "--------------------------------------------------------"
        echo ">>> Starting Pipeline for: [$exp] vs [$out]"
        echo "--------------------------------------------------------"

        # === 步骤 1: Discovery (这步最耗时，会利用多核) ===
        echo "[Master] Launching Phase 1..."
        bash 01_run_phase1_discovery.sh "$exp" "$out"
        
        # 检查 01 的退出状态 (可选，取决于你的脚本是否正确设置 exit code)
        if [ $? -ne 0 ]; then
            echo "[Master] Error in Phase 1 for $exp - $out. Skipping..."
            continue
        fi

        # === 步骤 2: Process Results ===
        echo "[Master] Launching Data Processing..."
        bash 02_process_and_prep_phase2.sh "$exp" "$out"

        # === 步骤 3: Replication (利用多核) ===
        echo "[Master] Launching Phase 2..."
        bash 03_run_phase2_replication.sh "$exp" "$out"

        echo ">>> Finished Pipeline for: [$exp] vs [$out]"
        echo "--------------------------------------------------------"
        
        # 可选：每做完一对，等待几秒让系统 I/O 缓冲冷却
        sleep 5 
    done
done

END_TIME=$(date)
echo "========================================================"
echo "All Tasks Completed!"
echo "Start: $START_TIME"
echo "End:   $END_TIME"
echo "========================================================"
