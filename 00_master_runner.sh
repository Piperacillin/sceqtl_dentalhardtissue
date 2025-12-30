#!/bin/bash
#bash 00_master_runner.sh 2>&1 | tee master_runner.log

# ================= 配置 =================
# 1. 引入配置文件 (核心修改)
source ./aa_discovery_config.sh


START_TIME=$(date)
echo "Master Pipeline Started at: $START_TIME"

PHASE1_ERR_COUNT=0
for exp in "${EXPOSURE_LIST[@]}"; do
    for out in "${OUTCOME_LIST[@]}"; do
        
        echo "--------------------------------------------------------"
        echo ">>> [Master] Processing: $exp vs $out"
        
        # === Phase 1: Discovery (由 Parallel 智能控制内存) ===
        bash 01_run_phase1_discovery.sh "$exp" "$out"
        EXIT_CODE=$?  # 立即捕获退出码
        
        # 2. 修改判断逻辑：仅计数和报讯，但不跳过 (Remove continue)
        if [ $EXIT_CODE -ne 0 ]; then
            ((PHASE1_ERR_COUNT++))
            echo "!!! [Master] Warning: Phase 1 returned error code $EXIT_CODE."
            echo "!!! [Master] Total Phase 1 Errors so far: $PHASE1_ERR_COUNT"
            echo "!!! [Master] Proceeding to Phase 2/3..."
        fi

        # === Phase 2 & 3 ===
        # 这些步骤通常没那么吃内存，或者可以复用同样的 Parallel 逻辑
        echo ">>> [Master] Processing Phase 2..."
        bash 02_process_and_prep_phase2.sh "$exp" "$out"

        echo ">>> [Master] Running Phase 2 Replication..."
        bash 03_run_phase2_replication.sh "$exp" "$out"

        sleep 3
    done
done

echo "========================================================"
echo "All Tasks Completed. Duration: $START_TIME - $(date)"

