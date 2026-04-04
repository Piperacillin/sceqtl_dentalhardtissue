#!/bin/bash
# 00.2_master_phase4_bulk.sh
# Phase 4 (Bulk eQTL) 验证主控脚本

# 1. 引入配置
if [ -f "./aa_discovery_config.sh" ]; then
    source ./aa_discovery_config.sh
else
    echo "Error: aa_discovery_config.sh not found!"
    exit 1
fi

RESULT_DIR="result"
mkdir -p "$RESULT_DIR"

# 2. 初始化清理
echo ">>> [Init] Cleaning old Phase 4 results..."
for exp in "${EXPOSURE_PHASE4_LIST[@]}"; do
    rm -f "$RESULT_DIR/${exp}_Phase4_Bulk_All_Results.txt"
done

START_TIME=$(date)
echo ">>> Phase 4 (Bulk) Pipeline Started at: $START_TIME"

# 3. 主循环
for exp in "${EXPOSURE_PHASE4_LIST[@]}"; do
    
    FINAL_EXP_FILE="$RESULT_DIR/${exp}_Phase4_Bulk_All_Results.txt"
    
    for out in "${OUTCOME_LIST[@]}"; do
        echo "--------------------------------------------------------"
        echo ">>> Processing: $exp vs $out"
        
        # [修改点] 调用 Phase 4 专用 Runner         
        bash 06_run_phase4_bulk.sh "$exp" "$out"
        
        # 实时合并
        CURRENT_RES="${exp}_${out}/${exp}_${out}_Phase4_Bulk_Result.txt"
        
        if [ -f "$CURRENT_RES" ]; then
            if [ ! -f "$FINAL_EXP_FILE" ]; then
                cat "$CURRENT_RES" > "$FINAL_EXP_FILE"
                echo "    [Create] Summary file created."
            else
                tail -n +2 "$CURRENT_RES" >> "$FINAL_EXP_FILE"
                echo "    [Append] Results added."
            fi
        fi
    done
    
    # 简单统计
    if [ -f "$FINAL_EXP_FILE" ]; then
        echo ">>> Finished $exp. Stats:"
        awk -F"\t" 'NR>1 {count[$2]++} END {for (c in count) print "   " c ": " count[c]}' "$FINAL_EXP_FILE" | sort
    fi
done

echo "========================================================"
echo "All Phase 4 Tasks Completed."
