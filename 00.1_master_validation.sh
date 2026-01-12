#!/bin/bash
# 00.1_master_validation.sh
# 验证阶段主控脚本

# ==============================================================================
# 1. 配置与环境
# ==============================================================================

# 1. 引入配置文件 (包含 EXPOSURE_VAL_LIST 和 OUTCOME_LIST)
if [ -f "./aa_discovery_config.sh" ]; then
    source ./aa_discovery_config.sh
else
    echo "Error: aa_discovery_config.sh not found!"
    exit 1
fi

# 2. 设置并发
# export GLOBAL_MAX_JOBS=30  #aa_discovery_config.sh 已经设置
RESULT_DIR="result"
mkdir -p "$RESULT_DIR"

# ==============================================================================
# 3. 初始化清理
# ==============================================================================
# 3. [清理] 删除 result 文件夹中旧的验证汇总文件，防止追加重复
echo ">>> Cleaning old validation results in $RESULT_DIR..."
for exp in "${EXPOSURE_VAL_LIST[@]}"; do
    rm -f "$RESULT_DIR/${exp}_Phase3_Validation_All_Results.txt"
done

START_TIME=$(date)
echo ">>> Validation Pipeline Started at: $START_TIME"

# ==============================================================================
# 4. 主循环
# ==============================================================================
for exp in "${EXPOSURE_VAL_LIST[@]}"; do
    
    # 定义该 Exposure 的最终大表路径
    FINAL_EXP_FILE="$RESULT_DIR/${exp}_Phase3_Validation_All_Results.txt"
    
    for out in "${OUTCOME_LIST[@]}"; do
        echo "--------------------------------------------------------"
        echo ">>> Processing: $exp vs $out"
        
        # 运行并行验证脚本
        bash 05_run_validation.sh "$exp" "$out"
        
        # --- 实时合并到 result 文件夹 ---
        # 读取刚刚生成的单个 Outcome 结果
        CURRENT_RES="${exp}_${out}/${exp}_${out}_Phase3_Validation_Result.txt"
        
        if [ -f "$CURRENT_RES" ]; then
            if [ ! -f "$FINAL_EXP_FILE" ]; then
                # 如果大表不存在，直接复制（包含表头）
                cat "$CURRENT_RES" > "$FINAL_EXP_FILE"
                echo "    [Create] Created summary file: $FINAL_EXP_FILE"
            else
                # 如果大表存在，跳过表头追加
                tail -n +2 "$CURRENT_RES" >> "$FINAL_EXP_FILE"
                echo "    [Append] Added results to summary."
            fi
        fi
        
    done
    
        # --- 统计展示 (可选) ---
    if [ -f "$FINAL_EXP_FILE" ]; then
        echo ">>> Finished $exp. Cell Type Stats:"
        # 简单的统计一下 cell_type (假设在第二列)
        # 如果有 act_time, cell_type 也是第二列; 如果没有, 也是第二列 (R脚本保证了顺序)
        awk -F"\t" 'NR>1 {count[$2]++} END {for (c in count) print "   " c ": " count[c]}' "$FINAL_EXP_FILE" | sort
    fi
    
done

echo "========================================================"
echo "All Validation Tasks Completed."
