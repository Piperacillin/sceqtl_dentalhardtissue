#!/bin/bash
# 04_merge_results_for_gen_cd4dynamic_dice_validation.sh


# 1. 引入配置文件 (动态读取 00 脚本同款配置)
if [ -f "./aa_discovery_config.sh" ]; then
    source ./aa_discovery_config.sh
else
    echo "Error: aa_discovery_config.sh not found!"
    exit 1
fi

RESULT_DIR="result"

# ==============================================================================
# 2. [新增] 全局清理模块 (Global Cleanup)
#    在开始任何合并前，强制删除旧的合成文件，防止内容重复堆积
# ==============================================================================
echo "========================================================"
echo ">>> [Init] Global Cleanup: Removing old merged files..."

# 2.1 清理 Shell 合并生成的中间文件
for exp in "${EXPOSURE_LIST[@]}"; do
    # 定义要清理的文件名
    FILE_1="$RESULT_DIR/${exp}_Phase1_All_Results.txt"
    FILE_2="$RESULT_DIR/${exp}_Phase1_Egger_Cochransq_Results.txt"
    FILE_3="$RESULT_DIR/${exp}_Phase1_sensitivity_test_Results.txt"
    FILE_4="$RESULT_DIR/${exp}_Phase2_Final_Report.csv"

    # 强制删除
    if [ -f "$FILE_1" ] || [ -f "$FILE_2" ] || [ -f "$FILE_3" ] || [ -f "$FILE_4" ]; then
        echo "  [Clean] Removing old files for exposure: $exp"
        rm -f "$FILE_1" "$FILE_2" "$FILE_3" "$FILE_4"
    fi
done

# # 2.2 清理 R 脚本生成的最终结果文件 (防止 R 脚本 append 模式导致的问题)
# R_OUT_1="$RESULT_DIR/1mscblood_CD4T_Dynamic_Intersect.csv"
# R_OUT_2="$RESULT_DIR/onek1k_DICE_Validation_Subset.csv"
# 
# if [ -f "$R_OUT_1" ] || [ -f "$R_OUT_2" ]; then
#     echo "  [Clean] Removing old R script output files..."
#     rm -f "$R_OUT_1" "$R_OUT_2"
# fi

echo ">>> Cleanup Completed. Starting fresh merge..."
echo "========================================================"



# ================= 函数定义 =================

# 1. 合并文件内容 (保持不变)
merge_file_content() {
    local src="$1"
    local tgt="$2"

    if [ -f "$src" ]; then
        if [ ! -f "$tgt" ]; then
            cat "$src" > "$tgt"
            echo "  [Create] Created $tgt from $src"
        else
            tail -n +2 "$src" >> "$tgt"
            echo "  [Append] Merged content from $src"
        fi
    else
        : # Silent skip
    fi
}

# 2. [核心修改] 后处理：智能拆分 exposure 列 (静态 vs 非静态)
process_split_columns() {
    local tgt="$1"
    if [ ! -f "$tgt" ]; then return; fi

    echo "  [Process] Intelligent splitting for: $(basename "$tgt")"

    local fs_val="\t"
    if [[ "$tgt" == *.csv ]]; then fs_val=","; fi

    awk -v FS="$fs_val" -v OFS="$fs_val" '
    BEGIN { col_idx=0 }
    
    # --- 处理表头 ---
    NR==1 {
        for(i=1; i<=NF; i++) {
            clean_header = $i
            gsub(/"/, "", clean_header)
            if(clean_header == "exposure") {
                col_idx = i
                break
            }
        }
        # 新增三列: gene_symbol, cell_type, act_time
        print $0, "gene_symbol", "cell_type", "act_time"
    }
    
    # --- 处理数据行 ---
    NR>1 {
        if(col_idx > 0) {
            val = $col_idx
            gsub(/"/, "", val) # 去除引号
            
            # 1. 寻找第一个下划线的位置
            idx1 = index(val, "_")
            
            if(idx1 > 0) {
                # [Gene] 是第一个下划线之前的部分
                gene = substr(val, 1, idx1 - 1)
                
                # 获取第一个下划线之后的所有内容 (temp_rest)
                temp_rest = substr(val, idx1 + 1)
                
                # 2. 在剩余部分中寻找第二个下划线
                idx2 = index(temp_rest, "_")
                
                if(idx2 > 0) {
                    # === 非静态 (Dynamic): Gene_Cell_Time ===
                    # [Cell] 是第一个和第二个下划线中间的部分
                    cell = substr(temp_rest, 1, idx2 - 1)
                    # [Time] 是第二个下划线之后的部分
                    time_val = substr(temp_rest, idx2 + 1)
                } else {
                    # === 静态 (Static): Gene_Cell ===
                    # [Cell] 是第一个下划线之后的所有内容
                    cell = temp_rest
                    # [Time] 占位符
                    time_val = "None"
                }
            } else {
                # 异常情况: 没有下划线
                gene = val
                cell = "NA"
                time_val = "NA"
            }
            
            # 输出原行 + 3个新列
            print $0, gene, cell, time_val
        } else {
            # 没找到 exposure 列
            print $0, "NA", "NA", "NA"
        }
    }
    ' "$tgt" > "${tgt}.tmp" && mv "${tgt}.tmp" "$tgt"
}

# 3. 统计 Cell Type (简单更新，不影响主逻辑)
print_cell_type_stats() {
    local tgt="$1"
    if [ ! -f "$tgt" ]; then return; fi

    echo "  ------------------------------------------------"
    echo "  [Stats] Cell Type Distribution for: $(basename "$tgt")"
    
    local fs_val="\t"
    if [[ "$tgt" == *.csv ]]; then fs_val=","; fi

    # 统计倒数第二列 (因为倒数第一列现在是 act_time)
    # gene_symbol(NF-2), cell_type(NF-1), act_time(NF)
    awk -v FS="$fs_val" '
    NR>1 {
        cell_val = $(NF-1) 
        gsub(/"/, "", cell_val)
        count[cell_val]++ 
    }
    END {
        printf "    %-30s %s\n", "CELL_TYPE", "COUNT"
        print "    ----------------------------------------"
        for (c in count) {
            printf "    %-30s %d\n", c, count[c]
        }
    }
    ' "$tgt" | sort
    echo ""
}

# ================= 主逻辑 =================

echo ">>> Starting Intelligent Merge Pipeline..."
mkdir -p "$RESULT_DIR"

for exp in "${EXPOSURE_LIST[@]}"; do
    echo "========================================================"
    echo "Processing Exposure Dataset: $exp"
    
    TGT_P1_ALL="$RESULT_DIR/${exp}_Phase1_All_Results.txt"
    TGT_P1_HET="$RESULT_DIR/${exp}_Phase1_Egger_Cochransq_Results.txt"
    TGT_P1_SEN="$RESULT_DIR/${exp}_Phase1_sensitivity_test_Results.txt"
    TGT_P2_FIN="$RESULT_DIR/${exp}_Phase2_Final_Report.csv"

    for out in "${OUTCOME_LIST[@]}"; do
        CURRENT_DIR="${exp}_${out}"
        if [ -d "$CURRENT_DIR" ]; then
            merge_file_content "$CURRENT_DIR/Phase1_All_Results.txt" "$TGT_P1_ALL"
            merge_file_content "$CURRENT_DIR/Phase1_Egger_Cochransq_Results.txt" "$TGT_P1_HET"
            merge_file_content "$CURRENT_DIR/Phase1_sensitivity test_Results.txt" "$TGT_P1_SEN"
            
            SRC_P2="$CURRENT_DIR/${exp}_${out}_Phase2_Final_Report.csv"
            merge_file_content "$SRC_P2" "$TGT_P2_FIN"
        fi
    done

    echo "  >>> Post-processing: Intelligent splitting..."
    process_split_columns "$TGT_P1_ALL"
    process_split_columns "$TGT_P1_HET"
    process_split_columns "$TGT_P1_SEN"
    process_split_columns "$TGT_P2_FIN"
    
    if [ -f "$TGT_P1_ALL" ]; then
        print_cell_type_stats "$TGT_P1_ALL"
    fi
    
    echo "Dataset $exp completed."
done

# --- 阶段 2: [新增] 调用 R 脚本进行 CD4 动态验证与 DICE 映射 ---
echo "========================================================"
echo ">>> [Extension] Running R script for CD4 Dynamic & DICE Validation..."

R_SCRIPT="03filter_and_mapping_cd4dynamic_dice_validation.R"

if [ -f "$R_SCRIPT" ]; then
    # 检查是否有 R 环境 (可选)
    if command -v Rscript >/dev/null 2>&1; then
        echo "  Found Rscript. Executing $R_SCRIPT..."
        
        # 执行 R 脚本
        Rscript "$R_SCRIPT"
        EXIT_CODE=$?
        
        if [ $EXIT_CODE -eq 0 ]; then
            echo "  [Success] R script executed successfully."
            echo "  Check '$RESULT_DIR' for new qs files."
        else
            echo "  [Error] R script failed with exit code $EXIT_CODE."
            # exit 1  # 如果希望 R 失败则整个脚本报错，取消此行注释
        fi
    else
        echo "  [Error] Rscript command not found! Please install R."
    fi
else
    echo "  [Warning] R script '$R_SCRIPT' not found in current directory."
    echo "  Skipping validation subset generation."
fi

echo "========================================================"
echo "All tasks completed."