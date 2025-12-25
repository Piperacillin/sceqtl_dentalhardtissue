#!/bin/bash

# 使用 pgrep 获取进程 ID
pids=$(pgrep -f '/usr/lib/R/bin/exec/R')

# 检查是否有任何匹配的进程
if [ -z "$pids" ]; then
    echo "No matching processes found."
else
    # 逐个终止每个进程
    for pid in $pids; do
        echo "Terminating process ID $pid"
        kill -9 $pid
    done
    echo "All matching processes have been terminated."
fi
