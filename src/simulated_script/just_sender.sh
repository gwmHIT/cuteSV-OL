#!/bin/bash
src_dir=$1
tgt_dir=$2
send_rate=$3

mkdir -p "$tgt_dir"
split_files=($src_dir/*.fq.gz)

# 模拟“发送”过程，将拆分的文件按照速率移动到指定目录
for file in "${split_files[@]}"; do
    # 模拟文件发送，可以使用 mv 或 cp 将文件移到目标文件夹
    cp "$file" "$tgt_dir"
    current_time=$(date "+%Y-%m-%d %H:%M:%S")
    echo "file $file done time：$current_time"
    
    # 等待指定的发送速率（控制速率）
    sleep $send_rate
done
echo "ALL DONE"