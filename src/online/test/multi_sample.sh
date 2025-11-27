#!/bin/bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "用法: $0 <sample_dir> <reference.fa> <output_dir> <work_dir>" >&2
    exit 1
fi

SAMPLE_DIR="$1"
REFERENCE="$2"
OUTPUT_DIR="$3"
WORK_DIR="$4"

# 创建输出和工作目录
mkdir -p "$OUTPUT_DIR" "$WORK_DIR"

echo "[`date`] 样本目录:   $SAMPLE_DIR"
echo "[`date`] 参考基因组: $REFERENCE"
echo "[`date`] 结果目录:   $OUTPUT_DIR"
echo "[`date`] 工作目录:   $WORK_DIR"

start_time=$(date +%s)

# 找到所有 BAM 文件（如果你只想要 *.sorted.bam，可以改成 *.sorted.bam）
shopt -s nullglob
bam_files=( "$SAMPLE_DIR"/*.bam )
shopt -u nullglob

if [ "${#bam_files[@]}" -eq 0 ]; then
    echo "错误: 在 $SAMPLE_DIR 中没有找到任何 .bam 文件" >&2
    exit 1
fi

# 逐个样本运行 cuteSV
for bam in "${bam_files[@]}"; do
    sample_name=$(basename "$bam")
    sample_name="${sample_name%.bam}"      # 去掉后缀 .bam

    sample_work_dir="${WORK_DIR}/${sample_name}"
    mkdir -p "$sample_work_dir"

    out_vcf="${OUTPUT_DIR}/${sample_name}.vcf"

    echo "[`date`] 开始处理样本: $sample_name"
    echo "         输入 BAM: $bam"
    echo "         输出 VCF: $out_vcf"
    echo "         工作目录: $sample_work_dir"

    # 核心命令：结构变异检测
    cuteSV "$bam" "$REFERENCE" "$out_vcf" "$sample_work_dir"

    echo "[`date`] 样本 $sample_name 处理完成"
done

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo "--------------------------------------------"
echo "总耗时: ${elapsed} 秒 (~$(printf '%.2f' "$(echo "$elapsed/60" | bc -l)") 分钟)"

if [ "$elapsed" -le 3600 ]; then
    echo "程序在 1 小时内结束。"
else
    echo "程序运行超过 1 小时。"
fi
echo "--------------------------------------------"

# 检查是否生成了每个样本对应的 VCF 文件
missing=0
for bam in "${bam_files[@]}"; do
    sample_name=$(basename "$bam")
    sample_name="${sample_name%.bam}"
    out_vcf="${OUTPUT_DIR}/${sample_name}.vcf"

    if [ -s "$out_vcf" ]; then
        echo "[OK] 找到 VCF: $out_vcf"
    else
        echo "[缺失] 未找到或文件为空: $out_vcf"
        missing=1
    fi
done

if [ "$missing" -eq 0 ]; then
    echo "所有样本的结构变异检测结果（VCF）均已在 $OUTPUT_DIR 中生成。"
else
    echo "存在缺失或为空的 VCF 文件，请检查日志和输入数据。"
    exit 1
fi
