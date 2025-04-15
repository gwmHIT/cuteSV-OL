def check_vcf_format(vcf_path):
    with open(vcf_path, 'r') as f:
        for i, line in enumerate(f, 1):
            if line.startswith('#'):
                continue  # 跳过注释头部
            fields = line.strip().split('\t')
            if len(fields) < 8:
                print(f"[行 {i}] 字段数不足：{line.strip()}")
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = fields[:8]

            # 检查 pos 是否为整数
            if pos == '':
                print(f"[行 {i}] pos 为空：{line.strip()}")
                continue
            try:
                int(pos)
            except ValueError:
                print(f"[行 {i}] pos 非法：{pos}")
                continue

            # 从 info 字段中提取 SVLEN
            svlen = None
            for item in info.split(';'):
                if item.startswith('SVLEN='):
                    svlen_str = item.replace('SVLEN=', '')
                    if svlen_str == '':
                        print(f"[行 {i}] SVLEN 为空：{line.strip()}")
                        break
                    try:
                        svlen = int(svlen_str)
                    except ValueError:
                        print(f"[行 {i}] SVLEN 非法：{svlen_str}")
                    break
            else:
                print(f"[行 {i}] 没有 SVLEN 字段：{line.strip()}")


# 使用方法
check_vcf_format("/home/user/guoweimin/data/cutesv_ol/experiment_last/hg38_detect_rate/work_dir/e2t1_output/41.6_output.vcf")
