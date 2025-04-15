#!/usr/bin/env python3
import sys
import bisect

def parse_info(info):
    """
    解析 VCF INFO 字段，将其转换为字典
    """
    info_dict = {}
    for field in info.split(';'):
        if '=' in field:
            key, value = field.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[field] = True
    return info_dict

def load_highfreq_file(highfreq_file, af_threshold, mode):
    """
    读取高频变异 VCF 文件，或者自定义SV的 TXT 文件

    返回字典： {染色体: [(POS, ID, SVTYPE, SVLEN), ...]}，列表按POS排序。
    """
    total_hight_variants = 0
    highfreq = {}

    with open(highfreq_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                    continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            id = parts[2]
            info_field = parts[7]
            info_dict = parse_info(info_field)
                
            svtype = info_dict.get('SVTYPE', '')
            try:
                svlen = int(info_dict.get('SVLEN', '0'))
            except ValueError:
                svlen = 0
            if mode == 2:
                af = info_dict.get('AF', '')
                if float(af) < af_threshold:
                    continue            
            total_hight_variants += 1
            if chrom not in highfreq:
                highfreq[chrom] = []
            if mode == 1:
                highfreq[chrom].append((int(pos), id, svtype, int(svlen)))
            else:
                highfreq[chrom].append((int(pos), id, svtype, int(svlen), af))
    # else:
    #     with open(highfreq_file, 'r') as f:
    #         for line in f:
    #             parts = line.strip().split('\t')
    #             chrom = parts[0]
    #             pos = parts[1]
    #             id = parts[2]
    #             svtype = parts[3]
    #             svlen = parts[4]
    #             if chrom not in highfreq:
    #                 highfreq[chrom] = []
    #             highfreq[chrom].append((int(pos), id, svtype, int(svlen)))
    #             total_hight_variants += 1
    # 对每个染色体的记录按位置排序
    for chrom in highfreq:
        highfreq[chrom].sort(key=lambda x: x[0])
    return total_hight_variants,highfreq

def compare_vcf_highfreq_mapping(vcf_file, highfreq_file, output_file, af_threshold, mode,
                                 pos_tolerance=1000, length_lower_ratio=0.9, length_upper_ratio=1.1):
    # 加载高频变异数据
    total_hight_variants,highfreq_data = load_highfreq_file(highfreq_file, af_threshold, mode)
    # 对于每个染色体构建位置列表（用于二分查找）
    highfreq_positions = {}
    for chrom, records in highfreq_data.items():
        positions = [rec[0] for rec in records]
        highfreq_positions[chrom] = positions
    
    high_set = set()
    with open(vcf_file, 'r') as fin, open(output_file, 'w') as fout:
        # 写入映射关系文件的表头
        fout.write("CHROM\tInput_VCF_ID\tInput_VCF_POS\tInput_VCF_SVLEN\tHighfreq_ID\tPOS\tSVLEN\tAF\n")
        for line in fin:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            # VCF 标准列：CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, ...
            if len(parts) < 8:
                continue
            chrom = parts[0]
            if not chrom.startswith("chr"):
                    chrom = 'chr' + chrom
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            input_id = parts[2]
            info_field = parts[7]
            info_dict = parse_info(info_field)
            svtype = info_dict.get('SVTYPE', '')
            try:
                svlen = int(info_dict.get('SVLEN', '0'))
            except ValueError:
                svlen = 0
            abs_svlen = abs(svlen)
            
            # 若当前染色体在高频数据中存在候选
            if chrom in highfreq_data:
                pos_list = highfreq_positions[chrom]
                records_list = highfreq_data[chrom]
                # 利用二分查找确定候选记录的起始位置
                left_index = bisect.bisect_left(pos_list, pos - pos_tolerance)
                for i in range(left_index, len(records_list)):
                    if mode == 2:
                        candidate_pos, candidate_id, candidate_svtype, candidate_svlen, af = records_list[i]
                    else:
                        candidate_pos, candidate_id, candidate_svtype, candidate_svlen = records_list[i]
                    if candidate_pos > pos + pos_tolerance:
                        break
                    if svtype != candidate_svtype:
                        continue
                    candidate_abs_svlen = abs(candidate_svlen)
                    lower_bound = length_lower_ratio * candidate_abs_svlen
                    upper_bound = length_upper_ratio * candidate_abs_svlen
                    if lower_bound <= abs_svlen <= upper_bound:
                        if candidate_id == ".":
                            candidate_id = f'{chrom}_{candidate_pos}_{candidate_svtype}_{candidate_svlen}'
                        high_set.add(candidate_id)
                        if mode == 2:
                            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, input_id, pos, svlen, candidate_id, candidate_pos, candidate_svlen, af))
                        else:
                            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, input_id, pos, svlen, candidate_id, candidate_pos, candidate_svlen))
                        # mapping_lines += 1
            # if candidate_found:
                # variants_with_match += 1

        detection_rate = (len(high_set) / total_hight_variants * 100) if total_hight_variants > 0 else 0
        print("检出数: {}, 总高频变异数: {}, 检出率: {:.2f}%\n".format(len(high_set), total_hight_variants, detection_rate))
        fout.write("检出数: {}, 总高频变异数: {}, 检出率: {:.2f}%\n".format(len(high_set), total_hight_variants, detection_rate))
    return detection_rate

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("用法: {} <输入_vcf文件> <高频变异_vcf文件> <输出_mapping_txt文件> <高频变异频率>".format(sys.argv[0]))
        sys.exit(1)
    input_vcf = sys.argv[1]
    highfreq_txt = sys.argv[2]
    output_mapping = sys.argv[3]
    af_threshold = float(sys.argv[4])
    compare_vcf_highfreq_mapping(input_vcf, highfreq_txt, output_mapping, af_threshold, 1)
