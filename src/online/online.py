import argparse
import logging
from pathlib import Path
import os
import multiprocessing
import subprocess
import sys
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import json
import datetime
from online.compare_model import compare_vcf_highfreq_mapping
import glob
import gzip
import pkg_resources

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="cuteSV_ONLINE", 
		formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # **************Parameters of input******************
    parser.add_argument('fastq_dir', 
		type = str, 
		help = "The fastq folder monitored by cuteSV-OL.")
    parser.add_argument("reference",  
		type = str, 
		help ="The reference genome in fasta format.")
    parser.add_argument("work_dir", 
		type = str, 
		help ="Work diretory for cuteSV-OL.")
    parser.add_argument('output_vcf', 
		type = str, 
		help = "The vcf folder where cuteSV-OL outputs real-time test results to.")
    parser.add_argument("--threads", 
		type = int, 
        default = 4,
		help ="Number of threads to use.")
    parser.add_argument("--mmi_path", 
		type = str, 
		help ="Minimizer index for the reference in minimap2.",
        default='')
    parser.add_argument('--platform', 
		type = str, 
		help = argparse.SUPPRESS,
        default = "map-ont")
    parser.add_argument('--monitor_fade', 
		type = int, 
		help = "Monitor will close if no new files are detected after monitor_fade second.",
        default = 600)
    parser.add_argument('--high_freq_file', 
		type = str, 
		help = "high frequence SV file or user-defined recall set[vcf]",
        default = "")
    parser.add_argument('--target_set', 
		type = str, 
		help = "high frequence SV file or user-defined recall set[vcf]",
        default = "")
    parser.add_argument('--user_defined',
        action = 'store_true',
        help = 'The recall set[vcf] is user-defined')
    parser.add_argument('--sv_freq', 
		type = float, 
		help = "Target SV frequence for detection",
        default = 1.0)
    parser.add_argument('--pctsize', 
		type = float, 
		help = "Min pct allele size similarity",
        default = 0.9)
    parser.add_argument('--ref_dist', 
		type = float, 
		help = "Max reference location distance",
        default = 1000)
    parser.add_argument('--recall_file', 
		type = str, 
		help = "candidate SV mapping to the high frequence SV",
        default = "")
    parser.add_argument('--target_rate', 
		type = float, 
		help = "stop sequency if the detected rate is higher than target_rate",
        default = 100.0)
    parser.add_argument('--batch_interval', 
		type = int, 
		help = "Real-time results are generated every batch_interval batches",
        default = 4)
    args = parser.parse_args(argv)
    return args


def clean_fastq_inplace(filepath):
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        lines = f.readlines()

    cleaned = []
    i = 0
    while i < len(lines):
        if lines[i].startswith('@') and i + 3 < len(lines):
            header = lines[i].strip()
            seq = lines[i + 1].strip()
            plus = lines[i + 2].strip()
            qual = lines[i + 3].strip()
            if plus.startswith('+') and len(seq) == len(qual):
                cleaned.extend([header + '\n', seq + '\n', plus + '\n', qual + '\n'])
                i += 4
            else:
                i += 1  # skip to next potential read
        else:
            i += 1
    with gzip.open(filepath, 'wt', encoding='utf-8') as f:
        f.writelines(cleaned)


def handle_fault_one(fq_path, bam_path):
    path_1 = bam_path
    path_2 = f"{bam_path}.bai"
    path_3 = f"{bam_path}.coverage.chr.stat.gz"
    path_4 = f'{bam_path}.coverage.chr.stat'
    if os.path.exists(path_1) and os.path.isfile(path_1):
        os.remove(path_1)
    if os.path.exists(path_2) and os.path.isfile(path_2):
        os.remove(path_2)
    if os.path.exists(path_3) and os.path.isfile(path_3):
        os.remove(path_3)
    if os.path.exists(path_4) and os.path.isfile(path_4):
        os.remove(path_4)        
    clean_fastq_inplace(fq_path)


def handle_fault_two(signatures_file,bam_name):
    delete_signature_files(signatures_file,bam_name)


def handle_fault_three():
    return


def minimap2_step(work_dir, fq_path, mmi_path, platform, bam_path, thread, pandepth_path):
    while True:
        try:
            command = f'minimap2 -t {thread} -ax {platform} {mmi_path} {fq_path} | samtools sort -o {bam_path};'
            command_add = f'samtools index {bam_path};{pandepth_path} -i {bam_path} -o {bam_path}.coverage -t 16;gzip -d {bam_path}.coverage.chr.stat.gz'
            subprocess.run(command + command_add, shell=True, check=True)
            break
        
        except subprocess.CalledProcessError as e:
            logging.info(f"命令执行失败，错误码: {e.returncode}")
            handle_fault_one(fq_path,bam_path)


def generate_mmi(reference_path, mmi_path):
    command = f'minimap2 -d {mmi_path} {reference_path}'
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.info(f"命令执行失败，错误码: {e.returncode}")


def cutesv_extract_sigs(work_dir, fa_path, bam_path, bam_name, thread, task_dir):
    while True:
        command = f'cuteSV --input {bam_path} --reference {fa_path} --work_dir {work_dir} --bam_name {bam_name} --threads {thread} --mode 1'
        try:
            subprocess.run(command, shell=True, check=True)
            break
        except subprocess.CalledProcessError as e:
            signatures_file = f"{work_dir}/signatures"
            handle_fault_two(signatures_file,bam_name)
    with open(f'{bam_path}.coverage.chr.stat', 'r') as file:
        last_line = file.readlines()[-1]

    # 提取 MeanDepth
    mean_depth_str = last_line.split('\t')[-1]  # 通过制表符分割，获取最后一列
    mean_depth = float(mean_depth_str.split(': ')[-1])  # 去除 'MeanDepth: '，保留数值并转为浮点数
    with open(f'{task_dir}coverage_list.txt', 'a') as file:
        file.write(f"{mean_depth}\n")


def cutesv_combine_cluster(work_dir, vcf_output, thread, reference, high_freq_file, sv_freq, recall_file, user_defined, pctsize, ref_dist):
    file_path = f'{work_dir}coverage_list.txt'
    total_sum = 0.0
    with open(file_path, 'r') as file:
        for line in file:   
            try:
                total_sum += float(line.strip())  
            except ValueError:
                print(f"无法转换这一行: {line.strip()}")
    if total_sum <= 0.1: #深度太低直接返回
        if high_freq_file == "":
            return None
        else:
            return 0
    elif 1 <= total_sum <= 10:
        min_support = 2
    elif 10 < total_sum <= 15:
        min_support = 3
    elif 15 < total_sum <= 20:
        min_support = 4
    else:
        min_support = 5
    cutesv_work_dir = work_dir + "cutesv_work_dir"
    detect_rate = 0
    vcf_path = f'{vcf_output}{total_sum:.1f}_output.vcf'
    command = f'cuteSV --retain_work_dir --write_old_sigs --genotype --output {vcf_path} --reference {reference} --work_dir {cutesv_work_dir} --threads {thread} --min_support {min_support} --mode 2'
    while True:
        try:
            subprocess.run(command, shell=True, check=True)
            break
        except subprocess.CalledProcessError as e:
            handle_fault_three()
    if high_freq_file != "":
        length_lower_ratio = pctsize
        length_upper_ratio = 1 + (1 - pctsize)
        if user_defined is True:
            detect_rate = compare_vcf_highfreq_mapping(vcf_path,high_freq_file,recall_file,sv_freq,1, ref_dist, length_lower_ratio, length_upper_ratio)
        else:
            detect_rate = compare_vcf_highfreq_mapping(vcf_path,high_freq_file,recall_file,sv_freq,2, ref_dist, length_lower_ratio, length_upper_ratio)
        with open(f'{work_dir}depth_performance_rate.txt', 'a') as file:
            file.write(f"{total_sum},{detect_rate}\n")
        return detect_rate
    else:
        return None


def worker(task_queue, fa_path, work_dir, fq_dir, platform, mmi_path, per_thread, shared_value, 
           batch_interval, high_freq_file, output_vcf, user_defined, pctsize, ref_dist, sv_freq, 
           recall_file, target_rate, pandepth_path):
    counter = 0
    do_flag = False
    while True:
        task = task_queue.get()
        counter += 1
        if counter == batch_interval:
            counter = 0
            do_flag = True
        if task is None:
            # 接收到终止信号，退出
            break
        with open(f'{work_dir}debug.txt','w') as f:
            f.write(task)
        name, ext = os.path.splitext(task)
        fq_path = fq_dir + task
        bam_path = work_dir + "bam/" + name + ".bam"
        cutesv_work_dir = work_dir + "cutesv_work_dir/"
        minimap2_step(work_dir, fq_path, mmi_path, platform, bam_path, per_thread, pandepth_path)
        cutesv_extract_sigs(cutesv_work_dir, fa_path, bam_path, name, per_thread, work_dir)
        if do_flag == True:
            if high_freq_file != "":
                detect_rate = cutesv_combine_cluster(work_dir, output_vcf, per_thread, fa_path, high_freq_file, sv_freq, recall_file, user_defined, pctsize, ref_dist)
                if detect_rate >= target_rate:
                    shared_value.value = True  # 设置停止标志
            else:
                cutesv_combine_cluster(work_dir, output_vcf, per_thread, fa_path, high_freq_file, sv_freq, recall_file, user_defined, pctsize, ref_dist)
            do_flag = False
        result = f"{task}"
        with open(f"{work_dir}finished.txt", 'a', encoding='utf-8') as file:
            file.write(result + '\n')


def arrange_task(fastq_dir, output_txt, finished_path):
    dir_path = Path(fastq_dir)

    if not dir_path.exists():
        logging.info(f"Error:dir {fastq_dir} is not exist。")
        return

    if not dir_path.is_dir():
        logging.info(f"Error:dir {fastq_dir} is not a directory")
        return

    # 使用 glob 方法筛选 .fastq 文件，忽略大小写
    fastq_files = [f.name for f in dir_path.iterdir() if f.is_file() and (f.suffix.lower() == '.fastq' or f.suffix.lower() == '.fq' or f.suffix.lower() == '.gz')]
    unfinished_files = []
    if not fastq_files:
        return

    try:
        # 将文件名写入输出文本文件
        with open(output_txt, 'w') as fout, open(finished_path,'r') as fin:
            finish_names = [line.strip() for line in fin if line.strip()]
            for filename in fastq_files:
                if filename not in finish_names:
                    unfinished_files.append(filename)
                fout.write(filename + '\n')
        # logging.info(f"成功将 {len(fastq_files)} 个 .fastq 文件名写入 {output_txt}。")
    except IOError as e:
        logging.info(f"写入文件时发生错误：{e}")
    return unfinished_files


def is_file_complete(file_path, check_interval=2, retries=1):
    """通过检测文件大小是否稳定来判断文件是否写入完成。
    
    :param file_path: 文件路径
    :param check_interval: 每次检查之间的等待时间（秒）
    :param retries: 连续检测相同大小的次数，达到该次数则认为文件写入完成
    :return: True 表示文件稳定，否则 False
    """
    previous_size = -1
    stable_count = 0

    while stable_count < retries:
        try:
            current_size = os.path.getsize(file_path)
        except OSError:
            # 文件可能暂时不可访问
            current_size = -1

        if current_size == previous_size and current_size != -1:
            stable_count += 1
        else:
            stable_count = 0
            previous_size = current_size

        time.sleep(check_interval)
    return True


class FQFileHandler(FileSystemEventHandler):
    def __init__(self, process_param1, process_param2):
        """
        :param process_param1: 主进程传入的参数1
        :param process_param2: 主进程传入的参数2
        """
        super().__init__()
        self.last_event_time = time.time()  # 记录上次事件时间
        self.task_queue = process_param1
        self.task_list_path = process_param2

    def on_created(self, event):
        # 只处理文件，并且判断扩展名为 .fq
        if not event.is_directory:
            # 判断文件扩展名是否为 fq, fastq, fq.gz, fastq.gz
            valid_extensions = ('.fq', '.fastq', '.fq.gz', '.fastq.gz')
            if event.src_path.endswith(valid_extensions):
                # 等待文件写入完成（检测文件大小是否稳定）
                if is_file_complete(event.src_path):
                    self.last_event_time = time.time()
                    file_name = os.path.basename(event.src_path)
                    self.task_queue.put(file_name)
                    with open(self.task_list_path, 'a') as file:
                        file.write(file_name + '\n')



def delete_signature_files(temp_dir, bam_name):
    """
    在 temp_dir 目录下查找所有以 bam_name 开头，后面跟任意字符，
    最后以 .pickle 结尾的文件，并删除这些文件。
    
    参数:
    - temp_dir: 临时目录路径（字符串）。
    - bam_name: 文件名前缀（字符串）。
    """
    # 构造文件搜索模式
    pattern = os.path.join(temp_dir, bam_name + "*" + ".pickle")
    
    # 查找匹配的文件列表
    files = glob.glob(pattern)
    
    # 遍历并删除匹配的文件
    for file_path in files:
        try:
            os.remove(file_path)
            logging.info(f"已删除文件: {file_path}")
        except Exception as e:
            logging.error(f"删除文件 {file_path} 时出错: {e}")



# def check_condition_periodically(precision_target,recall_target,f1_target,gt_target,interval):
#     global stop_flag
#     while not stop_flag:
#         time.sleep(interval)  # 每300秒检查一次
#         precision,recall,f1,gt = cutesv_combine_cluster
#         if precision >= precision_target and recall >= recall_target and f1 >= f1_target and gt >= gt_target:
#             stop_flag = True  # 设置停止标志


def main_function():
    args = parseArgs(sys.argv[1:])
    if not args.work_dir.endswith('/'):
        args.work_dir += '/'
    task_queue = multiprocessing.Queue()
    args.high_freq_file = args.target_set
    if not os.path.exists(f'{args.work_dir}debug.txt'):
        os.mkdir("%sbam"%args.work_dir)
        os.mkdir("%scutesv_work_dir"%args.work_dir)
        os.mkdir("%scutesv_work_dir/signatures"%args.work_dir)
        os.mkdir(f'{args.output_vcf}')
        with open(f'{args.work_dir}coverage_list.txt', 'w') as file:
            pass    
        with open(f'{args.work_dir}depth_performance_rate.txt', 'w') as file:
            pass  
        with open(f'{args.work_dir}debug.txt','w') as f:
            pass
        with open(f'{args.work_dir}finished.txt','w') as f:
            pass      
    else:
        with open(f'{args.work_dir}debug.txt','r') as f:
            bam_name = f.readline()
        signatures_file = f"{args.work_dir}cutesv_work_dir/signatures"
        delete_signature_files(signatures_file, bam_name)
        name, ext = os.path.splitext(bam_name)
        bam_path = args.work_dir + "bam/" + name + ".bam"
        path_1 = bam_path
        path_2 = f"{bam_path}.bai"
        path_3 = f"{bam_path}.coverage.chr.stat.gz"
        path_4 = f'{bam_path}.coverage.chr.stat'
        if os.path.exists(path_1) and os.path.isfile(path_1):
            os.remove(path_1)
        if os.path.exists(path_2) and os.path.isfile(path_2):
            os.remove(path_2)
        if os.path.exists(path_3) and os.path.isfile(path_3):
            os.remove(path_3)
        if os.path.exists(path_4) and os.path.isfile(path_4):
            os.remove(path_4)
    mmi_path = args.mmi_path
    if mmi_path == '':
        mmi_path = f'{args.work_dir}ref.mmi'
        generate_mmi(args.reference, mmi_path)
    pandepth_path = pkg_resources.resource_filename("online", "bin/pandepth")
    shared_value = multiprocessing.Value("b",False)
    if args.recall_file == "":
        args.recall_file = f'{args.work_dir}recall_file.txt'
    p = multiprocessing.Process(target=worker, args=(task_queue, 
                                                        args.reference, 
                                                        args.work_dir, 
                                                        args.fastq_dir, 
                                                        args.platform, 
                                                        mmi_path, 
                                                        args.threads, 
                                                        shared_value, 
                                                        args.batch_interval,
                                                        args.high_freq_file,
                                                        args.output_vcf,
                                                        args.user_defined,
                                                        args.pctsize,
                                                        args.ref_dist,
                                                        args.sv_freq, 
                                                        args.recall_file,
                                                        args.target_rate,
                                                        pandepth_path))
    p.daemon = True
    p.start()

    finished_path = args.work_dir + "/finished.txt"
    task_list_path = args.work_dir + "/all_task.txt"
    tasks = arrange_task(args.fastq_dir, task_list_path, finished_path)
    
    if tasks is not None:
        for task in tasks:
                task_queue.put(task)

    if not os.path.exists(args.fastq_dir):
        raise FileNotFoundError("[Errno 2] No such directory: '%s'"%args.fastq_dir)
    event_handler = FQFileHandler(task_queue, task_list_path)
    observer = Observer()
    observer.schedule(event_handler, args.fastq_dir, recursive=False)
    observer.start()
    time_out = False
    interval = 150
    while True:  # 防止对空的sigs文件夹进行聚类
        time.sleep(interval)  # 
        if time_out == False:
            if time.time() - event_handler.last_event_time > args.monitor_fade or shared_value == True:
                time_out = True
                with open(task_list_path,"a") as f:
                    f.write(f"end observer at {time.time()}\n")
                observer.stop()
                observer.join()
                task_queue.put(None)
        if time_out and task_queue.empty():
            break

    p.join()
    if shared_value != True:
        cutesv_combine_cluster(args.work_dir, args.output_vcf, args.threads, args.reference,  
                                args.high_freq_file, args.sv_freq, args.recall_file, args.user_defined, args.pctsize, args.ref_dist)



