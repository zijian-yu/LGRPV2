#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_pfam.py
# @Time :2023/6/3 0:29
# python3 run_pfam.py
import glob
import subprocess
import concurrent.futures
import time
import os

# pfam_scan.pl -fasta Lpu.pep -dir /www/tools/PfamScan/pfamdata -outfile Lpu.pfam -cpu 1

def run_pfam(input_file, CPU):
    file_name = os.path.basename(input_file).split(".")[0]
    cmd = ["pfam_scan.pl", "-fasta", input_file, "-dir", "/www/tools/PfamScan/pfamdata", "-outfile",  f"{file_name}.pfam", "-cpu", str(CPU)]
    print("Running command: ", " ".join(cmd))
    subprocess.run(cmd)


if __name__ == '__main__':
    start_time = time.time()
    pep_file = glob.glob("input/*.pep")  # （可改） - 文件的路径和文件后缀
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:  # 改 - 线程池最大线程数
        for file in pep_file:
            executor.submit(run_pfam, file, 6)  # 改 - 运行单个文件所占的CPU
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")
