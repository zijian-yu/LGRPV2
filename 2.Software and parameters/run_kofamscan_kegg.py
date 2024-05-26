#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_kofamscan_kegg.py
# @Time :2024/5/25 14:20
# python3 run_kofamscan_kegg.py
import os
import glob
import time
import subprocess
import concurrent.futures

# kofamscan -f mapper -p /bio/yzj/tools/kofam_scan-1.3.0/profiles -k /bio/yzj/tools/kofam_scan-1.3.0/ko_list --cpu 2 --tmp-dir Adu.tmp -o Adu.kofamscan3.txt Adu.pep

def run_pfam(input_file, CPU):
    file_name = os.path.basename(input_file).split(".")[0]
    cmd = ["kofamscan", "-f", "mapper", "-p", "/bio/yzj/tools/kofam_scan-1.3.0/profiles", "-k", "/bio/yzj/tools/kofam_scan-1.3.0/ko_list", "--cpu", str(CPU), "--tmp-dir", f"{file_name}.tmp", "-o", f"{file_name}.kofamscan3.txt", input_file]
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
