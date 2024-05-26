#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :iTAK_run.py
# @Time :2023/8/31 17:34
# python3 iTAK_run.py
# 如果蛋白文件的序列长度小于50（猜的）会报错。
import glob
import subprocess
import concurrent.futures
import time
import os
import shutil


# 1.run iTAK
def run_pfam(input_file, CPU):
    file_name = os.path.basename(input_file).split(".")[0]
    cmd = ["iTAK.pl", input_file]
    print("Running command: ", " ".join(cmd))
    subprocess.run(cmd)


# 1.结果移动到output文件夹
def move_files():
    """
    这个函数将文件从'input_pep'文件夹中的文件夹移动到'output'文件夹。
    """
    os.makedirs("output", exist_ok=True)
    input_files = glob.glob("input_pep/*")
    for file in input_files:
        if os.path.isdir(file):
            shutil.move(file, os.path.join("output", os.path.basename(file)))


if __name__ == '__main__':
    start_time = time.time()
    pep_file = glob.glob("input_pep/*.pep")  # （可改） - 文件的路径和文件后缀
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:  # 改 - 线程池最大线程数
        for file in pep_file:
            executor.submit(run_pfam, file, 6)  # 改 - 运行单个文件所占的CPU
    end_time = time.time()
    
    move_files()
    print(f"Time taken: {end_time - start_time} seconds")


