#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author : yuzijian1010@163.com
# @FileName : Batch_MSA_Tree.py
# @Time : 2024/1/11 22:33
# @Last time : 2024/1/11 22:33
# python3 Batch_MSA_Tree.py
# 使用的是 mafft, trimal, iqtree2 和软件
# 实现批量的多序列比对和构树程序
import os
import glob
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor


def run_analysis(input_file):
    folder_name = os.path.basename(input_file).split(".")[0]
    shutil.rmtree(folder_name, ignore_errors=True)
    os.makedirs(folder_name)

    shutil.copy(input_file, folder_name)
    subprocess.run(f"mafft --auto --thread 10 {folder_name}/{os.path.basename(input_file)} > {folder_name}/{os.path.basename(input_file)}.fas", shell=True)  # MAFFT
    subprocess.run(f"trimal -in {folder_name}/{os.path.basename(input_file)}.fas -out {folder_name}/{os.path.basename(input_file)}.trimal.aln -gt 0.5", shell=True)  # TrimAl
    subprocess.run(f"iqtree2 -s {folder_name}/{os.path.basename(input_file)}.trimal.aln -m MFP -bb 1000 -nt 10 -abayes", shell=True)  # IQTree


if __name__ == "__main__":
    input_files = glob.glob("input_pep/*")  # 改 - 输入文件夹
    with ThreadPoolExecutor(max_workers=4) as executor:  # 改 - 最大线程数
        executor.map(run_analysis, input_files)

