#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_hummer.py
# @Time :2023/12/17 18:08
# @Last time :2023/12/17 18:08
# python3 run_hummer.py
import os
import glob
import time
import shutil
import subprocess
import concurrent.futures

def run_hmmsearch(pep_file, hmm_model, CPU):
    pep_name = os.path.basename(pep_file).split(".")[0]
    model_name = os.path.basename(hmm_model).split(".")[0]
    output_file = f"{pep_name}.{model_name}.out"
    cmd = ["hmmsearch", "-o", f"output_result/{output_file}", "--noali", "-E", "1e-3", hmm_model, pep_file]
    print("Running command:", " ".join(cmd))
    subprocess.run(cmd)


if __name__ == '__main__':

    shutil.rmtree("output_result", ignore_errors=True)
    os.makedirs("output_result")

    start_time = time.time()
    pep_files = glob.glob("input_pep/*.pep")  # 改 - 输入文件夹的 pep 文件
    hmm_models = glob.glob("input_hmmer_model/*.hmm")  # 改 - hmm 模型
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:  # 改 - 最大线程数
        for pep_file in pep_files:
            for hmm_model in hmm_models:
                executor.submit(run_hmmsearch, pep_file, hmm_model, 6)  # 改 - CPU 数字
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")


