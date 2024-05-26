#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_interproscan.py
# @Time :2023/6/3 0:29
# @Last time :2023/12/19 22:03
# python3 run_interproscan.py
import os
import time
import glob
import shutil
import subprocess
import concurrent.futures


def run_interproscan(input_file, CPU):
    cmd = ["interproscan.sh", "-i", input_file, "-b", "result/", "-goterms", "-iprlookup", "-pa", "-dp", "-f", "tsv", "-cpu", str(CPU)]
    print("Running command: ", " ".join(cmd))
    subprocess.run(cmd)


if __name__ == '__main__':

    shutil.rmtree("result", ignore_errors=True)
    os.makedirs("result")

    start_time = time.time()
    pep_file = glob.glob("input/*.pep")  # （可改） - 文件的路径和文件后缀
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:  # 改 - 线程池最大线程数
        for file in pep_file:
            executor.submit(run_interproscan, file, 10)  # 改 - 运行单个文件所占的CPU
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")

# InterProScan的TSV结果文件中的每一列都包含不同的注释信息，以下是常见的一些列及其含义：
# 1. **Protein Accession**: 这一列通常包含蛋白质的访问号，如UniProtKB/TrEMBL号码或其他标识符。
# 2. **MD5-Hash**: 这一列包含与蛋白质序列的MD5哈希值相关的信息，用于唯一标识蛋白质。
# 3. **Sequence Length**: 这一列包含蛋白质序列的长度。
# 4. **Analysis Accession**: 这一列可能包含与分析或数据库的访问号相关的信息，通常用于跟踪分析的来源或版本。
# 5. **Signature Accession**: 这一列包含与InterPro数据库中的特定签名或域的访问号相关的信息。
# 6. **Signature Description**: 这一列包含与InterPro签名或域的描述性信息，通常会解释该签名或域的功能。
# 7. **Start Position**: 这一列包含匹配的签名或域在蛋白质序列中的起始位置。
# 8. **End Position**: 这一列包含匹配的签名或域在蛋白质序列中的结束位置。
# 9. **e-value**: 这一列通常包含匹配的签名或域的e-value。
# 10. **Status**: 这一列可能包含有关匹配状态的信息，例如"T"表示正常匹配，"U"表示未知。
# 11. **Date**: 这一列包含匹配或注释的日期信息。
# 12. **InterPro Accession**: 这一列包含与InterPro数据库中的相关信息或注释的访问号。
# 13. **InterPro Description**: 这一列包含与InterPro数据库中相关信息或注释的描述性信息，通常会解释注释的功能。
# 14. **GO Terms**: 这一列通常包含与Gene Ontology (GO)注释相关的信息，包括GO号。
# 15. **Pathway**: 这一列可能包含与生物学途径或通路相关的信息。
# 16. **Protein Data Source**: 这一列包含蛋白质数据的来源，例如UniProtKB。
