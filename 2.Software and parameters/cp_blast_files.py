#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_dupGen_finder.py
# @Time :2024/3/29 0:12
# python3 cp_blast_files.py
# 从源路径，复制需要的blast文件到目标路径
import glob
import os
import shutil


def copy_files_by_pattern(source_folder, destination_folder):
    file_paths = glob.glob(os.path.join(source_folder, '*'))
    for file_path in file_paths:
        file_name = os.path.basename(file_path)
        first_part = file_name.split(".")[0]
        prefix = first_part.split("_")[0]
        suffix = first_part.split("_")[-1]
        if prefix == suffix or prefix == "Vvi" or suffix == "Vvi":
            destination_file_path = os.path.join(destination_folder, file_name)
            shutil.copy(file_path, destination_file_path)


if __name__ == "__main__":
    # 调用函数并传入路径参数
    source_folder = '/wangjp/yuzijian/F_plantform/diamond_run/Bva_Seca/output_pep/blastp'  # 改 - 源路径文件
    destination_folder = '/bio/yzj/tools/DupGen_finder-master/input_gff_blast'  # 改 - 目标路径
    copy_files_by_pattern(source_folder, destination_folder)
    print("over!")
