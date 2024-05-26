#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :ks_circos_web_gff_lens.py
# @Time :2024/5/20 15:42
# @Last time: 2024/5/20 15:42
# python3 ks_circos_web_gff_lens.py
import os
import glob
import shutil
import pandas as pd


# 1. 生成新的四列 gff 文件
def gff_files_deal(input_folder, output_folder):
    input_files = glob.glob(os.path.join(input_folder, "*.gff"))
    for input_file in input_files:
        file_name = os.path.basename(input_file)
        sp_name = os.path.basename(input_file).split(".")[0]
        df = pd.read_csv(input_file, sep='\t', header=None, names=["Column1", "Column2", "Column3", "Column4", "Column5", "Column6", "Column7"])
        df["Column1"] = sp_name + df["Column1"].astype(str)
        output_file = os.path.join(output_folder, file_name)
        selected_columns = df[["Column1", "Column6", "Column2", "Column3"]]
        selected_columns.to_csv(output_file, sep='\t', index=False, header=False)

# 2. 生成两列 lens 文件
def lens_files_deal(input_folder, output_folder):
    gff_files = glob.glob(os.path.join(input_folder, '*.gff'))
    for file_path in gff_files:
        file_name = os.path.basename(file_path).split(".")[0] + ".lens"
        df = pd.read_csv(file_path, sep='\t', header=None)
        df = df.groupby(0, sort=False).apply(lambda x: x.iloc[[-1], [0, 3]])
        output_file_path = os.path.join(output_folder, file_name)
        df.to_csv(output_file_path, sep='\t', header=False, index=False)


if __name__ == "__main__":
    input_gff = "input_gff"  # 改 - 输入文件夹的名字
    output_gff = "output_gff"  # 改 - 输出文件夹的名字
    output_lens = "output_lens"
    shutil.rmtree(output_gff, ignore_errors=True)
    os.makedirs(output_gff)
    shutil.rmtree(output_lens, ignore_errors=True)
    os.makedirs(output_lens)

    # 1. 生成新的四列 gff 文件
    gff_files_deal(input_gff, output_gff)
    # 2. 生成两列 lens 文件
    lens_files_deal(output_gff, output_lens)


