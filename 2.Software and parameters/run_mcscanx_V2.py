#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_mcscanx_V2.py
# @Time :2024/5/11 12:03
# @Last time: 2024/5/11 12:03
# python3 run_mcscanx_V2.py
import os
import time
import glob
import shutil
import itertools
import subprocess
import pandas as pd
import concurrent.futures


class MCScanX:
    def __init__(self):
        self.gff_files = glob.glob("blastp/*.gff")
        self.name = self.get_name()

    def get_name(self):
        sp_name = []
        for line in self.gff_files:
            sp_name.append(shutil.os.path.basename(line).split('.')[0])
        new = list(itertools.permutations(sp_name, 2))
        return new

    def deal_blast(self, one, two):
        cat_command = f"cat blastp/{one}_{one}.blast blastp/{one}_{two}.blast blastp/{two}_{two}.blast blastp/{two}_{one}.blast > out_gff_blast/{one}_{two}.blast"
        subprocess.run(cat_command, shell=True)
        print(f"Blast file for {one}_{two} cat finish!!!")


    def deal_gff(self, one, two):
        sp_brr_one = one
        sp_brr_two = two
        # 处理 one 的 gff 文件
        df_one = pd.read_csv(f"blastp/{sp_brr_one}.new.gff", sep='\t', header=None)
        df_one = df_one.iloc[:, [0, 5, 1, 2]]
        df_one[0] = sp_brr_one + df_one[0].astype(str)
        df_one.to_csv(f"out_gff_blast/{sp_brr_one}.gff", sep='\t', header=None, index=None)
        print(f"gff file {one}.gff deal 4 col finish!")
        # 处理 two 的 gff 文件
        df_two = pd.read_csv(f"blastp/{sp_brr_two}.new.gff", sep='\t', header=None)
        df_two = df_two.iloc[:, [0, 5, 1, 2]]
        df_two[0] = sp_brr_two + df_two[0].astype(str)
        df_two.to_csv(f"out_gff_blast/{sp_brr_two}.gff", sep='\t', header=None, index=None)
        print(f"gff file {two}.gff deal 4 col finish!")

        gff_one = pd.read_csv(f"out_gff_blast/{one}.gff", sep='\t', header=None)
        gff_two = pd.read_csv(f"out_gff_blast/{two}.gff", sep='\t', header=None)
        new_gff = pd.concat([gff_one, gff_two], ignore_index=True)
        new_gff.to_csv(f"out_gff_blast/{one}_{two}.gff", sep='\t', header=None, index=None)
        print(f"gff files {one}_{two}.gff deal MCScanX finish!")

    def run_mcscanx(self, one, two):
        print(one, two)
        self.deal_gff(one, two)
        self.deal_blast(one, two)
        
        shutil.copy(f"out_gff_blast/{one}_{two}.gff", ".")
        shutil.copy(f"out_gff_blast/{one}_{two}.blast", ".")

        subprocess.run(f"MCScanX {one}_{two}", shell=True)

        # synvisio-共线性要求文件（网站用）
        shutil.copy(f"{one}_{two}.collinearity", "output_synvisio")
        shutil.move(f"output_synvisio/{one}_{two}.collinearity", f"output_synvisio/{one}_{two}_collinear.collinearity")
        shutil.copy(f"{one}_{two}.gff", "output_synvisio")
        shutil.move(f"output_synvisio/{one}_{two}.gff", f"output_synvisio/{one}_{two}_coordinate.gff")

        shutil.move(f"{one}_{two}.tandem", "output")
        shutil.move(f"{one}_{two}.collinearity", "output")
        shutil.move(f"{one}_{two}.html", "output")

        os.remove(f"{one}_{two}.gff")
        os.remove(f"{one}_{two}.blast")
        os.remove(f"out_gff_blast/{one}_{two}.gff")
        os.remove(f"out_gff_blast/{one}_{two}.blast")


if __name__ == '__main__':
    start_time = time.time()
    mc = MCScanX()
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(mc.run_mcscanx, one, two) for one, two in mc.name]
        concurrent.futures.wait(futures)
    end_time = time.time()
    print(f"Total running time: {end_time - start_time} seconds")

