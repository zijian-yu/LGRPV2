#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_dupGen_finder.py
# @Time :2023/6/15 3:04
# python3 run_dupGen_finder.py

"""
此程序需要在指定路径下运行：/bio/yzj/tools/DupGen_finder-master

或者安装 DupGen_finder软件, 并在 DupGen_finder 文件夹中创建 input_gff_blast, output_gff_blast. output_geneDuplication 三个文件夹,。此程序放到 DupGen_finder 路径下，然后运行。

安装 DupGen_finder 的方式
    git clone https://github.com/qiao-xin/DupGen_finder.git
    cd DupGen_finder
    make
"""
import glob
import os
import concurrent.futures
import re
import shutil
import subprocess


class GeneLeixingGff:
    def __init__(self):
        self.gff_files = glob.glob(os.path.join('input_gff_blast', '*.gff'))
        self.blast_files = glob.glob(os.path.join('input_gff_blast', '*.blast'))
        self.cat_blasts = glob.glob("output_gff_blast/*.blast")

    def deal_input_gff(self, gff_file):
        file_name = os.path.basename(gff_file)
        with open(gff_file, 'r') as f, open(f"{file_name.split('.')[0]}.gff", 'w') as new_gff:
            for line_new in f:
                new_line = line_new.split()
                one = re.findall(r'[a-zA-Z]+', new_line[5].split('g')[0])[0] + "-" + "Chr" + new_line[0]
                new_gff.write(f"{one}\t{new_line[5]}\t{new_line[1]}\t{new_line[2]}\n")
            new_gff.close()

    def cat_gff(self, file_gff):
        new_file_name = f"{file_gff.split('.')[0]}_Vvi.gff"
        with open(file_gff, 'r') as f, open("Vvi.gff", 'r') as vvi, open(new_file_name, 'w') as new_gff:
            new_gff.write(f.read())
            new_gff.write(vvi.read())
        shutil.move(f"{file_gff.split('.')[0]}_Vvi.gff", "output_gff_blast")
        shutil.move(f"{file_gff.split('.')[0]}.gff", "output_gff_blast")

    def cat_blast(self, blast_file):
        file_name = os.path.basename(blast_file)
        if file_name.split(".")[0].split('_')[0] == file_name.split(".")[0].split('_')[1]:
            with open(blast_file, 'r') as f, open('input_gff_blast/Vvi_Vvi.blast', 'r') as vvi, open(f"output_gff_blast/{file_name.split('_')[0]}_Vvi.blast", 'w') as new_blast:
                new_blast.write(f.read())
                new_blast.write(vvi.read())

                species_name = file_name.split(".")[0].split('_')[0]
                shutil.copy(blast_file, f"{species_name}.blast")
                shutil.move(f"{species_name}.blast", "output_gff_blast")

    def run_dupGen_finder(self, cat_blasts, mode):
        name_one = os.path.basename(cat_blasts).split('_')[0]
        [shutil.copy(file, ".") for file in glob.glob(os.path.join("output_gff_blast", f"{name_one}*"))]
        if mode:
            cmd = f"DupGen_finder.pl -i output_gff_blast -t {name_one.split('.')[0]} -c Vvi -o output_geneDuplication"
        else:
            cmd = f"DupGen_finder-unique.pl -i output_gff_blast -t {name_one.split('.')[0]} -c Vvi -o output_geneDuplication"
        print(cmd)
        subprocess.run(cmd, shell=True)
        [os.remove(file) for file in glob.glob(f"{name_one}*")]

    def main(self):
        shutil.rmtree("output_gff_blast", ignore_errors=True)
        os.makedirs("output_gff_blast", exist_ok=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:  # （可改）
            results = [executor.submit(self.deal_input_gff, file_gff) for file_gff in self.gff_files]

        new_1_gff_files = glob.glob("*.gff")
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:  # （可改）
            results = [executor.submit(self.cat_gff, file_gff) for file_gff in new_1_gff_files]

        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:  # （可改）
            results = [executor.submit(self.cat_blast, blast_file) for blast_file in self.blast_files]

        shutil.rmtree("output_geneDuplication", ignore_errors=True)
        os.makedirs("output_geneDuplication", exist_ok=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:  # （可改）
            # 改 - mode=True：一般模式；mode=False：严格模式
            results = [executor.submit(self.run_dupGen_finder, blast, mode=True) for blast in self.cat_blasts]


if __name__ == '__main__':
    gene_leixing_gff = GeneLeixingGff()
    gene_leixing_gff.main()
