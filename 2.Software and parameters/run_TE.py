#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @File    ：run_TE.py
# @Author  ：yuzijian1010@163.com
# @Date    ：2023/5/31 10:43 
# @Last time: 2023/6/19 0:21
# python3.6 run_TE.py
"""
基因组文件：sp_name + .fasta。例如（Lpu.fasta）
    "Mtr", "Adu" 基因组文件例子
"""
import os
import glob
import time
import shutil
import subprocess
import concurrent.futures
from Bio import SeqIO


class Run_TE:
    def __init__(self, species_name):
        self.species_name = species_name

    def opt_DeepTE(self):
        opt_dict = {}
        with open("opt_DeepTE.txt", 'r') as f:
            for line in f:
                new_line = line.split()
                opt_dict[new_line[0]] = new_line[1]

        with open(f"new-{self.species_name}-families.fa", 'w') as new_file:
            for line in SeqIO.parse(f"{self.species_name}-families.fa", 'fasta'):
                line_id_gai = line.id.split("#")[1]
                new_id = opt_dict[line.id].replace("_", "/")
                new_id = line.id.split("#")[0] + "#" + new_id
                print(line_id_gai, line.id, new_id)
                new_file.write(">" + str(new_id) + '\n' + str(line.seq) + '\n')

    def process_species(self):
        os.chdir(f"{self.species_name}")

        subprocess.run(["BuildDatabase", "-name", self.species_name, f"{self.species_name}.fasta"])
        subprocess.run(["/bio/yzj/tools/RepeatModeler-2.0.3/RepeatModeler", "-database", self.species_name, "-pa", "8"])  # 改 - CPU 数字
        subprocess.run(["/bio/yzj/tools/DeepTE-master/DeepTE.py", "-o", ".", "-d", "working_dir", "-i", f"{self.species_name}-families.fa", "-m_dir", "/bio/yzj/tools/DeepTE-master/Plants_model", "-sp", "P"])
        self.opt_DeepTE()
        subprocess.run(["/bio/yzj/tools/RepeatMasker_software/RepeatMasker/RepeatMasker", "-e", "ncbi", "-dir", "../TE_result/", "-pa", "25", "-qq", "-lib", f"new-{self.species_name}-families.fa", f"{self.species_name}.fasta"])  # 改 - CPU 数字

        os.chdir("..")


if __name__ == '__main__':
    start_time = time.time()
    species_list = [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join("input_genome", '*.fasta'))]
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:  # 改 - 最大线程数
        for species in species_list:
            os.makedirs(species, exist_ok=True)
            shutil.copy(os.path.join("input_genome", f"{species}.fasta"), os.path.join(species, f"{species}.fasta"))
            deep_te = Run_TE(species)
            executor.submit(deep_te.process_species)
    end_time = time.time()
    print(f"The program run time is: {end_time - start_time} 秒")

