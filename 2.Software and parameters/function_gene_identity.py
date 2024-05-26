#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :function_gene_identity.py
# @Time :2023/11/14 12:00
# @Last time: 2023/11/20 12:02
# python3 function_gene_identity.py
# 与拟南芥的家族基因pep文件进行blastp比对，然后根据物种与拟南芥的blastp比对结果，提取物种的家族信息。
import os
import glob
import time
import shutil
import itertools
import subprocess
from Bio import SeqIO
from multiprocessing import Pool
from collections import defaultdict


# 1 与拟南芥进行 blastp 比对
class RunBioSoftware:
    def __init__(self, spec, repeat, align, noalign, flag, num_process, thread):
        self.species = spec
        self.repeat = repeat
        self.subs = align
        self.nosubs = noalign
        self.flag = flag
        self.num_process = num_process
        self.thread = thread
        self.main()

    def run_blast(self, pair):
        spec1, spec2 = pair.split('_')
        p = subprocess.call(f'blastp -query {spec1}.pep -out {pair}.blast -db {spec2}.db'
                            f' -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads '
                            f'{self.thread}', shell=True)
        if p != 0:
            print(f'Warning!!! {p}')

    def run_diamond(self, pair):
        spec1, spec2 = pair.split('_')
        p = subprocess.call(f'diamond blastp -d {spec2}.nr -q {spec1}.pep -e 1e-5 -f 6 '
                            f'-p {self.thread} -o {pair}.diamond -k 50', shell=True)
        if p != 0:
            print(f'Warning!!! {p}')

    def get_spec(self):
        spec = []
        for pair in self.subs:
            spec.extend(pair.split('_'))
        return list(set(spec))

    def get_pair(self):
        if self.repeat:
            pairs = list(itertools.product(self.species, repeat=2))
        else:
            pairs = list(itertools.combinations(self.species, 2))
        return ['_'.join(pa) for pa in pairs]

    def remove_sub(self):
        new_sub = []
        for pair in self.subs:
            if pair not in self.nosubs:
                new_sub.append(pair)
        self.subs = new_sub

    def main(self):
        p = Pool(self.num_process)
        if self.subs:
            self.species = self.get_spec()
        else:
            self.subs = self.get_pair()

        if self.nosubs:
            self.remove_sub()

        if self.flag == 'blast':
            for s in self.species:
                subprocess.call(f'makeblastdb -in {s}.pep -dbtype prot -out {s}.db', shell=True)
            for pair in self.subs:
                timestr = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                print(f'run ------------------------------{pair}--{timestr}')
                p.apply_async(self.run_blast, args=(pair,))
        else:
            for s in self.species:
                subprocess.call(f'diamond makedb --in {s}.pep --db {s}.nr', shell=True)
            for pair in self.subs:
                timestr = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                print(f'run ------------------------------{pair}--{timestr}')
                p.apply_async(self.run_diamond, args=(pair,))

        p.close()
        p.join()


# 2 将 blast 比对结果，移动到 input_blastp 文件夹中
def move_blast_results(input_folder):
    shutil.rmtree(input_folder, ignore_errors=True)
    os.makedirs(input_folder)
    blast_files = [file for file in os.listdir() if file.endswith('.blast')]
    for file in blast_files:
        shutil.move(file, os.path.join(input_folder, file))


# 3 根据物种与拟南芥的blastp比对结果，提取物种的家族信息。
class GeneProcessor_Ath_function:
    def __init__(self, input_blastp, output_1_sp_GF_id, output_2_GF_id_number, gene_dict_file):
        self.input_blastp = input_blastp
        self.output_1_sp_GF_id = output_1_sp_GF_id
        self.output_2_GF_id_number = output_2_GF_id_number
        self.gene_dict_result = self.read_gene_file(gene_dict_file)

    def read_gene_file(self, file_path):
        gene_dict = {}
        with open(file_path, "r") as gene_file:
            for line in gene_file:
                gene_id, gene_value = map(str.strip, line.split("\t"))
                gene_dict[gene_value.upper()] = gene_id
        return gene_dict

    def process_blastp_file(self, input_file_path, output_folder):
        file_name = f"{os.path.basename(input_file_path).split('_')[0]}.GF.txt"
        output_file_path = os.path.join(output_folder, file_name)
        seen_pairs = set()
        with open(input_file_path, "r") as blastp_file, open(output_file_path, 'w') as output_file:
            for line in blastp_file:
                fields = line.strip().split("\t")
                gene_id = fields[1].upper()
                if gene_id in self.gene_dict_result:
                    pair = (fields[0], self.gene_dict_result[gene_id])
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        output_file.write(f"{fields[0]}\t{self.gene_dict_result[gene_id]}\n")

    def process_output_go_files(self, input_folder):
        file_list = glob.glob(os.path.join(input_folder, "*"))
        for file_path in file_list:
            go_dict = defaultdict(list)
            with open(file_path, 'r') as file:
                for line in file:
                    gene_id, go_terms = map(str.strip, line.split('\t'))
                    go_dict[go_terms].append(gene_id)
            output_file_path = os.path.join(self.output_2_GF_id_number, os.path.basename(file_path))
            with open(output_file_path, 'w') as output_file:
                for go_term, gene_ids in go_dict.items():
                    output_file.write(f"{go_term}\t{len(gene_ids)}\t{','.join(gene_ids)}\n")

    def process_genes(self):
        # 1 根据物种与拟南芥的blastp比对结果，提取物种的家族信息
        shutil.rmtree(self.output_1_sp_GF_id, ignore_errors=True)
        os.makedirs(self.output_1_sp_GF_id)

        for file_path in glob.glob(f"{self.input_blastp}/*"):
            self.process_blastp_file(file_path, self.output_1_sp_GF_id)

        # 2 统计为基因家族对应物种基因数目和基因id
        shutil.rmtree(self.output_2_GF_id_number, ignore_errors=True)
        os.makedirs(self.output_2_GF_id_number)

        self.process_output_go_files(self.output_1_sp_GF_id)


if __name__ == '__main__':
    start_time = time.time()
    # 1 与拟南芥进行 blastp 比对
    process = 10  # 改 - 线程数
    thread = 12  # 改 - blast占CPU数
    repeat = True
    software = 'blast'  # 改 - blast or diamond
    species = []
    sub = []
    nosub = []
    add_species = ["AthGF"]  # 拟南芥的pep文件名称
    species = [fn.split('.')[0] for fn in glob.glob('*.pep')]
    if add_species:
        sub = list(set([f'{spec}_{na}' for na in add_species for spec in species]))
    RunBioSoftware(species, repeat, sub, nosub, software, process, thread)
    print("1. blastp run over!")

    # 2 将 blast 比对结果，移动到 input_blastp 文件夹中
    move_blast_results("input_blastp")   # (可改) - 输入

    # 3 根据物种与拟南芥的blastp比对结果，提取物种的家族信息。
    gene_processor = GeneProcessor_Ath_function(
        input_blastp="input_blastp",  #（可改）- blast文件
        output_1_sp_GF_id="output_1_sp_GF_id",  #（可改）- 基因id对应家族文件
        output_2_GF_id_number="output_2_GF_id_number",  #（可改）- 家族对应物种的基因id文件
        gene_dict_file="Arh_GF_id.txt"  #（可改）- 拟南芥家族对应的基因id文件
    )
    gene_processor.process_genes()
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Functional family gene identification completed! Total execution time: {elapsed_time} seconds")


