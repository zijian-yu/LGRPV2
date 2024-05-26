#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :cds_genome_dot.py
# @Time :2023/9/25 21:55
# @Last time :2023/9/25 21:55
# python3 cds_genome_dot.py
# 需要的文件 1.基因组文件。2.参考物种的gff文件。3.绘制点图的cds、gff、lens文件
import os
import re
import glob
import subprocess
import pandas as pd
from Bio import SeqIO
from functools import reduce
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor


# 1. 基因组建库，参考物种进行blastn比对
class CDSBlast_one:
    def __init__(self, genome_file, species_name, cds_file):
        self.genome_file = genome_file
        self.species_name = species_name
        self.cds_file = cds_file

    def create_blastdb(self):
        db_name = self.species_name + ".db"
        db_command = f"makeblastdb -in {self.genome_file} -dbtype nucl -out {db_name}"
        subprocess.run(db_command, shell=True, check=True)

    def run_blastn(self):
        db_name = self.species_name + ".db"
        output_file = f"{os.path.splitext(self.cds_file)[0]}_{self.species_name}.blastn"
        blastn_command = f"blastn -query {self.cds_file} -db {db_name} -out {output_file} -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -num_threads 8"  # Modify blastn parameters as needed
        subprocess.run(blastn_command, shell=True, check=True)

    def process_cds_file(self):
        self.create_blastdb()
        self.run_blastn()


# 2. 根据 cds（参考物种）_genome（基因组）.blastn 文件、 genome 文件和 参考物种的gff文件，生成 基因组的cds、lens和gff文件。
class GeneAnalysis:
    def __init__(self, genome_file, blast_file, gff_file):
        self.genome_file = genome_file
        self.blast_file = blast_file
        self.gff_file = gff_file
        self.tin_genome_dict = {}
        self.pvu_gff_dict = {}
        self.blast_df = None

    def read_genome(self):
        for record in SeqIO.parse(self.genome_file, "fasta"):
            self.tin_genome_dict[record.id] = str(record.seq)

    def read_gff(self, gff_file):
        with open(gff_file, "r") as file:
            self.pvu_gff_dict = {line.split("\t")[5]: line.split("\t")[3] for line in file}

    def process_blast(self):
        self.blast_df = pd.read_csv(self.blast_file, delimiter='\t', header=None)
        quer_name = os.path.basename(self.blast_file).split(".")[0].split("_")[1]
        self.blast_df.columns = ["QueryID", "Chromosome", "PercentIdentity", "AlignmentLength", "Mismatches", "GapOpens", "QueryStart", "QueryEnd", "SubjectStart", "SubjectEnd", "EValue", "BitScore"]

        self.blast_df = self.blast_df.sort_values(by=["Chromosome", "SubjectStart"])
        new_gene_ids = []
        current_chromosome = ""
        current_gene_count = 0
        for _, row in self.blast_df.iterrows():
            if row["Chromosome"] == current_chromosome:
                current_gene_count += 1
            else:
                current_chromosome = row["Chromosome"]
                current_gene_count = 1

            # chromosome_number = current_chromosome.split('_')[2].zfill(2)
            chromosome_number = current_chromosome.split(".")[0]
            gene_count_str = str(current_gene_count).zfill(5)
            new_gene_id = f"{quer_name}{chromosome_number}g{gene_count_str}"
            new_gene_ids.append(new_gene_id)

        self.blast_df.insert(1, "NewGeneID", new_gene_ids)
        self.blast_df.insert(3, "Chain", self.blast_df["QueryID"].map(self.pvu_gff_dict))
        # self.blast_df.insert(0, "Chr", self.blast_df["Chromosome"].str.split("_").str[2])
        # self.blast_df.insert(0, "Chr", self.blast_df["Chromosome"].str.split("_").str[0].str[-2:])
        # self.blast_df.insert(0, "Chr", self.blast_df["Chromosome"].str.split(".").str[0].str[-2:])
        # self.blast_df.insert(4, "Gene_number", self.blast_df["NewGeneID"].str.split("g").str[2])
        self.blast_df.insert(0, "Chr", self.blast_df["Chromosome"].str.split(".").str[0].str[-2:])
        self.blast_df.insert(4, "Gene_number", self.blast_df["NewGeneID"].str.split("g").str[1].astype(int))

        self.blast_df = self.blast_df.drop_duplicates(subset=["SubjectEnd", "SubjectStart"])
        self.blast_df[['Chr', 'SubjectStart', 'SubjectEnd', 'Chain', 'Chromosome', 'NewGeneID', 'Gene_number']].to_csv(
            self.gff_file, sep='\t', index=False, header=False)  # 改 - 生成 gff 文件

        # lens_df = self.blast_df.groupby("Chr")["Gene_number"].max().reset_index()
        # lens_df.to_csv(f"{quer_name}.lens", sep="\t", index=False, header=False)  # 改 - 生成 lens 文件

        lens_df = self.blast_df.groupby("Chr")["Gene_number"].max().reset_index()
        lens_df = lens_df.sort_values(by="Gene_number", ascending=False)
        top_30_rows = lens_df.head(30)
        top_30_rows.to_csv(f"{quer_name}.lens", sep="\t", index=False, header=False)  # 生成 lens 文件


        cds_file = open(f"{quer_name}.cds", "w")  # 改 - 生成 cds 文件
        for _, row in self.blast_df.iterrows():
            chromosome = row["Chromosome"]
            subject_start = row["SubjectStart"]
            subject_end = row["SubjectEnd"]

            if subject_start > subject_end:
                cds_sequence = self.tin_genome_dict[chromosome][subject_end - 1:subject_start]
            else:
                cds_sequence = self.tin_genome_dict[chromosome][subject_start - 1:subject_end]
            cds_file.write(f">{row['NewGeneID']}\n{cds_sequence}\n")


# 3.所有cds 与 提取出来的基因组cds进行blastn
class CDSBlast:
    def __init__(self, species_name):
        self.species_name = species_name

    def create_blastdb(self, cds_file):
        species = os.path.splitext(cds_file)[0]
        db_name = species + ".db"
        db_command = f"makeblastdb -in {cds_file} -dbtype nucl -out {db_name}"
        self.run_command(db_command)

    def run_command(self, command):
        subprocess.run(command, shell=True, check=True)

    def run_blastn(self, query_file, db_name, output_file):
        blastn_command = f"blastn -query {query_file} -db {db_name} -out {output_file} -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -num_threads 6"
        self.run_command(blastn_command)

    def process_cds_files(self):
        cds_files = glob.glob("*.cds")
        # Create BLAST databases for all species
        with ProcessPoolExecutor(max_workers=20) as executor:
            executor.map(self.create_blastdb, cds_files)

        with ProcessPoolExecutor(max_workers=15) as executor:
            for cds_file in cds_files:
                species = os.path.splitext(cds_file)[0]
                if species != self.species_name:
                    # Perform BLASTn for A vs B
                    output_file = f"{species}_{self.species_name}.blastn"
                    executor.submit(self.run_blastn, cds_file, f"{self.species_name}.db", output_file)
                    # Perform reciprocal BLASTn for B vs A
                    reciprocal_output_file = f"{self.species_name}_{species}.blastn"
                    executor.submit(self.run_blastn, f"{self.species_name}.cds", f"{species}.db", reciprocal_output_file)


# 4. 绘制点图
class DotPlot:
    def __init__(self, num_process, suffix, bed_path, blast_path, result_path, name_change):
        self.num_process = num_process
        self.suffix = suffix
        self.name_change = name_change
        self.bed = os.path.join('.', bed_path)
        self.blast_path = os.path.join('.', blast_path)
        self.result_path = os.path.join('.', result_path)

        self.colors = ['red', 'blue', 'white']
        # self.colors = ['orangered', 'blue', 'gray']
        self.hitnum = 5
        self.align = dict(family='Times New Roman', horizontalalignment="center",
                          verticalalignment="center", weight='semibold')

    @staticmethod
    def check_file(file_path):
        n = os.path.isfile(file_path)
        file_name = os.path.basename(file_path)
        if not n:
            print(f'{file_name}---File not exist!!!')

    @staticmethod
    def get_spec_chr(spec_chr):
        spec_chr = re.sub('^\\D+', '', spec_chr)
        spec_chr = re.sub('^0', '', spec_chr)
        return spec_chr

    def gene_location(self, spec, gl):
        len_pos, gff_pos = 1, 6
        chr_dict, chr_lens, loc_gene, n = {}, {}, {}, 0
        self.check_file(f'{self.bed}/{spec}.lens')
        for li in open(f'{self.bed}/{spec}.lens'):
            lis = li.strip().split()
            spec_chr = self.get_spec_chr(lis[0])
            chr_lens[spec_chr] = float(lis[len_pos])
            chr_dict[spec_chr] = float(n)
            n += float(lis[len_pos])
        total_lens = reduce(lambda x, y: int(x) + int(y), chr_lens.values())
        step = gl / total_lens
        self.check_file(f'{self.bed}/{spec}.new.gff')
        for li in open(f'{self.bed}/{spec}.new.gff'):
            lis = li.strip().split()
            spec_chr = self.get_spec_chr(lis[0])
            if spec_chr not in chr_dict:
                continue
            loc = (chr_dict[spec_chr] + float(lis[gff_pos])) * step
            loc_gene[lis[5]] = loc
        return loc_gene, step, chr_lens

    def getnewblast(self, blast, score, evalue, repnum, loc_1, loc_2):
        newblast = {}
        for li in open(blast):
            lis = li.strip().split()
            if not all((float(lis[11]) >= score, float(lis[10]) < evalue, lis[0] != lis[1])):
                continue
            if not all((lis[0] in loc_1, lis[1] in loc_2)):
                continue
            if lis[0] in newblast and lis[1] in newblast[lis[0]]:
                continue
            if lis[0] in newblast and len(newblast[lis[0]]) < repnum:
                newblast[lis[0]].append(lis[1])
            else:
                newblast[lis[0]] = [lis[1]]
        return newblast

    def pair_positon(self, blast, loc1, loc2):
        pos1, pos2, newcolor = [], [], []
        gl_start1, gl_start2 = 11 / 12, 1 / 12
        for k, v in blast.items():
            for i in range(len(v)):
                if i == 0:
                    color = self.colors[0]
                elif i <= self.hitnum:
                    color = self.colors[1]
                else:
                    color = self.colors[2]
                pos1.append(gl_start2 + loc2[v[i]])
                pos2.append(gl_start1 - loc1[k])
                newcolor.append(color)
        return pos1, pos2, newcolor

    def plot_line(self, x, y):
        plt.plot(x, y, linestyle='-', color='black', linewidth=0.25)
        plt.plot(x, y, linestyle='-', color='black', linewidth=0.75, alpha=0.5)

    def plot_chr1(self, lens, gl, gl2, step, mark, name):
        gl_start, n, start_x = 11 / 12, 0, 1 / 12
        mark_y = 17 / 240
        for k in lens.keys():
            n += lens[k]
            mark_new = str(mark) + str(k)
            x = gl_start - float(n) * step
            mark_x = x + 0.5 * lens[k] * step
            self.plot_line([start_x, start_x + gl2], [x, x])
            plt.text(mark_y, mark_x, mark_new, color='black', fontsize=14, rotation=0, **self.align, style='normal')
        self.plot_line([start_x, start_x + gl2], [gl_start, gl_start])
        plt.text(mark_y - 0.03, 0.5 * (2 * gl_start - gl), name, color='black', fontsize=18, rotation=90,
                 **self.align, fontstyle='italic')

    def plot_chr2(self, lens, gl, gl2, step, mark, name):
        gl_start, n, start_x = 1 / 12, 0, 11 / 12
        mark_y = 223 / 240
        for k in lens.keys():
            n += lens[k]
            mark_new = str(mark) + str(k)
            x = gl_start + float(n) * step
            mark_x = x - 0.5 * lens[k] * step
            self.plot_line([x, x], [start_x, start_x - gl2])
            plt.text(mark_x, mark_y, mark_new, color='black', fontsize=14, rotation=0, **self.align, style='normal')
        self.plot_line([gl_start, gl_start], [start_x, start_x - gl2])
        plt.text(0.5 * (2 * gl_start + gl), mark_y + 0.03, name, color='black', fontsize=18, rotation=0,
                 **self.align, fontstyle='italic')

    def plotfig(self, spec1, spec2):
        plt.figure(figsize=(8, 8), dpi=300)
        root = plt.axes([0, 0, 1, 1])
        gl1, gl2 = 5 / 6, 5 / 6
        gene_loc_1, step1, chr1_lens = self.gene_location(spec1, gl1)
        gene_loc_2, step2, chr2_lens = self.gene_location(spec2, gl2)
        self.plot_chr1(chr1_lens, gl1, gl2, step1, '', self.name_change[spec1])
        self.plot_chr2(chr2_lens, gl1, gl2, step2, '', self.name_change[spec2])
        score, evalue, repnum = 100, 1e-5, 20
        blast = f'{self.blast_path}/{spec1}_{spec2}.{self.suffix}'
        blast = self.getnewblast(blast, score, evalue, repnum, gene_loc_1, gene_loc_2)
        x, y, colors = self.pair_positon(blast, gene_loc_1, gene_loc_2)
        plt.scatter(x, y, s=0.5, c=colors, alpha=0.5, edgecolors=None, linewidths=0, marker='o')
        root.set_xlim(0, 1)
        root.set_ylim(0, 1)
        root.set_axis_off()
        plt.savefig(f'{self.result_path}/{spec1}_{spec2}.png')

    def main(self):
        blasts = glob.glob(f'{self.blast_path}/*.{self.suffix}')
        with ProcessPoolExecutor(max_workers=self.num_process) as executor:
            futures = []
            for file in blasts:
                fna = os.path.basename(file)
                spec1, spec2 = fna.split('.')[0].split('_')
                futures.append(executor.submit(self.plotfig, spec1, spec2))
            for future in futures:
                future.result()


if __name__ == "__main__":
	
    genome_file = "Fam_genome.fasta"
    species_name = "Fam"
    refer_sp_name= "Oeu"  # 改 - 参考物种的简称

    # 1. 基因组建库，参考物种进行blastn比对
    cds_file = f"{refer_sp_name}.cds"
    cds_blast = CDSBlast_one(genome_file, species_name, cds_file)
    cds_blast.process_cds_file()

    # 2. 生成 基因组的cds、lens和gff文件。
    gene_analysis = GeneAnalysis(genome_file, f"{refer_sp_name}_{species_name}.blastn", f"{species_name}.new.gff")
    gene_analysis.read_genome()
    gene_analysis.read_gff(f"{refer_sp_name}.new.gff")
    gene_analysis.process_blast()

    # 3. 所有cds 与 提取出来的基因组cds进行blastn
    cds_blast = CDSBlast(species_name)
    cds_blast.process_cds_files()

    # 4. 绘制点图
    num_process = 15  # (可改) - 进程 
    blast_suffix = 'blastn'  # blast的后缀,也可以是diamond
    bed_path = ''  # 文件路径
    blast_path = ''
    result_path = ''
    # 豆科 -- 改
    name_change = {
        "Adu": "Arachis duranensis",
        "Aed": "Amphicarpaea edgeworthii",
        "Aev": "Aeschynomene evenia",
        "Ahy": "Arachis hypogaea",
        "Aip": "Arachis ipaensis",
        "Amo": "Arachis monticola",
        "Apr": "Abrus precatorius",
        "Car": "Cicer arietinum",
        "Cca": "Cajanus cajan",
        "Gma": "Glycine max",
        "Dod": "Dalbergia odorifera",
        "Gso": "Glycine soja",
        "Lal": "Lupinus albus",
        "Lan": "Lupinus angustifolius",
        "Lja": "Lotus japonicus",
        "Mal": "Melilotus albus",
        "Mpo": "Medicago polymorpha",
        "Mru": "Medicago ruthenica",
        "Msa": "Medicago sativa",
        "Mtr": "Medicago truncatula",
        "Pac": "Phaseolus acutifolius",
        "Plu": "Phaseolus lunatus",
        "Psa": "Pisum sativum",
        "Pvu": "Phaseolus vulgaris",
        "Ssu": "Spatholobus suberectus",
        "Sto": "Senna tora",
        "Tpr": "Trifolium pratense",
        "Tsu": "Trifolium subterranum",
        "Van": "Vigna angularis",
        "Vra": "Vigna radiata",
        "Vun": "Vigna unguiculata",
        "Vvi": "Vitis vinifera",
        "Bva": "Bauhinia variegata",
        "Tin":"Tin",
        "Amna":"Ammopiptanthus nanus",
        "Daol":"Dalbergia oliveri",
        "Cire":"Cicer reticulatum",
        "Abme":"Abrus melanospermus",
        "Cech":"Cercis chinensis",
        "Sesi":"Senna siamea",
        "Asmo":"Astragalus mongholicus",
        "Daco":"Dalbergia cochinchinensis",
        "Dacu":"Dalbergia cultrata",
        "Ptma":"Pterocarpus macrocarpus",
        "Ptsa":"Pterocarpus santalinus",
        "Sigl":"Sindora glabra",
        "Pral":"Prosopis alba",
        "Vapa":"Vachellia pachyceras",
        "Ceca":"Cercis canadensis",
        "Mipu":"Mimosa pudica",
        "Chfa":"Chamaecrista fasciculata",
        "Nisc":"Nissolia schottii",
        "Mupr":"Mucuna pruriens",
        "Meru":"Medicago ruthenica",
        "Latu":"Lathyrus tuberosus",
        "Maun":"Macrotyloma uniflorum",
        "Glur":"Glycyrrhiza uralensis",
        "Glla":"Glycine latifolia",
        "Oxoc":"Oxytropis ochrocephala",
        "Popi":"Pongamia pinnata",
        "Sofl":"Sophora flavescens",
        "Troc":"Trifolium occidentale",
        "Visa":"Vicia sativa",
        "Vihi":"Vigna hirtella",
        "Vitr":"Vigna trinervia",
        "Vire":"Vigna reflexopilosa",
        "Acpy":"Acacia pycnantha",
        "Faal":"Faidherbia albida",
        "Fam": "Fraxinus americana",
        "Sob": "Syringa oblata",
        "Fpe": "Fraxinus pennsylvanica",
        "Oeu": "Olea europaea",
        "Jsa": "Jasminum sambac",
    }
    p = DotPlot(num_process, blast_suffix, bed_path, blast_path, result_path, name_change)
    p.main()

