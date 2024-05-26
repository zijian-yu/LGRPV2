# -*- coding: utf-8 -*-
# @Author: wjq
# @Date:   2021-07-12 15:16:31
# @Last Modified by:   wjq
# @Last Modified time: 2022-02-15 19:11:19
"""
run blast and diamond

"""
import glob
import time
import itertools
import subprocess
from multiprocessing import Pool


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
                            f' -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -num_threads '
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
            # repeat
            pairs = list(itertools.product(self.species, repeat=2))
        else:
            # no repeat
            pairs = list(itertools.combinations(self.species, 2))
        return ['_'.join(pa) for pa in pairs]
        # for pa in pairs:
            # yield '_'.join(pa)

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


if __name__ == '__main__':

    process = 10
    thread = 8
    repeat = True  # 是否出现重复的 Os_At At_Os
    software = 'blast'  # 1.改 blast \ diamond
    species = []
    sub = []  # 可以为空
    nosub = []  # 可以为空

    # add_species = ["Tpr"]
    add_species = []

    # 当物种过多时
    files = glob.glob('*.pep')
    for fn in files:
        species.append(fn.split('.')[0])


    if add_species:
        for na in add_species:
            for spec in species:
                sub.append(f'{na}_{spec}')
                sub.append(f'{spec}_{na}')

        sub = list(set(sub))
    RunBioSoftware(species, repeat, sub, nosub, software, process, thread)


