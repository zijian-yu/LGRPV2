# -*- coding: utf-8 -*-
# @Author: wjq
# @Date:   2021-07-12 15:16:31
# @Last Modified by:   wjq
# @Last Modified time: 2021-09-22 23:24:55
"""
calculata ks
"""
import os
import glob
import subprocess
from multiprocessing import Pool


class RunBioSoftware:
    def __init__(self, num_process):
        self.num_process = num_process
        self.main()

    @staticmethod
    def check_file(file_path):
        n = os.path.isfile(file_path)
        # file_name = os.path.basename(file_path)
        # if not n:
            # print(f'{file_name}---File not exist!!!')
        return n
        # assert n, "File not exist!!!"

    def run_ks(self, path, blk, cds1, cds2):
        print('start!!!')
        spec1, spec2 = path.split('_')
        ks_sf = f'{spec1}_{spec2}.ks.txt'
        if not os.path.isdir(path):
            os.mkdir(path)
        if spec1 == spec2:
            cds_f = f'{path}/{spec1}.cds'
            os.system(f'cp {cds1} {cds_f}')
        else:
            cds_f = f'{path}/{spec1}_{spec2}.fasta'
            os.system(f'cat {cds1} {cds2} > {cds_f}')

        p = subprocess.call(f'perl calculate.Ks.pl {path} {spec1} {blk} {cds_f} {ks_sf}', shell=True)
        if p != 0:
            print(f'Warning!!! {p}')
        os.system(f'rm -rf {path}')

    def main(self):
        p = Pool(self.num_process)
        files = glob.glob(f'./blk/*.block.rr.txt')
        for fn in files:
            if 'partially' in fn or 'less' in fn:
                continue
            name = os.path.basename(fn).split('.')[0]
            spec1, spec2 = name.split('_')
            n1 = self.check_file(f'./cds//{spec1}.cds')
            n2 = self.check_file(f'./cds/{spec2}.cds')
            if not n1:
                 cds1 = f'./cds/{spec1}.cds.fasta'
            else:
                cds1 = f'./cds/{spec1}.cds'
            if not n2:
                cds2 = f'./cds/{spec2}.cds.fasta'
            else:
                cds2 = f'./cds/{spec2}.cds'
            # ks_sf = f'{spec1}_{spec2}.ks.txt'
            # if not os.path.isdir(name):
            #     os.mkdir(name)
            # if spec1 == spec2:
            #     cds_f = f'{name}/{spec1}.cds'
            #     os.system(f'cp ./cds/{spec1}/{spec1}.cds {cds_f}')
            #     p.apply_async(self.run_ks, args=((name, spec1, fn, cds_f, ks_sf)))
            # else:
            #     cds_f = f'{name}/{spec1}_{spec2}.fasta'
            #     os.system(f'cat ./cds/{spec1}/{spec1}.cds ./cds/{spec2}/{spec2}.cds > {cds_f}')
            p.apply_async(self.run_ks, args=((name, fn, cds1, cds2)))

        p.close()
        p.join()


if __name__ == '__main__':

    process = 12  # default
    RunBioSoftware(process)


