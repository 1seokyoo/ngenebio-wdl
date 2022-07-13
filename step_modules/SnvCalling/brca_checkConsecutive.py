#!/usr/bin/env python

import re
import sys
import os
import argparse
from os.path import join

class consecutive_variant(): 
    def __init__(self, input_vcf, output_vcf, sample_name, output_dir, input_bam, merge_db, variant_caller, reference, vt, min_alt_fraction):
        self._input_vcf = input_vcf # vcf_f
        self._output_vcf = output_vcf
        self._tmp_vcf = join(output_dir, '{0}_consecutive_tmp.vcf'.format(sample_name) )
        self._input_bam = input_bam # bam
        self._merge_db = merge_db # merge_db
        self._variant_caller = variant_caller # caller
        self._reference = reference # hg19
        self._vt = vt # vt
        self._min_alt_fraction = min_alt_fraction # minaltfrac


    def read_vcf(self, vcf_f):
        tmp_dic = {"header":[], "var":{}, "var_list":[]}
        with open(vcf_f) as vcf:
            for line in vcf:
                if '#' in line.split()[0]:
                    tmp_dic["header"].append(line)
                    continue
                key_tmp = line.split()
                key_tmp.remove(key_tmp[2])
                key = ".".join(key_tmp[0:4])
                tmp_dic["var"][key] = line.split()
                tmp_dic["var_list"].append(key)                                
        return tmp_dic


    def read_DB(self, txt):
        tmp_dic = {}
        with open(txt) as db:
            for line in db:
                key = ".".join(line.split()[0:7])
                tmp_dic[key] = line.split()
        return tmp_dic


    def merge_vcf(self, vcf_dic, DB_dic):
        tmp_dic = vcf_dic.copy()
        for key in DB_dic:
            total_var = key.split(".")
            var_1 = ".".join(total_var[0:4])
            var_2 = total_var[0] + "." + ".".join(total_var[4:])
            if var_1 in vcf_dic["var"] and var_2 in vcf_dic["var"]:
                wl = vcf_dic["var"][var_1][:]
                wl[0] = DB_dic[key][0]
                wl[1] = DB_dic[key][7]
                wl[3] = DB_dic[key][8]
                wl[4] = DB_dic[key][9]
                wl[-3] = vcf_dic["var"][var_1][-3].split(";CIGAR=")[0] + ";CIGAR=" + DB_dic[key][10] + ";DP=" + vcf_dic["var"][var_1][-3].split(";DP=")[1]
                del tmp_dic["var"][var_2]
                tmp_dic["var_list"].remove(var_2)
                tmp_dic["var"][var_1] = wl              
        return tmp_dic


    def check_var(self, vcf_dic, caller, reference, tmp_vcf, bam, minaltfrac):
        var_d = {}
        for var in vcf_dic["var_list"]:
            chrom = var.split(".")[0]
            pos = var.split(".")[1]
            ref = var.split(".")[2]
            var_key = chrom + "." + pos + "." + ref
            if var_key not in var_d:
                var_d[var_key] = [var]
            else:
                var_d[var_key].append(var)
                var_region = chrom + ":" + str(int(pos)-1) + "-" + str( (int(pos)-1) + len(ref)+1 )
                recall_dic = self.run_caller(caller, reference, var_region, tmp_vcf, bam, minaltfrac)
                for rm_var in var_d[var_key]:
                    del vcf_dic["var"][rm_var]
                    vcf_dic["var_list"].remove(rm_var)
                for new_var in recall_dic["var_list"]:
                    vcf_dic["var_list"].append(new_var)
                    vcf_dic["var"][new_var] = recall_dic["var"][new_var]
        return vcf_dic


    def second_check(self, vcf_dic, caller, reference, tmp_vcf, bam, minaltfrac):
        var_list = vcf_dic["var_list"]
        for var in var_list:
            chrom = var.split(".")[0]
            pos = var.split(".")[1]
            ref = var.split(".")[2]
            var_info = vcf_dic["var"][var]
            cigar = var_info[7].split(";CIGAR=")[1].split(";")[0]
            var_type = var_info[7].split(";TYPE=")[1].split(";")[0]
            if var_type == "complex":
                if len(var.split(".")[2]) == len(var.split(".")[3]) and "X" in cigar:
                    if len(cigar.split("X")) >= 3 and cigar.split("X")[-1] == "":
                        var_region = chrom + ":" + str(int(pos)-1) + "-" + str( (int(pos)-1) + len(ref)+1 )
                        recall_dic = self.run_caller(caller, reference, var_region, tmp_vcf, bam, minaltfrac)
                        del vcf_dic["var"][var]
                        vcf_dic["var_list"].remove(var)
                        for new_var in recall_dic["var_list"]:
                            if new_var not in vcf_dic["var_list"]:
                                vcf_dic["var_list"].append(new_var)
                            vcf_dic["var"][new_var] = recall_dic["var"][new_var]
        return vcf_dic


    def run_caller(self, caller, reference, var_region, tmp_vcf, bam, minaltfrac):
        cmd = '{0} -f {1} -r {2} --min-coverage 20 -F {3} --max-complex-gap 0  {4} > {5}'.format(caller, reference, var_region, minaltfrac, bam, tmp_vcf)
        os.system(cmd)
        tmp_dic = self.read_vcf(tmp_vcf)
        cmd = 'rm {0}'.format(tmp_vcf)
        os.system(cmd)
        return tmp_dic


    def write_vcf(self, vcf_dic, tmp_vcf):
        with open(tmp_vcf,'w') as new_f:
            for line in vcf_dic["header"]:
                new_f.write(line)
            for var in vcf_dic["var_list"]:
                wl = "\t".join(vcf_dic["var"][var])
                new_f.write(wl+"\n")


    def sort_vcf(self, tmp_vcf, output_vcf, vt):
        cmd = '{0} sort -o {1} {2}'.format(vt, output_vcf, tmp_vcf)
        os.system(cmd)

        cmd = 'rm {0}'.format(tmp_vcf)
        os.system(cmd)


    def main(self):
        ###run function
        self._vcf_dic = self.read_vcf(self._input_vcf)
        self._DB_dic_m = self.read_DB(self._merge_db)
        self._vcf_dic_merged = self.merge_vcf(self._vcf_dic, self._DB_dic_m)
        self._vcf_dic_1 = self.check_var(self._vcf_dic_merged, self._variant_caller, self._reference, self._tmp_vcf, self._input_bam, self._min_alt_fraction)
        self._vcf_dic_2 = self.second_check(self._vcf_dic_1, self._variant_caller, self._reference, self._tmp_vcf, self._input_bam, self._min_alt_fraction)
        self.write_vcf(self._vcf_dic_2, self._tmp_vcf)
        self.sort_vcf(self._tmp_vcf, self._output_vcf, self._vt)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Edit Consecutive Variant in BRCAaccuTest Pipeline')
    parser.add_argument('-i', '--input_vcf', type=str, help='')
    parser.add_argument('-o', '--output_vcf', type=str, help='')
    parser.add_argument('-s', '--sample_name', type=str, help='')
    parser.add_argument('-d', '--output_dir', type=str, help='')
    parser.add_argument('-b', '--input_bam', type=str, help='')
    parser.add_argument('-m', '--merge_db', type=str, help='')
    parser.add_argument('-c', '--variant_caller', type=str, help='')
    parser.add_argument('-r', '--reference', type=str, help='')
    parser.add_argument('-t', '--vt', type=str, help='')
    parser.add_argument('-f', '--min_alt_fraction', type=str, help='')

    args = parser.parse_args()

    my_class = consecutive_variant(args.input_vcf, args.output_vcf, args.sample_name, args.output_dir, args.input_bam, args.merge_db, args.variant_caller, args.reference, args.vt, args.min_alt_fraction)
    my_class.main()


