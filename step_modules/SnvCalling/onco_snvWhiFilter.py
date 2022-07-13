#!/usr/bin/python

import sys
import argparse
import vcf

class whitelist_filter():

    def __init__(self, input_name, output_vcf, output_txt, target):

        self._input_name = input_name
        self._output_vcf = output_vcf
        self._output_txt = output_txt

        self._target = target
        
        self._variant_dic = {}
        self._depth_dic = {}


    def read_vcf(self, input_data):

        for record in input_data:
            
            input_chr = record.CHROM
            input_pos = str(record.POS)
            input_ref = record.REF
            input_alt = str(record.ALT).strip().replace("[","").replace("]","").replace(" ","")
            input_end = str(record.INFO['END'])

            if not input_alt == 'None':
                input_key = '-'.join([input_chr, input_pos, input_ref, input_alt])           
                self._variant_dic[input_key] = record

            input_range = range(int(input_pos), int(input_end) + 1)
            for pos in input_range:
                total_key = '-'.join([input_chr, str(pos)])
                self._depth_dic[total_key] = record.INFO['DP']
                
                # +1 Postion
                total_key = '-'.join([input_chr, str(pos +1)])
                self._depth_dic[total_key] = record.INFO['DP']
                

    def run(self):

        # Read VCF File
        input_file = open(self._input_name, 'r')
        input_reader = vcf.Reader(input_file)

        # Write VCF File
        output_vcf_file = open(self._output_vcf, 'w')
        output_vcf_writer = vcf.Writer(output_vcf_file, input_reader)

        # Write Txt File
        output_txt_file = open(self._output_txt, 'w')

        self.read_vcf(input_reader)

        with open(self._target, 'r') as input_data:
            for input_line in input_data:
                input_list = input_line.strip().split('\t')

                input_chr = input_list[0]
                input_pos = input_list[1]
                input_ref = input_list[2]
                input_alt = input_list[3]
                input_gene = input_list[4].split('_')[0]
                input_aaChange = input_list[4].split('_')[1]

                input_key = '-'.join([input_chr, input_pos, input_ref, input_alt])
                if input_key in self._variant_dic:
                    input_depth = str(self._variant_dic[input_key].INFO['DP'])
                    input_vaf = str(format(float(self._variant_dic[input_key].INFO['AF'][0]) * 100, ".2f"))
                    output_vcf_writer.write_record(self._variant_dic[input_key])
                else:
                    input_key = '-'.join([input_chr, input_pos])
                    input_depth = str(self._depth_dic[input_key])
                    input_vaf = '0'

                output_list = [input_chr, input_pos, input_ref, input_alt, input_gene, input_aaChange, input_depth, input_vaf]
                output_txt_file.write('\t'.join(output_list)+'\n')

        # VCF File Close
        input_file.close()
        output_vcf_file.close()
        output_txt_file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='VCF Format Filtering')
    parser.add_argument('-i', '--input', type=str, help='VCF input file name', default='None')
    parser.add_argument('-o', '--output_vcf', type=str, help='VCF output file name', default='None')
    parser.add_argument('-u', '--output_txt', type=str, help='TXT output file name', default='None')
    parser.add_argument('-t', '--target', type=str, help='Whitelist target file name', default='None')
    args = parser.parse_args()

    ngb_whitelist_filter = whitelist_filter(args.input, args.output_vcf, args.output_txt, args.target)
    ngb_whitelist_filter.run()
