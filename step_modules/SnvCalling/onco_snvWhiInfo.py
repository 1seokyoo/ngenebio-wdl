#!/usr/bin/python

import sys
import argparse
import vcf

class whitelist_info():

    def __init__(self, input_vcf_name, input_txt_name, output_name):

        self._input_vcf_name = input_vcf_name
        self._input_txt_name = input_txt_name
        self._output_name = output_name
        
        self._variant_dic = {}


    def read_vcf(self, input_data):

        for record in input_data:
            
            input_chr = record.CHROM
            input_pos = str(record.POS)
            input_ref = record.REF
            input_alt = str(record.ALT).strip().replace("[","").replace("]","").replace(" ","")
            input_end = str(record.INFO['END'])

            if record.FILTER == []:
                input_filter = "PASS"
            else:
                input_filter = "REJECT"

            input_key = '-'.join([input_chr, input_pos, input_ref, input_alt])
            self._variant_dic[input_key] = input_filter
            

    def run(self):

        # Read VCF File
        input_vcf_file = open(self._input_vcf_name, 'r')
        input_reader = vcf.Reader(input_vcf_file)
        self.read_vcf(input_reader)

        # Write Txt File
        output_file = open(self._output_name, 'w')

        with open(self._input_txt_name, 'r') as input_data:
            for input_line in input_data:
                input_list = input_line.strip().split('\t')

                input_chr = input_list[0]
                input_pos = input_list[1]
                input_ref = input_list[2]
                input_alt = input_list[3]

                input_key = '-'.join([input_chr, input_pos, input_ref, input_alt])
                if input_key in self._variant_dic:
                    input_filter = self._variant_dic[input_key]
                else:
                    input_filter = 'REJECT'

                output_list = input_list + [input_filter]
                output_file.write('\t'.join(output_list)+'\n')

        # VCF File Close
        input_vcf_file.close()
        output_file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='VCF Format Filtering')
    parser.add_argument('-i', '--input_vcf', type=str, help='VCF input file name', default='None')
    parser.add_argument('-j', '--input_txt', type=str, help='TXT input file name', default='None')
    parser.add_argument('-o', '--output', type=str, help='TXT output file name', default='None')
    args = parser.parse_args()

    ngb_whitelist_info = whitelist_info(args.input_vcf, args.input_txt, args.output)
    ngb_whitelist_info.run()
