#! /usr/bin/python
import optparse
import os
import sys
import subprocess
import argparse
import glob
import json


class sv_filter():

    def __init__(self, sample_name, input_dir, output_name, hrd_file, mmr_file, target_file, pon_file, disease_name):
    
        self._sample_name = sample_name
        self._input_dir = input_dir
        # self._input_dir = "{0}/output".format(input_dir)
        self._input_list = glob.glob('%s/output/*_rearrangement_svs.out' % self._input_dir)
        self._output_name = output_name
        self._hrd_file = hrd_file
        self._mmr_file = mmr_file
        self._target_file = target_file
        self._panel_of_normal = pon_file
        self._disease_name = disease_name


    def read_db(self, input_file):
        output_list = {}

        with open(input_file, 'r') as input_data:
            for line in input_data:
                input_line = line.strip().split('\t')

                if len(input_line) == 1:
                    output_list[input_line[0]] = ''
                else:
                    output_list[input_line[0]] = input_line[1]
                    
        return output_list


    def read_pon(self, input_file):
        output_list = {}

        with open(input_file, 'r') as input_data:
            for line in input_data:
                input_line = line.strip().split('\t')
                input_key = '{0}:{1}'.format(input_line[0], input_line[1])

                output_list[input_key] = input_line[1]
                    
        return output_list


    def read_target(self, input_file, input_disease):
        output_list = []

        if input_disease == 'CNS_Brain_Cancer_NOS':
            disease_list = [input_disease, 'Glioma', 'Malignant_Brain_Tumor', 'Glioblastoma']
        elif input_disease == 'Lung_Cancer_NOS':
            disease_list = [input_disease, 'Non_Small_Cell_Lung_Cancer', 'Small_Cell_Lung_Cancer']
        elif input_disease == 'Skin_Cancer_NOS':
            disease_list = [input_disease, 'Melanoma']
        elif input_disease == 'Soft_Tissue_Cancer_Sarcoma_NOS':
            disease_list = [input_disease, 'Gastrointestinal_Stromal_Tumor', 'Ewing_Sarcoma']
        elif input_disease == 'Uterine_Cancer_NOS':
            disease_list = [input_disease, 'Endometrial_Cancer']
        else:
            disease_list = [input_disease]
        

        with open(input_file, 'r') as input_data:
            for line in input_data:
                input_line = line.strip().split('\t')
                disease = input_line[0]
                gene = input_line[1]

                if disease == 'ALL':
                    output_list.append(gene)
                elif disease in disease_list:
                    output_list.append(gene)

        return output_list

    def check_db(self, input_gene, input_list):

        if input_list.has_key(input_gene):
            result = 'Yes'
        else:
            result = 'No'

        return result


    def check_blacklist(self, input_genes):

        gene_list = input_genes.split(',')

        if 'LINC00486' in gene_list:
            result = 'Yes'
        else:
            result = 'No'

        return result


    def write_dic(self, input_data, output_name):

        with open(output_name, 'w') as output_file:
            for key in sorted(input_data.keys(), key=str.lower):
                sv_json={}
                sv_json[key]=input_data[key]
                sv_json=json.dumps(sv_json, sort_keys=True)
                output_file.write(sv_json+"\n")        


    def write_txt(self, input_data, output_name):

        with open(output_name, 'w') as output_file:
            basic_hdr = ['Sample', 'Fusion Name', 'Target', 'Type', 'Discordant Read']
            left_hdr = ['Left Gene', 'Left Chromosome', 'Left Position', 'Left Cigar', 'Left Strand', 'Left Split Read', 'Left Coverage', 'Left HRD', 'Left MMR']
            right_hdr = ['Right Gene', 'Right Chromosome', 'Right Position', 'Right Cigar', 'Right Strand', 'Right Split Read', 'Right Coverage', 'Right HRD', 'Right MMR']
            output_file.write('\t'.join(basic_hdr + left_hdr + right_hdr) + '\n')
            
            for key in sorted(input_data.keys(), key=str.lower):
                input_dic = input_data[key]
                left_dic = input_dic['left']
                right_dic = input_dic['right']

                basic_value = [self._sample_name, input_dic['fusion_name'], input_dic['target'], input_dic['type'], str(input_dic['discordant_read'])]
                left_value = [left_dic['gene'], left_dic['chr'], left_dic['pos'], left_dic['cigar'], left_dic['strand'], str(left_dic['split_read']), str(left_dic['coverage']), left_dic['hrd'], left_dic['mmr']]
                right_value = [right_dic['gene'], right_dic['chr'], right_dic['pos'], right_dic['cigar'], right_dic['strand'], str(right_dic['split_read']), str(right_dic['coverage']), right_dic['hrd'], right_dic['mmr']]

                output_file.write('\t'.join(basic_value + left_value + right_value) + '\n')


    def output_dic(self, input_data, input_type):

        # left_dic = {'gene': '', 'gene_summary': '', 'gene_background': '', 'chr': '', 'pos': '', 'cigar': '', 'strand': '', 'overlap_segment': '', 'split_read': '', 'hrd': '', 'mmr': ''}
        # right_dic = {'gene': '', 'gene_summary': '', 'gene_background': '', 'chr': '', 'pos': '', 'cigar': '', 'strand': '', 'overlap_segment': '', 'split_read': '', 'hrd': '', 'mmr': ''}
        left_dic = {}
        right_dic = {}
        output_dic = {}

        # Gene
        left_gene = input_data[0].split(',')[0].upper()
        right_gene = input_data[0].split(',')[1].upper()
        Fusion_name = '{0}--{1}'.format(left_gene, right_gene)

        output_dic['fusion_name'] = Fusion_name
        left_dic['gene'] = left_gene
        right_dic['gene'] = right_gene

        # Type
        output_dic['type'] = input_type

        # Reference Version
        output_dic['refver'] = 'hg19'

        # Target
        if left_gene in self._target_list:
            output_dic['target'] = 'Y'
        elif right_gene in self._target_list:
            output_dic['target'] = 'Y'
        else:
            output_dic['target'] = 'N'

        # Genomic coordinate
        left_dic['chr'], left_dic['pos'], right_dic['chr'], right_dic['pos'] = input_data[1].replace(',',':').split(':')
        output_key = input_data[1].replace(',','-')

        # Cigar String
        left_dic['cigar'], right_dic['cigar'] = input_data[2].split(',')

        # Mismatches
        left_dic['mismatch'], right_dic['mismatch'] = input_data[3].split(',')

        # Strand
        left_dic['strand'], right_dic['strand'] = input_data[4].split(',')

        # Breakpoint
        left_dic['breakpoint'] = ':'.join([left_dic['chr'], left_dic['pos'], left_dic['strand']])
        right_dic['breakpoint'] = ':'.join([right_dic['chr'], right_dic['pos'], right_dic['strand']])

        # Overlap segment legnth
        left_dic['overlap_segment'], right_dic['overlap_segment'] = input_data[5].split(',')

        # Split read
        left_split_read, right_split_read = input_data[7].split(',')
        left_dic['split_read'] = int(left_split_read)
        right_dic['split_read'] = int(right_split_read)

        # K-mer
        output_dic['k_mer'] = int(input_data[8])

        # Discordant read
        output_dic['discordant_read'] = int(input_data[9])

        # Breakpoint Coverage
        left_coverage, right_coverage = input_data[10].split(',')
        left_dic['coverage'] = int(left_coverage)
        right_dic['coverage'] = int(right_coverage)

        # Contig ID
        output_dic['contig_id'] = input_data[11]

        # Contig Seq
        output_dic['contig_seq'] = input_data[12]

        # HRD
        left_dic['hrd'] = self.check_db(left_gene, self._hrd_list)
        right_dic['hrd'] = self.check_db(right_gene, self._hrd_list)

        # MMR
        left_dic['mmr'] = self.check_db(left_gene, self._mmr_list)
        right_dic['mmr'] = self.check_db(right_gene, self._mmr_list)

        output_dic['left'] = left_dic
        output_dic['right'] = right_dic
        
        return output_dic, output_key


    def run(self):

        self._hrd_list = self.read_db(self._hrd_file)
        self._mmr_list = self.read_db(self._mmr_file)
        self._pon_list = self.read_pon(self._panel_of_normal)
        self._target_list = self.read_target(self._target_file, self._disease_name)

        final_list = []
        final_dic = {}

        if len(self._input_list) == 0:
            pass

        else:
            for input_file in self._input_list:

                with open(input_file, 'r') as input_data:
                    for variant_line in input_data:
                        
                        # Exclude Header
                        if variant_line[0:5] == 'genes':
                            continue
                            
                        variant_data = variant_line.strip().split('\t')
                        gene_list = variant_data[0]
                        left_breakpoint = variant_data[1].split(',')[0]
                        right_breakpoint = variant_data[1].split(',')[1]
                        variant_type = variant_data[6]

                        if not len(gene_list.split(',')) == 2:
                            pass

                        elif self.check_blacklist(gene_list) == 'Yes':
                            pass

                        elif self.check_db(left_breakpoint, self._pon_list) == 'Yes':
                            pass

                        elif self.check_db(right_breakpoint, self._pon_list) == 'Yes':
                            pass

                        elif not variant_type == 'rearrangement':
                            pass

                        else:
                            variant_dic, variant_key = self.output_dic(variant_data, 'Rearrangement')

                            final_dic[variant_key] = variant_dic

        self.write_txt(final_dic, self._output_name)
        # self.write_dic(final_dic, self._output_name)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='BreaKmer Output Filtering')
    parser.add_argument('-s', '--sample_name', type=str, help='Sample Name')
    parser.add_argument('-i', '--input_dir', type=str, help='BreaKmer Output Directory')
    parser.add_argument('-o', '--output', type=str, help='Output file name')
    parser.add_argument('-r', '--hrd_file', type=str, help='HRD list file')
    parser.add_argument('-m', '--mmr_file', type=str, help='MMR list file')
    parser.add_argument('-t', '--target_file', type=str, help='Target list file')  
    parser.add_argument('-p', '--pon', type=str, help='Panel of Normal file')    
    parser.add_argument('-d', '--disease', type=str, help='Sample\'s disease information')
    args = parser.parse_args()

    onco_sv_filter = sv_filter(args.sample_name, args.input_dir, args.output, args.hrd_file, args.mmr_file, args.target_file, args.pon, args.disease)
    onco_sv_filter.run()
