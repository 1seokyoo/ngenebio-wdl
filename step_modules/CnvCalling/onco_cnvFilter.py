#! /usr/bin/python
import sys
import argparse

class cnv_filter():

    def __init__(self, input_call, input_raw, input_target, output_name, sample_name, disease_name, hrd_file, mmr_file, ):
    
        self._input_name = input_call
        self._output_name = output_name
        self._sample_name = sample_name
        self._disease_name = disease_name
        self._hrd_list = self.read_db(hrd_file)
        self._mmr_list = self.read_db(mmr_file)

        # Set target, raw cns, blacklist
        self.set_target(input_target)
        self.set_raw_cns(input_raw)
        self._black_db = ['TRG', 'TRB', 'TRA', 'IGK', 'IGH', 'IGL','10p','18p', '21p', 'xp']

        # 1p 19q Deletion
        self._1p19q_disease = ['Cancer_Unknown_Primary_NOS', 'CNS_Brain_Cancer_NOS', 'Malignant_Brain_Tumor', 'Glioma']
        self._1p_del = 0
        self._19q_del = 0
    

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


    def check_db(self, input_gene, input_list):

        if input_list.has_key(input_gene):
            result = 'Yes'
        else:
            result = 'No'

        return result


    def copy_num(self, copy_number):

        if copy_number >= 5:
            val = 'copy_gain'
        elif copy_number <= 0:
            val = 'copy_loss'
        else:
            val = 'normal'

        return val


    def check_1p_del(self, input_line):
        input_list = input_line.strip().split('\t')
        start_pos = float(input_list[1])
        end_pos = float(input_list[2])
        chr1_p_size = float(125000000)
        
        if end_pos > chr1_p_size:
            pass
        elif (float(end_pos) - float(start_pos)) / chr1_p_size < 0.95:
            pass
        else:
            self._1p_del = 1


    def check_19q_del(self, input_line):
        input_list = input_line.strip().split('\t')
        start_pos = float(input_list[1])
        end_pos = float(input_list[2])
        chr19_q_size = float(32628983)
       
        if start_pos < float(27732282):
            pass
        elif (float(end_pos) - float(start_pos)) / chr19_q_size < 0.95:
            pass
        else:
            self._19q_del = 1


    def set_raw_cns(self, input_file):

        self._raw_cns = {}

        with open(input_file, 'r') as input_data:

            for input_line in input_data:
                input_list = input_line.strip().split('\t')

                if input_list[0] == 'chromosome':
                    continue

                input_key = '-'.join(input_list[0:3])
                self._raw_cns[input_key] = input_list[4]


    def set_target(self, input_file):

        self._target_info = {}

        with open(input_file, 'r') as input_data:

            for input_line in input_data:
                input_list = input_line.strip().split('\t')

                if input_list[0] == 'chromosome':
                    continue

                key = input_list[3] # gene
                value = '-'.join(input_list[0:3]) # chr-start-end
                if not input_list[3] in self._target_info:
                    self._target_info[key] = value
                else:
                    self._target_info[key] = ','.join([self._target_info[key], value])


    def check_target(self, input_gene, input_key):
        input_chr, input_start, input_end = input_key.split('-')

        if not input_gene in self._target_info:
            return 'None'
        else:
            gene_list = self._target_info[input_gene].split(',')

            output_chr = gene_list[0].split('-')[0]
            output_start = int(gene_list[0].split('-')[1])
            output_end = int(gene_list[-1].split('-')[2])            

            if output_start < int(input_start):
                output_start = int(input_start)
            if int(input_end) < output_end:
                output_end = int(input_end)

            for gene_value in gene_list:
                chromosome = gene_value.split('-')[0]
                start = int(gene_value.split('-')[1])
                end = int(gene_value.split('-')[2])

                if chromosome != input_chr: 
                    continue


                if end < int(input_start):
                    continue
                elif start > int(input_end):
                    continue
                else:
                    if int(input_start) <= start and start <= output_start:
                        output_start = start
                    elif start <= int(input_start) and output_start <= int(input_start):
                        output_start = input_start

                    if int(input_end) <= end:
                        output_end = int(input_end)
                    elif end <= int(input_end):
                        output_end = end

            output_txt = '-'.join([output_chr, str(output_start), str(output_end)])
            return output_txt


    def write_txt(self, input_data, output_name):

        # output_header = 'Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\tGene\tCopy_Number\tLog2_Ratio\n'
        output_header = 'Sample\tHRD\tMMR\tCnv\tChromosome\tStart Position\tEnd Position\tGene\tLog2_Ratio\tCopy_Number\n'

        with open(output_name, 'w') as output_file:
            output_file.write(output_header)

            for data in input_data:
                input_list = data.split('\t')

                output_hrd = self.check_db(input_list[6], self._hrd_list)
                output_mmr = self.check_db(input_list[6], self._mmr_list)
                output_txt = '\t'.join([input_list[0], output_hrd, output_mmr, input_list[7],
                                        input_list[1], input_list[2], input_list[3],
                                        input_list[6], input_list[8], input_list[7]])

                output_file.write(output_txt+'\n')


    def run(self):

        final_list = []
        cnv_count = 0

        with open(self._input_name, 'r') as input_data:

            for input_line in input_data:

                variant_data = input_line.strip().split('\t')
                
                # Exclude header
                if variant_data[0] == 'chromosome':
                    continue
               
                variant_key = '-'.join(variant_data[0:3])
                gene_list = variant_data[3].split(',')
                copy_number = float(variant_data[5])
                copy_value = self.copy_num(copy_number)
                log2_value = self._raw_cns[variant_key]


                # Variant filtering
                if gene_list == ['-']:
                    continue
                elif copy_value == 'normal':
                    continue
                elif variant_data[0] == 'chrX':
                    continue 

                for gene in gene_list:

                    if gene in self._black_db:
                        continue

                    gene_pos = self.check_target(gene, variant_key).split('-')
                    
                    output_list = [self._sample_name] + gene_pos + ['0', copy_value, gene, str(copy_number), log2_value]
                
                    output_txt = '\t'.join(output_list)

                    final_list.append(output_txt)

                # Check 1p deletion
                if variant_data[0] == 'chr1' and copy_value == 'copy_loss':
                    self.check_1p_del(input_line)

                # Check 19q deletion
                if variant_data[0] == 'chr19' and copy_value == 'copy_loss':
                    self.check_19q_del(input_line)

        if not self._disease_name in self._1p19q_disease:
            pass
        
        elif self._1p_del == 1 and self._19q_del == 1:
            output_list = [self._sample_name, 'chr1', '10500', '121352354', '0', 'copy_loss', '1p/19q', '-1', '1.0']
        
            output_txt = '\t'.join(output_list)

            final_list.append(output_txt)            

        self.write_txt(final_list, self._output_name)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CNVKit Call Output Filtering')
    parser.add_argument('-i', '--input_call', type=str, help='CNVkit Call CNS file')
    parser.add_argument('-j', '--input_raw', type=str, help='CNVkit Batch CNS file')
    parser.add_argument('-t', '--target', type=str, help='CNVkit targetcoverage file')
    parser.add_argument('-o', '--output', type=str, help='Output file name')
    parser.add_argument('-s', '--sample', type=str, help='Sample name')
    parser.add_argument('-d', '--disease', type=str, help='Sample\'s disease information')
    parser.add_argument('-r', '--hrd_file', type=str, help='HRD list file')
    parser.add_argument('-m', '--mmr_file', type=str, help='MMR list file')
    args = parser.parse_args()

    onco_cnv_filter = cnv_filter(args.input_call, args.input_raw, args.target, args.output, args.sample, args.disease, args.hrd_file, args.mmr_file)
    onco_cnv_filter.run()


