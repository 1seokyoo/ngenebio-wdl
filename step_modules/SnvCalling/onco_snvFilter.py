#!/usr/bin/python

import sys
import argparse
import vcf
import json

class filter_format():

    def __init__(self, input_name, output_name, db_list, vcf_format, cut_off, clinvar_filter, remove_flag, recover_file):

        self._input_name = input_name
        self._output_name = output_name
        self._db_list = db_list
        self._vcf_format = vcf_format
        self._cut_off = cut_off
        self._clinvar_filter = clinvar_filter
        self._remove_flag = remove_flag
        
        self._recover_dic = None
        with open(recover_file) as sj:
            self._recover_dic = json.load(sj)


    def check_recover(self, key, info):

        output_list = []
        if key in self._recover_dic:
            if info in self._recover_dic[key]:
                output_list = str(self._recover_dic[key][info]).split('|')

        return output_list

    def check_database(self, key, record, db_list):
        
        db_list = db_list.split(',')
        false_db = []
        recover_db = []
        recover_list = self.check_recover(key, 'FP_DB')

        for filter_db in db_list:
            if filter_db in recover_list:
                recover_db.append(filter_db)
            elif filter_db in record.INFO:
                record.FILTER = 'REJECT'
                false_db.append(filter_db)

        if 'ExAC_EAS_AF' in db_list and 'ExAC_EAS_AC' in record.INFO:
            input_totalCount = float(record.INFO["ExAC_EAS_AN"])
            input_altCount = float(record.INFO["ExAC_EAS_AC"][0])
            input_alleleFrequency = 0
            
            if not input_totalCount == 0:
                input_alleleFrequency = round((input_altCount / input_totalCount), 4)
                record.INFO['ExAC_EAS_AF'] = str(input_alleleFrequency)

            if input_alleleFrequency >= 0.01:
                record.FILTER = 'REJECT'
                false_db.append('ExAC_EAS_AF')

        return record, recover_db, false_db


    def check_format(self, record, vcf_format, cut_off):

        false_value = 'NONE'
        if not record.samples[0][vcf_format] >= int(cut_off):
            record.FILTER = 'REJECT'
            false_value = '{0}<{1}'.format(vcf_format, str(cut_off))

        return record, false_value


    def check_clinvar(self, record):

        clinvar_info = 'None'

        if not 'CLNSIG' in record.INFO:
            pass
        else:
            input_clnsig = str(record.INFO["CLNSIG"]).replace("]","").replace("[","").replace("'","").replace('"',"").replace(" ","").strip()
            input_review = str(record.INFO["CLNREVSTAT"]).replace("]","").replace("[","").replace("'","").replace('"',"").replace(" ","").strip()

            if input_review == 'criteria_provided,_multiple_submitters,_no_conflicts' and 'enign' in input_clnsig:
                clinvar_info = '2:{0}'.format(input_clnsig)
                record.FILTER = 'REJECT'
            elif input_review == 'reviewed_by_expert_panel' and 'enign' in input_clnsig:
                clinvar_info = '3:{0}'.format(input_clnsig)
                record.FILTER = 'REJECT'

        return record, clinvar_info


    def run(self):

        # Read VCF File
        self._input_file = open(self._input_name, 'r')
        self._vcf_reader = vcf.Reader(self._input_file)

        # Write VCF File
        self._output_file = open(self._output_name, 'w')
        self._vcf_writer = vcf.Writer(self._output_file, self._vcf_reader)


        # Filtering
        for record in self._vcf_reader:
            ALT = str(record.ALT).replace("[","").replace("]","").replace(" ","").strip()
            key = '{0}-{1}-{2}-{3}'.format(record.CHROM, record.POS, record.REF, ALT)

            false_list = []

            # 1st VarDic Filter Column
            if record.FILTER == []: # FILTER == PASS
                record.FILTER = 'PASS'
            else:
                record.INFO['FP_FILTER'] = record.FILTER
                record.FILTER = 'REJECT'

            # 2nd Filter VCF Population DB
            recover_db = []
            false_db = []
            if not self._db_list == 'None':
                record, recover_db, false_db = self.check_database(key, record, self._db_list)
            if not recover_db == []:           
                record.INFO['TP_RECOVER'] = '|'.join(recover_db)
            if not false_db == []:           
                record.INFO['FP_DB'] = '|'.join(false_db)

            # 3rd ClinVar Filter
            clinvar_info = 'None'
            if self._clinvar_filter == 'Yes':
                record, clinvar_info = self.check_clinvar(record)
            if not clinvar_info == 'None':
                record.INFO['FP_CLINVAR'] = clinvar_info

            # 4th VCF Format Filter
            false_value = 'NONE'
            if not self._vcf_format == 'None':
                record, false_value = self.check_format(record, self._vcf_format, self._cut_off)
                
            if not false_value == 'NONE':
                record.INFO['FP_FORMAT'] = false_value

            self._vcf_writer.write_record(record)

        # VCF File Close
        self._input_file.close()
        self._output_file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='VCF Format Filtering')
    parser.add_argument('--input', type=str, help='VCF input file name', default='None')
    parser.add_argument('--output', type=str, help='VCF output file name', default='None')
    parser.add_argument('--database', type=str, help='Filtering database list (list delimeter is ",")', default='None')
    parser.add_argument('--format', type=str, help='Filtering VCF format column', default='None')
    parser.add_argument('--cutoff', type=str, help='Format cut-off value', default='None')
    parser.add_argument('--clinvar', type=str, help='ClinVar benign filter (Yes or No)', default='None')
    parser.add_argument('--remove', type=str, help='Variant filtering out in VCF file (Yes or No)', default='No')
    parser.add_argument('--recover', type=str, help='Variant recovery list at manual review', default='No')
    args = parser.parse_args()

    ngb_filter_format = filter_format(args.input, args.output, args.database, args.format, args.cutoff, args.clinvar, args.remove, args.recover)
    ngb_filter_format.run()
