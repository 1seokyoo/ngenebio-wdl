#!/usr/bin/python

import os
import sys
import argparse
import vcf
from difflib import ndiff

class filter_format():

    def __init__(self, sample_name, normal_name, whitelist_name, hrd_name, mutect_name, output_name, tmp_name):

        self._sample_name = sample_name
        self._normal_vcf = normal_name
        self._whitelist_vcf = whitelist_name
        self._hrd_vcf = hrd_name
        self._mutect_vcf = mutect_name
        self._output_name = output_name
        self._tmp_name = tmp_name
        

    def init_db(self, input_file):

        input_data = vcf.Reader(filename=input_file)

        output_db = {}

        for record in input_data:

            record = self.set_format(record)
            ALT = str(record.ALT).replace("[","").replace("]","").replace(" ","").strip()
            key = '{0}-{1}-{2}-{3}'.format(record.CHROM, record.POS, record.REF, ALT)

            if record.FILTER == [] or record.FILTER == 'PASS':
                record.FILTER = 'PASS'
                info_dic = {'TP_VCF':'NOR', 'FP_VCF': ''}
                if 'TP_RECOVER' in record.INFO:
                    info_dic['TP_RECOVER'] = 'NOR|{0}'.format(record.INFO['TP_RECOVER'][0])

            else:
                info_dic = {'TP_VCF':'', 'FP_VCF': 'NOR'}

                if 'TP_RECOVER' in record.INFO:
                    info_dic['TP_RECOVER'] = 'NOR|{0}'.format(record.INFO['TP_RECOVER'][0])
                if 'FP_FILTER' in record.INFO:
                    info_dic['FP_FILTER'] = 'NOR|{0}'.format(record.INFO['FP_FILTER'][0])
                if 'FP_DB' in record.INFO:
                    info_dic['FP_DB'] = 'NOR|{0}'.format(record.INFO['FP_DB'][0])
                if 'FP_FORMAT' in record.INFO:
                    info_dic['FP_FORMAT'] = 'NOR|{0}'.format(record.INFO['FP_FORMAT'][0])

            record.INFO = info_dic
            
            output_db[key] = record

        return output_db


    def set_format(self, input_data):

        input_info = input_data.INFO
        input_sample = input_data.samples[0]
        
        input_data.FORMAT = 'GT:DP:AD:F1R2:F2R1'
        GT_value = '0/0'
        AD_value = '0,0'
        DP_value = '0' 

        GT_value = input_sample['GT']
        AD_value = input_sample['AD']
        if 'DP' in input_info: # Mutect2
            DP_value = input_info['DP']        
        else: # VarDict
            DP_value = input_sample['DP']

        if 'VARBIAS' in input_info and 'REFBIAS' in input_info: # VarDict
            ref_count = input_info['REFBIAS'].split(':')[0]
            alt_count = input_info['VARBIAS'].split(':')[0]
            F1R2_value = '{0},{1}'.format(ref_count, alt_count)
        else: # Mutect2
            F1R2_value = input_sample['F1R2']
        
        if 'VARBIAS' in input_info and 'REFBIAS' in input_info: # VarDict
            ref_count = input_info['REFBIAS'].split(':')[1]
            alt_count = input_info['VARBIAS'].split(':')[1]
            F2R1_value = '{0},{1}'.format(ref_count, alt_count)
        else: # Mutect2
            F2R1_value = input_sample['F2R1']
        
        format_list = tuple(['GT', 'DP', 'AD', 'F1R2', 'F2R1'])
        output_data = vcf.model.make_calldata_tuple(format_list)(GT=GT_value, DP=DP_value, AD=AD_value, \
                                                                 F1R2=F1R2_value, F2R1=F2R1_value)
        
        input_data.samples[0].sample = self._sample_name
        input_data.samples[0].data = output_data

        return input_data


    def merge_format(self, input_record, caller):

        output_dic = {}
        if caller == 'VDT':
            output_dic['VDT_DP'] = input_record.samples[0]['DP']
            output_dic['VDT_AD'] = input_record.samples[0]['AD']
        elif caller == 'WHI':
            output_dic['WHI_DP'] = input_record.samples[0]['DP']
            output_dic['WHI_AD'] = input_record.samples[0]['AD']
        elif caller == 'HRD':
            output_dic['HRD_DP'] = input_record.samples[0]['DP']
            output_dic['HRD_AD'] = input_record.samples[0]['AD']
        elif caller == 'MUT':
            output_dic['MUT_DP'] = input_record.INFO['DP']
            output_dic['MUT_AD'] = input_record.samples[0]['AD']
            
        return output_dic


    def merge_info(self, input_dic, pre_db, now_db, tag, value):

        output_list = []

        if value in pre_db:

            output_list.append(pre_db[value])

        if value in now_db:
            output_list.append('{0}|{1}'.format(tag, now_db[value][0]))

        if not output_list == []:
            input_dic[value] = ','.join(output_list)

        return input_dic


    def merge_db(self, input_file, variant_db, variant_tag):

        if os.path.getsize(input_file) == 0:
            pass
        
        else:
            input_data = vcf.Reader(filename=input_file)

            for record in input_data:

                ALT = str(record.ALT).replace("[","").replace("]","").replace(" ","").strip()
                key = '{0}-{1}-{2}-{3}'.format(record.CHROM, record.POS, record.REF, ALT)

                ## The variant exist in DB
                if key in variant_db:

                    if record.FILTER == []: # This variant is PASS
                        info_dic = {'FP_VCF': variant_db[key].INFO['FP_VCF']}
                        if variant_db[key].INFO['TP_VCF'] == '': # Old and New variant are PASS
                            info_dic['TP_VCF'] = variant_tag
                        else:
                            info_dic['TP_VCF'] = '|'.join([variant_db[key].INFO['TP_VCF'], variant_tag])
                    else: 
                        info_dic = {'TP_VCF': variant_db[key].INFO['TP_VCF']}
                        if variant_db[key].INFO['FP_VCF'] == '':
                            info_dic['FP_VCF'] = variant_tag
                        else:
                            info_dic['FP_VCF'] = '|'.join([variant_db[key].INFO['FP_VCF'], variant_tag])
                        
                    info_dic = self.merge_info(info_dic, variant_db[key].INFO, record.INFO, variant_tag, 'TP_RECOVER')
                    info_dic = self.merge_info(info_dic, variant_db[key].INFO, record.INFO, variant_tag, 'FP_FILTER')
                    info_dic = self.merge_info(info_dic, variant_db[key].INFO, record.INFO, variant_tag, 'FP_DB')
                    info_dic = self.merge_info(info_dic, variant_db[key].INFO, record.INFO, variant_tag, 'FP_FORMAT')
                    info_dic = self.merge_info(info_dic, variant_db[key].INFO, record.INFO, variant_tag, 'FP_CONSEQUENCE')
    
                    if variant_db[key].FILTER == 'PASS': # Old variant is PASS
                        input_info = self.merge_format(record, variant_tag)
                        info_dic.update(input_info)
                        record = variant_db[key]
                        


                    elif record.FILTER == []: # Old variant is REJECT but new variant is PASS
                        record.FILTER = 'PASS'
                        record = self.set_format(record)
                        if variant_tag == 'MUT':
                            record.FILTER = 'Mutect2'

                    else: # Both variants are REJECT
                        record = variant_db[key]              
                    
                else: # The variant does not exist in DB
                    if record.FILTER == []: # This variant is PASS
                        info_dic = {'TP_VCF':variant_tag, 'FP_VCF': ''}                      
                        if variant_tag == 'MUT':
                            record.FILTER = 'Mutect2'
                        else:                    
                            record.FILTER = 'PASS'
                    else:
                        info_dic = {'TP_VCF':'', 'FP_VCF': variant_tag}

                    record = self.set_format(record)

                    if 'FP_FILTER' in record.INFO:
                        info_dic['FP_FILTER'] = '{0}|{1}'.format(variant_tag, record.INFO['FP_FILTER'][0])
                    if 'FP_DB' in record.INFO:
                        info_dic['FP_DB'] = '{0}|{1}'.format(variant_tag, record.INFO['FP_DB'][0])
                    if 'FP_FORMAT' in record.INFO:
                        info_dic['FP_FORMAT'] = '{0}|{1}'.format(variant_tag, record.INFO['FP_FORMAT'][0])

                record.INFO = info_dic
                variant_db[key] = record

        return variant_db


    def set_type(self, input_db):

        key_list = input_db.keys()
        key_list.sort()

        for key in key_list:
            
            var_ref = input_db[key].REF
            var_alt = str(input_db[key].ALT).replace("[","").replace("]","").replace(" ","").strip()
        
            if len(var_ref) == 1 and len(var_alt) == 1:
                var_type = 'SNV'
            elif len(var_ref) == 1 and len(var_alt) > 1:
                var_type = 'Ins'
            elif len(var_ref) > 1 and len(var_alt) == 1:
                var_type = 'Del'
            elif len(var_ref) > 1 and len(var_alt) > 1:
                var_type = 'Complex'
                                    
            input_db[key].INFO['TYPE'] = var_type

            if input_db[key].INFO['TP_VCF'] == '':
                input_db[key].INFO.pop('TP_VCF')
            if input_db[key].INFO['FP_VCF'] == '':
                input_db[key].INFO.pop('FP_VCF')

        return input_db


    def write_vcf(self, variant_db, output_file, tmp_file):

        ## VCF header File
        vcf_format = vcf.Reader(filename=tmp_file)

        ## Set vcf writing
        vcf_output = open(output_file, 'w')
        vcf_writer = vcf.Writer(vcf_output, vcf_format)

        ## Write variant data
        key_list = variant_db.keys()
        key_list.sort()

        for key in key_list:
            vcf_writer.write_record(variant_db[key])

        vcf_output.close()


    def run(self):

        ## Set variant database used by Normal vcf File
        output_db = self.init_db(self._normal_vcf)

        ## Merge Whitelist vcf File
        output_db = self.merge_db(self._whitelist_vcf, output_db, 'WHI')

        ## Merge HRD vcf File
        output_db = self.merge_db(self._hrd_vcf, output_db, 'HRD')
        
        ## Merge Mutect vcf File
        output_db = self.merge_db(self._mutect_vcf, output_db, 'MUT')

        ## Set Variant Type
        output_db = self.set_type(output_db)

        ## Write VCF
        self.write_vcf(output_db, self._output_name, self._tmp_name)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Merge vcf File')
    parser.add_argument('-s', '--sample', type=str, help='Sample name')
    parser.add_argument('-n', '--normal', type=str, help='Normal vcf file name')
    parser.add_argument('-w', '--whitelist', type=str, help='Whitelist vcf file name')
    parser.add_argument('-r', '--hrd', type=str, help='HRD vcf file name')
    parser.add_argument('-m', '--mutect', type=str, help='Mutect2 vcf file name')
    parser.add_argument('-o', '--output', type=str, help='VCF output file name')
    parser.add_argument('-t', '--temp', type=str, help='VCF format file name')
    args = parser.parse_args()

    ngb_merge_vcf = filter_format(args.sample, args.normal, args.whitelist, args.hrd, args.mutect, args.output, args.temp)
    ngb_merge_vcf.run()
