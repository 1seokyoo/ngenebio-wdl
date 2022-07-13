import sys
import pysam
import re

def make_amp_dic(amp_bed):
    amp_info_dic = { "BRCA1":{}, "BRCA2":{}, "check_amp_f":{}, "check_amp_r":{} }
    amp_name_list_br1 = []
    amp_name_list_br2 = []
    with open(amp_bed) as bed:
        for line in bed:
            if line.split()[0] == "chr13":
                amp_info_dic["BRCA2"][line.split()[3]] = [line.split()[0],line.split()[1],line.split()[2],str(int(line.split()[2])-int(line.split()[1])),line.split()[3]]
                amp_name_list_br2.append(line.split()[3])
            else:    
                amp_info_dic["BRCA1"][line.split()[3]] = [line.split()[0],line.split()[1],line.split()[2],str(int(line.split()[2])-int(line.split()[1])),line.split()[3]]
                amp_name_list_br1.append(line.split()[3])
            amp_info_dic["check_amp_f"][str(line.split()[1])] = line.split()[3]
            amp_info_dic["check_amp_r"][str(line.split()[2])] = line.split()[3]

    return amp_info_dic, amp_name_list_br1, amp_name_list_br2

def find_amp(read, chrom, amp_dic, amp_name_list, amp_count_dic):
    if chrom == 0:
        gene = "BRCA2"
    else:
        gene = "BRCA1"
    if read.rname == chrom and read.pos >= int(amp_dic[gene][amp_name_list[0]][1]) and read.pos <= int(amp_dic[gene][amp_name_list[-1]][2]):
        if str(read.pos) in amp_dic["check_amp_f"]:
            if amp_dic["check_amp_f"][str(read.pos)] in amp_count_dic:
                amp_count_dic[amp_dic["check_amp_f"][str(read.pos)]] += 1
            else:
                amp_count_dic[amp_dic["check_amp_f"][str(read.pos)]] = 1
        
        elif str(read.reference_end) in amp_dic["check_amp_r"]:
            if amp_dic["check_amp_r"][str(read.reference_end)] in amp_count_dic:
                amp_count_dic[amp_dic["check_amp_r"][str(read.reference_end)]] += 1
            else:
                amp_count_dic[amp_dic["check_amp_r"][str(read.reference_end)]] = 1
        else:
            cigar = read.cigarstring
            cigar_d = re.compile('([0-9]+)D')
            del_base = sum([ int(x) for x in cigar_d.findall(cigar) ])
            for amp in amp_name_list:
                if read.pos - (del_base/2) <= int(amp_dic[gene][amp][1]) and int(amp_dic[gene][amp][2]) <= read.reference_end + (del_base/2):
                    if amp in amp_count_dic:
                        amp_count_dic[amp] += 1
                    else:
                        amp_count_dic[amp] = 1
                    break

    return amp_count_dic


def read_bam(bam_file, amp_info_dic, amp_name_list_br1, amp_name_list_br2):
    bam = pysam.AlignmentFile(bam_file,"rb")

    total = 0
    num = 0
    amplicon_dic = {}

    for read in bam.fetch(until_eof=True):
        total += 1
        if read.is_secondary: continue
        if read.is_paired: continue
        if read.is_unmapped: continue
        if read.rname != 0 and read.rname != 1: continue 
        if read.rname == 0:
            amp_name_list = amp_name_list_br2
        else:
            amp_name_list = amp_name_list_br1
        amplicon_dic = find_amp(read, read.rname, amp_info_dic, amp_name_list, amplicon_dic) 
     
    return amplicon_dic    

amp_info_dic, amp_name_list_br1, amp_name_list_br2 = make_amp_dic(sys.argv[1])
amplicon_dic = read_bam(sys.argv[2], amp_info_dic, amp_name_list_br1, amp_name_list_br2)

total_amplicon_count = sum(amplicon_dic.values())

amplicon_dic_300x = {}
amp_300x_total = 0
amp_name_list = amp_name_list_br2 + amp_name_list_br1
for i in amplicon_dic:
    convert_300x = round(int(amplicon_dic[i]) / float(total_amplicon_count) * len(amp_name_list) * 300 , 0)
    amplicon_dic_300x[i] = convert_300x
    amp_300x_total += convert_300x

for i in amp_name_list:
    if i not in amplicon_dic:
        amplicon_dic[i] = 0
    if i not in amplicon_dic_300x:
        amplicon_dic_300x[i] = 0

with open(sys.argv[3], 'w') as output_file:
    output_file.write('\t'.join(["amplicon_name", "raw_count", "300x_count"]) + '\n')
    
    for i in amp_name_list:
        output_file.write('\t'.join([str(i), str(amplicon_dic[i]), str(amplicon_dic_300x[i])]) + '\n')
 
