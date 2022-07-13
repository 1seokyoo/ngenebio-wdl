#!/usr/bin/env python
from __future__ import print_function

__author__="Changbum Hong"
__email__="cb.hong@ngenebio.com"
__modname__="clip_finder.py"
__created__="2016_06_03"
__company__="ngenebio"

import re
import sys
import argparse
import logging
from sys import argv

"""
usage: python /NGENEBIO/workflow/ngb_workflow/utils/amplicon_check.py amplicon=/NGENEBIO/workflow_out/200/assay_reference/Qiagen.BRCA.bed input_sam=/NGENEBIO/workflow_out/200/data/basecall/alignment/HCT-15-1_final.sam sample_name=HCT-15-1 
=/NGENEBIO/workflow_out/200/data/stat
"""

"""
sam,bam cigar parser
"""
cigarPattern = '([0-9]+[MIDNSHP])'
cigarSearch = re.compile(cigarPattern)
atomicCigarPattern = '([0-9]+)([MIDNSHP])'
atomicCigarSearch = re.compile(atomicCigarPattern)
softclipPattern = '(^[0-9]+[S])'
softclipSearch = re.compile(softclipPattern)

def main(args):

   if len(sys.argv) < 1:
      raise Exception('Error running Analysis pipeline: Invalid arguments ({0})'.format(str(sys.argv)))

   for arg in sys.argv:
      try:
         name, val = arg.split("=")
         if name == "amplicon":
            inFile = val
         elif name == "input_sam":
            samFile = val
         elif name == "sample_name":
            sample_name = val
         elif name == "output_bed":
            output_bed = val
         elif name == "output_txt":
            output_txt = val
      except:
         pass

   """
   parse the amplicon(primer) inforamtion
   """
   warn_bed=open(output_bed,'w')
   clip_stat=open(output_txt,'w')
   data = open(inFile, 'r')
   amp={}
   amp2={}
   match_am={}
   k=[]
   for line in data:
      ampchr=line.split('\t')[0]
      ampstart=line.split('\t')[1]
      ampend=line.split('\t')[2]
      ampname=line.split('\t')[3]
      ampkey='%s_%s'%(ampchr,ampstart)
      amp[ampkey]=[ampname]
      match_am[ampname]={}


      amp2chr=line.split('\t')[0]
      amp2start=line.split('\t')[1]
      amp2end=line.split('\t')[2]
      amp2name=line.split('\t')[3]
      amp2key='%s_%s_%s'%(amp2chr,amp2start,amp2end)
      amp2[ampname]=[amp2key]

   """
   amp2 {'BRCA1.CDS.9.line.6.chr17.41243401.41246927_tile_18': ['chr17_41245531_41245751']}
   """

   qsPos=0
   qePos=0
   qLen=0
   sam = open(samFile,'r')
   tot_read=0
   tot_matchnum=0
   tot_front_soft=0
   tot_end_soft=0
   tot_both_soft=0
   mat_amplicons=[]
   debug="no"
   warn=[]
   soft_clip_info={}
   #read sam file
   for line in sam:
      tot_read=tot_read+1
      samreadname=line.split("\t")[0]
      samchr=line.split("\t")[2]
      samstart=line.split("\t")[3]
      samcigr=line.split("\t")[5]
      samend=0
      samlen=0

      if True:
      #if samcigr.find("H") != -1 or samcigr.find("S") != -1:
         readCigarOps=extractCigarOps(samcigr)
         if debug=="yes":
            print(">>>>>>>>>>>>%s\t%s-%s<<<<<<<<<<<<"%(samreadname,samchr,samstart))
            print(samcigr)
            for a in readCigarOps:
               print("%s\t%s"%(a.op,a.length))
         readQueryPos = calcQueryPosFromCigar(readCigarOps)
         if debug=="yes":
            print("read-position(start-end-len-clipsize): %s %s %s %s \tsoft-clipping(front-end-both): %s %s %s"%(readQueryPos.qsPos,readQueryPos.qePos,readQueryPos.qLen,readQueryPos.softLen,readQueryPos.f_soft,readQueryPos.e_soft,readQueryPos.b_soft))

         #count soft clipped read
         if readQueryPos.f_soft >=1:
            tot_front_soft = tot_front_soft +1
         elif readQueryPos.e_soft >=1:
            tot_end_soft = tot_end_soft +1
         if readQueryPos.b_soft >=1:
            tot_fornt_soft = tot_front_soft -1
            tot_end_soft = tot_end_soft -1
            tot_both_soft = tot_both_soft +1
         if debug=="yes":
            print("%s %s %s"%(tot_front_soft,tot_end_soft,tot_both_soft))
         orgsamstart=samstart
         samstart=int(orgsamstart)-readQueryPos.softLen-1
         samend=readQueryPos.qLen+int(orgsamstart)-readQueryPos.softLen-1
         samlen=readQueryPos.qLen

         samkey='%s_%s'%(samchr,samstart)
         samkey_plus_1='%s_%s'%(samchr,int(samstart)+1)
         samkey_minus_1='%s_%s'%(samchr,int(samstart)-1)

         if debug=="yes":
            print("genomic location: %s:%s-%s"%(samchr,samstart,samend))

         matchnum_=0
         hhh=0   
         for j in range(50):
            #print j
            orgsamstart = int(samstart)


            samstart_plus = int(samstart)
            samstart_plus = int(samstart_plus)+j
            samkey='%s_%s'%(samchr,samstart_plus)

            matchnum=0
            if amp.get(samkey):
               amp1=''.join(amp[samkey])
               ###print("+ Matched Amplicon %d: %s %s"%(matchnum,amp1,amp2[amp1]))
               mampstart = str(amp2[amp1]).split("_")[1]
               mampend = str(amp2[amp1]).split("_")[2].strip('\']')
               diffstart = int(mampstart)-int(orgsamstart)
               diffend = int(mampend)-int(samend)
               ###print("diff start: %d"%(diffstart))
               ###print("diff end: %d"%(diffend))
               mat_amplicons.append(str(amp1))
               matchnum_ = 1

            samstart_minus = int(samstart)-1
            samstart_minus = int(samstart_minus)-j
            samkey='%s_%s'%(samchr,samstart_minus)
            if amp.get(samkey):
               amp1=''.join(amp[samkey])
               if debug=="yes":
                  print("Matched Amplicon %d: %s %s"%(matchnum,amp1,amp2[amp1]))
               mampstart = str(amp2[amp1]).split("_")[1]
               mampend = str(amp2[amp1]).split("_")[2].strip('\']')
               diffstart = int(mampstart)-int(orgsamstart)
               diffend = int(mampend)-int(samend)
               ###print("diff start: %d"%(diffstart))
               ###print("diff end: %d"%(diffend))
               mat_amplicons.append(str(amp1))
               matchnum_ = 1
               if str(amp1) in soft_clip_info:
                  update_clip=soft_clip_info[str(amp1)]
                  f_soft=update_clip.split("-")[0]
                  e_soft=update_clip.split("-")[1]
                  b_soft=update_clip.split("-")[2]
                  l_f_soft=int(f_soft)+readQueryPos.f_soft
                  l_e_soft=int(e_soft)+readQueryPos.e_soft
                  l_b_soft=int(b_soft)+readQueryPos.b_soft
                  soft_clip_info[str(amp1)]="%s-%s-%s"%(l_f_soft,l_e_soft,l_b_soft)
               else:
                  soft_clip_info[str(amp1)]="%s-%s-%s"%(readQueryPos.f_soft,readQueryPos.e_soft,readQueryPos.b_soft)

         tot_matchnum=tot_matchnum+matchnum_
         if debug=="yes":
            #print(mat_amplicons)
            print(soft_clip_info)

   # find read in amplicon position
   clip_stat.write("Amplicon\tCount\tPercentage\tFront Clip\tEnd Clip\tBoth Clip\tChr\tStart\tEnd\n")
   logging.warn("Amplicon\tCount\tPercentage\tFront Clip\tEnd Clip\tBoth Clip\tChr\tStart\tEnd\n")
   for k in match_am.keys():
      amp=str(amp2[k]).replace("]","").replace("[","").replace("'","")
      amp_chr=amp.split("_")[0]
      amp_start=amp.split("_")[1]
      amp_end=amp.split("_")[2]
      amp_count=str(mat_amplicons.count(k)).replace("\n","")

      if k in soft_clip_info:
         soft_info=soft_clip_info[k]
         front_soft=soft_info.split("-")[0]
         end_soft=soft_info.split("-")[1]
         both_soft=soft_info.split("-")[2]
      else:
         front_soft=0
         end_soft=0
         both_soft=0
      if amp_count != "0":
         percentage_clip=(float(front_soft)+float(end_soft)+float(both_soft))/float(amp_count)*100
      else:
         percentage_clip=0
      if percentage_clip >= 40:
         extend_start=int(amp_start)-100
         extend_end=int(amp_end)+100         
         warn_bed.write("%s\t%s\t%s\t%s\t+\t%s\n"%(amp_chr,extend_start,extend_end,k,percentage_clip))
         warn.append("%s\t%s\t%d%%\t%s\t%s\t%s\t%s\t%d\t%d\n"%(k,amp_count,round(percentage_clip,2),front_soft,end_soft,both_soft,amp_chr,int(amp_start),int(amp_end)))
      clip_stat.write("%s\t%s\t%d%%\t%s\t%s\t%s\t%s\t%d\t%d\n"%(k,amp_count,round(percentage_clip,2),front_soft,end_soft,both_soft,amp_chr,int(amp_start),int(amp_end)))
      logging.warn("%s\t%s\t%d%%\t%s\t%s\t%s\t%s\t%d\t%d\n"%(k,amp_count,round(percentage_clip,2),front_soft,end_soft,both_soft,amp_chr,int(amp_start),int(amp_end)))
   print("=======================================================")
   print("total_reads\t%s\t0\t0\t0"%(tot_read))
   print("total_matched_reads\t%s\t0\t0\t0"%(tot_matchnum))
   percentage_front_soft=float(tot_front_soft)/float(tot_read)*100
   percentage_end_soft=float(tot_end_soft)/float(tot_read)*100
   percentage_both_soft=float(tot_both_soft)/float(tot_read)*100
   print("front soft clipping read: %s %f%%"%(tot_front_soft,percentage_front_soft))
   print("end soft clipping read: %s %f%%"%(tot_end_soft,percentage_end_soft))
   print("both soft clipping read: %s %f%%"%(tot_both_soft,percentage_both_soft))
   print(warn)
   warn_bed.close()
   clip_stat.close()

class queryPos (object):
    """
    struct to store the start and end positions of query CIGAR operations
    """
    def __init__(self, qsPos, qePos, qLen, softLen, f_soft, e_soft, b_soft):
        self.qsPos = int(qsPos)
        self.qePos = int(qePos)
        self.qLen  = int(qLen)
	self.softLen = int(softLen)
        self.f_soft = int(f_soft)
        self.e_soft = int(e_soft)
        self.b_soft = int(b_soft)

class cigarOp (object):
    """
    sturct to store a discrete CIGAR operations
    """
    def __init__(self, opLength, op):
        self.length = int(opLength)
        self.op     = op

def extractCigarOps(cigar):
        if (cigar == "*"):
                cigarOps = []
        else:
           cigarOpStrings = cigarSearch.findall(cigar)
           cigarOps = []
           for opString in cigarOpStrings:
              cigarOpList = atomicCigarSearch.findall(opString)
              # "struct" for the op and it's length
              cigar = cigarOp(cigarOpList[0][0], cigarOpList[0][1])
              # add to the list of cigarOps
              cigarOps.append(cigar)
        return(cigarOps)

def calcQueryPosFromCigar(cigarOps):
        qsPos = 0
        qePos = 0
        qLen  = 0
        softLen = 0
        f_soft = 0
        e_soft = 0
        b_soft = 0 
        # if first op is a H, need to shift start position
        # the opPosition counter sees if the for loop is looking at the first index of the cigar object
        opPosition = 0
        for cigar in cigarOps:
                #soft-clipping in the front read
                if opPosition == 0 and (cigar.op == 'H' or cigar.op == 'S'):
                        qsPos += cigar.length
                        qePos += cigar.length
                        qLen  += cigar.length
			softLen += cigar.length
                        #soft clip more than 5bp
                        if cigar.length >= 1:
                           f_soft += 1
                elif opPosition > 0 and (cigar.op == 'H' or cigar.op == 'S'):
                        qLen  += cigar.length
                        #soft clip more than 5bp
                        if cigar.length >= 1:
			   e_soft += 1
                elif cigar.op == 'M' or cigar.op == 'I':
                        qePos += cigar.length
                        qLen  += cigar.length
                        opPosition += 1
                if f_soft >= 1 and e_soft >= 1:
                   b_soft += 1
        d = queryPos(qsPos, qePos, qLen,softLen,f_soft,e_soft,b_soft);
        return d

if __name__=="__main__":
   main(argv)
