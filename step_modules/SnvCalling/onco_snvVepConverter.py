#!/usr/bin/python

import sys
import argparse
import vcf
import pandas
from vcf.parser import _vcf_metadata_parser as vcf_parser

class vep_converter():

    def __init__(self, input_name, output_name, transcript_file):

        self._input_name = input_name
        self._output_name = output_name
        
        self._transcript_dic = {}
        with open(transcript_file, 'r') as transcript_data:
            for input_line in transcript_data:
                
                # Remove header
                if input_line[0] == '#':
                    pass
                else:
                    input_list = input_line.strip().split('\t')
                    ensembl = input_list[0]

                    if len(input_list) == 2:
                        refseq = 'NONE'
                    elif input_list[2] == '':
                        refseq = 'NONE'
                    else:
                        refseq = input_list[2]

                    self._transcript_dic[ensembl] = refseq
            
        # Consequence condition
        self._exonic_region = ["Frame_Shift_Del", "Frame_Shift_Ins",   
                               "In_Frame_Del", "In_Frame_Ins",      
                               "Missense_Mutation", 
                               "Nonsense_Mutation", 
                               "Nonstop_Mutation",  
                               "Splice_Site"]

        self._conseq_dict = {'transcript_ablation' : ['Splice_Site', 1], # A feature ablation whereby the deleted region includes a transcript feature
                             'exon_loss_variant' : ['Splice_Site', 1], # A sequence variant whereby an exon is lost from the transcript
                             'splice_donor_variant' : ['Splice_Site', 2], # A splice variant that changes the 2 base region at the 5' end of an intron
                             'splice_acceptor_variant' : ['Splice_Site', 2], # A splice variant that changes the 2 base region at the 3' end of an intron
                             'stop_gained' : ['Nonsense_Mutation', 3], # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
                             'frameshift_variant' : ['Frame_Shift', 3], # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
                             'stop_lost' : ['Nonstop_Mutation', 3], # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
                             'start_lost' : ['Translation_Start_Site', 4], # A codon variant that changes at least one base of the canonical start codon
                             'initiator_codon_variant' : ['Translation_Start_Site', 4], # A codon variant that changes at least one base of the first codon of a transcript
                             'disruptive_inframe_insertion' : ['In_Frame_Ins', 5], # An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon
                             'disruptive_inframe_deletion' : ['In_Frame_Del', 5], # An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon
                             'inframe_insertion' : ['In_Frame_Ins', 5], # An inframe non synonymous variant that inserts bases into the coding sequence
                             'inframe_deletion' : ['In_Frame_Del', 5], # An inframe non synonymous variant that deletes bases from the coding sequence
                             'protein_altering_variant' : ['Protein_Alter', 5], # A sequence variant which is predicted to change the protein encoded in the coding sequence
                             'missense_variant' : ['Missense_Mutation', 6], # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
                             'conservative_missense_variant' : ['Missense_Mutation', 6], # A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
                             'rare_amino_acid_variant' : ['Missense_Mutation', 6], # A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
                             'transcript_amplification' : ['Intron', 7], # A feature amplification of a region containing a transcript
                             'splice_region_variant' : ['Splice_Region', 8], # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
                             'stop_retained_variant' : ['Silent', 9], # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
                             'synonymous_variant' : ['Silent', 9], # A sequence variant where there is no resulting change to the encoded amino acid
                             'incomplete_terminal_codon_variant' : ['Silent', 10], # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
                             'coding_sequence_variant' : ['Missense_Mutation', 11], # A sequence variant that changes the coding sequence
                             'mature_miRNA_variant' : ['RNA', 11], # A transcript variant located with the sequence of the mature miRNA
                             'exon_variant' : ['RNA', 11], # A sequence variant that changes exon sequence
                             '5_prime_UTR_variant' : ['5\'UTR', 12], # A UTR variant of the 5' UTR
                             '5_prime_UTR_premature_start_codon_gain_variant' : ['5\'UTR', 12], # snpEff-specific effect, creating a start codon in 5' UTR
                             '3_prime_UTR_variant' : ['3\'UTR', 12], # A UTR variant of the 3' UTR
                             'non_coding_exon_variant' : ['RNA', 13], # A sequence variant that changes non-coding exon sequence
                             'non_coding_transcript_exon_variant' : ['RNA', 13], # snpEff-specific synonym for non_coding_exon_variant
                             'non_coding_transcript_variant' : ['RNA', 14], # A transcript variant of a non coding RNA gene
                             'nc_transcript_variant' : ['RNA', 14], # A transcript variant of a non coding RNA gene (older alias for non_coding_transcript_variant)
                             'intron_variant' : ['Intron', 14], # A transcript variant occurring within an intron
                             'intragenic_variant' : ['Intron', 14], # A variant that occurs within a gene but falls outside of all transcript 
                             'INTRAGENIC' : ['Intron', 14], # snpEff-specific synonym of intragenic_variant
                             'NMD_transcript_variant' : ['Silent', 15], # A variant in a transcript that is the target of NMD
                             'upstream_gene_variant' : ['5\'Flank', 16], # A sequence variant located 5' of a gene
                             'downstream_gene_variant' : ['3\'Flank', 16], # A sequence variant located 3' of a gene
                             'TFBS_ablation' : ['Targeted_Region', 17], # A feature ablation whereby the deleted region includes a transcription factor binding site
                             'TFBS_amplification' : ['Targeted_Region', 17], # A feature amplification of a region containing a transcription factor binding site
                             'TF_binding_site_variant' : ['Intergenic', 17], # A sequence variant located within a transcription factor binding site
                             'regulatory_region_ablation' : ['Targeted_Region', 17], # A feature ablation whereby the deleted region includes a regulatory region
                             'regulatory_region_amplification' : ['Targeted_Region', 17], # A feature amplification of a region containing a regulatory region
                             'regulatory_region_variant' : ['Intergenic', 17], # A sequence variant located within a regulatory region
                             'regulatory_region' : ['Intergenic', 17], # snpEff-specific effect that should really be regulatory_region_variant
                             'feature_elongation' : ['Targeted_Region', 18], # A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
                             'feature_truncation' : ['Targeted_Region', 18], # A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
                             'intergenic_variant' : ['Intergenic', 19], # A sequence variant located in the intergenic region, between genes
                             'intergenic_region' : ['Intergenic', 19], # snpEff-specific effect that should really be intergenic_variant
                             'ETC' : ['None', 20] # none 
                             }

        self._biotype_dict = {'protein_coding' : 1, # Contains an open reading frame (ORF)
                              'LRG_gene' : 2, # Gene in a "Locus Reference Genomic" region known to have disease-related sequence variations
                              'IG_C_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                              'IG_D_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                              'IG_J_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                              'IG_LV_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                              'IG_V_gene' : 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
                              'TR_C_gene' : 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
                              'TR_D_gene' : 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
                              'TR_J_gene' : 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
                              'TR_V_gene' : 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
                              'miRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'snRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'snoRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'ribozyme' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'tRNA' : 3, #Added by Y. Boursin
                              'sRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'scaRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'rRNA' : 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'lincRNA' : 3, # Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily conserved, intergenic regions
                              'bidirectional_promoter_lncrna' : 3, # A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand
                              'bidirectional_promoter_lncRNA' : 3, # A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand
                              'known_ncrna' : 4,
                              'vaultRNA' : 4, # Short non coding RNA genes that form part of the vault ribonucleoprotein complex
                              'macro_lncRNA' : 4, # unspliced lncRNAs that are several kb in size
                              'Mt_tRNA' : 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'Mt_rRNA' : 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'antisense' : 5, # Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand
                              'antisense_RNA' : 5, # Alias for antisense (Y. Boursin)
                              'sense_intronic' : 5, # Long non-coding transcript in introns of a coding gene that does not overlap any exons
                              'sense_overlapping' : 5, # Long non-coding transcript that contains a coding gene in its intron on the same strand
                              '3prime_overlapping_ncrna' : 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
                              '3prime_overlapping_ncRNA' : 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
                              'misc_RNA' : 5, # Non-coding RNA predicted using sequences from RFAM and miRBase
                              'non_coding' : 5, # Transcript which is known from the literature to not be protein coding
                              'regulatory_region' : 6, # A region of sequence that is involved in the control of a biological process
                              'disrupted_domain' : 6, # Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain
                              'processed_transcript' : 6, # Doesn't contain an ORF
                              'TEC' : 6, # To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies
                              'TF_binding_site' : 7, # A region of a nucleotide molecule that binds a Transcription Factor or Transcription Factor complex
                              'CTCF_binding_site' :7, # A transcription factor binding site with consensus sequence CCGCGNGGNGGCAG, bound by CCCTF-binding factor
                              'promoter_flanking_region' : 7, # A region immediately adjacent to a promoter which may or may not contain transcription factor binding sites
                              'enhancer' : 7, # A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter
                              'promoter' : 7, # A regulatory_region composed of the TSS(s) and binding sites for TF_complexes of the basal transcription machinery
                              'open_chromatin_region' : 7, # A DNA sequence that in the normal state of the chromosome corresponds to an unfolded, un-complexed stretch of double-stranded DNA
                              'retained_intron' : 7, # Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants
                              'nonsense_mediated_decay' : 7, # If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD
                              'non_stop_decay' : 7, # Transcripts that have polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation
                              'ambiguous_orf' : 7, # Transcript believed to be protein coding, but with more than one possible open reading frame
                              'pseudogene' : 8, # Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following
                              'processed_pseudogene' : 8, # Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome
                              'polymorphic_pseudogene' : 8, # Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated
                              'retrotransposed' : 8, # Pseudogene owing to a reverse transcribed and re-inserted sequence
                              'translated_processed_pseudogene' : 8, # Pseudogenes that have mass spec data suggesting that they are also translated
                              'translated_unprocessed_pseudogene' : 8, # Pseudogenes that have mass spec data suggesting that they are also translated
                              'transcribed_processed_pseudogene' : 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
                              'transcribed_unprocessed_pseudogene' : 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
                              'transcribed_unitary_pseudogene' : 8, #Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
                              'unitary_pseudogene' : 8, # A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species
                              'unprocessed_pseudogene' : 8, # Pseudogene that can contain introns since produced by gene duplication
                              'Mt_tRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                              'tRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                              'snoRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                              'snRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                              'scRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                              'rRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                              'misc_RNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                              'miRNA_pseudogene' : 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
                              'IG_C_pseudogene' : 8, # Inactivated immunoglobulin gene
                              'IG_D_pseudogene' : 8, # Inactivated immunoglobulin gene
                              'IG_J_pseudogene' : 8, # Inactivated immunoglobulin gene
                              'IG_V_pseudogene' : 8, # Inactivated immunoglobulin gene
                              'TR_J_pseudogene' : 8, # Inactivated immunoglobulin gene
                              'TR_V_pseudogene' : 8, # Inactivated immunoglobulin gene
                              'artifact' : 9, # Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
                              'ETC' : 10 # none
                              }


    def pick_cons(self, input_list):

        best_cons = 'intergenic_variant'
        best_score = 19
        
        if input_list == []:
            pass

        else:
            for input_cons in input_list:
                if input_cons in self._conseq_dict:
                    one_cons, one_score = self._conseq_dict[input_cons]
                else:
                    one_cons, one_score = self._conseq_dict['ETC']

                if one_cons == 'Frame_Shift':
                    if self._var_type == 'Del':
                        one_cons = 'Frame_Shift_Del'
                    elif self._var_type == 'Ins':
                        one_cons = 'Frame_Shift_Ins'
                elif one_cons == 'Protein_Alter':
                    if self._var_inframe == 'FALSE' and self._var_type == 'Del':
                        one_cons = 'Frame_Shift_Del'
                    elif self._var_inframe == 'FALSE' and self._var_type == 'Ins':
                        one_cons = 'Frame_Shift_Ins'
                    elif self._var_inframe == 'TRUE' and self._var_type == 'Del':
                        one_cons = 'In_Frame_Del'
                    elif self._var_inframe == 'TURE' and self._var_type == 'Ins':
                        one_cons = 'In_Frame_Ins'

                if one_score < best_score:
                    best_cons = one_cons
                    best_score = one_score

        return best_cons, best_score


    def filter_cons(self, input_record, input_cons):

        input_filter = input_record.FILTER
        input_alt = str(input_record.ALT).replace("[","").replace("]","").replace(" ","").strip()

        if input_cons in self._exonic_region:
            pass

        elif input_record.CHROM == 'chr7' and input_record.POS > 116411708 and input_record.POS <= 116411903: # Met Exon 14 Skipping
            pass

        elif input_record.CHROM == 'chr5' and input_record.POS == 1295228 and input_alt == 'A': # TERT Promoter - C228T
            pass

        elif input_record.CHROM == 'chr5' and input_record.POS == 1295250 and input_alt == 'A': # TERT Promoter - C250T
            pass

        else:
            input_record.INFO['FP_CONSEQUENCE'] = input_cons
            input_record.FILTER = 'REJECT'

        return input_record


    def sort_csq(self, input_list):

        column_names = ['#', 'gene', 'biotype', 'consequence', 'length', 'transcript', 'vep']
        raw_df = pandas.DataFrame(columns=column_names)

        for num in range(len(input_list)):
            csq_list = str(input_list[num]).strip().split("|")
            csq_gene = csq_list[3].strip()
            csq_enst = csq_list[6].strip()
            if csq_list[33] == '':
                csq_refseq = 'NONE'
            else:
                csq_refseq = csq_list[33].strip()

            # Check biotype score
            csq_biotype = csq_list[7].strip()
            if csq_biotype in self._biotype_dict:
                bio_score = self._biotype_dict[csq_biotype]
            else:
                bio_score = self._biotype_dict['ETC']

            # Convert consequence
            csq_cons = csq_list[1].strip().split('&')
            cons_value, cons_score = self.pick_cons(csq_cons)

            # Check transcript length
            if not csq_list[12] == '':
                tr_length = csq_list[12].strip().split("/")[1]
            else:
                tr_length = 0

            # Set custom transcript to default transcript
            if csq_enst in self._transcript_dic:
                transcript_value = 0
                csq_list[33] = self._transcript_dic[csq_enst]
            else:
                transcript_value = 1

            # Set VEP canonical transcript to default transcript
            if csq_list[26] == 'YES':
                vep_value = 0
            else:
                vep_value = 1

            # Convert CSQ data
            csq_list[1] = cons_value
            csq_list[26] = 'NO'
            input_list[num] = "|".join(csq_list)

            # Set dataframe
            raw_df.loc[num] = [num, csq_gene, bio_score, cons_score, tr_length, transcript_value, vep_value]

        # Sort by bio_score, cons_score, tr_length        
        output_df = raw_df.sort_values(['biotype', 'consequence', 'length'])

        return input_list, output_df

    
    def select_transcript(self, input_list, input_df):

        # Set the highest priority gene symbol
        priority_gene = input_df.iloc[0]['gene']
        same_gene_df = input_df.loc[input_df['gene'] == priority_gene]
        other_gene_df = input_df.loc[input_df['gene'] != priority_gene]

        # Same the highest priority gene
        ## Check custom canonical transcript
        if not same_gene_df.loc[same_gene_df['transcript'] == 0].empty:
            default_num = same_gene_df.loc[same_gene_df['transcript'] == 0]['#'].values[0]

        ## Check vep canonical transcript
        elif not same_gene_df.loc[same_gene_df['vep'] == 0].empty:
            default_num = same_gene_df.loc[same_gene_df['vep'] == 0]['#'].values[0]

        # Other the highest priority gene
        ## Check custom canonical transcript
        elif not other_gene_df.loc[other_gene_df['transcript'] == 0].empty:
            default_num = other_gene_df.loc[other_gene_df['transcript'] == 0]['#'].values[0]

        ## Check vep canonical transcript
        elif not other_gene_df.loc[other_gene_df['vep'] == 0].empty:
            default_num = other_gene_df.loc[other_gene_df['vep'] == 0]['#'].values[0]
        
        else:
            default_num = 0
    
        csq_list = str(input_list[default_num]).strip().split("|")
        default_cons = csq_list[1]
        csq_list[26] = 'YES'
        input_list[default_num] = "|".join(csq_list)

        return input_list, default_cons


    def run(self):

        # Read VCF File
        self._input_file = open(self._input_name, 'r')
        self._vcf_reader = vcf.Reader(self._input_file)

        # Write VCF File
        self._output_file = open(self._output_name, 'w')
        self._vcf_writer = vcf.Writer(self._output_file, self._vcf_reader)

        # Filtering
        for record in self._vcf_reader:

            variant_name="GRCh37-%s-%s-%s-%s"%(record.CHROM, record.POS, record.REF, str(record.ALT).replace("[","").replace("]","").replace(" ","").strip()) 

            # Check Inframe variant
            ref_length = len(record.REF)
            alt_length = len(str(record.ALT).replace("[","").replace("]","").replace(" ","").strip())
            length_diff = abs (ref_length - alt_length) % 3
            if length_diff == 0:
                self._var_inframe = 'FALSE' # Not inframe variant
            else:
                self._var_inframe = 'TRUE' # Inframe variant

            # Check variant type
            if ref_length == 1 and alt_length == 1:
                self._var_type = 'SNV'
            elif ref_length > alt_length:
                self._var_type = 'Del'
            elif ref_length < alt_length:
                self._var_type = 'Ins'
            elif ref_length > alt_length:
                self._var_type = 'Del'
            elif ref_length < alt_length:
                self._var_type = 'Ins'
            else:
                self._var_type = 'Complex'

            # CSQ parsing (VEP)

            if 'CSQ'in record.INFO:
                vep_list = []
                vep_csq = record.INFO["CSQ"]

                # Sort CSQ
                vep_csq, csq_df = self.sort_csq(vep_csq)

                # Select default transcript
                vep_csq, vep_cons = self.select_transcript(vep_csq, csq_df)

                record.INFO["CSQ"] = ",".join(vep_csq)

            else:
                vep_cons = 'intergenic_variant'

            record.INFO['NGB'] = vep_cons
            
            record = self.filter_cons(record, vep_cons)

            self._vcf_writer.write_record(record)

        # VCF File Close
        self._input_file.close()
        self._output_file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='VCF Format Filtering')
    parser.add_argument('--input', type=str, help='VCF input file name', default='None')
    parser.add_argument('--output', type=str, help='VCF output file name', default='None')
    parser.add_argument('--custom_enst', type=str, help='Custom canonical transcript name', default='/NGENEBIO/workflow/assay_reference/isoform_overrides_uniprot')
    args = parser.parse_args()

    ngb_vep_convert = vep_converter(args.input, args.output, args.custom_enst)
    ngb_vep_convert.run()
