import pandas as pd
import argparse
import warnings
warnings.filterwarnings('ignore')

class make_cnv_file() :

    def __init__(self, input_file, domain_file, output_name, loss_cutoff, gain_cutoff, exon_threshold) :
        
        self.input = pd.read_table(input_file, sep='\t')
        self.domain = pd.read_table(domain_file, sep='\t', header=None)
        self.output_name = output_name
        self.Lcut = float(loss_cutoff)
        self.Gcut = float(gain_cutoff)
        self.Ecut = float(exon_threshold)

        #cnv terms
        self.normal = 'Normal'
        self.loss = 'Loss'
        self.gain = 'Gain'

        #column list
        self.geneCol = list(self.input.columns)[0]
        self.exonCol = list(self.input.columns)[1]
        self.normal_range_minCol = list(self.input.columns)[3]
        self.normal_range_maxCol = list(self.input.columns)[4]
        self.Sample_ratioCol = list(self.input.columns)[5]
        self.domain_col = [self.geneCol, self.exonCol, 'Domain']
        self.output_col = [self.geneCol, self.exonCol, 'CNV', 'Domain', 'Total Amplicon', 'Loss', 'Normal', 'Gain']

    def make_joinCol(self, df1, col1, col2) :
        
        df1['joinCol'] = df1[col1] + '_' + df1[col2]

    def make_cnvCol(self, df) :

        #calculate threshold
        df['Loss threshold'] = df[self.normal_range_minCol] - self.Lcut
        df['Gain threshold'] = df[self.normal_range_maxCol] + self.Gcut

        #determine CNV
        df['CNV'] = self.normal
        df.loc[df[self.Sample_ratioCol] < df['Loss threshold'], 'CNV'] = self.loss
        df.loc[df[self.Sample_ratioCol] > df['Gain threshold'], 'CNV'] = self.gain

    def determine_CNV(self, input, domain_df, output) :

        #sort by gene, exon
        exon_str = input[input[self.exonCol].str.contains("MLPA") | input[self.exonCol].str.contains("Promoter")]
        exon_str = list(exon_str['joinCol'].drop_duplicates())
        join_df = input[~input[self.exonCol].str.contains("MLPA") & ~input[self.exonCol].str.contains("Promoter")]
        join_df = join_df.astype({self.exonCol : int})
        join_df = join_df.sort_values(by = [self.geneCol, self.exonCol])

        #make gene+exon list
        join_list = list(join_df['joinCol'].drop_duplicates())
        if "BRCA1_Promoter" in exon_str:
            join_list.insert(0, "BRCA1_Promoter")
        if "BRCA2_5'MLPA" in exon_str:
            index = join_list.index('BRCA2_1')
            join_list.insert(index, "BRCA2_5'MLPA")
        if "BRCA2_3'MLPA" in exon_str:
            join_list.append("BRCA2_3'MLPA")

        for comb in join_list:

            tmp = input[input['joinCol'] == comb]
            gene = comb.split('_')[0]
            exon = comb.split('_')[1]
            total = float(tmp.shape[0])
            normal = tmp[tmp['CNV'] == self.normal].shape[0]
            loss = float(tmp[tmp['CNV'] == self.loss].shape[0])
            gain = float(tmp[tmp['CNV'] == self.gain].shape[0])

            domain = domain_df[domain_df['joinCol'] == comb]['Domain']

            if domain.shape[0] == 0 :
                print('Check your domain file again.')
                domain = '.'
            else :
                domain = domain.values[0]

            loss_fraction = loss/total
            gain_fraction = gain/total

            if loss_fraction >= self.Ecut :
                cnv = self.loss
            elif gain_fraction >= self.Ecut :
                cnv = self.gain
            else :
                cnv = self.normal
            
            output = output.append(pd.Series([gene, exon, cnv, domain, int(total), int(loss), int(normal), int(gain)], index = self.output_df.columns), ignore_index=True)
        return(output)

    def run(self) :

        #fillna & make col name of domain df
        self.domain = self.domain.fillna('.')
        self.domain.columns = self.domain_col

        self.make_joinCol(self.input, self.geneCol, self.exonCol)
        self.make_joinCol(self.domain, self.geneCol, self.exonCol)

        self.make_cnvCol(self.input)

        #make output file
        self.output_df = pd.DataFrame(columns = self.output_col)
        self.output_df = self.determine_CNV(self.input, self.domain, self.output_df)

        #save output file
        self.output_df.to_csv(self.output_name, sep='\t', header=True, index=False)

def main() :

    parser = argparse.ArgumentParser(description='Hi')

    parser.add_argument('-I', '--input', action='store', dest='input', default=False, help='CNV input file')
    parser.add_argument('-D', '--domain', action='store', dest='domain', default=False, help='domain file')
    parser.add_argument('-O', '--output', action='store', dest='output', default=False, help='CNV output file')
    parser.add_argument('-L', '--loss', action='store', dest='Lcut', default=0.2, help='Loss Cutoff (default = 0.2)')
    parser.add_argument('-G', '--gain', action='store', dest='Gcut', default=0.25, help='Gain Cutoff (default = 0.25)')
    parser.add_argument('-T', '--exon', action='store', dest='Ecut', default=0.6, help='Exon threshold (default = 0.6)')

    args = parser.parse_args()
    input = args.input
    domain = args.domain
    output = args.output
    Lcut = args.Lcut
    Gcut = args.Gcut
    Ecut = args.Ecut

    cnv_def = make_cnv_file(input, domain, output, Lcut, Gcut, Ecut)
    cnv_def.run()

if __name__ == '__main__':
    main()
