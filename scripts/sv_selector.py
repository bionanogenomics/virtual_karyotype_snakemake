import argparse
import pandas as pd
import os

def select_svs(input, output, sample_n):
    """
    Function groups SV table by Type and chromosome and selects 1 per group to be simulated
    """
    sv_table = pd.read_table(input)
    sv_table['Group'] = sv_table['Group'].astype(str)
    sample_svs = sv_table[sv_table['Group'] == sample_n].iloc[:,:-1]
    # sampled_svs = sv_table.groupby(['Type','chrom1']).sample(1)
    sample_svs.to_csv(output,index=False,sep='\t')
        
def main():
    parser = argparse.ArgumentParser(
        """Function accepts an SV table with columns: ['Type', 'chrom1', 'chrom2', 'breakpoint1', 'breakpoint2', 'SVsize', 'Orientation']
        groups by both Type and chrom1 and selects one SV per group to be simulated"""
    )
    parser.add_argument('--input', type=str, help="")
    parser.add_argument('--output', type=str, help="")
    parser.add_argument('--sample', type=str, help="")

    # parse command line arguments
    args = parser.parse_args()
    print(args)
    input = args.input
    output = args.output
    sample_n = args.sample

    select_svs(input=input, output=output, sample_n=sample_n)

if __name__ == "__main__":
    main()