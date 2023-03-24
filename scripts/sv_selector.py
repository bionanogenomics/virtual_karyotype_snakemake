import argparse
import pandas as pd
import os

def select_svs(input, output):
    """
    Function groups SV table by Type and chromosome and selects 1 per group to be simulated
    """
    sv_table = pd.read_table(input)
    sampled_svs = sv_table.groupby(['Type','chrom1']).sample(1)
    sampled_svs.to_csv(output,index=False,sep='\t')
        
def main():
    parser = argparse.ArgumentParser(
        """Function accepts an SV table with columns: ['Type', 'chrom1', 'chrom2', 'breakpoint1', 'breakpoint2', 'SVsize', 'Orientation']
        groups by both Type and chrom1 and selects one SV per group to be simulated"""
    )
    parser.add_argument('--input', type=str, help="")
    parser.add_argument('--output', type=str, help="")

    # parse command line arguments
    args = parser.parse_args()
    print(args)
    input = args.input
    output = args.output

    select_svs(input=input, output=output)

if __name__ == "__main__":
    main()