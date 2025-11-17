#! /usr/bin/env python

# import python modules
import pandas as pd
from datetime import date
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import sys
import argparse
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(description="Format aggregate results into a FASTA file with metadata headers.")
    parser.add_argument("--aggregate", required=True, help="Input TSV file with aggregate results.")
    parser.add_argument("--sample_name", required=False, help="sample name.")
    parser.add_argument("--out_fn", required=False, help="")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    df = pd.read_csv(args.aggregate, sep="\t", 
                     dtype= {'lineages': str, 
                             'abundances': str})
    
    # clean up aggregate file
    df = df.rename(columns = {'Unnamed: 0': 'file_name'})
    df['sample_name'] = df['file_name'].str.split('.').str[0]

    # reformat to long file format
    sample_name = df.sample_name[0]
    coverage = df.coverage[0]
    resid = df.resid[0]
    lineages = df.lineages[0].split(',')
    abundances = df.abundances[0].split(',')
    lineage_dict = dict(zip(lineages, abundances))

    columns = ['sample_name', 'lineage', 'abundance', 'coverage', 'resid']
    new_df = pd.DataFrame(columns=columns)

    for lineage, abundance in lineage_dict.items():
        new_row = [sample_name, lineage, abundance, coverage, resid]
        new_df.loc[len(df)] = new_row
    
    new_df.to_csv(args.out_fn, index=False)

if __name__ == "__main__":
    main()