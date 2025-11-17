#!/usr/bin/env python

import sys
import re
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description='Calculate alignment metrics and percent coverage.')
    parser.add_argument('--sample_name',  help='sample name')
    parser.add_argument('--project_name', help='Name of the project')
    parser.add_argument('--consensus_fasta', help='consensus FASTA file')
    parser.add_argument('--reference', help='Reference FASTA file')
    parser.add_argument('--out_fn', help = 'Out filename')
    parser.add_argument('--fastqc_clean_summary', help=' fastqc clean summary results')
    parser.add_argument('--fastqc_raw_summary', help='fastqc raw summary results')
    parser.add_argument('--samtools_coverages',  help='samtools coverage.txt file')
    return parser.parse_args()

def calculate_percent_coverage(reference, consensus, sample_name):
    ref_record = SeqIO.read(reference, "fasta")
    ref_length = len(ref_record.seq)

    seq_record = SeqIO.read(consensus, "fasta")
    seq_length = len(seq_record.seq)

    percent_coverage = round((seq_length / ref_length) * 100, 2)

    df = pd.DataFrame()
    df['sample_name'] = [str(sample_name)]
    df['percent_coverage'] = [percent_coverage]

    return df

def read_alignment_metrics(samtools_coverage, sample_name):
    d = pd.read_csv(samtools_coverage, sep="\t")
    
    df = pd.DataFrame()
    df['sample_name'] = [str(sample_name)]
    df['reads_mapped'] = d.numreads
    df['mean_depth'] = d.meandepth
    df['mean_baseq'] = d.meanbaseq
    df['mean_mapq'] = d.meanmapq

    return df

def read_fastqc_summary(fastqc_clean_summary, fastqc_raw_summary, sample_name):
    df_clean = pd.read_csv(fastqc_clean_summary, dtype = {'sample_name': object})
    if df_clean.r1_total_reads[0] == df_clean.r2_total_reads[0]:
        df_clean['paired_reads'] = df_clean.r1_total_reads
    for column in df_clean:
        if column != 'sample_name':
            df_clean = df_clean.rename(columns={column: f"clean_{column}"})
    
    df_raw = pd.read_csv(fastqc_raw_summary, dtype = {'sample_name': object})   
    if df_raw.r1_total_reads[0] == df_raw.r2_total_reads[0]:
        df_raw['paired_reads'] = df_raw.r1_total_reads
    for column in df_raw:   
        if column != 'sample_name':
            df_raw = df_raw.rename(columns={column: f"raw_{column}"})
    
    df = pd.merge(df_clean, df_raw, on='sample_name')
    
    return df


def main():
    args = get_args()
    print('Input arguments:')
    for arg in vars(args):
        print(f'{arg}: {getattr(args, arg)}')
    print('-------------------')

    percent_coverage_df = calculate_percent_coverage(args.reference, args.consensus_fasta, args.sample_name)
    alignment_metrics_df = read_alignment_metrics(args.samtools_coverages, args.sample_name)
    fastqc_summary_df = read_fastqc_summary(args.fastqc_clean_summary, args.fastqc_raw_summary, args.sample_name)

    print('percent_coverage_df:')
    print(percent_coverage_df)
    print('-------------------')

    print('alignment_metrics_df:')
    print(alignment_metrics_df)
    print('-------------------')

    print('fastqc_summary_df:')
    print(fastqc_summary_df)
    print('-------------------')

    final_df = percent_coverage_df.merge(alignment_metrics_df, on='sample_name')
    final_df = final_df.merge(fastqc_summary_df, on='sample_name')
    
    print('final_df before calculations:')
    print(final_df)
    print('-------------------')

    final_df['percent_reads_mapped'] = round((final_df['reads_mapped'] / final_df['clean_paired_reads']) * 100, 2)   
    final_df['project_name'] = args.project_name

    # order columns
    column_order = ['sample_name','project_name', 'percent_coverage', 'mean_depth', 'reads_mapped', 'percent_reads_mapped', 'mean_baseq', 'mean_mapq'] + \
                   [col for col in final_df.columns if col.startswith('clean_')] + \
                   [col for col in final_df.columns if col.startswith('raw_')]
    
    final_df = final_df[column_order]

    final_df.to_csv(args.out_fn, index=False)

if __name__ == "__main__":
    main()