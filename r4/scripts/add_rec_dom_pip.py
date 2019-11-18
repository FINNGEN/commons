#!/usr/bin/env python3

import sys
import os
import glob
import datetime
import pandas as pd
import gzip
import math

def run():

    df = pd.read_csv(sys.argv[1], sep='\t').fillna('NA')
    df['pheno_variant'] = df['pheno'] + '|' + df['variant']
    df = df[['pheno_variant', 'pheno','phenoname','variant','rsid','pval','beta','sebeta','variant_category','most_severe','gene_most_severe','INFO','AF','AC_Hom','enrichment_nfsee','grch37_locus']]
    pv_dict = {pv: {} for pv in df['pheno_variant']}
    for pv in pv_dict:
        pv_dict[pv]['pval_recessive'] = 'NA'
        pv_dict[pv]['pval_dominant'] = 'NA'
        pv_dict[pv]['pip'] = 'NA'

    for i,t in enumerate(['recessive', 'dominant']):
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\tadding ' + t)
        with gzip.open(sys.argv[i+2], 'rt') as f:
            header = {h:j for j,h in enumerate(f.readline().rstrip().split('\t'))}
            for line in f:
                line = line.strip().split('\t')
                variant = line[header['variant']].replace('chr', '').replace('X', '23').replace('_', ':')
                pv = line[header['pheno']] + '|' + variant
                if pv in pv_dict and line[header['p.value']] != 'NA':
                    pv_dict[pv]['pval_' + t] = float(line[header['p.value']])

    files = glob.glob(sys.argv[4])
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\tadding pip from ' + str(len(files)) + ' files')
    for filename in files:
        pheno = os.path.basename(filename).split('.')[0]
        with open(filename) as f:
            header = {h:i for i,h in enumerate(f.readline().rstrip().split('\t'))}
            for line in f:
                line = line.strip().split('\t')
                variant = line[header['rsid']].replace('chr', '').replace('X', '23').replace('_', ':')
                pv = pheno + '|' + variant
                if pv in pv_dict and int(line[header['cs']]) > 0:
                    pv_dict[pv]['pip'] = float(line[header['prob']])

    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\tadding columns to data frame')
    df['pval_recessive'] = df.apply(lambda row: pv_dict[row['pheno_variant']]['pval_recessive'], axis=1)
    df['pval_dominant'] = df.apply(lambda row: pv_dict[row['pheno_variant']]['pval_dominant'], axis=1)
    df['rec_dom_log_ratio'] = df.apply(lambda row: math.log10(row['pval_recessive']/row['pval_dominant']) if row['pval_recessive'] != "NA" and row['pval_dominant'] != "NA" else "NA", axis=1)
    df['pip'] = df.apply(lambda row: pv_dict[row['pheno_variant']]['pip'], axis=1)
    df.drop('pheno_variant', axis=1, inplace=True)
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\twriting data frame')
    df.to_csv('coding_web.txt', sep='\t', index=False)

if __name__ == '__main__':
    run()
    
