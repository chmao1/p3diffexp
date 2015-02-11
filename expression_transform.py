#!/usr/bin/python

import argparse
import pandas as pd
import json
import sys

#Input
#1. PATRIC (Gene Matrix || Gene List) in csv, tsv, xls, or  xlsx formats
#2. (optional) PATRIC (Metadata template) in csv, tsv, xls, or xlsx formats
#3. trasformation metadata in json string with the following:

#{source_id_type:"refseq_locus_tag || alt_locus_tag || feature_id", 
#data_type: "Transcriptomics || Proteomics || Phenomics", 
#experiment_title: "User input", experiment_description: "User input",
#organism name: "user input", pubmed_id: "user_input"}

#Sample Output
#experiment.json
#{"origFileName":"filename","geneMapped":4886,"samples":8,"geneTotal":4985,"cdate":"2013-01-28 13:40:47","desc":"user input","organism":"some org","owner":"user name","title":"user input","pmid":"user input","expid":"whatever","collectionType":"ExpressionExperiment","genesMissed":99,"mdate":"2013-01-28 13:40:47"}

#expression.json
#{"expression":[{"log_ratio":"0.912","na_feature_id":"36731006","exp_locus_tag":"VBISalEnt101322_0001","pid":"8f2e7338-9f04-4ba5-9fe2-5365c857d57fS0","z_score":"-0.23331085637221843"}]

#mapping.json
#{"mapping":{"unmapped_list":[{"exp_locus_tag":"VBISalEnt101322_pg001"}],"unmapped_ids":99,"mapped_list":[{"na_feature_id":"36731006","exp_locus_tag":"VBISalEnt101322_0001"}],"mapped_ids":4886}}

#sample.json
#{"sample":[{"sig_log_ratio":2675,"expmean":"1.258","sampleUserGivenId":"LB_stat_AerobicM9_stat_aerobic","expname":"LB_stat_AerobicM9_stat_aerobic","pid":"8f2e7338-9f04-4ba5-9fe2-5365c857d57fS0","genes":4429,"sig_z_score":139,"expstddev":"1.483"}]}

def gene_list_to_matrix(cur_table):
    comparisons=set(cur_table['Comparison ID'])
    genes=set(cur_table['Gene ID'])
    result=pd.DataFrame(index=sorted(list(genes)), columns=sorted(list(comparisons)))
    gene_pos=cur_table.columns.get_loc('Gene ID')
    comparison_pos=cur_table.columns.get_loc('Comparison ID')
    ratio_pos=cur_table.columns.get_loc('Log Ratio')
    for row in cur_table.iterrows():
        gene_id=row[-1][gene_pos]
        comp=row[-1][comparison_pos]
        ratio=row[-1][ratio_pos]
        result[comp][gene_id]=ratio
    return result
        

def main():
    matrix_columns=['Gene ID']
    list_columns=['Gene ID', 'Comparison ID', 'Log Ratio']
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', required=True, help='comparisons file')
    parser.add_argument('--xformat', required=True, help='format of comparisons', choices=['csv', 'tsv', 'xls', 'xlsx'])
    parser.add_argument('--xsetup', required=True, help='setup for comparisons file', choices=['gene_matrix', 'gene_list'])
    parser.add_argument('-m', help='PATRIC metadata template')
    parser.add_argument('--mformat', help='format of PATRIC metadata template', choices=['csv', 'tsv', 'xls', 'xlsx'])
    parser.add_argument('-u', required=True, help='json string from user input')
    args = parser.parse_args()
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)

    if (args.m!=None or args.mformat!=None) and (args.m==None or args.mformat==None):
        sys.stderr.write("Expression transformation: (file,format) pair must be given\n")
        sys.exit(2)

    #read comparisons file
    comparisons_table=None
    if args.xformat == 'csv':
        comparisons_table=pd.read_csv(args.x, header=0)
    if args.xformat == 'tsv':
        comparisons_table=pd.read_table(args.x, header=0)
    if args.xformat == 'xls' or args.xformat == 'xlsx':
        comparisons_table=pd.io.excel.read_excel(args.x, 0, index_col=None, na_values=['NA'])

    check_columns=None
    if args.xsetup == 'gene_matrix':
        check_columns=matrix_columns
    else:
        check_columns=list_columns
    columns_ok = True
    for i in check_columns:
        columns_ok=columns_ok and i in comparisons_table.columns
    if not columns_ok:
            sys.stderr.write("Missing appropriate column names in "+args.xsetup+"\n")
            sys.exit(2)

    #parse user provided data
    form_data=None
    try:
        form_data=json.loads(args.u)
    except:
        sys.stderr.write("Failed to parse user provided form data "+args.u+"\n")
        raise

    

    
     
if __name__ == "__main__":
    main()
