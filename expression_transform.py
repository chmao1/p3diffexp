#!/usr/bin/python

import argparse
import pandas as pd
import json
import sys
import numpy as np
import requests
import os

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


def pretty_print_POST(req):
    """
    At this point it is completely built and ready
    to be fired; it is "prepared".

    However pay attention at the formatting used in
    this function because it is programmed to be pretty
    printed and may differ from the actual request.
    """
    print('{}\n{}\n{}\n\n{}'.format(
        '-----------START-----------',
        req.method + ' ' + req.url,
        '\n'.join('{}: {}'.format(k, v) for k, v in req.headers.items()),
        req.body,
    ))

#convert gene list format to gene matrix
#there is definitely a more efficient conversion than this...
def gene_list_to_matrix(cur_table):
    comparisons=set(cur_table['Comparison ID'])
    genes=set(cur_table['Gene ID'])
    result=pd.DataFrame(index=sorted(list(genes)), columns=sorted(list(comparisons)))
    result['Gene ID']=result.index
    gene_pos=cur_table.columns.get_loc('Gene ID')
    comparison_pos=cur_table.columns.get_loc('Comparison ID')
    ratio_pos=cur_table.columns.get_loc('Log Ratio')
    for row in cur_table.iterrows():
        gene_id=row[-1][gene_pos]
        comp=row[-1][comparison_pos]
        ratio=row[-1][ratio_pos]
        result[comp][gene_id]=ratio
    return result

def list_to_mapping_table(cur_table):
    genes=set(cur_table['Gene ID'])
    result=pd.DataFrame(index=sorted(list(genes)))
    result['Gene ID']=result.index
    return result


#mapping.json
#{"mapping":{"unmapped_list":[{"exp_locus_tag":"VBISalEnt101322_pg001"}],"unmapped_ids":99,"mapped_list":[{"na_feature_id":"36731006","exp_locus_tag":"VBISalEnt101322_0001"}],"mapped_ids":4886}}

#creates mapping.json for results
def create_mapping_file(output_path, mapping_table, form_data):
    mapping_dict={"mapping":{"unmapped_list":[],"unmapped_ids":0,"mapped_list":[],"mapped_ids":0}}
    mapping_dict['mapping']['unmapped_list']=mapping_table[mapping_table.isnull().any(axis=1)][['Gene ID']].rename(columns={'Gene ID': 'exp_locus_tag'}).to_dict(outtype='records')
    mapping_dict['mapping']['mapped_list']=mapping_table[mapping_table.notnull().all(axis=1)].rename(columns={'Gene ID': 'exp_locus_tag', 'Map ID': "feature_id"}).to_dict(outtype='records')
    mapping_dict['mapping']['unmapped_ids']=len(mapping_dict['mapping']['unmapped_list'])
    mapping_dict['mapping']['mapped_ids']=len(mapping_dict['mapping']['mapped_list'])
    output_file=os.path.join(output_path, 'mapping.json')
    out_handle=open(output_file, 'w')
    json.dump(mapping_dict, out_handle)
    out_handle.close()

    #mapped_list=[{form_data["source_id_type"]: i["Map ID"], "exp_locus_tag":i['Gene ID']} for i in mapping_table[mapping_table.notnull().any(axis=1)]]
    #mapped_list=[{form_data["source_id_type"]: i["Map ID"], "exp_locus_tag":i["Gene ID"]} for i in mapping_table.query('Gene ID != @np.nan')]

#convert gene matrix format to gene list
#there is definitely a more efficient conversion than this...
def gene_matrix_to_list(cur_table):
    result=pd.melt(cur_table, id_vars=['Gene ID'], var_name='Comparison ID', value_name='Log Ratio')
    return result

def place_ids(query_results,cur_table,form_data):
    for d in query_results.json()['response']['docs']:
        source_id=None
        target_id=None
        if form_data["source_id_type"] in d:
            source_id=d[form_data["source_id_type"]]
        if 'feature_id' in d:
            target_id=d['feature_id']
        if source_id and target_id:
            cur_table["Map ID"][source_id]=target_id

def make_map_query(id_list, form_data, server_setup, chunk_size):
    current_query={'q':form_data["source_id_type"]+":("+" OR ".join(id_list)+")"}
    current_query["fl"]="feature_id,"+form_data["source_id_type"]
    current_query["rows"]=str(chunk_size)
    current_query["wt"]="json"
    #headers = {'Content-Type': 'application/solrquery+x-www-form-urlencoded'}
    headers = {'Content-Type': 'application/x-www-form-urlencoded; charset=utf-8'}
    req = requests.Request('POST', server_setup["data_api"], headers=headers, data=current_query)
    prepared = req.prepare()
    #pretty_print_POST(prepared)
    s = requests.Session()
    response=s.send(prepared)
    if not response.ok:
        sys.stderr.write("Mapping API not responding. Please try again later.\n")
        #sys.exit(2)
    return response

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

def map_gene_ids(cur_table, form_data, server_setup):
    cur_table["Map ID"]=np.nan
    chunk_size=10000
    for i in chunker(cur_table['Gene ID'], chunk_size):
        mapping_results=make_map_query(i, form_data, server_setup, chunk_size)
        place_ids(mapping_results, cur_table, form_data)
        

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
    parser.add_argument('-s', required=True, help='server setup JSON file')
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

    #convert gene matrix to list
    if args.xsetup == 'gene_matrix':
        comparisons_table=gene_matrix_to_list(comparisons_table)

    #parse user provided data
    form_data=None
    try:
        form_data=json.loads(args.u)
    except:
        sys.stderr.write("Failed to parse user provided form data "+args.u+"\n")
        raise
    try:
        server_setup=json.loads(args.s)
    except:
        sys.stderr.write("Failed to parse server data "+args.s+"\n")
        raise
    mapping_table=list_to_mapping_table(comparisons_table)
    #map gene ids
    map_gene_ids(mapping_table, form_data, server_setup)
    comparisons_table=comparisons_table.merge(mapping_table, how='left', on="Gene ID")
    create_mapping_file('./', mapping_table, form_data)

    

    
     
if __name__ == "__main__":
    main()
