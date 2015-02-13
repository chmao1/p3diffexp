#!/usr/bin/python

import argparse
import pandas as pd
import json
import sys
import numpy as np
import requests
import os
import uuid
from scipy import stats


#Input
#1. metadata in json with the following:

"""
{xfile:"comparisons file",
xformat:"csv || tsv || xls ||  xlsx",
xsetup:"gene_matrix || gene_list",
source_id_type:"refseq_locus_tag || alt_locus_tag || feature_id", 
data_type: "Transcriptomics || Proteomics || Phenomics", 
title: "User input", 
description: "User input",
organism: "user input", 
pmid: "user_input",
output_path: "path",
"metadata_template":"file",
"metadata_format":"csv || tsv || xls ||  xlsx"}
"""
#2. server info for the data api
"""
{"data_api":"url"}
"""

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

#{source_id_type:"refseq_locus_tag || alt_locus_tag || feature_id", 
#data_type: "Transcriptomics || Proteomics || Phenomics", 
#experiment_title: "User input", experiment_description: "User input",
#organism name: "user input", pubmed_id: "user_input"}

#Sample Output
#experiment.json
#{"origFileName":"filename","geneMapped":4886,"samples":8,"geneTotal":4985,"cdate":"2013-01-28 13:40:47","desc":"user input","organism":"some org","owner":"user name","title":"user input","pmid":"user input","expid":"whatever","collectionType":"ExpressionExperiment","genesMissed":99,"mdate":"2013-01-28 13:40:47"}
def create_experiment_file(output_path, mapping_dict, sample_dict, expression_dict, form_data, experiment_id):
    experiment_dict={"geneMapped":mapping_dict["mapping"]["mapped_ids"],"samples":len(sample_dict['sample']),"geneTotal":mapping_dict["mapping"]["mapped_ids"]+mapping_dict["mapping"]["unmapped_ids"],"desc":(form_data.get('description','')),"organism":form_data.get('organism',''),"title":form_data.get("title",""),"pmid":form_data.get("pmid",""),"expid":experiment_id,"collectionType":"ExpressionExperiment","genesMissed":mapping_dict["mapping"]["unmapped_ids"]}
    output_file=os.path.join(output_path, 'experiment.json')
    out_handle=open(output_file, 'w')
    json.dump(experiment_dict, out_handle)
    out_handle.close()
    return experiment_dict
    

#expression.json
#{"expression":[{"log_ratio":"0.912","na_feature_id":"36731006","exp_locus_tag":"VBISalEnt101322_0001","pid":"8f2e7338-9f04-4ba5-9fe2-5365c857d57fS0","z_score":"-0.23331085637221843"}]
 
#sample.json
#{"sample":[{"sig_log_ratio":2675,"expmean":"1.258","sampleUserGivenId":"LB_stat_AerobicM9_stat_aerobic","expname":"LB_stat_AerobicM9_stat_aerobic","pid":"8f2e7338-9f04-4ba5-9fe2-5365c857d57fS0","genes":4429,"sig_z_score":139,"expstddev":"1.483"}]}

def create_comparison_files(output_path, comparisons_table, form_data, experiment_id, sig_z, sig_log):
    #create dicts for json
    sample_dict={'sample':[]}
    expression_dict={'expression':[]}
    #create stats table for sample.json
    grouped=comparisons_table.groupby(["Comparison ID"])
    sample_stats=grouped.agg([np.mean, np.std])['Log Ratio']
    sample_stats=sample_stats.rename(columns={'mean':'expmean','std':'expstddev'})
    sample_stats["genes"]=grouped.count()["Gene ID"]
    sample_stats["pid"]=[str(experiment_id)+"S"+str(i) for i in range(0,len(sample_stats))]
    sample_stats["sampleUserGivenId"]=sample_stats.index
    sample_stats["expname"]=sample_stats.index
    #get zscore and significance columns
    comparisons_table["z_score"]=grouped.transform(stats.zscore)
    comparisons_table["sig_z"]=comparisons_table["z_score"] >= sig_z
    comparisons_table["sig_log"]=comparisons_table["Log Ratio"] >= sig_log
    #store counts in stats
    sample_stats["sig_z_score"]=comparisons_table.groupby(["Comparison ID","sig_z"]).count()['sig_z'].unstack()[True]
    sample_stats["sig_log_ratio"]=comparisons_table.groupby(["Comparison ID","sig_log"]).count()['sig_log'].unstack()[True]
    sample_stats["sig_log_ratio"]=sample_stats["sig_log_ratio"].fillna(0).astype('int64')
    sample_stats["sig_z_score"]=sample_stats["sig_z_score"].fillna(0).astype('int64')
    #set pid's for expression.json
    comparisons_table=comparisons_table.rename(columns={'Comparison ID':'sampleUserGivenId','Gene ID': 'exp_locus_tag', 'Map ID': "feature_id","Log Ratio":'log_ratio'})
    comparisons_table=comparisons_table.merge(sample_stats[["pid","sampleUserGivenId"]], how="left", on="sampleUserGivenId")
    #pull in metadata spreadsheet if provided
    if 'metadata_template' in form_data and 'metadata_format' in form_data and form_data['metadata_template']:
        meta_table=None
        meta_cols=["Comparison ID","Title","PubMed","Accession","Organism","Strain","Gene Modification","Experiment Condition","Time Point"]
        try:
            if form_data['metadata_format'] == 'xls' or form_data['metadata_format'] == 'xlsx':
                meta_table=pd.io.excel.read_excel(form_data['metadata_template'], 0, index_col=None)
            elif form_data['metadata_format'] == 'csv':        
                meta_table=pd.read_csv(form_data['metadata_template'], header=0)
            else:
                meta_table=pd.read_table(form_data['metadata_template'], header=0)
            meta_key="sampleUserGivenId"
            meta_table=meta_table[meta_cols].rename(columns={'Comparison ID':'sampleUserGivenId', 'Title':'expname', 'PubMed':'pubmed', 'Accession':'accession', 'Organism':'organism', 'Strain':'strain', 'Gene Modification':'mutant', 'Experiment Condition':'condition', 'Time Point':'timepoint'})
            to_add=meta_table.columns-sample_stats.columns
            meta_table=meta_table.set_index('sampleUserGivenId')
            sample_stats.update(meta_table)
            sample_stats=sample_stats.merge(meta_table[to_add], left_index=True, right_index=True)
        except:
            sys.stderr.write("failed to use user provide metadata template\n")
            pass
    #populate json dicts
    sample_stats=sample_stats.fillna("")
    sample_dict['sample']=sample_stats.to_dict(outtype='records')
    cols = [col for col in comparisons_table.columns if col not in ['sig_z', 'sig_log']]
    expression_dict['expression']=comparisons_table[cols].to_dict(outtype='records')
    output_file=os.path.join(output_path, 'sample.json')
    out_handle=open(output_file, 'w')
    json.dump(sample_dict, out_handle)
    out_handle.close()
    output_file=os.path.join(output_path, 'expression.json')
    out_handle=open(output_file, 'w')
    json.dump(expression_dict, out_handle)
    out_handle.close()
    return (sample_dict, expression_dict)
    
    

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
    return mapping_dict

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
    current_query={'q':form_data["source_id_type"]+":("+" OR ".join(id_list)+") AND annotation:PATRIC"}
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
    sig_z=2
    sig_log=1
    valid_formats=set(['csv', 'tsv', 'xls', 'xlsx'])
    valid_setups=set(['gene_matrix','gene_list'])
    
    req_info=['output_path', 'xfile','xformat','xsetup','source_id_type','data_type','experiment_title','experiment_description','organism'] 

    parser = argparse.ArgumentParser()
    userinfo = parser.add_mutually_exclusive_group(required=True)
    userinfo.add_argument('--ufile', help='json file from user input')
    userinfo.add_argument('--ustring', help='json string from user input')
    serverinfo = parser.add_mutually_exclusive_group(required=True)
    serverinfo.add_argument('--sfile', help='server setup JSON file')
    serverinfo.add_argument('--sstring', help='server setup JSON string')
    args = parser.parse_args()
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)

    #parse user form data
    form_data=None
    user_parse=None
    server_parse=None
    parse_server = json.loads if 'sstring' in args else json.load
        
    try:
        form_data = json.loads(args.ustring) if 'ustring' in args else json.load(args.ufile)
    except:
        sys.stderr.write("Failed to parse user provided form data \n")
        raise
    #parse setup data
    try:
        server_setup= json.loads(args.sstring) if 'sstring' in args else json.load(args.sfile)
    except:
        sys.stderr.write("Failed to parse server data\n")
        raise

    #make sure all required info present
    missing=[x not in form_data for x in req_info]
    if (any(missing)):
        sys.stderr.write("Missing required user input data: "+" ".join([req_info[i] for i in range(len(missing)) if missing[i]])+"\n")
        sys.exit(2)
    if ('metadata_template' in form_data or 'metadata_format' in form_data) and ('metadata_format' not in form_data or 'metadata_template' not in form_data):
        sys.stderr.write("Expression transformation: (file,format) pair must be given for metadata template\n")
        #sys.exit(2)

    #read comparisons file
    comparisons_table=None
    if form_data['xformat'] == 'csv':
        comparisons_table=pd.read_csv(form_data['xfile'], header=0)
    if form_data['xformat'] == 'tsv':
        comparisons_table=pd.read_table(form_data['xfile'], header=0)
    if form_data['xformat'] == 'xls' or form_data['xformat'] == 'xlsx':
        comparisons_table=pd.io.excel.read_excel(form_data['xfile'], 0, index_col=None)

    check_columns=None
    if form_data['xsetup'] == 'gene_matrix':
        check_columns=matrix_columns
    else:
        check_columns=list_columns
    columns_ok = True
    for i in check_columns:
        columns_ok=columns_ok and i in comparisons_table.columns
    if not columns_ok:
            sys.stderr.write("Missing appropriate column names in "+form_data['xfile']+"\n")
            sys.exit(2)

    output_path=form_data["output_path"]

    #convert gene matrix to list
    if form_data['xsetup'] == 'gene_matrix':
        comparisons_table=gene_matrix_to_list(comparisons_table)

    #limit log ratios
    comparisons_table.ix[comparisons_table["Log Ratio"] > 1000000, 'Log Ratio']=1000000
    comparisons_table.ix[comparisons_table["Log Ratio"] < -1000000, 'Log Ratio']=-1000000
    comparisons_table=comparisons_table.dropna()

    #map gene ids
    mapping_table=list_to_mapping_table(comparisons_table)
    map_gene_ids(mapping_table, form_data, server_setup)
    comparisons_table=comparisons_table.merge(mapping_table, how='left', on="Gene ID")

    #create json files to represent experiment
    experiment_id=str(uuid.uuid1())
    mapping_dict=create_mapping_file(output_path, mapping_table, form_data)
    (sample_dict, expression_dict) = create_comparison_files(output_path, comparisons_table, form_data, experiment_id, sig_z, sig_log)
    experiment_dict=create_experiment_file(output_path, mapping_dict, sample_dict, expression_dict, form_data, experiment_id)
    sys.stdout.write(json.dumps(experiment_dict)+"\n")
    
     
if __name__ == "__main__":
    main()
