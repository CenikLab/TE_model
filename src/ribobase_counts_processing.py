######author Yue Liu, email yliu5@utexas.edu 
######process data from raw reads counts talbe
######return processed counts, CPM, Quantile table
######removed non-polyA genes for paired data
######do not need to run under the ribobase environment
######run python ribobase_counts_processing.py -h for help
######ribo only:ribobase_counts_processing.py -i "ribo_input" -m "only"
######paired:ribobase_counts_processing.py -i "ribo_input" -r "rna_input" -m "paired"

from bioinfokit.analys import norm
from optparse import OptionParser
import numpy as np
import pandas as pd

####input paramater####
parser = OptionParser()
parser.add_option("-i", "--input_ribo_file", help='ribo raw file', dest='ribof')
parser.add_option("-r", "--input_rna_file", help='rna raw file', dest='rnaf')
parser.add_option("-c", "--cpm_cut_off", help='filter genes step1', dest='cpm_cutoff',default=int(1))
parser.add_option("-a", "--overall_cut_off", help='filter genes step2', dest='overall_cutoff',default=int(70))
parser.add_option("-m", "--mode", help='ribo only or paired', dest='mode')
parser.add_option('-o', '--dir', help='work directory', dest='workdir', default='.')
(options, args) = parser.parse_args()

options.cpm_cutoff = int(options.cpm_cutoff)
options.overall_cutoff = int(options.overall_cutoff)

####prepare file for ribo only data, counts, CPM and Quantile table
####normalized
def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    #df_log=np.log2(df_qn+1)
    return(df_qn)

def CPM_normalize(df):
    # now, normalize raw counts using CPM method 
    nm = norm()
    nm.cpm(df=df)
    # get CPM normalized dataframe
    cpm_df = nm.cpm_norm
    #df_log=np.log2(cpm_df+1)
    return cpm_df

def data_process(df):
    df_all = pd.read_csv(df,index_col=0)
    df_count=df_all.groupby(df_all.index).mean()
    df_count.columns=df_count.columns.str.rstrip("dedup")
    df_cpm=CPM_normalize(df_count)
    df_quantile=quantile_normalize(df_count)
    return df_count,df_cpm,df_quantile

#for propotional data (CLR, IQLR)check rare expressed genes with CPM
def dummy_gene_df(df,cpm_cutoff=1,overall_cutoff=70):
    row_cut_off = int(overall_cutoff/100*len(df.columns))
    df_dummy = df[(df<cpm_cutoff).sum(axis='columns') > row_cut_off]
    dummy_gene=df_dummy.index.to_series().rename("dummy_gene_series")
    return dummy_gene

def combine_dummy_gene(dummy_gene,df):
    non_dummy=df[~df.index.isin(dummy_gene.index)]
    dummy_df=df[df.index.isin(dummy_gene.index)]
    dummy_result=pd.DataFrame(dummy_df.sum())
    dummy_result.rename(columns={0:'dummy_gene'}, inplace=True)
    frames=[non_dummy,dummy_result.T]
    df_count_dummy=pd.concat(frames)
    return df_count_dummy

def count_ranking_pearson(df):
    #df_tmp=df[:-1]
    df_rank = df.assign(**df.iloc[:, :].rank(axis = 1, ascending = False).astype(int))
    return(df_rank)

####prepare file for paired data, counts, CPM and Quantile table
if options.mode == "only":
    print("####preprocessing raw data####")
    ribo_count_temp,ribo_CPM_temp,ribo_Q_temp = data_process(options.ribof)
    print("####dymmy gene####")
    ribo_dummy_gene=dummy_gene_df(ribo_CPM_temp,options.cpm_cutoff,options.overall_cutoff)
    print("####generate final table####")
    ribo_count_dummy=combine_dummy_gene(ribo_dummy_gene,ribo_count_temp)
    ribo_CPM=ribo_CPM_temp[~ribo_CPM_temp.index.isin(ribo_dummy_gene.index)]
    ribo_Q=ribo_Q_temp[~ribo_Q_temp.index.isin(ribo_dummy_gene.index)]
    ribo_count_dummy.to_csv(options.workdir + "/ribo_only_count_dummy_%s.csv"%(options.overall_cutoff))
    ribo_CPM.to_csv(options.workdir + "/ribo_only_cpm_dummy_%s.csv"%(options.overall_cutoff))
    ribo_Q.to_csv(options.workdir + "/ribo_only_quantile_dummy_%s.csv"%(options.overall_cutoff))
    print("####generate spearman cor table####")
    ribo_spearman_pre = ribo_count_temp[~ribo_count_temp.index.isin(ribo_dummy_gene.index)]
    ribo_spearman=count_ranking_pearson(ribo_spearman_pre)
    ribo_spearman.to_csv(options.workdir + "/ribo_only_spearman_ranking_%s.csv"%(options.overall_cutoff))
    
    
elif options.mode =="paired":
    print("####preprocessing raw ribo data####")
    ribo_count_temp,ribo_CPM_temp,ribo_Q_temp = data_process(options.ribof)
    print("####preprocessing raw rna data####")
    rna_count_temp,rna_CPM_temp,rna_Q_temp = data_process(options.rnaf)
    print("####dummy gene####")
    ribo_dummy_gene=dummy_gene_df(ribo_CPM_temp,options.cpm_cutoff,options.overall_cutoff)
    rna_dummy_gene=dummy_gene_df(rna_CPM_temp,options.cpm_cutoff,options.overall_cutoff)
    all_dummy=pd.merge(ribo_dummy_gene,rna_dummy_gene,how="outer",left_index=True, right_index=True)
    print("####generate pre-final table####")
    ribo_count_dummy_wployA=combine_dummy_gene(all_dummy,ribo_count_temp)
    ribo_CPM_wployA=ribo_CPM_temp[~ribo_CPM_temp.index.isin(all_dummy.index)]
    ribo_Q_wployA=ribo_Q_temp[~ribo_Q_temp.index.isin(all_dummy.index)]
    rna_count_dummy_wployA=combine_dummy_gene(all_dummy,rna_count_temp)
    rna_CPM_wployA=rna_CPM_temp[~rna_CPM_temp.index.isin(all_dummy.index)]
    rna_Q_wployA=rna_Q_temp[~rna_Q_temp.index.isin(all_dummy.index)]
    print("####remove polyA genes####")
    polyA=pd.read_csv("./data/nonpolyA_gene.csv",index_col=0)
    ribo_count_dummy=ribo_count_dummy_wployA[~ribo_count_dummy_wployA.index.isin(polyA.index)]
    ribo_CPM=ribo_CPM_wployA[~ribo_CPM_wployA.index.isin(polyA.index)]
    ribo_Q=ribo_Q_wployA[~ribo_Q_wployA.index.isin(polyA.index)]
    rna_count_dummy=rna_count_dummy_wployA[~rna_count_dummy_wployA.index.isin(polyA.index)]
    rna_CPM=rna_CPM_wployA[~rna_CPM_wployA.index.isin(polyA.index)]
    rna_Q=rna_Q_wployA[~rna_Q_wployA.index.isin(polyA.index)]
    print("####generate final table####")
    ribo_count_dummy.to_csv(options.workdir + "/ribo_paired_count_dummy.csv")
    ribo_CPM.to_csv(options.workdir + "/ribo_paired_cpm_dummy_%s.csv"%(options.overall_cutoff))
    ribo_Q.to_csv(options.workdir + "/ribo_paired_quantile_dummy_%s.csv"%(options.overall_cutoff))
    rna_count_dummy.to_csv(options.workdir + "/rna_paired_count_dummy.csv")
    rna_CPM.to_csv(options.workdir + "/rna_paired_cpm_dummy_%s.csv"%(options.overall_cutoff))
    rna_Q.to_csv(options.workdir + "/rna_paired_quantile_dummy_%s.csv"%(options.overall_cutoff))

