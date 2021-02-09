# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 19:16:22 2021

@author: bwjh9
"""

import os
import pandas as pd
import itertools
from matplotlib import pyplot as plt
import numpy as np

def get_genes(file):
    """ retrieve identified genes from bedtools output """
    print('Getting genes from file:', file)
    
    columns = ['anchor_chr','anchor_start','anchor_end','fdrBL','ind','sample',
               'feature_chr','source','feature_type',
               'feature_start','feature_end',
               '?','strand','??','details']
    df = pd.read_csv(file, sep = '\t',header=None)

    df.columns=columns
    df = df[df.feature_type == 'gene'] #check if feature is gene
    df[['gene','name','biotype','description','gene_id',
        'logic_name','version']] = df['details'].str.split(';',expand=True) #expand hg38 gene description
    
    df['gene'] = df['gene'].str[8:]
    df['name'] = df['name'].str[5:]
    df['biotype'] = df['biotype'].str[8:]
    df['description'] = df['description'].str[12:]
    df['description'] = df['description'].str.split('[',0)
    df.loc[:, 'description'] = df.description.map(lambda x: x[0][:-1])
    df = df.reset_index(drop=True)
    
    df = df[['anchor_chr','anchor_start','anchor_end','fdrBL','ind','sample',
             'feature_chr','feature_start','feature_end','strand',
             'gene','name','biotype','description']]
    
    #return a new file genes_mapped_xxx.csv
    df.to_csv('genes_'+file.split('.')[0]+'.csv', sep = '\t',header=True)
    

def merge_ind(df):
    """ merge indices for loops associated with a given gene into a single row """
    print('Merging index... ', end='')
    for name in set(df.name):
        df.loc[df['name'] == name, 'ind'] =  '; '.join(df.loc[df['name'] == name].ind.to_list())
        df.loc[df['name'] == name, 'fdrBL'] =  '; '.join(df.loc[df['name'] == name].fdrBL.to_list())       
        
    print('Completed')
    return df

def collate(file1, file2):
    """ retrieve information from output file 1 and 2 """
    print('Files:', file1, file2, "...", end='')
    df1 = pd.read_csv(file1,sep='\t')
    df2 = pd.read_csv(file2,sep='\t')
    df = pd.concat([df1,df2])
    df['fdrBL'] = df['fdrBL'].astype(str)
    
    retaincolumns = ['gene', 'name', 'biotype','description',
                     'feature_chr', 'feature_start', 'feature_end', 'strand', 'ind', 'fdrBL']
    
    
    df=df[df.columns[1:]].reset_index(drop=True)
    df = (merge_ind(df)[retaincolumns]).drop_duplicates(ignore_index=True)
    print('Completed')
    return df
    
def cancerdata(df, SANGERfile):
    """ append cancer-linked gene information. FOS and JUN data can be ignored at present """
    print('Appending cancer-linked information... ', end='')
    onco_columns = ['gene', 'name', 'biotype', 'description', 'feature_chr',
                'feature_start', 'feature_end', 'strand',
                'Name', 'Role in Cancer','ind','fdrBL']
    
    cosmicdb = pd.read_csv(SANGERfile,sep='\t')
    cosmicdb.columns = ['name']+cosmicdb.columns[1:].to_list()
    
    df=df.merge(cosmicdb, on='name',how='outer').dropna(subset=['gene'])[onco_columns].sort_values(by=['Role in Cancer','feature_chr','feature_start']).reset_index(drop=True)
    df=df.drop_duplicates(ignore_index=True)
    print('Completed')
    return df

def exportdata(df,output_filename):
    print('Exporting...', end='')
    df.to_csv("query_" + output_filename + ".txt",sep='\t',index=False)
    print('Completed')
    
def pipeline(file1, file2, output_filename, mode='add_data', SANGERfile=''):
    df = collate(file1, file2)
    if mode == 'add_data':
        df = cancerdata(df, SANGERfile)
    exportdata(df, output_filename)

#%% Full pipeline from bedtools output

for file in os.listdir():
    if file.split('.')[-1] == "csv" and file.split('_')[0] == "mapped":
        get_genes(file)

files = [file for file in os.listdir() if file.split('.')[-1] == "csv" and file.split('_')[0] == "genes"]
filepairs = [(i[0], i[1]) for i in itertools.combinations(files,2) if (i[0].split('.')[0].split('_')[:-1] == i[1].split('.')[0].split('_')[:-1])]

for file1, file2 in filepairs:
    pipeline(file1=file1, file2=file2,
             output_filename = file1.split('.')[0][:-2],
             mode = 'add_data',
             SANGERfile = 'COSMICDB.tsv')