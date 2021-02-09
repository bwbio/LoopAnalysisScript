# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 23:05:47 2021

@author: bwjh9
"""
import os
import pandas as pd
import itertools
from tqdm import tqdm
import pickle
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt

def loops_add_index(): 
    """ add index to loop.bedpe files in same directory """
    statdf = pd.DataFrame(columns = ['sample','loops','sharedloops'],index=range(6))
    i=0
    for file in os.listdir():
        if file.split('.')[-1] == 'bedpe':
            df = pd.read_csv(file, sep = '\t')
            if ('ind' not in df.columns) or ('sample' not in df.columns):
                print("Adding loop index to %s." %file, len(df), "loops found")
                df['ind'] = [file.split('_')[0]+"_"+str(i) for i in range(len(df))]
                df['sample'] = file.split('_')[0]
                df.to_csv(file,sep = '\t',index=False)
            else:
                print(file, "already has loop index.", len(df), "loops found")
            
            statdf.loc[i]['sample'] = file.split('_')[0]
            statdf.loc[i]['loops'] = len(df)
            i+=1
            
    print('Loop indexing complete!')
    print()

def loop_blacklist(file1, file2, sharedLoops, blacklistLoops):
    """ blacklist shared loops """
    df1 = pd.read_csv(file1, sep = '\t')[['chr1','x1','x2','chr2','y1','y2','fdrBL','ind']]
    df2 = pd.read_csv(file2, sep = '\t')[['chr1','x1','x2','chr2','y1','y2','fdrBL','ind']]
    print("Comparing",file1.split('.')[0],'and',file2.split('.')[0])

    chromosomes = list(set(df1.chr1.tolist()))
    chromosomes.sort()
    
    sharedLoopsSample = 0 #counter for loops shared between file1, file2
    
    for chromosome in chromosomes:

        sharedLoopsChr = [] #container for shared loop within chromosome
        
        df1chr = df1.loc[df1['chr1'] == chromosome]
        df2chr = df2.loc[df2['chr1'] == chromosome]
        print('Chromosome %s...' % chromosome, end="")
        
        chromosomeLength = max(df1chr.y2.values.tolist()+df2chr.y2.values.tolist())
        k=150 #default number of partitions to split chromosome into (larger = faster runtime)
        partitions = [(int(i/k*chromosomeLength),int((i+1)/k*chromosomeLength)) for i in range(k)]
        
        for partition in partitions:
            
            L1, L2 = partition
            df1part = df1chr.loc[(df1chr['y2']<=L2+1000) & (df1chr['y2']>=L1+1000)] #internal check to prevent missed loops
            df2part = df2chr.loc[(df2chr['y2']<=(L2+100000)) & (df2chr['y1']>=(L1-100000))] #external check to prevent missed comparisons
        
            for i1, loop1 in df1part.iterrows():
                for i2, loop2 in df2part.iterrows():
                    if min(loop1.x2,loop2.x2)-max(loop1.x1,loop2.x1) > 0:
                        if min(loop1.y2,loop2.y2)-max(loop1.y1,loop2.y1) > 0:
                            sharedLoopsChr.append((loop1.ind,loop2.ind))
                            blacklistLoops.append(loop1.ind)
                            blacklistLoops.append(loop2.ind)
                            sharedLoopsSample +=1

        sharedLoopsChr = list(set(sharedLoopsChr))
        print('',len(sharedLoopsChr), 'shared loops')
        sharedLoops.extend(sharedLoopsChr)
    print(sharedLoopsSample,'shared loops found between cross-class')
    print()

def loop_whitelist(file1, file2, file3, sameclassLoops, whitelistLoops, blacklistLoops):
    """ whitelist common loops between 3 sampless """
    df1 = pd.read_csv(file1, sep = '\t')[['chr1','x1','x2','chr2','y1','y2','fdrBL','ind']]
    df2 = pd.read_csv(file2, sep = '\t')[['chr1','x1','x2','chr2','y1','y2','fdrBL','ind']]
    df3 = pd.read_csv(file3, sep = '\t')[['chr1','x1','x2','chr2','y1','y2','fdrBL','ind']]
    df1 = df1[~df1['ind'].isin(blacklistLoops)]
    df2 = df2[~df2['ind'].isin(blacklistLoops)]
    df3 = df3[~df3['ind'].isin(blacklistLoops)]
    print("Comparing",file1.split('.')[0],'and',file2.split('.')[0],'and',file3.split('.')[0])

    chromosomes = list(set(df1.chr1.tolist()))
    chromosomes.sort()
    
    sharedLoopsSample = 0 #counter for loops shared between file1, file2, file3
    
    for chromosome in chromosomes: #container for shared loop within chromosome

        sharedLoopsChr = [] #container for shared loop within chromosome
        
        df1chr = df1.loc[df1['chr1'] == chromosome]
        df2chr = df2.loc[df2['chr1'] == chromosome]
        df3chr = df3.loc[df3['chr1'] == chromosome]
        print('Chromosome %s...' % chromosome, end="")
        
        chromosomeLength = max(df1chr.y2.values.tolist()+df2chr.y2.values.tolist()+df3chr.y2.values.tolist())
        k=150 #default number of partitions to split chromosome into (larger = faster runtime)
        partitions = [(int(i/k*chromosomeLength),int((i+1)/k*chromosomeLength)) for i in range(k)]
        
        for partition in partitions:
            
            L1, L2 = partition
            df1part = df1chr.loc[(df1chr['y2']<=L2+1000) & (df1chr['y2']>=L1+1000)]  #internal check to prevent missed loops
            df2part = df2chr.loc[(df2chr['y2']<=(L2+100000)) & (df2chr['y1']>=(L1-100000))] #external check to prevent missed comparisons
            df3part = df3chr.loc[(df3chr['y2']<=(L2+100000)) & (df3chr['y1']>=(L1-100000))] #external check to prevent missed comparisons
            
            for i1, loop1 in df1part.iterrows():
                for i2, loop2 in df2part.iterrows():
                    if min(loop1.x2,loop2.x2)-max(loop1.x1,loop2.x1) > 0: #anchor 1 overlap
                        if min(loop1.y2,loop2.y2)-max(loop1.y1,loop2.y1) > 0: #anchor 2 overlap
                            for i3, loop3 in df3part.iterrows():
                                if min(loop1.x2,loop3.x2)-max(loop1.x1,loop3.x1) > 0: #anchor 1 overlap
                                    if min(loop1.y2,loop3.y2)-max(loop1.y1,loop3.y1) > 0: #anchor 2 overlap
                                        sharedLoopsChr.append((loop1.ind,loop2.ind,loop3.ind))
                                        whitelistLoops.append(loop1.ind)
                                        whitelistLoops.append(loop2.ind)
                                        whitelistLoops.append(loop3.ind)
                                        sharedLoopsSample +=1

        sharedLoopsChr = list(set(sharedLoopsChr))
        print('',len(sharedLoopsChr), 'shared loops')
        sameclassLoops.extend(sharedLoopsChr)
    print(sharedLoopsSample,'shared loops found between cross-class')
    print()

#%% Loop indexing 

loops_add_index() # add index to .bedpe files in same directory

#%% Loop blacklisting

sharedLoops = [] #container for pairs of shared loops
blacklistLoops = [] #container for shared loop indices

print('Blacklisting cross-class shared loops')
for file1,file2 in itertools.combinations([i for i in os.listdir() if i.split('.')[-1] == "bedpe"],2):
    if file1[2] != file2[2]:
        loop_blacklist(file1,file2,sharedLoops,blacklistLoops)
        
sharedLoops = sorted(sharedLoops, key=lambda x: (x[0], x[1]))
print(len(sharedLoops),' total pairs of shared blacklisted loops found')
print(len(blacklistLoops),' total shared blacklisted loops found')
print("Shared loop blacklisting complete")
print()

#%% Loop whitelisting

sameclassLoops = [] #container for 3-tuple of common loops
whitelistLoops = [] #container for common loop indices
print('Whitelisting same-class shared loops')
for file1,file2,file3 in itertools.permutations([i for i in os.listdir() if i.split('.')[-1] == "bedpe"],3):       
    if file1[2] == file2[2] == file3[2]: 
        loop_whitelist(file1, file2, file3, sameclassLoops, whitelistLoops, blacklistLoops)

sameclassLoops = sorted(sameclassLoops, key=lambda x: (x[0], x[1]))
print(len(sameclassLoops),' total shared whitelisted loops found')
print("Shared loop whitelisting complete")
print()
print()

#save blacklists and whitelists
with open('blacklistLoopsPart.data', 'wb') as filehandle:
    pickle.dump(blacklistLoops, filehandle)
with open('whitelistLoopsPart.data', 'wb') as filehandle:
    pickle.dump(whitelistLoops, filehandle)    
with open('differentclassLoopsPart.data', 'wb') as filehandle:
    pickle.dump(sharedLoops, filehandle)  
with open('sameclassLoopsPart.data', 'wb') as filehandle:
    pickle.dump(sameclassLoops, filehandle)
   
#%% Pre-processing 1

AML_multiSharedLoops = [] #container for graph components
Knee_multiSharedLoops = [] #container for graph components

whitelistLoops = list(set(whitelistLoops))
for loop in whitelistLoops:
    if loop[:3] == 'AML':
        tempList = []
        for sharedLoop in sameclassLoops:
            if loop in sharedLoop:
                tempList = tempList + list(sharedLoop)  #get vertices of all shared loops
        if tempList != []:
            AML_multiSharedLoops.append(set(tempList))  #get graph containing all vertices (shared loops)
    elif loop[:3] == 'Kne':
        tempList = []
        for sharedLoop in sameclassLoops:
            if loop in sharedLoop:
                tempList = tempList + list(sharedLoop)  #get vertices of all shared loops
        if tempList != []:
            Knee_multiSharedLoops.append(set(tempList)) #get graph containing all vertices (shared loops)


#%% Pre-processing 2

AML_allSampleSharedLoops = [] #container for 3-sample (COMMON) loop graphs
Knee_allSampleSharedLoops = [] #container for 3-sample (COMMON) loop graphs

for i in AML_multiSharedLoops:
    loops = set([loop[:5] for loop in i]) #check how many samples exist in loop graph
    if len(loops) >= 3: #ensure all 3 AML samples are present
        AML_allSampleSharedLoops.append(i)       
for i in Knee_multiSharedLoops:
    loops = set([loop[:6] for loop in i]) #check how many samples exist in loop graph
    if len(loops) >= 3: #ensure all 3 Knee samples are present
        Knee_allSampleSharedLoops.append(i)      

AML_commonLoops = [] #container for individual loop indices
Knee_commonLoops = [] #container for individual loop indices

for i in AML_allSampleSharedLoops:
    AML_commonLoops = AML_commonLoops + list(i)     
for i in Knee_allSampleSharedLoops:
    Knee_commonLoops = Knee_commonLoops + list(i)   
    
#remove duplicate loop indices
AML_commonLoops = set(AML_commonLoops) 
Knee_commonLoops = set(Knee_commonLoops)

#%% Save individual anchor-bed files by stratification


columns = ['chr1','x1','x2','chr2','y1','y2','fdrBL','ind','sample']
AML_masterdf = pd.DataFrame(columns=columns) #Common AML loops
Knee_masterdf = pd.DataFrame(columns=columns) #Common Knee loops
BL_masterdf = pd.DataFrame(columns=columns) #AML-Knee shared loops
Any_AML_masterdf = pd.DataFrame(columns=columns) #Any AML loops
Any_Knee_masterdf = pd.DataFrame(columns=columns) #Any Knee loops
Middle_AML_masterdf = pd.DataFrame(columns=columns) #Any AML loops - Common AML loops
Middle_Knee_masterdf = pd.DataFrame(columns=columns)  #Any Knee loops - Common Knee loops
AML_samples = ('AML28','AML29','AML30')
Knee_samples = ('Knee47','Knee49','Knee50')

for file in os.listdir():
    if file.split('.')[-1] == 'bedpe':
        print("Retreiving common loops from", file)
        df = pd.read_csv(file, sep = '\t')[['chr1','x1','x2','chr2','y1','y2','fdrBL','ind','sample']]
        
        AML_df = df[df['ind'].isin(AML_commonLoops)]
        AML_masterdf = AML_masterdf.append(AML_df)
        
        Knee_df = df[df['ind'].isin(Knee_commonLoops)]
        Knee_masterdf = Knee_masterdf.append(Knee_df)
        
        BL_df = df[df['ind'].isin(blacklistLoops)]
        BL_masterdf = BL_masterdf.append(BL_df)
        
        Any_AML_df = df[~df['ind'].isin(blacklistLoops)][df['sample'].isin(AML_samples)]
        Any_AML_masterdf = Any_AML_masterdf.append(Any_AML_df)
        
        Any_Knee_df = df[~df['ind'].isin(blacklistLoops)][df['sample'].isin(Knee_samples)]
        Any_Knee_masterdf = Any_Knee_masterdf.append(Any_Knee_df)
    
        Middle_AML_df = df[~df['ind'].isin(blacklistLoops)][~df['ind'].isin(AML_commonLoops)][df['sample'].isin(AML_samples)]
        Middle_AML_masterdf = Middle_AML_masterdf.append(Middle_AML_df)
        
        Middle_Knee_df = df[~df['ind'].isin(blacklistLoops)][~df['ind'].isin(Knee_commonLoops)][df['sample'].isin(Knee_samples)]
        Middle_Knee_masterdf = Middle_Knee_masterdf.append(Middle_Knee_df)


print("Generating BED copies of AML common loops")
path = "commonLoops_AML"
anchor1= AML_masterdf[['chr1','x1','x2','fdrBL','ind','sample']]
anchor2 = AML_masterdf[['chr2','y1','y2','fdrBL','ind','sample']]
anchor1name = path+"_1.bed" 
anchor2name = path+"_2.bed"
anchor1.to_csv(anchor1name,sep = '\t',index=False,header=False)
anchor2.to_csv(anchor2name,sep = '\t',index=False,header=False)
AML_masterdf.to_csv(path+'.bedpe.new',sep = '\t',index=False,header=True)

print("Generating BED copies of Knee common loops")
path = "commonLoops_Knee"
anchor1= Knee_masterdf[['chr1','x1','x2','fdrBL','ind','sample']]
anchor2 = Knee_masterdf[['chr2','y1','y2','fdrBL','ind','sample']]
anchor1name = path+"_1.bed" 
anchor2name = path+"_2.bed"
anchor1.to_csv(anchor1name,sep = '\t',index=False,header=False)
anchor2.to_csv(anchor2name,sep = '\t',index=False,header=False)
Knee_masterdf.to_csv(path+'.bedpe.new',sep = '\t',index=False,header=True)

print("Generating BED copies of blackisted shared loops")
path = "blacklistLoops"
anchor1= BL_masterdf[['chr1','x1','x2','fdrBL','ind','sample']]
anchor2 = BL_masterdf[['chr2','y1','y2','fdrBL','ind','sample']]
anchor1name = path+"_1.bed" 
anchor2name = path+"_2.bed"
anchor1.to_csv(anchor1name,sep = '\t',index=False,header=False)
anchor2.to_csv(anchor2name,sep = '\t',index=False,header=False)
BL_masterdf.to_csv(path+'.bedpe.new',sep = '\t',index=False,header=True)

print("Generating BED copies of Any AML common loops")
path = "anyLoops_AML"
anchor1= Any_AML_masterdf[['chr1','x1','x2','fdrBL','ind','sample']]
anchor2 = Any_AML_masterdf[['chr2','y1','y2','fdrBL','ind','sample']]
anchor1name = path+"_1.bed" 
anchor2name = path+"_2.bed"
anchor1.to_csv(anchor1name,sep = '\t',index=False,header=False)
anchor2.to_csv(anchor2name,sep = '\t',index=False,header=False)
Any_AML_masterdf.to_csv(path+'.bedpe.new',sep = '\t',index=False,header=True)


print("Generating BED copies of Any Knee shared loops")
path = "anyLoops_Knee"
anchor1= Any_Knee_masterdf[['chr1','x1','x2','fdrBL','ind','sample']]
anchor2 = Any_Knee_masterdf[['chr2','y1','y2','fdrBL','ind','sample']]
anchor1name = path+"_1.bed" 
anchor2name = path+"_2.bed"
anchor1.to_csv(anchor1name,sep = '\t',index=False,header=False)
anchor2.to_csv(anchor2name,sep = '\t',index=False,header=False)
Any_Knee_masterdf.to_csv(path+'.bedpe.new',sep = '\t',index=False,header=True)

print("Generating BED copies of Middle AML shared loops")
path = "middleLoops_AML"
anchor1= Middle_AML_masterdf[['chr1','x1','x2','fdrBL','ind','sample']]
anchor2 = Middle_AML_masterdf[['chr2','y1','y2','fdrBL','ind','sample']]
anchor1name = path+"_1.bed" 
anchor2name = path+"_2.bed"
anchor1.to_csv(anchor1name,sep = '\t',index=False,header=False)
anchor2.to_csv(anchor2name,sep = '\t',index=False,header=False)
Middle_AML_masterdf.to_csv(path+'.bedpe.new',sep = '\t',index=False,header=True)

print("Generating BED copies of Middle Knee shared loops")
path = "middleLoops_Knee"
anchor1= Middle_Knee_masterdf[['chr1','x1','x2','fdrBL','ind','sample']]
anchor2 = Middle_Knee_masterdf[['chr2','y1','y2','fdrBL','ind','sample']]
anchor1name = path+"_1.bed" 
anchor2name = path+"_2.bed"
anchor1.to_csv(anchor1name,sep = '\t',index=False,header=False)
anchor2.to_csv(anchor2name,sep = '\t',index=False,header=False)
Middle_Knee_masterdf.to_csv(path+'.bedpe.new',sep = '\t',index=False,header=True)

#%% #Generate linux command line output for bedtools

print("To identify genes associated with these loop files,")
print("please enter the following into the linux command line:")
print("When done, run the get_genes.py script")
print()
print()

for file in os.listdir():
    if file.split('.')[-1] == "bed":
        path='mapped_'+ file.split('.')[0]+'.csv'
        print('bedtools window -a', file, '-b hg38.gff -w 15000 >', path)
