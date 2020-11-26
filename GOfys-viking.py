# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:58:41 2019

@author: jr1214
"""

import sys
import os
import csv
import pandas as pd 
import numpy as np

import time
start_time = time.time()


###Guide to file inputs
#GOI_frame = sys.argv[1]
#EXPRMTRX_frame = sys.argv[2]
#Tig2GO_fields = sys.argv[3]
   



#Step 1 Creates frame from GOI input file.  
GOI_frame = pd.read_csv(sys.argv[1],delimiter="\t", header=None,names = ['GO_ID','GO_DESC'],encoding = 'unicode_escape')

#Step 1 Creates frame from TMM.matrix input file.
EXPRMTRX_frame = pd.read_csv(sys.argv[2],header=None,delimiter="\t")

#Step 1 Creates frame to input desired fields from trinotate summary report. Mainly tig number and ontology
Tig2GO_fields = ['#gene_id','gene_ontology_blast','gene_ontology_pfam']

#Creates frame from only the desired columns                 
Tig2GO_frame = pd.read_csv(sys.argv[3],delimiter="\t",usecols=Tig2GO_fields)



path = 'Subclusters/'  

###only need to run once if the files have been adapted for the script
for cluster in os.listdir(path):
    src=open(path + cluster,"r") #fixed it?
    fline="tig\t"    #Prepending string
    oline=src.readlines()
    oline.insert(0,fline)
    src.close()
    src=open(path + cluster,"w")
    src.writelines(oline)
    src.close()

###creating the cluster summary frame column names.
cluster_frames = []   
for cluster in os.listdir(path):    
    cluster_name = cluster[:12]
    cluster_frames.append(cluster_name)   
Cluster_Frame = pd.DataFrame(columns=[cluster_frames])

Text_Summary = []
Blast_summary_frame = pd.DataFrame()
Pfam_summary_frame = pd.DataFrame()
cluster_lengths = len(cluster_frames)

#For each GO term of interest
for index,item in enumerate(GOI_frame["GO_ID"]):
    #create empty dataframes from each variable in GO_ID
    
    pfam_frame = pd.DataFrame()
    blast_frame = pd.DataFrame()
    blast_sums =pd.DataFrame()
    pfam_sums =pd.DataFrame()
    cluster_GOI_counter = [0] * cluster_lengths
    cluster_GOI_frame = []  

    #Populates each frame with the rows in the tig2go frame if it contains the GO of interest    
    enriched_blast = Tig2GO_frame[Tig2GO_frame['gene_ontology_blast'].str.contains(item)]  
    enriched_pfam = Tig2GO_frame[Tig2GO_frame['gene_ontology_pfam'].str.contains(item)]
    #enriched_blast.drop_duplicates(subset=["#gene_id"], keep='first', inplace = True)
    #enriched_pfam.drop_duplicates(subset=["#gene_id"], keep='first', inplace = True)                                           
   
    #Creates a text report indicating how many blast and pfam matches are found for each respective GO in the ENTRIE dataset
    GOI_Sum = ((str(item) +' term codes for ' + GOI_frame.loc[index,'GO_DESC'] +' metabolism and contains '+
           str(len(enriched_blast)) + ' blast matches and ' + str(len(enriched_pfam)) +' pfam matches'))
    Text_Summary += (GOI_Sum,)
    
    
    #iterates through the contigs that are annotated with the GO of interest
    for blast_tig in enriched_blast["#gene_id"]:
       # print (blast_tig)
       enriched_Btig = EXPRMTRX_frame[EXPRMTRX_frame[0].str.match(blast_tig, case = False, na=False)]
       
       #creates a frame for the GO matches over timepoints
       blast_frame = blast_frame.append(enriched_Btig, ignore_index=True)
       
       #calculates sum of the columns
       blast_sums = blast_frame.append(blast_frame.select_dtypes(pd.np.number).sum().rename('Total'))
       
       # For each cluster count how many tigs are in there that match the GO of interst
       for clustidx,cluster in enumerate(os.listdir(path)):  
           cluster_tig_counter = [0]
           Cluster_data = pd.read_csv(path + cluster, delimiter="\t")
           tig_counts = Cluster_data.tig.str.contains(blast_tig).sum()
           cluster_tig_counter += tig_counts
           cluster_GOI_counter[clustidx] += int(cluster_tig_counter)
    
    Cluster_Frame.loc[GOI_frame.loc[index,'GO_DESC']] = cluster_GOI_counter
    

    
    Blast_summary_frame = Blast_summary_frame.append(blast_frame.select_dtypes(pd.np.number).sum().rename(GOI_frame.loc[index,'GO_DESC']))
    print ('Blast Results for ' + item + ' ' + GOI_frame.loc[index,'GO_DESC'])
    print (blast_sums)
           
    
    for pfam_tig in enriched_pfam["#gene_id"]:
       # print (blast_tig)
       enriched_Ptig = EXPRMTRX_frame[EXPRMTRX_frame[0].str.match(pfam_tig, case = False, na=False)]
       pfam_frame = pfam_frame.append(enriched_Ptig, ignore_index=True)
       #print (enriched_tig)
       pfam_sums = pfam_frame.append(pfam_frame.select_dtypes(pd.np.number).sum().rename('Total'))
    Pfam_summary_frame = Pfam_summary_frame.append(pfam_frame.select_dtypes(pd.np.number).sum().rename(GOI_frame.loc[index,'GO_DESC']))
    print ('Pfam Results for ' + item + ' ' + GOI_frame.loc[index,'GO_DESC'])
    print (pfam_sums)
    
#drops rows in frame with no matches
#Blast_summary_frame = Blast_summary_frame.dropna(how='all')
#Pfam_summary_frame = Pfam_summary_frame.dropna(how='all')
#Cluster_Frame = Cluster_Frame[(Cluster_Frame.T != 0).any()]


#pulls the timepoints into the summary frame
Blast_summary_frame.columns = [EXPRMTRX_frame.loc[0,1:]]
Pfam_summary_frame.columns = [EXPRMTRX_frame.loc[0,1:]]

print ('Blast results summmmry')
print (Blast_summary_frame)
print ('Pfam results summmmry')
print (Pfam_summary_frame)
for i in Text_Summary:   
    
    print (i)   
    
Cluster_Frame.to_csv(sys.argv[1]+r'_Cluster_Summary.csv')
Blast_summary_frame.to_csv(sys.argv[1]+r'_Blast_Summary.csv')
Pfam_summary_frame.to_csv(sys.argv[1]+r'_Pfam_Summary.csv')

Text_Summary_File = open(sys.argv[1]+"_Text_Summary.txt","w+")
for i in Text_Summary:
   Text_Summary_File.write(i+ "\n")
Text_Summary_File.close()     
   

    
print ("My program took", time.time() - start_time, "to run")  
    

                        
                        
                        
              
