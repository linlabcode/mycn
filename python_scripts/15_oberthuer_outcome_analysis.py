#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2016 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for mycn code

#See README for additional information on downloading and installing dependencies

#Requires linlab pipeline set of utilities 
#

#Requires bamliquidator
#

#Requires 

#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/home/cl6/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
import pickle
from collections import defaultdict
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'mycn'
genome ='hg19'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/mycn_resub/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figureCodeFolder = '%sfigureCode/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory
genomeDirectory = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ATAC-Seq
atac_dataFile = '%sdata_tables/ATAC_TABLE.txt' % (projectFolder)

#ChIP-Seq
be2c_dataFile = '%sdata_tables/BE2C_TABLE.txt' % (projectFolder)
mm1s_dataFile = '%sdata_tables/MM1S_TABLE.txt' % (projectFolder)
u87_dataFile = '%sdata_tables/U87_TABLE.txt' % (projectFolder)
nb_all_chip_dataFile = '%sdata_tables/NB_ALL.txt' % (projectFolder)
p4936_young_dataFile = '%sdata_tables/P493-6_YOUNG_TABLE.txt' % (projectFolder)
sclc_dataFile = '%sdata_tables/SCLC_DATA_TABLE.txt' % (projectFolder)
shep21_dataFile = '%sdata_tables/SHEP21_TABLE.txt' % (projectFolder)
shep_on_dataFile = '%sdata_tables/SHEP_ON_TABLE.txt' % (projectFolder)

chip_data_list = [be2c_dataFile,mm1s_dataFile,nb_all_chip_dataFile,p4936_young_dataFile,sclc_dataFile,shep21_dataFile,shep_on_dataFile,u87_dataFile]


#note: all mouse analysis of THMYCN tumors are in a separate script

#CHIP-RX
shep21_chiprx_dataFile = '%sdata_tables/SHEP21_CHIPRX_TABLE.txt' % (projectFolder)

#RNA-Seq
be2c_rna_drug_dataFile = '%sdata_tables/BE2C_RNA_DRUG_TABLE.txt' % (projectFolder)
be2c_rna_twist_dataFile = '%sdata_tables/BE2C_RNA_TWIST_TABLE.txt' % (projectFolder)
shep21_rna_dataFile = '%sdata_tables/SHEP21_DOX_RNA_TABLE.txt' % (projectFolder)

rna_data_list = [be2c_rna_drug_dataFile,be2c_rna_twist_dataFile,shep21_rna_dataFile]

all_data_list = [atac_dataFile,be2c_dataFile,mm1s_dataFile,nb_all_chip_dataFile,p4936_young_dataFile,sclc_dataFile,shep21_dataFile,shep_on_dataFile,u87_dataFile,shep21_chiprx_dataFile,be2c_rna_drug_dataFile,be2c_rna_twist_dataFile,shep21_rna_dataFile]

#==========================================================================
#======================GLOBAL PATHS FOR OBERTHUER DATA=====================
#==========================================================================

ob_s1_path = '%soberthuer_outcome/oberthuer_table_s1.txt' % (projectFolder)

#these are the array annotations
array_1_path = '%soberthuer_outcome/A-MEXP-1746.adf.txt' % (projectFolder)
array_2_path = '%soberthuer_outcome/A-MEXP-1747.adf.txt' % (projectFolder)

#this is the expression data
exp_1_path = '%soberthuer_outcome/251271410_processed_filled_to_ADF.txt' % (projectFolder)
exp_2_path = '%soberthuer_outcome/251496110_processed_filled_to_ADF.txt' % (projectFolder)


nb_mycn_conserved_path = '%senhancerPromoter/NB_MYCN_CONSERVED/NB_MYCN_CONSERVED_TOP_5000_ORDERED.txt' % (projectFolder)


#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for MYCN project')

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I. FORMATTING OBERTHUER DATA====================')
    print('#======================================================================')
    print('\n\n')

    #want to turn it into a table with the following columns
    #common gene name and then expression across patients using the array barcode identifier
    #e.g. 251271410154

    #first need to fix the table s1 if it's been butchered by excel
    fix_table_s1(ob_s1_path)

    #now we need to make a gene dict that goes from array id to gene name
    #make probe to gene dict
    probe_gene_dict = make_probe_to_gene_dict(annotFile,array_1_path,array_2_path)

    #now we need to make the exp table
    exp_table_path = '%soberthuer_outcome/oberthuer_expression_formatted.txt' % (projectFolder)
    patient_table_path = '%soberthuer_outcome/oberthuer_patient_formatted.txt' % (projectFolder)
    exp_table_path,patient_table_path = make_exp_table(probe_gene_dict,exp_1_path,exp_2_path,ob_s1_path,exp_table_path,patient_table_path)



    print('\n\n')
    print('#======================================================================')
    print('#=================II. MAKING ENHANCER/PROMOTER SIGS====================')
    print('#======================================================================')
    print('\n\n')



#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~FIXING TABLE S1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def fix_table_s1(ob_s1_path):
    
    '''
    fixes formatting of table s1
    '''

    s1 = open(ob_s1_path,'r')

    lines = s1.readlines()
    if len(lines) == 1:
        lines = lines[0].split('\r')
        s1_table = [line.rstrip().split('\t') for line in lines]
        utils.unParseTable(s1_table,ob_s1_path,'\t')

    return ob_s1_path


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~PROBE TO GENE DICT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def make_probe_to_gene_dict(annotFile,array_1_path,array_2_path):
    
    '''
    keyed by probe ID w/ gene as value
    '''
    #see if it already exists
    pickle_path = '%soberthuer_outcome/probe_dict.pkl' % (projectFolder)
    if utils.checkOutput(pickle_path,0,0):
        print('loading previously made probe dict at %s' % (pickle_path))
        probe_gene_dict =pickle.load(open( pickle_path, "rb" ))
        return probe_gene_dict
    
    #we want to intersect refseq common names w/ the array
    startDict = utils.makeStartDict(annotFile)

    ref_name_list = utils.uniquify([startDict[refID]['name'] for refID in startDict.keys()])
    probe_gene_dict = {}
    
    array_1 = utils.parseTable(array_1_path,'\t')
    array_2 = utils.parseTable(array_2_path,'\t')
    ticker = 0
    for line in array_1 + array_2:
        if len(line) < 5:
            continue
        ticker +=1
        probe_id = line[4]
        name = line[-1]
        # print(probe_id)
        # print(name)
        # if ticker== 10:
        #     sys.exit()
        # print(line)

        if ref_name_list.count(name) > 0:
            probe_gene_dict[probe_id] = name

    pickle.dump(probe_gene_dict,open(pickle_path,'wb'))
    return probe_gene_dict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAKE EXP TABLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
def make_exp_table(probe_gene_dict,exp_1_path,exp_2_path,ob_s1_path,exp_table_path,patient_table_path):
    
    '''
    making gene expression table
    first use an intermediary dictionary
    exp_dict[gene][patient] = list <- in case multiple exp values present at probe or patient level
    '''

    #check if this has already been done
    if utils.checkOutput(exp_table_path,0,0) and utils.checkOutput(patient_table_path,0,0):
        print('loading premade expression table at %s and formatted patient table at %s' % (exp_table_path,patient_table_path))
        return exp_table_path,patient_table_path

    #first get a master list of patients
    print('formatting patient data')
    table_s1 = utils.parseTable(ob_s1_path,'\t')
    patient_list = utils.uniquify([line[0] for line in table_s1[1:]])
    #make a dict w/ relevant patient info
    patient_dict = {}
    for line in table_s1[1:]:
        patient_line =[line[3]] + line[8:12]
        if line[12] == 'Amp':
            patient_line.append(1)
        else:
            patient_line.append(0)
        patient_dict[line[0]] = patient_line

    gene_list = utils.uniquify([probe_gene_dict[probe_id] for probe_id in probe_gene_dict.keys()])
    gene_list.sort()
    exp_dict = {}
    for gene in gene_list:
        exp_dict[gene] = {}
        for patient in patient_list:
            exp_dict[gene][patient] = []

    #next load up both expression tables
    exp_1 = utils.parseTable(exp_1_path,'\t')
    exp_2 = utils.parseTable(exp_2_path,'\t')

    exp_1_header = exp_1[0]
    exp_1_cols =[i for i in range(1,len(exp_1_header)) if patient_list.count(exp_1_header[i].split('_')[1]) > 0]
    exp_1_patients =  [exp_1_header[i].split('_')[1] for i in range(1,len(exp_1_header)) if patient_list.count(exp_1_header[i].split('_')[1]) > 0]

    
    exp_2_header = exp_2[0]
    exp_2_cols =[i for i in range(1,len(exp_2_header)) if patient_list.count(exp_2_header[i].split('_')[1]) > 0]
    exp_2_patients =  [exp_2_header[i].split('_')[1] for i in range(1,len(exp_2_header)) if patient_list.count(exp_2_header[i].split('_')[1]) > 0]

    print('loading expression dataset %s' % exp_1_path)
    for line in exp_1:
        if probe_gene_dict.has_key(line[0]):
            for i in range(len(exp_1_cols)):
                patient_name = exp_1_patients[i]
                col = exp_1_cols[i]
                expression = float(line[col])
                gene_name = probe_gene_dict[line[0]]
                exp_dict[gene_name][patient_name].append(expression)


    print('loading expression dataset %s' % exp_2_path)
    for line in exp_2:
        if probe_gene_dict.has_key(line[0]):
            for i in range(len(exp_2_cols)):
                patient_name = exp_2_patients[i]
                col = exp_2_cols[i]
                expression = float(line[col])
                gene_name = probe_gene_dict[line[0]]
                exp_dict[gene_name][patient_name].append(expression)
                
    print('making gene expression table')
    exp_table =[['GENE_NAME'] + exp_1_patients + exp_2_patients]
    for gene in gene_list:
        if min([len(exp_dict[gene][patient]) for patient in exp_1_patients+exp_2_patients]) == 0:
            print(gene)
            continue
        exp_line = [gene]+[numpy.mean(exp_dict[gene][patient]) for patient in exp_1_patients + exp_2_patients]
        exp_table.append(exp_line)

    patient_table = [['PATIENT','STAGE','EFS_D','EFS_STATUS','OS_D','OS_STATUS','MYCN_STATUS']]
    for patient in exp_1_patients + exp_2_patients:
        patient_line = [patient] + patient_dict[patient]
        patient_table.append(patient_line)

    utils.unParseTable(exp_table,exp_table_path,'\t')
    utils.unParseTable(patient_table,patient_table_path,'\t')

    return exp_table_path,patient_table_path

        

    


#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
