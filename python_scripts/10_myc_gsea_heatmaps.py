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


#Main method run script for generating heatmaps of enriched pathways from gsea analysis

#See README for additional information on downloading and installing dependencies

#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/home/chazlin/src/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'mycn'
genome ='hg19'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/grail/projects/mycn_resub/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

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
maskFolder = '%smasks/' % (projectFolder)
#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory
genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder,maskFolder]

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
nb_all_chip_dataFile = '%sdata_tables/NB_ALL.txt' % (projectFolder)
p4936_young_dataFile = '%sdata_tables/P493-6_YOUNG_TABLE.txt' % (projectFolder)
sclc_dataFile = '%sdata_tables/SCLC_DATA_TABLE.txt' % (projectFolder)
shep21_dataFile = '%sdata_tables/SHEP21_TABLE.txt' % (projectFolder)
shep_on_dataFile = '%sdata_tables/SHEP_ON_TABLE.txt' % (projectFolder)

chip_data_list = [be2c_dataFile,mm1s_dataFile,nb_all_chip_dataFile,p4936_young_dataFile,sclc_dataFile,shep21_dataFile,shep_on_dataFile]
#note: all mouse analysis of THMYCN tumors are in a separate script

#CHIP-RX
shep21_chiprx_dataFile = '%sdata_tables/SHEP21_CHIPRX_TABLE.txt' % (projectFolder)

#RNA-Seq
be2c_rna_drug_dataFile = '%sdata_tables/BE2C_RNA_DRUG_TABLE.txt' % (projectFolder)
be2c_rna_twist_dataFile = '%sdata_tables/BE2C_RNA_TWIST_TABLE.txt' % (projectFolder)
shep21_rna_dataFile = '%sdata_tables/SHEP21_DOX_RNA_TABLE.txt' % (projectFolder)


#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for MYCN project')

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I, LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for ChIP-Seq
    pipeline_dfci.summary(nb_all_chip_dataFile)



    print('\n\n')
    print('#======================================================================')
    print('#========================II. MAKING NES TABLE==========================')
    print('#======================================================================')
    print('\n\n')

    #at a given fdr cutoff, grab the NES pathways  
    
    nes_folder = utils.formatFolder('%snes_tables/' % (projectFolder),True)

    #for top 2k
    nes_path_list = ['%senhancerPromoter/NB_MYCN_CONSERVED/NB_MYCN_CONSERVED_top_5000_nes.txt' % (projectFolder),
                     '%senhancerPromoter/H2171_MYC_REGIONS_H2171_MYC/H2171_MYC_REGIONS_H2171_MYC_top_5000_nes.txt' % (projectFolder),
                     '%senhancerPromoter/MM1S_MYC_REGIONS_MM1S_MYC_DMSO/MM1S_MYC_REGIONS_MM1S_MYC_DMSO_top_5000_nes.txt' % (projectFolder),
                     '%senhancerPromoter/P493-6_T24_MYC_REGIONS_P493-6_T24_MYC/P493-6_T24_MYC_REGIONS_P493-6_T24_MYC_top_5000_nes.txt' % (projectFolder),
                     ]

    names_list = ['NB_MYCN_CONSERVED','H2171','MM1S','P493-6_T24']
    output_path = '%sMYC_HIGH_NES.txt' % (nes_folder)
    makeNESTable(nes_path_list,names_list,output_path)
                     




#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~MAKING THE COMBINED NES TABLE~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeNESTable(nes_path_list,names_list,output =''):

    '''
    combines the GSEA NES output from the enhancerPromoter analysis
    creates a table of all represented gene sets
    '''

    if len(nes_path_list) != len(names_list):
        print('please provide the same number of nes table paths and sample names')
        sys.exit()

    #nested dictionaries gsea['pathway']['name][nes,fdr]
    gsea_dict = defaultdict(dict)

    pathway_list = []
    #iterate once to just get all potential pathways
    for i in range(len(nes_path_list)):

        nes_path = nes_path_list[i]
        nes = utils.parseTable(nes_path,'\t')
        for line in nes[1:]:
            pathway_list.append(line[0])


    pathway_list = utils.uniquify(pathway_list)
    pathway_list.sort()

    #now blank the dictionary w defaul NES of 0 and FDR of 1
    for pathway in pathway_list:
        for name in names_list:
            gsea_dict[pathway][name] = [0.0,1.0]

    #now loop again to fill out properly
    for i in range(len(nes_path_list)):

        nes_path = nes_path_list[i]
        name = names_list[i]
        nes = utils.parseTable(nes_path,'\t')
        
        for line in nes[1:]:
            if line[2] == 'NA':
                continue
            try:
                nes_vector = [float(line[2]),float(line[3])]
            except ValueError:
                print(line)
                print(nes_path)
                print(name)
                sys.exit()
            pathway = line[0]
            gsea_dict[pathway][name] = nes_vector


    #set up the output table
    header = ['PATHWAY']
    for name in names_list:
        header += ['%s_NES' % (name),'%s_FDR' % (name)]

    nes_table = [header]
    for pathway in pathway_list:
        nes_line = [pathway]
        for name in names_list:
            nes_line += gsea_dict[pathway][name]
        nes_table.append(nes_line)

    if len(output) != 0:
        utils.unParseTable(nes_table,output,'\t')
    else:
        return nes_table


        



#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
