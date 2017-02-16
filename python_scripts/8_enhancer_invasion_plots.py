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

#8_enhancer_invasion_plots.py
#Main method run script for quantifying mycn fold changes at tss vs. distal regions

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
    #these are the datasets we will use
    pipeline_dfci.summary(shep_on_dataFile)
    pipeline_dfci.summary(shep21_dataFile)
    pipeline_dfci.summary(shep21_chiprx_dataFile)



    print('\n\n')
    print('#======================================================================')
    print('#=========================II. MAKE BOXPLOTS============================')
    print('#======================================================================')
    print('\n\n')

    #here we will wrap boxplots for each set of analysis

    region_prefix = 'SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE' #this is used to find the peak tables
    set_name = 'SHEP_MYCN' # this is the defacto title for the datasets
    scale_table_path = ''
    wrapInvasionBox(shep_on_dataFile,region_prefix,set_name,names_list = [],top=5000,scale_path = scale_table_path)

    region_prefix = 'SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE' #this is used to find the peak tables
    set_name = 'SHEP_MYCN_NOSPIKE' # this is the defacto title for the datasets
    scale_table_path = ''
    myc_list = ['SHEP21_0HR_MYCN_NOSPIKE','SHEP21_2HR_MYCN_NOSPIKE','SHEP21_24HR_MYCN_NOSPIKE']
    wrapInvasionBox(shep21_dataFile,region_prefix,set_name,names_list = myc_list,top=5000,scale_path = scale_table_path)

    region_prefix = 'SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE' #this is used to find the peak tables
    set_name = 'SHEP_MYCN_RX_NO_SCALE' # this is the defacto title for the datasets
    scale_table_path = ''
    myc_list = ['SHEP21_0HR_MYCN_RX','SHEP21_2HR_MYCN_RX','SHEP21_24HR_MYCN_RX']
    wrapInvasionBox(shep21_chiprx_dataFile,region_prefix,set_name,names_list = myc_list,top=5000,scale_path = scale_table_path)

    region_prefix = 'SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE' #this is used to find the peak tables
    set_name = 'SHEP_MYCN_RX' # this is the defacto title for the datasets
    scale_table_path = '%sHG19_SHEP21_CHIPRX_SCALE_FACTORS.txt' % (tableFolder)
    myc_list = ['SHEP21_0HR_MYCN_RX','SHEP21_2HR_MYCN_RX','SHEP21_24HR_MYCN_RX']
    wrapInvasionBox(shep21_chiprx_dataFile,region_prefix,set_name,names_list = myc_list,top=5000,scale_path = scale_table_path)




#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~WRAPPING ENHANCER INVASION BOXPLOT R CODE~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def wrapInvasionBox(data_file,region_prefix,set_name,names_list = [],top=5000,scale_path =''):



    '''
    wrapper for the enhancer invasion boxplots
    '''

    
    invasion_script = '%sr_scripts/7_enhancer_invasion_plots.R' % (projectFolder)

    #set the scale path default
    if len(scale_path) == 0:
        scale_path = 'NONE'

    dataDict = pipeline_dfci.loadDataTable(data_file)
    if len(names_list) == 0:
        names_list = [name for name in dataDict.keys() if name.count('MYC') > 0]
        names_list.sort()
    
    print('running enhancer invasion analysis on:')
    print(names_list)
    
    print('anchoring analysis on dataset: %s' % (names_list[0]))
    
    #need to get paths of the three peak tables
    #assumes formatting and naming conventions of the enhancerPromoter folder (see 5_nb_enhancer_promoter.py)
    peak_0_path = '%senhancerPromoter/%s_%s/%s_%s_PEAK_TABLE.txt' % (projectFolder,region_prefix,names_list[0],region_prefix,names_list[0])
    peak_1_path = '%senhancerPromoter/%s_%s/%s_%s_PEAK_TABLE.txt' % (projectFolder,region_prefix,names_list[1],region_prefix,names_list[1])
    peak_2_path = '%senhancerPromoter/%s_%s/%s_%s_PEAK_TABLE.txt' % (projectFolder,region_prefix,names_list[2],region_prefix,names_list[2])

    analysis_name = '%s_%s' % (region_prefix,set_name)
    print(analysis_name)

    sample_string = ','.join(names_list)
    print(sample_string)

    r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s' % (invasion_script,peak_0_path,peak_1_path,peak_2_path,analysis_name,sample_string,top,projectFolder,scale_path)
    
    print(r_cmd)
    os.system(r_cmd)


#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
