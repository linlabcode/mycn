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


#Main method run script for post alignment processing of RNA-Seq data for NB
#SHEP21 shutdown experiment, SHEP21  TWIST1 shRNA experiment
#JQ1 and CD532 in BE2C

#See README for additional information on downloading and installing dependencies

#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/home/cl6/src/pipeline/'

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
maskFolder = '%smasks/' % (projectFolder)
rnaFolder = '%srna_seq/' % (projectFolder)

#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#gft file for RNA-Seq
gtfFile = '%sgtf/genes_ercc.gtf' % (projectFolder)

#genomeDirectory
genomeDirectory = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder,maskFolder,rnaFolder]

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

    #for RNA-Seq
    pipeline_dfci.summary(shep21_rna_dataFile)
    pipeline_dfci.summary(be2c_rna_drug_dataFile)
    pipeline_dfci.summary(be2c_rna_twist_dataFile)

    
    print('\n\n')
    print('#======================================================================')
    print('#===================II, RUNNING LINE PLOT SCRIPTS======================')
    print('#======================================================================')
    print('\n\n')

    #make the folder to store output figures

    utils.formatFolder('%sfigures/6_rna_line_plots/' % (projectFolder),True)
    #we have 3 RNA-Seq datasets

    #first is shep21 at the mycn conserved regions w/ the replicate dropped
    #and at shep21 defined regions
    wrap_shep21()
    wrap_be2c_jq1()

    wrap_be2c_cd532()
    #now we need to run it w/ the replicate removed
    wrap_be2c_shTwist()

#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~SHEP21 RNA-SEQ LINE PLOTS~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def wrap_shep21():

    '''
    wraps the mycn shutdown timecourse in the original dataset w/o replicate removed 
    at the regions defined by the nb conserved dataset
    '''

    #these arguments are hard wired
    rna_r_script_path = '%sr_scripts/6_rna_line_plots.R' % (projectFolder)

    group_list = ['SHEP21_0HR',
                  'SHEP21_2HR',
                  'SHEP21_4HR',
                  'SHEP21_6HR',
                  'SHEP21_8HR',
                  'SHEP21_16HR',
                  'SHEP21_24HR'
                  ]

    group_string = ','.join(group_list)
    x_vector = [0,2,4,6,8,16,24]
    x_string = ','.join([str(x) for x in x_vector])
    
    color_string = '200,0,0'

    #for shep21 nb mycn conserved regions
    analysis_name = 'SHEP21_NB_MYCN_CONSERVED'
    gene_table_path = '%senhancerPromoter/NB_MYCN_CONSERVED/NB_MYCN_CONSERVED_GENE_TABLE.txt' % (projectFolder)
    exp_norm_path = '%srna_seq/shep21_cufflinks_no_rep2/SHEP21_cuffnorm/output/SHEP21_all_fpkm_exprs_norm.txt' % (projectFolder)
    exp_raw_path = '%srna_seq/shep21_cufflinks_no_rep2/SHEP21_cuffnorm/output/SHEP21_all_fpkm_exprs_raw.txt' % (projectFolder)

    top_list = [5000,5000,'all']
    bin_list = [3,5,5]
    for i in range(3):
        nBins = bin_list[i]
        top_count = top_list[i]
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s' %(rna_r_script_path,gene_table_path,exp_norm_path,exp_raw_path,group_string,str(top_count),str(nBins),analysis_name,x_string,color_string,projectFolder)
        os.system(r_cmd)
        print(r_cmd)

    #for shep21 nb mycn conserved regions
    analysis_name = 'SHEP21_0HR_MYCN_NOSPIKE'
    gene_table_path = '%senhancerPromoter/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_SHEP21_0HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_SHEP21_0HR_MYCN_NOSPIKE_GENE_TABLE.txt' % (projectFolder)
    exp_norm_path = '%srna_seq/shep21_cufflinks_no_rep2/SHEP21_cuffnorm/output/SHEP21_all_fpkm_exprs_norm.txt' % (projectFolder)
    exp_raw_path = '%srna_seq/shep21_cufflinks_no_rep2/SHEP21_cuffnorm/output/SHEP21_all_fpkm_exprs_raw.txt' % (projectFolder)

    for i in range(3):
        nBins = bin_list[i]
        top_count = top_list[i]
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s' %(rna_r_script_path,gene_table_path,exp_norm_path,exp_raw_path,group_string,str(top_count),str(nBins),analysis_name,x_string,color_string,projectFolder)
        os.system(r_cmd)
        print(r_cmd)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~BE2C JQ1 RNA-SEQ LINE PLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def wrap_be2c_jq1():

    '''
    Wraps the be2c jq1 timecourse at the regions defined by the nb conserved dataset
    '''

    #these arguments are hard wired
    rna_r_script_path = '%sr_scripts/6_rna_line_plots.R' % (projectFolder)

    group_list = ['BE2C_DMSO_B',
                  'BE2C_JQ1_4HR',
                  'BE2C_JQ1_8HR',
                  'BE2C_JQ1_24HR',
                  ]
    group_string = ','.join(group_list)

    x_vector = [0,4,8,24]
    x_string = ','.join([str(x) for x in x_vector])
    
    color_string = '200,0,0'

    #for nb mycn conserved regions
    analysis_name = 'BE2C_JQ1_NB_MYCN_CONSERVED'
    gene_table_path = '%senhancerPromoter/NB_MYCN_CONSERVED/NB_MYCN_CONSERVED_GENE_TABLE.txt' % (projectFolder)
    exp_norm_path = '%srna_seq/be2c_drug_cufflinks/BE2C_DRUG_cuffnorm/output/BE2C_DRUG_all_fpkm_exprs_norm.txt' % (projectFolder)
    exp_raw_path = '%srna_seq/be2c_drug_cufflinks/BE2C_DRUG_cuffnorm/output/BE2C_DRUG_all_fpkm_exprs_raw.txt' % (projectFolder)

    top_list = [5000,5000,'all']
    bin_list = [3,5,5]
    for i in range(3):
        nBins = bin_list[i]
        top_count = top_list[i]
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s' %(rna_r_script_path,gene_table_path,exp_norm_path,exp_raw_path,group_string,str(top_count),str(nBins),analysis_name,x_string,color_string,projectFolder)
        os.system(r_cmd)
        print(r_cmd)

    #for be2c mycn regions
    analysis_name = 'BE2C_JQ1_BE2C_MYCN'
    gene_table_path = '%senhancerPromoter/BE2C_MYCN/BE2C_MYCN_GENE_TABLE.txt' % (projectFolder)
    exp_norm_path = '%srna_seq/be2c_drug_cufflinks/BE2C_DRUG_cuffnorm/output/BE2C_DRUG_all_fpkm_exprs_norm.txt' % (projectFolder)
    exp_raw_path = '%srna_seq/be2c_drug_cufflinks/BE2C_DRUG_cuffnorm/output/BE2C_DRUG_all_fpkm_exprs_raw.txt' % (projectFolder)

    for i in range(3):
        nBins = bin_list[i]
        top_count = top_list[i]
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s' %(rna_r_script_path,gene_table_path,exp_norm_path,exp_raw_path,group_string,str(top_count),str(nBins),analysis_name,x_string,color_string,projectFolder)
        os.system(r_cmd)
        print(r_cmd)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~BE2C CD532 RNA-SEQ LINE PLOTS~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def wrap_be2c_cd532():

    '''
    Wraps the be2c jq1 timecourse at the regions defined by the nb conserved dataset
    '''

    #these arguments are hard wired
    rna_r_script_path = '%sr_scripts/6_rna_line_plots.R' % (projectFolder)

    group_list = ['BE2C_DMSO_B',
                  'BE2C_CD532_4HR',
                  'BE2C_CD532_8HR',
                  'BE2C_CD532_24HR',
                  ]
    group_string = ','.join(group_list)

    x_vector = [0,4,8,24]
    x_string = ','.join([str(x) for x in x_vector])
    
    color_string = '200,0,0'

    #for nb mycn conserved regions
    analysis_name = 'BE2C_CD532_NB_MYCN_CONSERVED'
    gene_table_path = '%senhancerPromoter/NB_MYCN_CONSERVED/NB_MYCN_CONSERVED_GENE_TABLE.txt' % (projectFolder)
    exp_norm_path = '%srna_seq/be2c_drug_cufflinks/BE2C_DRUG_cuffnorm/output/BE2C_DRUG_all_fpkm_exprs_norm.txt' % (projectFolder)
    exp_raw_path = '%srna_seq/be2c_drug_cufflinks/BE2C_DRUG_cuffnorm/output/BE2C_DRUG_all_fpkm_exprs_raw.txt' % (projectFolder)

    top_list = [5000,5000,'all']
    bin_list = [3,5,5]
    for i in range(3):
        nBins = bin_list[i]
        top_count = top_list[i]
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s' %(rna_r_script_path,gene_table_path,exp_norm_path,exp_raw_path,group_string,str(top_count),str(nBins),analysis_name,x_string,color_string,projectFolder)
        os.system(r_cmd)
        print(r_cmd)

    #for be2c mycn regions
    analysis_name = 'BE2C_CD532_BE2C_MYCN'
    gene_table_path = '%senhancerPromoter/BE2C_MYCN/BE2C_MYCN_GENE_TABLE.txt' % (projectFolder)
    exp_norm_path = '%srna_seq/be2c_drug_cufflinks/BE2C_DRUG_cuffnorm/output/BE2C_DRUG_all_fpkm_exprs_norm.txt' % (projectFolder)
    exp_raw_path = '%srna_seq/be2c_drug_cufflinks/BE2C_DRUG_cuffnorm/output/BE2C_DRUG_all_fpkm_exprs_raw.txt' % (projectFolder)

    for i in range(3):
        nBins = bin_list[i]
        top_count = top_list[i]
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s' %(rna_r_script_path,gene_table_path,exp_norm_path,exp_raw_path,group_string,str(top_count),str(nBins),analysis_name,x_string,color_string,projectFolder)
        os.system(r_cmd)
        print(r_cmd)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~BE2C CD532 RNA-SEQ LINE PLOTS~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def wrap_be2c_shTwist():

    '''
    Wraps the be2c jq1 timecourse at the regions defined by the nb conserved dataset
    '''

    #these arguments are hard wired
    rna_r_script_path = '%sr_scripts/6_rna_line_plots.R' % (projectFolder)

    group_list = ['BE2C_shT_nodox',
                  'BE2C_shT_3HR',
                  'BE2C_shT_6HR',
                  'BE2C_shT_12HR',
                  'BE2C_shT_24HR',
                  'BE2C_shT_48HR',
                  ]
    group_string = ','.join(group_list)

    x_vector = [0,3,6,12,24,48]
    x_string = ','.join([str(x) for x in x_vector])
    
    color_string = '200,0,0'

    #for nb mycn conserved regions
    analysis_name = 'BE2C_TWIST_NB_MYCN_CONSERVED'
    gene_table_path = '%senhancerPromoter/NB_MYCN_CONSERVED/NB_MYCN_CONSERVED_GENE_TABLE.txt' % (projectFolder)
    exp_norm_path = '%srna_seq/be2c_twist_cufflinks/BE2C_TWIST_cuffnorm/output/BE2C_TWIST_all_fpkm_exprs_norm.txt' % (projectFolder)
    exp_raw_path = '%srna_seq/be2c_twist_cufflinks/BE2C_TWIST_cuffnorm/output/BE2C_TWIST_all_fpkm_exprs_raw.txt' % (projectFolder)


    top_list = [5000,5000,'all']
    bin_list = [3,5,5]
    for i in range(3):
        nBins = bin_list[i]
        top_count = top_list[i]
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s' %(rna_r_script_path,gene_table_path,exp_norm_path,exp_raw_path,group_string,str(top_count),str(nBins),analysis_name,x_string,color_string,projectFolder)
        os.system(r_cmd)
        print(r_cmd)

    #for be2c mycn regions
    analysis_name = 'BE2C_TWIST_BE2C_MYCN'
    gene_table_path = '%senhancerPromoter/BE2C_MYCN/BE2C_MYCN_GENE_TABLE.txt' % (projectFolder)
    exp_norm_path = '%srna_seq/be2c_twist_cufflinks/BE2C_TWIST_cuffnorm/output/BE2C_TWIST_all_fpkm_exprs_norm.txt' % (projectFolder)
    exp_raw_path = '%srna_seq/be2c_twist_cufflinks/BE2C_TWIST_cuffnorm/output/BE2C_TWIST_all_fpkm_exprs_raw.txt' % (projectFolder)

    for i in range(3):
        nBins = bin_list[i]
        top_count = top_list[i]
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s' %(rna_r_script_path,gene_table_path,exp_norm_path,exp_raw_path,group_string,str(top_count),str(nBins),analysis_name,x_string,color_string,projectFolder)
        os.system(r_cmd)
        print(r_cmd)




#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
