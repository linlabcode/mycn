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


#Main method run script for quantifying MYCN enhancer promoter binding

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
    print('#===================II. ENHANCER PROMOTER FOR ALL NB===================')
    print('#======================================================================')
    print('\n\n')

    # input_path = '%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder)
    # activity_path = '%sHG19_NB_H3K27AC_ACTIVE_UNION.txt' % (geneListFolder)
    # analysis_name = 'NB_MYCN_CONSERVED'
    # nb_enhancer_promoter_bash = wrap_enhancer_promoter(nb_all_chip_dataFile,input_path,activity_path,analysis_name)
    # os.system('bash %s' % (nb_enhancer_promoter_bash))


    print('\n\n')
    print('#======================================================================')
    print('#===============III. ENHANCER PROMOTER IN SHEP21 SYSTEM================')
    print('#======================================================================')
    print('\n\n')

    # #for SHEP21 nospike
    # mycn_list = ['SHEP21_0HR_MYCN_NOSPIKE','SHEP21_2HR_MYCN_NOSPIKE','SHEP21_24HR_MYCN_NOSPIKE']
    # for mycn_name in mycn_list:
    #     input_path = '%sSHEP21_0HR_MYCN_NOSPIKE_peaks.bed' % (macsEnrichedFolder)
    #     activity_path = '%sHG19_SHEP21_0HR_H3K27AC_NOSPIKE_ACTIVE.txt' % (geneListFolder)
    #     analysis_name = 'SHEP21_0HR_MYCN_NOSPIKE_REGIONS_%s' % (mycn_name)
    #     nb_enhancer_promoter_bash = wrap_enhancer_promoter(shep21_dataFile,input_path,activity_path,analysis_name,names_list = [mycn_name])
    #     os.system('bash %s' % (nb_enhancer_promoter_bash))

    # #for SHEP21 chiprx
    # mycn_list = ['SHEP21_0HR_MYCN_RX','SHEP21_2HR_MYCN_RX','SHEP21_24HR_MYCN_RX']
    # for mycn_name in mycn_list:
    #     input_path = '%sSHEP21_0HR_MYCN_NOSPIKE_peaks.bed' % (macsEnrichedFolder)
    #     activity_path = '%sHG19_SHEP21_0HR_H3K27AC_NOSPIKE_ACTIVE.txt' % (geneListFolder)
    #     analysis_name = 'SHEP21_0HR_MYCN_NOSPIKE_REGIONS_%s' % (mycn_name)
    #     nb_enhancer_promoter_bash = wrap_enhancer_promoter(shep21_chiprx_dataFile,input_path,activity_path,analysis_name,names_list = [mycn_name])
    #     os.system('bash %s' % (nb_enhancer_promoter_bash))


    print('\n\n')
    print('#======================================================================')
    print('#===============IV. ENHANCER PROMOTER IN SHEP ON SYSTEM================')
    print('#======================================================================')
    print('\n\n')

    #for SHEP21 on
    mycn_list = ['SHEP_0HR_MYCN','SHEP_2HR_MYCN','SHEP_6HR_MYCN']
    for mycn_name in mycn_list:
        input_path = '%sSHEP21_0HR_MYCN_NOSPIKE_peaks.bed' % (macsEnrichedFolder)
        activity_path = '%sHG19_SHEP21_0HR_H3K27AC_NOSPIKE_ACTIVE.txt' % (geneListFolder)
        analysis_name = 'SHEP21_0HR_MYCN_NOSPIKE_REGIONS_%s' % (mycn_name)
        nb_enhancer_promoter_bash = wrap_enhancer_promoter(shep_on_dataFile,input_path,activity_path,analysis_name,names_list = [mycn_name])
        os.system('bash %s' % (nb_enhancer_promoter_bash))

    #for SHEP21 on @ NB conserved regions
    mycn_list = ['SHEP_0HR_MYCN','SHEP_2HR_MYCN','SHEP_6HR_MYCN']
    for mycn_name in mycn_list:
        input_path = '%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder)
        activity_path = '%sHG19_SHEP21_0HR_H3K27AC_NOSPIKE_ACTIVE.txt' % (geneListFolder)
        analysis_name = 'NB_MYCN_CONSERVED_%s' % (mycn_name)
        nb_enhancer_promoter_bash = wrap_enhancer_promoter(shep_on_dataFile,input_path,activity_path,analysis_name,names_list = [mycn_name])
        os.system('bash %s' % (nb_enhancer_promoter_bash))






    print('\n\n')
    print('#======================================================================')
    print('#================V. ENHANCER PROMOTER IN INDIVIDUAL NB=================')
    print('#======================================================================')
    print('\n\n')

    # #for BE2C, KELLY, NGP
    # mycn_list = ['BE2C','KELLY','NGP']
    # for mycn_name in mycn_list:
    #     input_path = '%s%s_MYCN_peaks.bed' % (macsEnrichedFolder,mycn_name)
    #     activity_path = '%sHG19_%s_H3K27AC_ACTIVE.txt' % (geneListFolder,mycn_name)
    #     analysis_name = '%s_MYCN' % (mycn_name)
    #     nb_enhancer_promoter_bash = wrap_enhancer_promoter(nb_all_chip_dataFile,input_path,activity_path,analysis_name,names_list = ['%s_MYCN' % (mycn_name)])
    #     os.system('bash %s' % (nb_enhancer_promoter_bash))



    print('\n\n')
    print('#======================================================================')
    print('#============VI. ENHANCER PROMOTER ANALYSIS IN OTHER CANCERS===========')
    print('#======================================================================')
    print('\n\n')

    # #for p493-6, mm1s, h2171, and h128

    # #for p493-6
    # myc_list = ['P493-6_T0_MYC','P493-6_T1_MYC','P493-6_T24_MYC']
    # for myc_name in myc_list:
    #     input_path = '%sP493-6_T24_MYC_peaks.bed' % (macsEnrichedFolder)
    #     activity_path = '%sHG19_P493-6_T24_H3K27AC_ACTIVE.txt' % (geneListFolder)
    #     analysis_name = 'P493-6_T24_MYC_REGIONS_%s' % (myc_name)
    #     enhancer_promoter_bash = wrap_enhancer_promoter(p4936_young_dataFile,input_path,activity_path,analysis_name,names_list = [myc_name])
    #     os.system('bash %s' % (enhancer_promoter_bash))


    # #for sclc
    # myc_list = ['H128_MYC','H2171_MYC']
    # for myc_name in myc_list:
    #     input_path = '%sH2171_MYC_peaks.bed' % (macsEnrichedFolder)
    #     activity_path = '%sHG19_H2171_H3K27AC_ACTIVE.txt' % (geneListFolder)
    #     analysis_name = 'H2171_MYC_REGIONS_%s' % (myc_name)
    #     enhancer_promoter_bash = wrap_enhancer_promoter(sclc_dataFile,input_path,activity_path,analysis_name,names_list = [myc_name])
    #     os.system('bash %s' % (enhancer_promoter_bash))


    # #for MM
    # myc_list = ['MM1S_MYC_DMSO']
    # for myc_name in myc_list:
    #     input_path = '%sMM1S_MYC_DMSO_peaks.bed' % (macsEnrichedFolder)
    #     activity_path = '%sHG19_MM1S_H3K27AC_ACTIVE.txt' % (geneListFolder)
    #     analysis_name = 'MM1S_MYC_REGIONS_%s' % (myc_name)
    #     enhancer_promoter_bash = wrap_enhancer_promoter(mm1s_dataFile,input_path,activity_path,analysis_name,names_list = [myc_name])
    #     os.system('bash %s' % (enhancer_promoter_bash))





#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~RUNNING ENHANCER PROMOTER ANALYSIS FOR NB CONSERVED PEAKS~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def wrap_enhancer_promoter(dataFile,input_path,activity_path,analysis_name,names_list = []):

    '''
    runs enhancer promoter on everybody with the conserved regions and union of active genes
    '''
    
    #hard coded paths
    tads_path ='%shESC_domains_hg19.bed' %(bedFolder)

    #setting the output folder
    ep_folder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)
    


    dataDict = pipeline_dfci.loadDataTable(dataFile)
    if len(names_list) == 0:
        names_list = [name for name in dataDict.keys() if name.count('MYC') > 0]
        names_list.sort()

    bams_list = [dataDict[name]['bam'] for name in names_list]
    bams_string = ' '.join(bams_list)
    
    background_names = [dataDict[name]['background'] for name in names_list]
    background_list = [dataDict[background_name]['bam'] for background_name in background_names]
    background_string = ' '.join(background_list)


    ep_bash_path = '%s%s_enhancer_promoter.sh' % (ep_folder,analysis_name)
    ep_bash = open(ep_bash_path,'w')

    ep_bash.write('#!/usr/bin/bash\n\n\n')
    
    ep_bash.write('#enhancer promoter analysis for %s\n\n' % (analysis_name))
    
    python_cmd = 'python %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 2000\n\n' % (pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)
    
    ep_bash.write(python_cmd)

    python_cmd = 'python %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 5000\n\n' % (pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)
    
    ep_bash.write(python_cmd)
    ep_bash.close()
    
    return(ep_bash_path)


#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
