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
rnaFolder = '%srna_seq/' % (projectFolder)

#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#gft file for RNA-Seq
gtfFile = '%sgtf/genes_ercc.gtf' % (projectFolder)

#genomeDirectory
genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

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
    print('#====================II. PROCESSING RNA_SEQ BAMS=======================')
    print('#======================================================================')
    print('\n\n')


    shep21_bash_path = process_shep_rna(shep21_rna_dataFile,gtfFile)

    shep21_drop_rep_bash_path = process_shep_rna_drop_rep(shep21_rna_dataFile,gtfFile)

    be2c_drug_bash_path = process_be2c_drug_rna(be2c_rna_drug_dataFile,gtfFile)

    be2c_twist_drug_path = process_be2c_twist_rna(be2c_rna_twist_dataFile,gtfFile)




#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis

def makeCuffTable(dataFile,analysisName,gtfFile,cufflinksFolder,groupList=[],bashFileName = ''):

    '''
    call cuffquant on each bam individually
    and then string the cbx files into cuffnorm
    groupList = [['A_1','A_2'],['B_1','B_2']]
    '''

    def long_substr(data):
        '''
        helper function to find longest substring for group naming
        '''
        substr = ''
        if len(data) > 1 and len(data[0]) > 0:
            for i in range(len(data[0])):
                for j in range(len(data[0])-i+1):
                    if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                        substr = data[0][i:i+j]
        return substr
    
    dataDict = pipeline_dfci.loadDataTable(dataFile)

    #if no grouplist is given
    #run every dataset as a single group
    #for now assumes that every dataset given is RNA Seq
    if len(groupList) == 0:
        namesList = dataDict.keys()
        namesList.sort()
        groupList = [[x] for x in namesList]
        namesString = ','.join(namesList)

    else:
        #only a single name per group
        namesList =[]
        namesStringList = []
        groupTicker = 1
        for group in groupList:
            
            namesList+=group
            coreName = long_substr(group)
            if len(coreName) ==0:
                coreName = '%s_GROUP_%s' % (analysisName,groupTicker)
            else:
                if '-_.'.count(coreName[-1]) == 1:  #get rid of any separators for a core name
                    coreName = coreName[:-1]
            namesStringList.append(coreName)
            groupTicker+=1
        namesString = ','.join(namesStringList)
            
    cufflinksFolder = utils.formatFolder(cufflinksFolder,True)

    #let's do this in bashfile format
    if len(bashFileName) ==0:
        bashFileName = '%scuffquant.sh' % (cufflinksFolder)
        
    
    bashFile = open(bashFileName,'w')

    bashFile.write('#!/usr/bin/bash\n')

    bashFile.write('cd %s\n\n' % (cufflinksFolder))

    bashFile.write("echo 'making cuffquant folders'\n")

    for name in namesList:
        bashFile.write('mkdir %s\n' % (name))

    bashFile.write("\necho 'calling cuffquant'\n")

    cuffquantList = [] # create a list to store cuffquant .cxb outputs so we can check for completeness
    for name in namesList:
        bamFileName = dataDict[name]['bam']
        bashFile.write('cuffquant -p 4 -o %s%s/ %s %s --library-type fr-firststrand\n' % (cufflinksFolder,name,gtfFile,bamFileName))
        cuffquantList.append('%s%s/abundances.cxb' % (cufflinksFolder,name))


    cxbList = []
    for group in groupList:
        
        groupString = ','.join(['%s%s/abundances.cxb' % (cufflinksFolder,name) for name in group])
        cxbList.append(groupString)

    cxbString = ' '.join(cxbList)

    #set up the analysis output folders
    cuffnormFolder = utils.formatFolder('%s%s_cuffnorm' % (cufflinksFolder,analysisName),True)
    rOutputFolder = utils.formatFolder('%s%s_cuffnorm/output/' % (cufflinksFolder,analysisName),True)


    #now run the cuffnorm    
    bashFile.write("\necho 'running cuffnorm command'\n")

    
    cuffNormCmd = 'cuffnorm -p 4 -o %s%s_cuffnorm/ -L %s %s %s\n' % (cufflinksFolder,analysisName,namesString,gtfFile,cxbString)

    bashFile.write(cuffNormCmd + '\n')


    #now we'll want to pipe the output into the R script for RNA_Seq normalization
    geneFPKMFile = '%s%s_cuffnorm/genes.fpkm_table' % (cufflinksFolder,analysisName)


    
    rCmd = 'R --no-save %s %s %s %s TRUE < %snormalizeRNASeq.R\n' % (geneFPKMFile,rOutputFolder,analysisName,namesString,pipeline_dir)

    bashFile.write(rCmd)
    bashFile.close()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~SHEP21 RNA-SEQ PROCESSING~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def process_shep_rna(shep21_rna_dataFile,gtfFile):

    '''
    quantifies gene expression to the hg19 ucsc refseq genes_ercc.gtf
    which has the spike included
    '''

    analysisName = 'SHEP21'

    cufflinksFolder = utils.formatFolder('%sshep21_cufflinks/' % (rnaFolder),True)

    groupList = [['SHEP21_0HR_rep1','SHEP21_0HR_rep2','SHEP21_0HR_rep3'],
                 ['SHEP21_2HR_rep1','SHEP21_2HR_rep2','SHEP21_2HR_rep3'],
                 ['SHEP21_4HR_rep1','SHEP21_4HR_rep2','SHEP21_4HR_rep3'],
                 ['SHEP21_6HR_rep1','SHEP21_6HR_rep2','SHEP21_6HR_rep3'],
                 ['SHEP21_8HR_rep1','SHEP21_8HR_rep2','SHEP21_8HR_rep3'],
                 ['SHEP21_16HR_rep1','SHEP21_16HR_rep2','SHEP21_16HR_rep3'],
                 ['SHEP21_24HR_rep1','SHEP21_24HR_rep2','SHEP21_24HR_rep3'],
                 ]

    bashFileName = '%sshep21_rna_seq_cuff.sh' % (cufflinksFolder)
    makeCuffTable(shep21_rna_dataFile,analysisName,gtfFile,cufflinksFolder,groupList,bashFileName)

    return bashFileName

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~SHEP21 RNA-SEQ PROCESSING WITH DROPPED REPLICATE~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def process_shep_rna_drop_rep(shep21_rna_dataFile,gtfFile):

    '''
    quantifies gene expression to the hg19 ucsc refseq genes_ercc.gtf
    which has the spike included
    '''

    analysisName = 'SHEP21'

    cufflinksFolder = utils.formatFolder('%sshep21_cufflinks_no_rep2/' % (rnaFolder),True)

    groupList = [['SHEP21_0HR_rep1','SHEP21_0HR_rep3'],
                 ['SHEP21_2HR_rep1','SHEP21_2HR_rep2','SHEP21_2HR_rep3'],
                 ['SHEP21_4HR_rep1','SHEP21_4HR_rep2','SHEP21_4HR_rep3'],
                 ['SHEP21_6HR_rep1','SHEP21_6HR_rep2','SHEP21_6HR_rep3'],
                 ['SHEP21_8HR_rep1','SHEP21_8HR_rep2','SHEP21_8HR_rep3'],
                 ['SHEP21_16HR_rep1','SHEP21_16HR_rep2','SHEP21_16HR_rep3'],
                 ['SHEP21_24HR_rep1','SHEP21_24HR_rep2','SHEP21_24HR_rep3'],
                 ]

    bashFileName = '%sshep21_rna_seq_cuff_no_rep2.sh' % (cufflinksFolder)
    makeCuffTable(shep21_rna_dataFile,analysisName,gtfFile,cufflinksFolder,groupList,bashFileName)

    return bashFileName



#==========================================================================
#=========================MAKE CUFF TABLE BE2C=============================
#==========================================================================


def process_be2c_drug_rna(be2c_rna_drug_dataFile,gtfFile):

    analysisName = 'BE2C_DRUG'

    cufflinksFolder = utils.formatFolder('%sbe2c_drug_cufflinks/' % (rnaFolder),True)

    groupList = [['BE2C_DMSO_A1','BE2C_DMSO_A2','BE2C_DMSO_A3'],
                 ['BE2C_JQ1_4HR_1','BE2C_JQ1_4HR_2','BE2C_JQ1_4HR_3'],
                 ['BE2C_JQ1_8HR_1','BE2C_JQ1_8HR_2','BE2C_JQ1_8HR_3'],
                 ['BE2C_JQ1_24HR_1','BE2C_JQ1_24HR_2','BE2C_JQ1_24HR_3'],
                 ['BE2C_DMSO_B1','BE2C_DMSO_B2','BE2C_DMSO_B3'],
                 ['BE2C_CD532_4HR_1','BE2C_CD532_4HR_2','BE2C_CD532_4HR_3'],
                 ['BE2C_CD532_8HR_1','BE2C_CD532_8HR_2','BE2C_CD532_8HR_3'],
                 ['BE2C_CD532_24HR_1','BE2C_CD532_24HR_2','BE2C_CD532_24HR_3'],
                 ]

    bashFileName = '%sbe2c_drug_rna_seq_cuff.sh' % (cufflinksFolder)
    makeCuffTable(be2c_rna_drug_dataFile,analysisName,gtfFile,cufflinksFolder,groupList,bashFileName)
    
    return bashFileName

#==========================================================================
#====================MAKE CUFF TABLE BE2C shTWIST==========================
#==========================================================================

def process_be2c_twist_rna(be2c_rna_twist_dataFile,gtfFile):


    analysisName = 'BE2C_TWIST'

    cufflinksFolder = utils.formatFolder('%sbe2c_twist_cufflinks/' % (rnaFolder),True)

    groupList = [['BE2C_shT_nodox_rep1','BE2C_shT_nodox_rep2','BE2C_shT_nodox_rep3'],
                 ['BE2C_shT_3HR_rep1','BE2C_shT_3HR_rep2','BE2C_shT_3HR_rep3'],
                 ['BE2C_shT_6HR_rep1','BE2C_shT_6HR_rep2','BE2C_shT_6HR_rep3'],
                 ['BE2C_shT_12HR_rep1','BE2C_shT_12HR_rep2','BE2C_shT_12HR_rep3'],
                 ['BE2C_shT_24HR_rep1','BE2C_shT_24HR_rep2','BE2C_shT_24HR_rep3'],
                 ['BE2C_shT_48HR_rep1','BE2C_shT_48HR_rep2','BE2C_shT_48HR_rep3'],
                 ]

    bashFileName = '%sbe2c_twist_rna_seq_cuff.sh' % (cufflinksFolder)
    makeCuffTable(be2c_rna_twist_dataFile,analysisName,gtfFile,cufflinksFolder,groupList,bashFileName)

    return bashFileName

#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
