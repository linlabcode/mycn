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


#Main method run script for processing of SHEP MYCN DOX ON ChIP-Seq data

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
    pipeline_dfci.summary(shep_on_dataFile)


    print('\n\n')
    print('#======================================================================')
    print('#=========================II. MAP ENHANCERS============================')
    print('#======================================================================')
    print('\n\n')

    #for enhancers
    enhancer_bashFileName,enhancer_region_map_path,namesList = define_enhancer_landscape(projectFolder,pipeline_dir,shep_on_dataFile)
    print(enhancer_bashFileName)
    #runs only if no output detected
    if not utils.checkOutput(enhancer_region_map_path,0,0):
        print(enhancer_bashFileName)
        os.system('bash %s' % (enhancer_bashFileName))

    print('\n\n')
    print('#======================================================================')
    print('#=======================III. MAP MYC LANDSCAPE=========================')
    print('#======================================================================')
    print('\n\n')

    #for mycn
    myc_bashFileName,myc_region_map_path,namesList = define_myc_landscape(projectFolder,pipeline_dir,shep_on_dataFile)

    if not utils.checkOutput(myc_region_map_path,0,0):
        print(myc_bashFileName)
        os.system('bash %s' % (myc_bashFileName))

    print('\n\n')
    print('#======================================================================')
    print('#============IV. MYCN ENHANCER PROMOTER BINDING IN SHEP ON=============')
    print('#======================================================================')
    print('\n\n')




    #big questions
    #do we see enhancer invasion?
    
    #map enhancers at each time point
    
    #define myc landscape

    #quantify binding at the consensus myc landscape
    #do this using the enhancer promoter code
    #at the union of myc sites
    #or at the nb_mycn conserved sites




#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~DEFINING SHEP ON H3K27AC ENHANCER LANDSCAPE~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_enhancer_landscape(projectFolder,pipeline_dir,shep_on_dataFile):

    '''
    defines the NB enhancer baseline using H3K27ac chips from NGP, KELLY, BE2C, and SHEP21
    enhancers defined using auto optimized stitching of nearby regions
    w/ a 2.5kb tss exclusion
    uses the meta rose code and writes out a .sh file for reproducibility
    '''

    #For H3K27AC
    #with TSS exclusion and auto stitching

    dataDict = pipeline_dfci.loadDataTable(shep_on_dataFile)
    analysisName = 'SHEP_ON_H3K27AC'
    namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]

    bamFileList = [dataDict[name]['bam'] for name in namesList]
    bamString = string.join(bamFileList,',')

    controlBams = [dataDict[name]['background'] for name in namesList]
    controlFileList = [dataDict[name]['bam'] for name in controlBams]
    controlBamString = string.join(controlFileList,',')

    bedFileList = [macsEnrichedFolder + dataDict[name]['enrichedMacs'] for name in namesList]
    bedString = string.join(bedFileList,',')

    roseFolder = '%smeta_rose/' % (projectFolder)
    roseFolder = utils.formatFolder(roseFolder,True)

    outputFolder = '%s%s/' % (roseFolder,analysisName)
    bashFileName = '%s%s_meta_rose.sh' % (roseFolder,analysisName)

    bashFile = open(bashFileName,'w')
    bashFile.write('#!/usr/bin/bash\n\n')
    bashFile.write('cd %s\n' % (pipeline_dir))

    metaRoseCmd = 'python %sROSE2_META.py -g hg19 -i %s -r %s -c %s -o %s -n %s -t 2500 --mask %s' % (pipeline_dir,bedString,bamString,controlBamString,outputFolder,analysisName,maskFile)

    bashFile.write(metaRoseCmd + '\n')
    bashFile.close()

    region_map_path = '%s%s/%s_AllEnhancers.table.txt' % (roseFolder,analysisName,analysisName)
    return bashFileName,region_map_path,namesList



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~DEFINING NB MYCN BINDING LANDSCAPE~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_myc_landscape(projectFolder,pipeline_dir,shep_on_dataFile):

    '''
    defines the myc baseline in shep on system across the union of all time points
    uses the meta rose code and writes out a .sh file for reproducibility
    '''

    #For MYC baseline
    #no TSS exclusion and no stitching

    dataDict = pipeline_dfci.loadDataTable(shep_on_dataFile)
    analysisName = 'SHEP_ON_MYC'
    namesList = [name for name in dataDict.keys() if name.count('MYC') == 1]

    bamFileList = [dataDict[name]['bam'] for name in namesList]
    bamString = string.join(bamFileList,',')

    controlBams = [dataDict[name]['background'] for name in namesList]
    controlFileList = [dataDict[name]['bam'] for name in controlBams]
    controlBamString = string.join(controlFileList,',')

    bedFileList = [macsEnrichedFolder + dataDict[name]['enrichedMacs'] for name in namesList]
    bedString = string.join(bedFileList,',')

    roseFolder = '%smeta_rose/' % (projectFolder)
    roseFolder = utils.formatFolder(roseFolder,True)

    outputFolder = '%s%s/' % (roseFolder,analysisName)
    bashFileName = '%s%s_meta_rose.sh' % (roseFolder,analysisName)

    bashFile = open(bashFileName,'w')
    bashFile.write('#!/usr/bin/bash\n\n')
    bashFile.write('cd %s\n' % (pipeline_dir))

    metaRoseCmd = 'python %sROSE2_META.py -g hg19 -i %s -r %s -c %s -o %s -n %s -t 0 -s 0 --mask %s' % (pipeline_dir,bedString,bamString,controlBamString,outputFolder,analysisName,maskFile)

    bashFile.write(metaRoseCmd + '\n')
    bashFile.close()

    #this is the expeceted region map output
    region_map_path = '%s%s/%s_0KB_STITCHED_ENHANCER_REGION_MAP.txt' % (roseFolder,analysisName,analysisName)
    return bashFileName,region_map_path,namesList


#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
