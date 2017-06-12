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

pipeline_dir = '/storage/cylin/home/cl6/pipeline/'

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
#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory
genomeDirectory = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

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

    # #for enhancers
    # enhancer_bashFileName,enhancer_region_map_path,namesList = define_enhancer_landscape(projectFolder,pipeline_dir,shep_on_dataFile)
    # print(enhancer_bashFileName)
    # #runs only if no output detected
    # if not utils.checkOutput(enhancer_region_map_path,0,0):
    #     print(enhancer_bashFileName)
    #     os.system('bash %s' % (enhancer_bashFileName))


    # #in individual systems
    # bash_path = map_shep_enhancers(shep_on_dataFile)
    # #os.system('bash %s' % (bash_path))


    print('\n\n')
    print('#======================================================================')
    print('#=======================III. MAP MYC LANDSCAPE=========================')
    print('#======================================================================')
    print('\n\n')

    # #for mycn
    # myc_bashFileName,myc_region_map_path,namesList = define_myc_landscape(projectFolder,pipeline_dir,shep_on_dataFile)

    # if not utils.checkOutput(myc_region_map_path,0,0):
    #     print(myc_bashFileName)
    #     os.system('bash %s' % (myc_bashFileName))



    print('\n\n')
    print('#======================================================================')
    print('#===================IV. MAKING +/- 5KB MYCN GFFs=======================')
    print('#======================================================================')
    print('\n\n')

    #make_shep_on_mycn_landscape(shep_on_dataFile)    



    print('\n\n')
    print('#======================================================================')
    print('#================V. MAPPING MYCN GFFs FOR METAS AND HEATMAP============')
    print('#======================================================================')
    print('\n\n')
    

    #mapping at shep on defined regions and same regions from 3_shep21_chiprx_heatmap
    #those regions are defined in shep21 data
    #map_shep_for_heatmap(shep_on_dataFile)

    print('\n\n')
    print('#======================================================================')
    print('#==================VI. MAPPING MYCN GFFs FOR BOX PLOT==================')
    print('#======================================================================')
    print('\n\n')

    #mapping @ a 1 bin scale for the shep21 conserved mycn regions

    gffList = [
               '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-1kb_+1kb.gff' % (gffFolder),
               '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb.gff' % (gffFolder),
               '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb.gff' % (gffFolder),
               '%sHG19_SHEP_MYCN_CONSERVED_-1kb_+1kb.gff' % (gffFolder),
               '%sHG19_SHEP_MYCN_CONSERVED_PROMOTER_-1kb_+1kb.gff' % (gffFolder),
               '%sHG19_SHEP_MYCN_CONSERVED_ENHANCER_-1kb_+1kb.gff' % (gffFolder),
        ]

    #map_regions(shep_on_dataFile,gffList,names_list=[])


    print('\n\n')
    print('#======================================================================')
    print('#=====================VII. MAKING BOX PLOTS============================')
    print('#======================================================================')
    print('\n\n')


    set_name = 'SHEP_MYCN'
    gff_name = 'SHEP_MYCN_CONSERVED_ENHANCER_-1kb_+1kb'
    names_list = ['SHEP_0HR_MYCN','SHEP_2HR_MYCN','SHEP_6HR_MYCN']
    makeBoxPlot(shep_on_dataFile,set_name,gff_name,names_list)

    set_name = 'SHEP_H3K27AC'
    names_list = ['SHEP_0HR_H3K27AC','SHEP_2HR_H3K27AC','SHEP_6HR_H3K27AC']
    makeBoxPlot(shep_on_dataFile,set_name,gff_name,names_list)



    set_name = 'SHEP_MYCN'
    gff_name = 'SHEP_MYCN_CONSERVED_PROMOTER_-1kb_+1kb'
    names_list = ['SHEP_0HR_MYCN','SHEP_2HR_MYCN','SHEP_6HR_MYCN']
    makeBoxPlot(shep_on_dataFile,set_name,gff_name,names_list)

    set_name = 'SHEP_H3K27AC'
    names_list = ['SHEP_0HR_H3K27AC','SHEP_2HR_H3K27AC','SHEP_6HR_H3K27AC']
    makeBoxPlot(shep_on_dataFile,set_name,gff_name,names_list)



    print('\n\n')
    print('#======================================================================')
    print('#===============VII. MAKING HEATMAPS AND METAS ========================')
    print('#======================================================================')
    print('\n\n')

    # #==========================================
    # #for shep mycn at shep mycn defined sites
    # plot_name = 'SHEP_MYCN'
    # names_list = ['SHEP_0HR_MYCN','SHEP_2HR_MYCN','SHEP_6HR_MYCN']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP_MYCN_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP_MYCN_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP_MYCN_CONSERVED_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'red'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)


    #==========================================
    # #for shep mycn at shep mycn defined sites
    # plot_name = 'SHEP_H3K27AC'
    # names_list = ['SHEP_0HR_H3K27AC','SHEP_2HR_H3K27AC','SHEP_6HR_H3K27AC']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP_MYCN_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP_MYCN_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP_MYCN_CONSERVED_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'blue'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)









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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~DEFINING ACTIVE ENHANCERS IN INDIVIDUAL NB SYSTEMS~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def map_shep_enhancers(shep_on_dataFile):
    '''
    for enhancers in individual systems defined by k27ac
    '''
    dataDict = pipeline_dfci.loadDataTable(shep_on_dataFile)
    namesList = dataDict.keys()

    print(namesList)

    parentFolder = '%senhancer_rose' % (projectFolder)
    parentFolder = utils.formatFolder(parentFolder,True)

    bashFileName = '%senhancer_rose/shep_on_enhancer_rose.sh' %(projectFolder)

    namesList = ['SHEP_0HR_H3K27AC','SHEP_2HR_H3K27AC','SHEP_6HR_H3K27AC']

    pipeline_dfci.callRose2(shep_on_dataFile,macsEnrichedFolder,parentFolder,namesList,[],'',2500,'',bashFileName,maskFile)

    return bashFileName



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~DEFINING SHEP ON MYCN SITES THAT ARE CONSERVED IN NB~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def make_shep_on_mycn_landscape(shep_on_dataFile):

    '''
    finds mycn peaks in shep21 that are conserved in nb and segregates them into promoter or enhancer
    '''
    dataDict = pipeline_dfci.loadDataTable(shep_on_dataFile)


    print('LOADING SHEP ON MYCN SITES')
    #load all of the shep_on sites
    # shep_on_gff_path = '%smeta_rose/SHEP_ON_MYC/gff/HG19_SHEP_ON_MYC_ALL_-0_+0.gff' % (projectFolder)
    # shep_on_gff = utils.parseTable(shep_on_gff_path,'\t')

    shep_on_bed_path = '%sSHEP_6HR_MYCN_peaks.bed' % (macsEnrichedFolder)
    shep_on_bed = utils.parseTable(shep_on_bed_path,'\t')
    shep_on_gff = utils.bedToGFF(shep_on_bed)
    
    #now get the conserved NB MYCN regions
    nb_conserved_mycn_gff_file = '%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder)
    nb_conserved_mycn_collection = utils.gffToLocusCollection(nb_conserved_mycn_gff_file)

    print('LOADING SHEP ACTIVE ENHANCERS') 
    #make a collection of enhancers
    shep_enhancer_file = '%smeta_rose/SHEP_ON_H3K27AC/SHEP_ON_H3K27AC_AllEnhancers.table.txt' % (projectFolder)
    shep_enhancer_collection = utils.makeSECollection(shep_enhancer_file,'SHEP_H3K27AC')

    #now get the active promoters
    print('LOADING SHEP ACTIVE PROMOTERS')
    startDict = utils.makeStartDict(annotFile)
    shep_transcribed_file = '%sHG19_SHEP_ON_H3K27AC_ACTIVE.txt' % (geneListFolder)
    shep_transcribed_table = utils.parseTable(shep_transcribed_file,'\t')
    transcribedList = [line[1] for line in shep_transcribed_table]
    tssLoci = []
    for refID in transcribedList:
        tssLoci.append(utils.makeTSSLocus(refID,startDict,1000,1000))

    shep_tss_collection = utils.LocusCollection(tssLoci,50)

    #now initialize the 6 gffs we will need
    shep_mycn_gff = [] 
    shep_mycn_gff_5kb = []
    shep_mycn_gff_1kb = []

    shep_mycn_promoter_gff = []
    shep_mycn_promoter_gff_1kb = []
    shep_mycn_promoter_gff_5kb = []

    shep_mycn_enhancer_gff = []
    shep_mycn_enhancer_gff_1kb = []
    shep_mycn_enhancer_gff_5kb = []

    #and their respective file names
    shep_mycn_gff_file = '%sHG19_SHEP_MYCN_CONSERVED_-0_+0.gff' % (gffFolder)
    shep_mycn_gff_5kb_file = '%sHG19_SHEP_MYCN_CONSERVED_-5kb_+5kb.gff' % (gffFolder)
    shep_mycn_gff_1kb_file = '%sHG19_SHEP_MYCN_CONSERVED_-1kb_+1kb.gff' % (gffFolder)

    shep_mycn_promoter_gff_file = '%sHG19_SHEP_MYCN_CONSERVED_PROMOTER_-0_+0.gff' % (gffFolder)
    shep_mycn_promoter_gff_5kb_file = '%sHG19_SHEP_MYCN_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder)
    shep_mycn_promoter_gff_1kb_file = '%sHG19_SHEP_MYCN_CONSERVED_PROMOTER_-1kb_+1kb.gff' % (gffFolder)

    shep_mycn_enhancer_gff_file = '%sHG19_SHEP_MYCN_CONSERVED_ENHANCER_-0_+0.gff' % (gffFolder)
    shep_mycn_enhancer_gff_5kb_file = '%sHG19_SHEP_MYCN_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder)
    shep_mycn_enhancer_gff_1kb_file = '%sHG19_SHEP_MYCN_CONSERVED_ENHANCER_-1kb_+1kb.gff' % (gffFolder)

    print('ITERATING THROUGH SHEP MYCN PEAKS')

    ticker = 0
    enhancer = 0
    promoter = 0 

    other = 0
    for line in shep_on_gff:
        if ticker % 1000 == 0:
            print ticker
        ticker+=1
        peakID = '%s_%s' % ('SHEP_MYCN',str(ticker))

        lineLocus = utils.Locus(line[0],line[3],line[4],'.',peakID)

        if nb_conserved_mycn_collection.getOverlap(lineLocus):
            gffLine = [line[0],peakID,peakID,line[3],line[4],'','.','',peakID]
            peakCenter = (int(line[3]) + int(line[4]))/2
            gffLine_5kb = [line[0],peakID,peakID,peakCenter - 5000,peakCenter + 5000,'','.','',peakID]
            #the 1kb is not a center +/- but a flank
            gffLine_1kb = [line[0],peakID,peakID,int(line[3]) - 1000,int(line[4]) + 1000,'','.','',peakID]

            shep_mycn_gff.append(gffLine)
            shep_mycn_gff_5kb.append(gffLine_5kb)
            shep_mycn_gff_1kb.append(gffLine_1kb)

            #tss overlap should take precedence over enhancer overlap
            if shep_tss_collection.getOverlap(lineLocus,'both'):
                shep_mycn_promoter_gff.append(gffLine)
                shep_mycn_promoter_gff_5kb.append(gffLine_5kb)
                shep_mycn_promoter_gff_1kb.append(gffLine_1kb)
                promoter+=1
            #now check for enhancer overlap
            elif shep_enhancer_collection.getOverlap(lineLocus,'both'):
                shep_mycn_enhancer_gff.append(gffLine)
                shep_mycn_enhancer_gff_5kb.append(gffLine_5kb)
                shep_mycn_enhancer_gff_1kb.append(gffLine_1kb)
                enhancer+=1
            else:
                other+=1
    
    print('Of %s shep on mycn peaks' % (len(shep_on_gff)))
    print('%s are promoter' % (promoter))
    print('%s are enhancer' % (enhancer))
    print('%s are other' % (other))
    #now write out the gffs
    utils.unParseTable(shep_mycn_gff,shep_mycn_gff_file,'\t')
    utils.unParseTable(shep_mycn_gff_5kb,shep_mycn_gff_5kb_file,'\t')
    utils.unParseTable(shep_mycn_gff_1kb,shep_mycn_gff_1kb_file,'\t')

    utils.unParseTable(shep_mycn_promoter_gff,shep_mycn_promoter_gff_file,'\t')
    utils.unParseTable(shep_mycn_promoter_gff_5kb,shep_mycn_promoter_gff_5kb_file,'\t')
    utils.unParseTable(shep_mycn_promoter_gff_1kb,shep_mycn_promoter_gff_1kb_file,'\t')

    utils.unParseTable(shep_mycn_enhancer_gff,shep_mycn_enhancer_gff_file,'\t')
    utils.unParseTable(shep_mycn_enhancer_gff_5kb,shep_mycn_enhancer_gff_5kb_file,'\t')
    utils.unParseTable(shep_mycn_enhancer_gff_1kb,shep_mycn_enhancer_gff_1kb_file,'\t')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~MAPPING SIGNAL AT SHEP21 MYCN SITES FOR METAS~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def map_shep_for_heatmap(shep_on_dataFile):

    '''
    map to both chiprx and regular chip bams
    to make quantification easier, all bams read lengths extended to 200

    '''

    dataDict = pipeline_dfci.loadDataTable(shep_on_dataFile)

    #gff files

    shep_mycn_conserved_gff_5kb_file = '%sHG19_SHEP_MYCN_CONSERVED_-5kb_+5kb.gff' % (gffFolder)

    shep_mycn_conserved_promoter_gff_5kb_file = '%sHG19_SHEP_MYCN_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder)

    shep_mycn_conserved_enhancer_gff_5kb_file = '%sHG19_SHEP_MYCN_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder)

    #shep21 gff files
    shep21_mycn_conserved_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-5kb_+5kb.gff' % (gffFolder)

    shep21_mycn_conserved_promoter_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder)

    shep21_mycn_conserved_enhancer_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder)


    #setting the list of gff's to map
    gffList = [shep_mycn_conserved_gff_5kb_file,
               shep_mycn_conserved_promoter_gff_5kb_file,
               shep_mycn_conserved_enhancer_gff_5kb_file,
               shep21_mycn_conserved_gff_5kb_file,
               shep21_mycn_conserved_promoter_gff_5kb_file,
               shep21_mycn_conserved_enhancer_gff_5kb_file]
    cellTypeList = ['SHEP']
    mapList = [] # map everything

    #for the non spike in
    #note, this data is 75bp reads
    pipeline_dfci.mapBams(shep_on_dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,rpm=True,nameList = mapList,extension=125)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAPPING REGIONS FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def map_regions(dataFile,gffList,names_list=[]):

    '''
    making a normalized binding signal table at all regions
    '''

    #since each bam has different read lengths, important to carefully normalize quantification
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    dataFile_name = dataFile.split('/')[-1].split('.')[0]

    if len(names_list) == 0:
        names_list = dataDict.keys()
    names_list.sort()
    
    for name in names_list:
        bam = utils.Bam(dataDict[name]['bam'])
        read_length = bam.getReadLengths()[0]
        bam_extension = 200-read_length
        print('For dataset %s using an extension of %s' % (name,bam_extension))
        pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = [name],extension=bam_extension,rpm=True)

    #want a signal table of all datasets to each gff
    print('Writing signal tables for each gff:')
    for gffFile in gffList:
        gffName = gffFile.split('/')[-1].split('.')[0]
        signal_table_path = '%s%s_%s_SIGNAL.txt' % (signalFolder,gffName,dataFile_name)
        print(signal_table_path)
        pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = names_list,medianNorm=False,output =signal_table_path)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~WRAPPING BOXPLOT R CODE~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def makeBoxPlot(dataFile,set_name,gff_name,names_list=[]):
        
    '''
    wrapping the boxplot script
    '''

    boxplot_script_path = '%sr_scripts/4_chiprx_plots.R' % (projectFolder)
    scale_table_path = '%sHG19_SHEP21_CHIPRX_SCALE_FACTORS.txt' % (tableFolder)

    dataDict=  pipeline_dfci.loadDataTable(dataFile)
    dataFile_name = dataFile.split('/')[-1].split('.')[0]
    if len(names_list) == 0:
        names_list = [name for name in dataDict.keys() if name.count(set_name) > 0]
        names_list.sort()

    background_list = [ dataDict[name]['background'] for name in names_list]
    names_string = ','.join(names_list)
    background_string = ','.join(background_list)


    signal_table_path = '%sHG19_%s_%s_SIGNAL.txt' % (signalFolder,gff_name,dataFile_name)
    
    plot_name = '%s_%s' % (gff_name,set_name)
    r_cmd = 'Rscript %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder)
    print(r_cmd)
    
    os.system(r_cmd)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~WRAPPING HEATMAP AND META R CODE~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeHeatmap(names_list,gff_list,plot_name,plot_color):

    '''
    wrapper for the heatmap and meta R script
    '''
    meta_heat_script = '%sr_scripts/5_chiprx_heatmaps.R' % (projectFolder)
    scale_table_path = '%sHG19_SHEP21_CHIPRX_SCALE_FACTORS.txt' % (tableFolder)


    names_string = ','.join(names_list)
    
    for gff in gff_list:
        gffName = gff.split('/')[-1].split('.')[0]
        mapped_list = ['%s%s/%s_%s.gff' % (mappedFolder,gffName,gffName,name) for name in names_list]
        mapped_string = ','.join(mapped_list)
        
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s' % (meta_heat_script,mapped_string,scale_table_path,names_string,plot_color,gffName,plot_name,'TRUE',projectFolder)
        print(r_cmd)
        os.system(r_cmd)


#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
