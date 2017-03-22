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

#10_twist_dynamic_analysis.py
#Main method run script for quantifying which regions gain/lose twist upon MYCN shutdown
#and the MYCN levels at those regions

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
u87_dataFile = '%sdata_tables/U87_TABLE.txt' % (projectFolder)

chip_data_list = [be2c_dataFile,mm1s_dataFile,nb_all_chip_dataFile,p4936_young_dataFile,sclc_dataFile,shep21_dataFile,shep_on_dataFile,u87_dataFile]
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
    pipeline_dfci.summary(shep21_dataFile)




    print('\n\n')
    print('#======================================================================')
    print('#================II. RUNNING DIFFERENTIAL ROSE ANALYSIS================')
    print('#======================================================================')
    print('\n\n')

    #use the dynamic rose tools to first map twist1 binding sites
    #and then quantify
    
    name1= 'SHEP21_0HR_TWIST'
    name2= 'SHEP21_24HR_B_TWIST'
    analysis_name = 'SHEP21_TWIST1'
    rank_gff_path = wrapDRose(shep21_dataFile,name1,name2,analysis_name)


    print('\n\n')
    print('#======================================================================')
    print('#=================III. MAPPING MYCN DATA TO RANK GFF===================')
    print('#======================================================================')
    print('\n\n')

    #for shep21 nospike
    gffList = [rank_gff_path]
    dataDict = pipeline_dfci.loadDataTable(shep21_dataFile)
    names_list = [name for name in dataDict.keys() if name.count('MYCN') == 1 or name.count('INPUT') == 1 or name.count('TWIST') == 1 and name.count('rep2') == 0]
    print(names_list)
    #map_regions(shep21_dataFile,gffList,names_list)

    gffList = ['%smacsEnriched/SHEP21_0HR_TWIST_peaks.bed' % (projectFolder)]
    #map_regions(shep21_dataFile,gffList,names_list)

    #make a gff of twist and mycn sites at 0hr
    twist_collection =utils.importBoundRegion('%smacsEnriched/SHEP21_0HR_TWIST_peaks.bed' % (projectFolder),'SHEP21_0HR_TWIST')

    mycn_collection =utils.importBoundRegion('%smacsEnriched/SHEP21_0HR_MYCN_NOSPIKE_peaks.bed' % (projectFolder),'SHEP21_0HR_MYCN_NOSPIKE')

    all_loci = twist_collection.getLoci() + mycn_collection.getLoci()
    all_collection = utils.LocusCollection(all_loci,50)
    stitched_collection = all_collection.stitchCollection()

    stitched_loci = stitched_collection.getLoci()
    
    overlap_loci = []
    for locus in stitched_loci:
        if len(twist_collection.getOverlap(locus,'both')) > 0 and len(mycn_collection.getOverlap(locus,'both')) > 0:
            overlap_loci.append(locus)

    overlap_collection = utils.LocusCollection(overlap_loci,50)
    overlap_gff = utils.locusCollectionToGFF(overlap_collection)
    overlap_gff_path = '%sHG19_SHEP21_0HR_TWIST_MYCN_INTERSECTION_-0_+0.gff' % (gffFolder)
    utils.unParseTable(overlap_gff,overlap_gff_path,'\t')

    gffList = [overlap_gff_path]
    map_regions(shep21_dataFile,gffList,names_list)


    
    

#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~WRAPPING DYNAMIC ROSE ANALYSIS OF TWIST~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def wrapDRose(dataFile,name1,name2,analysis_name):

    '''
    wraps the delta rose analysis that will be done here using rose w/ 0 tss and 0 stitch
    '''

    #first call rose
    parentFolder = utils.formatFolder('%stwist1_rose/' % (projectFolder),True)

    #determine what the eventual output will look like
    enhancer_path_1 = '%s%s_ROSE/%s_peaks_AllEnhancers.table.txt' % (parentFolder,name1,name1)
    enhancer_path_2 = '%s%s_ROSE/%s_peaks_AllEnhancers.table.txt' % (parentFolder,name2,name2)
    
    if utils.checkOutput(enhancer_path_1,0.1,0.1) and utils.checkOutput(enhancer_path_2,0.1,0.1):
        print('Found ROSE2 output for %s and %s in %s' % (name1,name2,parentFolder))
    else:
        print('Running ROSE2 on %s and %s with -t 0 and -s 0 parameters')    
        bashFileName = '%s%s_rose.sh' % (parentFolder,analysis_name)
        pipeline_dfci.callRose2(dataFile,macsEnrichedFolder,parentFolder,[name1,name2],[],'',0,0,bashFileName,maskFile,True)
        
        #os.system('bash %s' % (bashFileName))

    #next run dynamic rose
    
    dynamicFolder = utils.formatFolder('%sdynamic_rose/' % (projectFolder),True)
    rose_folder_1 = '%s%s_ROSE/' % (parentFolder,name1)
    rose_folder_2 = '%s%s_ROSE/' % (parentFolder,name2)
    bashFileName = '%s%s_dynamic.sh' % (dynamicFolder,analysis_name)
    bashFile = open(bashFileName,'w')
    bashFile.write('#!/usr/bin/bash\n\n')
    
    bashFile.write('#dynamic rose on twist datasets for %s\n\n' % (analysis_name))

    dynamic_cmd = 'python %sdynamicEnhancer.py -g %s -d %s -n %s,%s -r %s,%s -o %s%s/ -a' % (pipeline_dir,genome,dataFile,name1,name2,rose_folder_1,rose_folder_2,dynamicFolder,analysis_name)
    bashFile.write(dynamic_cmd+'\n\n')
    bashFile.close()

    rank_path = '%s%s/output/%s_%s_%s_merged_MERGED_ENHANCERS_RANK_TABLE.txt' % (dynamicFolder,analysis_name,genome.upper(),name1,name2)

    print(rank_path)
    if not utils.checkOutput(rank_path,0.1,0.1):
        #only run if you can't find the terminal output
        print('Running dynamic rose from %s' % (bashFileName))
        os.system('bash %s' % (bashFileName))
    

    if  utils.checkOutput(rank_path,1,30):
        print('Found dynamic rose output at %s' % (rank_path))
    
        rank_table= utils.parseTable(rank_path,'\t')
        rank_gff = []
        for line in rank_table[1:]:
            gff_line = [line[1],line[0],'',line[2],line[3],'','.','',line[0]]
            rank_gff.append(gff_line)
            
            
        rank_gff_path = '%s%s_%s_RANK.gff' % (gffFolder,genome.upper(),analysis_name)
        print('writing rank table as a gff to %s' % (rank_gff_path))
        utils.unParseTable(rank_gff,rank_gff_path,'\t')
        return rank_gff_path
    else:
        print('Error: operation timed out. Cannot find expected dynamic output at %s' % (rank_path))
        sys.exit()
        
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



#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
