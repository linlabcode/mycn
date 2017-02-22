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


#Main method run script for processing of THMYCN data
#folder is changed to prevent any horrific genome clashing

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
genome ='mm9'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/grail/projects/mycn_resub/%s/thmycn/' % (projectName) #PATH TO YOUR PROJECT FOLDER
projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff_mm9/' % (projectFolder)
macsFolder = '%smacsFolder_mm9/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched_mm9/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched_mm9/' % (projectFolder)
mappedFolder = '%smappedFolder_mm9/' % (projectFolder)
wiggleFolder = '%swiggles_mm9/' % (projectFolder)
metaFolder = '%smeta_mm9/' % (projectFolder)
metaRoseFolder = '%smeta_rose_mm9/' % (projectFolder)
roseFolder = '%srose_mm9/' % (projectFolder)
fastaFolder = '%sfasta_mm9/' % (projectFolder)
bedFolder = '%sbed_mm9/' % (projectFolder)
figuresFolder = '%sfigures_mm9/' % (projectFolder)
geneListFolder = '%sgeneListFolder_mm9/' % (projectFolder)
bedFolder = '%sbeds_mm9/' % (projectFolder)
signalFolder = '%ssignalTables_mm9/' % (projectFolder)
tableFolder = '%stables_mm9/' % (projectFolder)

#mask Files


#genomeDirectory
genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ChIP-Seq
mouse_dataFile = '%sdata_tables_mm9/THMYCN_TABLE.txt' % (projectFolder)




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
    pipeline_dfci.summary(mouse_dataFile)



    print('\n\n')
    print('#======================================================================')
    print('#==========================II. CALLING MACS============================')
    print('#======================================================================')
    print('\n\n')

    #running peak finding using macs 1.4.2 on all chip datasets
    #this usually takes ~2-3 hours on a reasonably fast machine
    #a 3 hour time out on this entire operation is set
    #if peak calling takes longer than 3 hours, simply run the script again after completion
    #run_macs(mouse_dataFile)


    print('\n\n')
    print('#======================================================================')
    print('#=================II. DEFINING ACTIVE GENES IN MOUSE===================')
    print('#======================================================================')
    print('\n\n')

    
    #here we will identify active promoters in various contexts as those with 
    #an H3K27AC peak in the +/- 1kb tss region
    #UCSC refseq annotations are used for all genes

    #make_active_gene_lists(mouse_dataFile)


    print('\n\n')
    print('#======================================================================')
    print('#==================III. CALLING ROSE TO MAP ENHANCERS==================')
    print('#======================================================================')
    print('\n\n')

    #for SCG_H3K27AC
    analysisName = 'SCG_H3K27AC'
    namesList = ['SCG_H3K27Ac']
    bashFileName,region_map_path,namesList=define_enhancer_landscape(mouse_dataFile,analysisName,namesList)


    #for CG_H3K27AC
    analysisName = 'CG_H3K27AC'
    namesList = ['CG_H3K27Ac']
    bashFileName,region_map_path,namesList=define_enhancer_landscape(mouse_dataFile,analysisName,namesList)


    #for GANGLIA_H3K27AC
    analysisName = 'GANGLIA_H3K27AC'
    namesList = ['CG_H3K27Ac','SCG_H3K27Ac']
    bashFileName,region_map_path,namesList=define_enhancer_landscape(mouse_dataFile,analysisName,namesList)

    #for THMYCN
    analysisName = 'THMYCN_H3K27AC'
    namesList = ['THMYCN_139076_H3K27Ac','THMYCN_139423_H3K27Ac','THMYCN1_H3K27Ac']
    bashFileName,region_map_path,namesList=define_enhancer_landscape(mouse_dataFile,analysisName,namesList)


    print('\n\n')
    print('#======================================================================')
    print('#======================VII. MAPPING TO REGIONS=========================')
    print('#======================================================================')
    print('\n\n')

    # #mapping ctcf to ctcf regions
    # gffList = ['%sHG19_SHEP21_CTCF_RX_UNION_-0_+0.gff' % (gffFolder), '%sHG19_SHEP21_CTCF_RX_INTERSECT_-0_+0.gff' % (gffFolder)]
    # names_list = ['SHEP21_0HR_CTCF_RX','SHEP21_2HR_CTCF_RX','SHEP21_24HR_CTCF_RX','SHEP21_0HR_INPUT_RX_2','SHEP21_2HR_INPUT_RX_2','SHEP21_24HR_INPUT_RX_2']
    # map_regions(shep21_chiprx_dataFile,gffList,names_list)

    # #mapping h3k4me3 to h3k4me3 regions
    # gffList = ['%sHG19_SHEP21_H3K4ME3_RX_UNION_-0_+0.gff' % (gffFolder), '%sHG19_SHEP21_H3K4ME3_RX_INTERSECT_-0_+0.gff' % (gffFolder)]
    # names_list = ['SHEP21_0HR_H3K4ME3_RX','SHEP21_2HR_H3K4ME3_RX','SHEP21_24HR_H3K4ME3_RX','SHEP21_0HR_INPUT_RX_2','SHEP21_2HR_INPUT_RX_2','SHEP21_24HR_INPUT_RX_2']
    # map_regions(shep21_chiprx_dataFile,gffList,names_list)
    
    # #mapping everybody to active TSS locations
    # gffList = ['%sHG19_TSS_NB_H3K27AC_ACTIVE_UNION_-1000_+1000.gff' % (gffFolder)]
    # map_regions(shep21_chiprx_dataFile,gffList,names_list=[])


    # #mapping everybody to mycn peaks
    # gffList = ['%sHG19_TSS_NB_H3K27AC_ACTIVE_UNION_-1000_+1000.gff' % (gffFolder)]
    # gffList = ['%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder),
    #            '%sHG19_NB_MYCN_CONSERVED_-500_+500.gff' % (gffFolder),
    #            '%sHG19_NB_MYCN_CONSERVED_ENHANCER_-0_+0.gff' % (gffFolder),
    #            '%sHG19_NB_MYCN_CONSERVED_ENHANCER_-500_+500.gff' % (gffFolder),
    #            '%sHG19_NB_MYCN_CONSERVED_PROMOTER_-0_+0.gff' % (gffFolder),
    #            '%sHG19_NB_MYCN_CONSERVED_PROMOTER_-500_+500.gff' % (gffFolder),
    #            ]
    # map_regions(shep21_chiprx_dataFile,gffList,names_list=[])





#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RUNNING MACS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def run_macs(dataFile):
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    namesList.sort()
    print(namesList)
    pipeline_dfci.callMacs(dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9')
    os.chdir(projectFolder) # the silly call macs script has to change into the output dir
    #so this takes us back to the project folder

    #to check for completeness, we will try to find all of the peak files
    peak_calling_done = False
    while not peak_calling_done:
        dataDict = pipeline_dfci.loadDataTable(dataFile)
        namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
        for name in namesList:
            peak_path = '%s%s/%s_summits.bed' % (macsFolder,name,name)
            print('searching for %s' % (peak_path))
            if utils.checkOutput(peak_path,1,180):
                peak_calling_done =True
                print('found %s' % (peak_path))
                continue
            else:
                print('Error: peak calling timed out')
                sys.exit()
    
    #now format the macs output
    print('formatting macs output')
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink ='',useBackground=True)
    print('Finished running Macs 1.4.2')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~DEFINING ACTIVE GENES IN OTHER LINES~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#now we need to make active gene lists in MM1S, P493-6, SCLC, and the SHEP ON system
def make_active_gene_lists(mouse_dataFile):

    '''
    making tss gffs and defining active genes
    '''
    
    #first define gffs
    pipeline_dfci.makeGeneGFFs(annotFile,gffFolder,species=genome.upper())

    #first map to enriched for each dataset

    #for thmycn models
    dataDict = pipeline_dfci.loadDataTable(mouse_dataFile)
    setName = 'MOUSE_TSS_H3K27AC'
    gffList = ['%sMM9_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    cellTypeList = ['CG','SCG','THMYCN','THMYCN1','THMYCN2']
    namesList = [name for name in dataDict.keys() if name.upper().count('H3K27AC') == 1]
    print(namesList)

    pipeline_dfci.mapEnrichedToGFF(mouse_dataFile,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)


    #====================
    #now create the gene lists
    #this is for THMYCN
    mappedEnrichedFile = '%sMM9_TSS_ALL_-1000_+1000/MM9_TSS_ALL_-1000_+1000_MOUSE_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['CG_H3K27Ac']]
    output = '%sgeneListFolder_mm9/MM9_CG_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)

    setList = [['SCG_H3K27Ac']]
    output = '%sgeneListFolder_mm9/MM9_SCG_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)

    setList = [['SCG_H3K27Ac'],['CG_H3K27Ac']]
    output = '%sgeneListFolder_mm9/MM9_GANGLIA_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)
    

    setList = [['THMYCN_139076_H3K27Ac'],['THMYCN_139423_H3K27Ac'],['THMYCN1_H3K27Ac'],['THMYCN2_H3K27Ac']]
    output = '%sgeneListFolder_mm9/MM9_THMYCN_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~DEFINING NB H3K27AC ENHANCER LANDSCAPE~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_enhancer_landscape(mouse_dataFile,analysisName,namesList=[]):

    '''
    define enhancers using h3k27ac in the 3 datasets that look good:
    CG, SCG, THMYCN_139076 using regular ROSE2
    '''


    #For SCG baseline
    #no TSS exclusion and no stitching

    dataDict = pipeline_dfci.loadDataTable(mouse_dataFile)
    
    if len(namesList) == 0:
        namesList = [name for name in dataDict.keys() if name.upper().count('H3K27AC') == 1]

    bamFileList = [dataDict[name]['bam'] for name in namesList]
    bamString = string.join(bamFileList,',')

    controlBams = [dataDict[name]['background'] for name in namesList]
    controlFileList = [dataDict[name]['bam'] for name in controlBams]
    controlBamString = string.join(controlFileList,',')

    bedFileList = [macsEnrichedFolder + dataDict[name]['enrichedMacs'] for name in namesList]
    bedString = string.join(bedFileList,',')

    
    outputFolder = '%s%s/' % (metaRoseFolder,analysisName)
    bashFileName = '%s%s_meta_rose.sh' % (metaRoseFolder,analysisName)

    bashFile = open(bashFileName,'w')
    bashFile.write('#!/usr/bin/bash\n\n')
    bashFile.write('cd %s\n' % (pipeline_dir))

    metaRoseCmd = 'python %sROSE2_META.py -g hg19 -i %s -r %s -c %s -o %s -n %s' % (pipeline_dir,bedString,bamString,controlBamString,outputFolder,analysisName)

    bashFile.write(metaRoseCmd + '\n')
    bashFile.close()

    region_map_path = '%s%s/%s_AllEnhancers.table.txt' % (metaRoseFolder,analysisName,analysisName)


    #runs only if no output detected
    if not utils.checkOutput(region_map_path,0,0):
        print(bashFileName)
        os.system('bash %s' % (bashFileName))
    return bashFileName,region_map_path,namesList



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~DEFINING NB MYCN BINDING LANDSCAPE~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_mycn_landscape(projectFolder,pipeline_dir,nb_all_chip_dataFile):

    '''
    defines the MYCN baseline using MYCN chips from NGP, KELLY, BE2C, and SHEP21
    uses the meta rose code and writes out a .sh file for reproducibility
    '''

    #For MYCN baseline
    #no TSS exclusion and no stitching

    dataDict = pipeline_dfci.loadDataTable(nb_all_chip_dataFile)
    analysisName = 'NB_MYCN'
    namesList = [name for name in dataDict.keys() if name.count('MYCN') == 1]

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
#~~~~~~~~~~~~~~~~~~~~~~~~~MAKING REGION GFFS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeStitchedGFF(dataFile,set_name,names_list):

    '''
    makes a stitched gff and dumps it into the gff folder
    '''

    dataDict = pipeline_dfci.loadDataTable(dataFile)
    
    loci = []
    collection_dict = {}
    for name in names_list:
        print(name)
        macsEnrichedFile = '%s%s_peaks_filtered.bed' % (macsEnrichedFolder,name)
        collection = utils.importBoundRegion(macsEnrichedFile,name)
        collection_dict[name]=collection
        loci+= collection.getLoci()

    all_collection = utils.LocusCollection(loci,50)
    stitched_collection = all_collection.stitchCollection()

    gff = utils.locusCollectionToGFF(stitched_collection)

    out_path = '%sHG19_%s_UNION_-0_+0.gff' % (gffFolder,set_name)
    print(out_path)
    utils.unParseTable(gff,out_path,'\t')


    #now get the intersect gff
    print('getting intersection gff')
    stitched_loci = stitched_collection.getLoci()
    intersect_loci = []
    ticker = 0
    for locus in stitched_loci:
        if ticker%1000==0:
            print(ticker)
        ticker+=1
        overlap = True
        for name in names_list:
            if len(collection_dict[name].getOverlap(locus,'both')) == 0:
                overlap = False

        if overlap == True:
            intersect_loci.append(locus)

            
    print('identified %s stitched loci' % (len(stitched_loci)))
    print('identified %s intersect loci' % (len(intersect_loci)))

    intersect_collection = utils.LocusCollection(intersect_loci,50)

    intersect_gff = utils.locusCollectionToGFF(intersect_collection)

    intersect_path = '%sHG19_%s_INTERSECT_-0_+0.gff' % (gffFolder,set_name)
    print(intersect_path)
    utils.unParseTable(intersect_gff,intersect_path,'\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~MAKING ACTIVS TSS GFFS~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_tss_gff(gene_list_path,name):

    '''
    makes a +/- 1kb tss gff from the provided gene list
    '''

    tss_gff = utils.parseTable('%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder),'\t')

    gene_list_table = utils.parseTable(gene_list_path,'\t')

    gene_rows = [int(line[0]) - 1 for line in gene_list_table]

    row_gff =[tss_gff[row] for row in gene_rows]

    row_gff_path = '%sHG19_TSS_%s_-1000_+1000.gff' % (gffFolder,name)

    utils.unParseTable(row_gff,row_gff_path,'\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~MAKING MYCN PROMOTER AND ENHANCER REGIONS~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_mycn_gffs(mycn_stats_path,window=0):

    '''
    makes promoter and enhancer gffs from the mycn stats table
    with an appropriate flanking window
    returns paths
    '''

    window = int(window)

    mycn_stats_table = utils.parseTable(mycn_stats_path,'\t')
    
    enhancer_gff = []
    promoter_gff = []
    for line in mycn_stats_table[1:]:
        
        gff_line = [line[1],line[0],'',int(line[2]) - window,int(line[3]) + window,'','.','',line[0]]

        if int(line[5]) == 1:
            promoter_gff.append(gff_line)
        
        if int(line[6]) == 1:
            enhancer_gff.append(gff_line)

    enhancer_gff_path = '%sHG19_NB_MYCN_CONSERVED_ENHANCER_-%s_+%s.gff' % (gffFolder,window,window)
    promoter_gff_path = '%sHG19_NB_MYCN_CONSERVED_PROMOTER_-%s_+%s.gff' % (gffFolder,window,window)
    utils.unParseTable(promoter_gff,promoter_gff_path,'\t')
    utils.unParseTable(enhancer_gff,enhancer_gff_path,'\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAPPING CHIPRX TO REGIONS~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def map_regions(nb_all_chip_dataFile,gffList,names_list=[]):

    '''
    making a normalized binding signal table at all regions
    '''

    #since each bam has different read lengths, important to carefully normalize quantification
    dataDict = pipeline_dfci.loadDataTable(nb_all_chip_dataFile)

    if len(names_list) == 0:
        names_list = dataDict.keys()
    names_list.sort()

    
    for name in names_list:
        bam = utils.Bam(dataDict[name]['bam'])
        read_length = bam.getReadLengths()[0]
        bam_extension = 200-read_length
        print('For dataset %s using an extension of %s' % (name,bam_extension))
        pipeline_dfci.mapBamsBatch(nb_all_chip_dataFile,gffList,mappedFolder,overWrite =False,namesList = [name],extension=bam_extension,rpm=True)

        

    #want a signal table of all datasets to each gff
    print('Writing signal tables for each gff:')
    for gffFile in gffList:
        gffName = gffFile.split('/')[-1].split('.')[0]
        signal_table_path = '%s%s_SIGNAL.txt' % (signalFolder,gffName)
        print(signal_table_path)
        pipeline_dfci.makeSignalTable(nb_all_chip_dataFile,gffFile,mappedFolder,namesList = names_list,medianNorm=False,output =signal_table_path)



#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
