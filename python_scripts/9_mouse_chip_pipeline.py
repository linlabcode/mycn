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
genome ='mm9'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/mycn_resub/%s/thmycn/' % (projectName) #PATH TO YOUR PROJECT FOLDER

hg19_projectFolder = '/storage/cylin/grail/projects/mycn_resub/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

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

    # #for SCG_H3K27AC
    # analysisName = 'SCG_H3K27AC'
    # namesList = ['SCG_H3K27Ac']
    # bashFileName,region_map_path,namesList=define_enhancer_landscape(mouse_dataFile,analysisName,namesList)


    # #for CG_H3K27AC
    # analysisName = 'CG_H3K27AC'
    # namesList = ['CG_H3K27Ac']
    # bashFileName,region_map_path,namesList=define_enhancer_landscape(mouse_dataFile,analysisName,namesList)


    # #for GANGLIA_H3K27AC
    # analysisName = 'GANGLIA_H3K27AC'
    # namesList = ['CG_H3K27Ac','SCG_H3K27Ac']
    # bashFileName,region_map_path,namesList=define_enhancer_landscape(mouse_dataFile,analysisName,namesList)

    # #for THMYCN
    # analysisName = 'THMYCN_H3K27AC'
    # namesList = ['THMYCN_139076_H3K27Ac','THMYCN_139423_H3K27Ac','THMYCN1_H3K27Ac']
    # bashFileName,region_map_path,namesList=define_enhancer_landscape(mouse_dataFile,analysisName,namesList)

    print('\n\n')
    print('#======================================================================')
    print('#=================IV. LIFTING OVER NB CONSERVED REGIONS================')
    print('#======================================================================')
    print('\n\n')

    # #liftover a pair of gffs
    # #first convert to bed
    # nb_promoter_gff_path = '%sgff/HG19_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000.gff' % (hg19_projectFolder)
    # nb_enhancer_gff_path = '%sgff/HG19_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000.gff' % (hg19_projectFolder)

    # nb_promoter_bed_path ='%sbeds/HG19_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000.bed' % (hg19_projectFolder)
    # nb_enhancer_bed_path ='%sbeds/HG19_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000.bed' % (hg19_projectFolder)

    # nb_promoter_gff = utils.parseTable(nb_promoter_gff_path,'\t')
    # nb_enhancer_gff = utils.parseTable(nb_enhancer_gff_path,'\t')


    # utils.gffToBed(nb_promoter_gff,nb_promoter_bed_path)
    # utils.gffToBed(nb_enhancer_gff,nb_enhancer_bed_path)

    # print('converted NB conserved gffs to beds at %s and %s' % (nb_promoter_bed_path,nb_enhancer_bed_path))

    # #note, now you have to liftover manually to create beds
    # mm9_promoter_bed_path = '%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000.bed' % (bedFolder)
    # mm9_enhancer_bed_path = '%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000.bed' % (bedFolder)

    # mm9_promoter_gff_path = '%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000.gff' % (gffFolder)
    # mm9_enhancer_gff_path = '%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000.gff' % (gffFolder)
    
    # utils.bedToGFF(mm9_promoter_bed_path,mm9_promoter_gff_path)
    # utils.bedToGFF(mm9_enhancer_bed_path,mm9_enhancer_gff_path)

    # print('writing mm9 nb mycn sites to %s and %s' % (mm9_promoter_gff_path,mm9_enhancer_gff_path))


    print('\n\n')
    print('#======================================================================')
    print('#======================V. MAPPING ENRICHED TO GFFS=====================')
    print('#======================================================================')
    print('\n\n')

    # setName = 'THMYCN'
    # gffList = [mm9_promoter_gff_path,mm9_enhancer_gff_path]
    # cellTypeList = ['THMYCN1','THMYCN2','THMYCN','CG','SCG']
    # mapList = ['CG_H3K27Ac',
    #             'SCG_H3K27Ac',
    #             'THMYCN1_H3K27Ac',
    #             'THMYCN_139423_H3K27Ac',
    #             'THMYCN_139076_H3K27Ac',
    #             ]

    # #pipeline_dfci.mapEnrichedToGFF(mouse_dataFile,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,macs=True,namesList=mapList,useBackground=True)

    # #summarize info for venn diagrams for each

    # promoter_mapped_path = '%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000/MM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000_THMYCN.txt' % (mappedEnrichedFolder)
    # promoter_venn_path = '%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000_VENN.txt' % (tableFolder)
    # summarizeVenn(promoter_mapped_path,group_list = ['CG','THMYCN'],output=promoter_venn_path)


    # enhancer_mapped_path = '%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000/MM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000_THMYCN.txt' % (mappedEnrichedFolder)
    # enhancer_venn_path = '%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000_VENN.txt' % (tableFolder)
    # summarizeVenn(enhancer_mapped_path,group_list = ['CG','THMYCN'],output=enhancer_venn_path)


    print('\n\n')
    print('#======================================================================')
    print('#=====================VI. MAKING MYCN REGIONS GFF======================')
    print('#======================================================================')
    print('\n\n')

    dataDict = pipeline_dfci.loadDataTable(mouse_dataFile)
    names_list = ['THMYCN2_MYCN',
                  'THMYCN_139076_MYCN',
                  'THMYCN_139423_MYCN',
                  ]
    
    mycn_loci = []
    for name in names_list:
        mycn_collection = utils.importBoundRegion('%s%s' % (macsEnrichedFolder,dataDict[name]['enrichedMacs']),name)
        mycn_loci+=mycn_collection.getLoci()

    mycn_collection = utils.LocusCollection(mycn_loci,50)
    mycn_collection.stitchCollection()
    mycn_gff = utils.locusCollectionToGFF(mycn_collection)
    mycn_gff_path = '%sMM9_THMYCN_MYCN_-0_+0.gff' % (gffFolder)
    utils.unParseTable(mycn_gff,mycn_gff_path,'\t')

    #make collections
    promoter_collection = utils.gffToLocusCollection('%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000.gff' % (gffFolder))
    enhancer_collection = utils.gffToLocusCollection('%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000.gff' % (gffFolder))
    #make the overlap table
    overlap_table = [['PROMOTER','ENHANCER','NONE']]
    promoter_count = 0
    enhancer_count = 0
    none_count = 0
    for line in mycn_gff:
        locus = utils.Locus(line[0],int(line[3])-10000,int(line[4])+10000,'.')
        if enhancer_collection.getOverlap(locus,'both'):
            enhancer_count +=1
            continue

        if promoter_collection.getOverlap(locus,'both'):
            promoter_count +=1
        else:
            none_count +=1

    overlap_table.append([promoter_count,enhancer_count,none_count])
    overlap_table_path = '%sMM9_THMYCN_OVERLAP.txt' % (tableFolder)
    utils.unParseTable(overlap_table,overlap_table_path,'\t')
        

    

    print('\n\n')
    print('#======================================================================')
    print('#=====================VI. MAPPING GFFS FOR HEATMAP=====================')
    print('#======================================================================')
    print('\n\n')

    #map_for_heatmap(mouse_dataFile)


    print('\n\n')
    print('#======================================================================')
    print('#=====================VII. AVERAGING MAPPED SIGNAL=====================')
    print('#======================================================================')
    print('\n\n')


    # set_list = ['GANGLIA_H3K27AC','THMYCN_H3K27AC','THMYCN_MYCN']
    # set_names = [
    #     ['CG_H3K27Ac','SCG_H3K27Ac'],
    #     ['THMYCN1_H3K27Ac','THMYCN_139423_H3K27Ac','THMYCN_139076_H3K27Ac'],
    #     ['THMYCN2_MYCN','THMYCN_139076_MYCN','THMYCN_139423_MYCN']
    # ]
    # for i in range(len(set_list)):
    #     setName = set_list[i]
    #     names_list =set_names[i]
    #     print(setName)
    #     print(names_list)
    #     #for promoters
    #     mapped_list = ['%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000/MM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000_%s.gff' % (mappedFolder,name) for name in names_list]
    #     output_path = '%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000/MM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000_%s.gff' % (mappedFolder,setName)
    #     print(output_path)
    #     averagingMappedSignal(mapped_list,output_path,setName)
        
    #     #for enhancers
    #     mapped_list = ['%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000/MM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000_%s.gff' % (mappedFolder,name) for name in names_list]
    #     output_path = '%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000/MM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000_%s.gff' % (mappedFolder,setName)
    #     print(output_path)
    #     averagingMappedSignal(mapped_list,output_path,setName)


    print('\n\n')
    print('#======================================================================')
    print('#=====================VIII. MAKING HEATMAPS/METAS======================')
    print('#======================================================================')
    print('\n\n')


    # set_list = ['THMYCN_MYCN','GANGLIA_H3K27AC','THMYCN_H3K27AC']
    # gff_list = ['%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000.gff' % (gffFolder),
    #             '%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000.gff' % (gffFolder),
    #             ]
    # plot_name = 'THMYCN_MYCN_ORDERED_BLUE'
    # plot_color = 'blue'
    # makeHeatmap(set_list,gff_list,plot_name,plot_color)


    # set_list = ['THMYCN_MYCN','THMYCN_MYCN','THMYCN_MYCN']
    # gff_list = ['%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000.gff' % (gffFolder),
    #             '%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000.gff' % (gffFolder),
    #             ]
    # plot_name = 'THMYCN_MYCN_ORDERED_RED'
    # plot_color = 'red'
    # makeHeatmap(set_list,gff_list,plot_name,plot_color)





    # print('\n\n')
    # print('#======================================================================')
    # print('#========================IV. RUNNING DYNAMIC===========================')
    # print('#======================================================================')
    # print('\n\n')

    # dynamic_meta_folder = utils.formatFolder('%sdynamic_meta_mm9' %(projectFolder),True)
    # #wrapping dynamic
    # meta_rose_1 = '%sGANGLIA_H3K27AC/' %  (metaRoseFolder)
    # meta_rose_2 = '%sTHMYCN_H3K27AC/' %  (metaRoseFolder)
    # group1_names = ['CG_H3K27Ac','SCG_H3K27Ac']
    # group2_names = ['THMYCN1_H3K27Ac','THMYCN_139076_H3K27Ac','THMYCN_139423_H3K27Ac']
    # output_folder = utils.formatFolder('%sGANGLIA_THMYCN_SE/' % (dynamic_meta_folder),True)
    # name_1 = 'GANGLIA'
    # name_2 = 'THMYCN'
    # wrap_dynamic_meta(mouse_dataFile,meta_rose_1,meta_rose_2,output_folder,group1_names,group2_names,name_1,name_2)


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


    setList = [['SCG_H3K27Ac'],['CG_H3K27Ac'],['THMYCN_139076_H3K27Ac'],['THMYCN_139423_H3K27Ac'],['THMYCN1_H3K27Ac'],['THMYCN2_H3K27Ac']]
    output = '%sgeneListFolder_mm9/MM9_ALL_H3K27AC_ACTIVE.txt' % (projectFolder)
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

    metaRoseCmd = 'python %sROSE2_META.py -g mm9 -i %s -r %s -c %s -o %s -n %s' % (pipeline_dir,bedString,bamString,controlBamString,outputFolder,analysisName)

    bashFile.write(metaRoseCmd + '\n')
    bashFile.close()

    region_map_path = '%s%s/%s_AllEnhancers.table.txt' % (metaRoseFolder,analysisName,analysisName)


    #runs only if no output detected
    if not utils.checkOutput(region_map_path,0,0):
        print(bashFileName)
        os.system('bash %s' % (bashFileName))
    return bashFileName,region_map_path,namesList



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~SUMMARIZING MAPPED ENRICHED TO MAKE VENN DIAGRAMS~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def summarizeVenn(mapped_path,group_list = ['CG','THMYCN'],output=''):

    '''
    summarizes binary occupancy across group to make a venn diagram
    '''

    group_table = [['GFF_LINE','ID'] + group_list]
    
    mapped_table = utils.parseTable(mapped_path,'\t')
    
    group_cols = []
    for group in group_list:
        group_names = [name for name in mapped_table[0] if name.count(group) > 0]
        group_cols.append([mapped_table[0].index(name) for name in group_names])

    print(group_cols)
    for line in mapped_table[1:]:
        binary_vector = [] #a 1/0 vector to hold mapping by group
        for i in range(len(group_list)):
            cols = group_cols[i]
            signal = max([int(line[x]) for x in cols])
            binary_vector.append(signal)

        new_line = line[0:2] + binary_vector
        group_table.append(new_line)

    print(group_table[0:5])

    #now add up the stats
    #this part assumes only 2 groups for now otherwise gets combinatorially challenging
    #permute all possible binary combinations given the vector length
    binary_combinations=[[0],[1]]
    for i in range(len(group_list)-1):
        new_combinations = []
        for x in binary_combinations:
            print(x)
            x1 = list(x) + [1]
            x0 = list(x) + [0]
            new_combinations.append(x1)
            new_combinations.append(x0)


            binary_combinations = list(new_combinations)

    print(binary_combinations)
    count_table = [group_list + ['count']]
    for combo in binary_combinations:
        count = len([line for line in group_table[1:] if line[2:] == combo])

        count_table.append(combo + [count])
    print(count_table)
    if len(output) > 0:
        utils.unParseTable(count_table,output,'\t')
    else:
        return count_table

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~MAPPING SIGNAL AT SHEP21 MYCN SITES FOR METAS~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def map_for_heatmap(mouse_dataFile):

    '''
    to make quantification easier, all bams read lengths extended to 200

    '''

    dataDict = pipeline_dfci.loadDataTable(mouse_dataFile)

    #gff files
    nb_conserved_promoter_gff_5kb_file = '%sMM9_NB_MYCN_CONSERVED_PROMOTER_-5000_+5000.gff' % (gffFolder)
    nb_conserved_enhancer_gff_5kb_file = '%sMM9_NB_MYCN_CONSERVED_ENHANCER_-5000_+5000.gff' % (gffFolder)



    #setting the list of gff's to map
    gffList = [
        nb_conserved_promoter_gff_5kb_file,
        nb_conserved_enhancer_gff_5kb_file,
        ]
    cellTypeList = ['CG','SCG','THMYCN1','THMYCN2','THMYCN']
    mapList = ['CG_H3K27Ac',
               'SCG_H3K27Ac',
               'THMYCN1_H3K27Ac',
               'THMYCN_139423_H3K27Ac',
               'THMYCN_139076_H3K27Ac',
               'THMYCN2_MYCN',
               'THMYCN_139076_MYCN',
               'THMYCN_139423_MYCN',
                ]

    #for the non spike in
    #note, this data is 75bp reads
    pipeline_dfci.mapBams(mouse_dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,rpm=True,nameList = mapList,extension=125)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~AVERAGING MAPPED GFF SIGNAL~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def averagingMappedSignal(mapped_list,output_path,setName):

    '''
    averages signal across a set of mapped gffs and writes the new output
    '''

    #create a list containing all of the tables
    table_list = [utils.parseTable(mapped_list[i],'\t') for i in range(len(mapped_list))]

    #first set up the output header
    output_header = ['GENE_ID','locusLine']
    nCols = len(table_list[0][0]) - 2
    for n in range(nCols):
        output_header.append('bin_%s_%s' % (n+1,setName))
        
    output_table = [output_header]
    #now iterate through each row to set up the gene ID and locus line
    for i in range(1,len(table_list[0])):
        
        line = table_list[0][i]
        if len(line) > 2:
            output_table.append(line[0:2])

    #now run through the whole matrix in i,j notation and put average signal into the final matrix
    
    #iterate through rows
    row_ticker = 1
    for i in range(1,len(table_list[0])):
        line = table_list[0][i]
        if len(line) == 2:
            continue
        signal_vector = []
        #iterate through columns
        for j in range(2,len(table_list[0][0])):
            try:
                signal_vector = [float(table[i][j]) for table in table_list]
            except IndexError:
                print(i,j)
                print(table_list[0][i])
                print(table_list[1][i])


            signal =max(round(numpy.average(signal_vector),4),0)

            output_table[row_ticker].append(signal)
        row_ticker+=1

    print(len(table_list[0]))
    print(len(output_table))
    utils.unParseTable(output_table,output_path,'\t')
    return output_path


    





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~WRAPPING HEATMAP AND META R CODE~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeHeatmap(names_list,gff_list,plot_name,plot_color):

    '''
    wrapper for the heatmap and meta R script
    '''
    meta_heat_script = '%sr_scripts/5_chiprx_heatmaps.R' % (hg19_projectFolder)
    scale_table_path = '%stables/HG19_SHEP21_CHIPRX_SCALE_FACTORS.txt' % (hg19_projectFolder)
    figures_path = utils.formatFolder('%sfigures/' % (projectFolder),True)
    figures_path = utils.formatFolder('%sfigures/5_chiprx_heatmaps/' % (projectFolder),True)
    

    names_string = ','.join(names_list)
    
    for gff in gff_list:
        gffName = gff.split('/')[-1].split('.')[0]
        mapped_list = ['%s%s/%s_%s.gff' % (mappedFolder,gffName,gffName,name) for name in names_list]
        mapped_string = ','.join(mapped_list)
        
        r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s' % (meta_heat_script,mapped_string,scale_table_path,names_string,plot_color,gffName,plot_name,'TRUE',projectFolder)
        print(r_cmd)
        os.system(r_cmd)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~WRAPPING DYNAMIC ENHANCER META~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def wrap_dynamic_meta(mouse_dataFile,meta_rose_1,meta_rose_2,output_folder,group1_names,group2_names,name_1,name_2):
    
    '''
    wraps the dynamic meta enhancer analysis
    '''
    
    output_folder = utils.formatFolder(output_folder,True)
    group1_string = ','.join(group1_names)
    group2_string = ','.join(group2_names)

    bash_path = '%s%s_%s_dynamic.sh' % (output_folder,name_1,name_2)
    bash_file = open(bash_path,'w')
    bash_file.write('#!/usr/bin/bash\n\n\n')

    bash_file.write('cd %s\n\n' % (projectFolder))
    cmd = 'srun --mem 16000 python %sdynamicEnhancer_meta.py -g MM9 -d %s -r %s,%s -o %s --group1 %s --group2 %s --name1 %s --name2 %s' % (pipeline_dir,mouse_dataFile,meta_rose_1,meta_rose_2,output_folder,group1_string,group2_string,name_1,name_2)

    print(cmd)
    bash_file.write(cmd)
    
    bash_file.close()
    print(bash_path)

#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
