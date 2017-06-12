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
genome ='mm9'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/mycn_resub/%s/thmycn/' % (projectName) #PATH TO YOUR PROJECT FOLDER
hg19_projectFolder ='/storage/cylin/grail/projects/mycn_resub/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER
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
genePlotFolder = '%sgene_plot_mm9/' % (projectFolder)
#mask Files


#genomeDirectory
genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder,genePlotFolder]

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
    print('#============II. MAKING A BED OUT OF HG19 FIGURE REGIONS===============')
    print('#======================================================================')
    print('\n\n')

    hg19_gff_path = '%sgff/HG19_NB_FIGURE_GENES.gff' % (hg19_projectFolder)
    
    hg19_gff = utils.parseTable(hg19_gff_path,'\t')
    print(hg19_gff)

    hg19_bed = utils.gffToBed(hg19_gff)
    print(hg19_bed)
    hg19_bed_path = '%sbeds/HG19_NB_FIGURE_GENES.bed' % (hg19_projectFolder)
    utils.unParseTable(hg19_bed,hg19_bed_path,'\t')
    #need to manually lift this over to mm9
    #https://genome.ucsc.edu/cgi-bin/hgLiftOver

    mm9_bed_path = '%sMM9_NB_FIGURE_GENES_LIFTOVER.bed' % (bedFolder)
    mm9_gff_path = '%sMM9_NB_FIGURE_GENES_LIFTOVER.gff' % (gffFolder)
    mm9_gff = utils.bedToGFF(mm9_bed_path)

    #now add some additional manual regions

    added_gff_regions = [
        ['chr12','TWIST1_ENHANCER','TWIST1_ENHANCER',34639818,34656263,'','-','','TWIST1_ENHANCER'],
        ['chr11','NPM1_PROMOTER_2','NPM1_PROMOTER_2',33049820,33065883,'','+','','NPM1_PROMOTER_2'],
        ['chr6','GATA2_ENHANCER','GATA2_ENHANCER',88135802,88159867,'','+','','GATA2_ENHANCER'],
        ['chr7','PHOX2A','PHOX2A',108964211,108974610,'','+','','PHOX2A'],
        ['chr15','LET7B','LET7B',85497440,85538754,'','+','','LET7B',],
        ['chr10','LIN28B','LIN28B',45161233,45217227,'','-','','LIN28B'],
        ]

    mm9_gff_full = mm9_gff+added_gff_regions

    utils.unParseTable(mm9_gff_full,mm9_gff_path,'\t')


    print('\n\n')
    print('#======================================================================')
    print('#=======================III. PLOTTING DATA IN MOUSE====================')
    print('#======================================================================')
    print('\n\n')
    
    #plot mouse regions
    plot_mouse_genes(mouse_dataFile,mm9_gff_path)


    
#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING FOR SHEP ON SYSTEM~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_mouse_genes(mouse_dataFile,mouse_figure_gff_path):

    '''
    plots all varieties and iterations of tracks @ lifted over mouse regions
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sTHMYCN/' % (genePlotFolder),True)
    plot_prefix = 'MM9_NB_FIGURE_GENES_LIFTOVER'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(mouse_dataFile)
    names_list = dataDict.keys()
    #initial check for consistency of read lengths
    # for name in names_list:
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,mouse_dataFile,bam_extension))
    # sys.exit()

    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (mouse_dataFile,bam_extension))
    
    #first do individuals
    for plot_group in ['_MYCN','H3K27AC']:
        plotList = [name for name in dataDict.keys() if name.upper().count(plot_group) > 0]
        print(plotList)
        if plot_group == '_MYCN':
            plotName = '%s_THMYCN%s' % (plot_prefix,plot_group)
        else:
            plotName = '%s_THMYCN_%s' % (plot_prefix,plot_group)

        print(plotName)
        pipeline_dfci.callBatchPlot(mouse_dataFile,mouse_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    #now as metas
    #we only have 3 good k27ac and 3 good mycn datasets
    plotList = ['CG_H3K27Ac',
                'SCG_H3K27Ac',
                'THMYCN1_H3K27Ac',
                'THMYCN_139423_H3K27Ac',
                'THMYCN_139076_H3K27Ac',
                'THMYCN2_MYCN',
                'THMYCN_139076_MYCN',
                'THMYCN_139423_MYCN',
                ]
    groupString = 'CG_,SCG,H3K27AC,H3K27AC,H3K27AC,MYCN,MYCN,MYCN'

    plotName = '%s_THMYCN_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(mouse_dataFile,mouse_figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed ='',plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_THMYCN_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(mouse_dataFile,mouse_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed ='',plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')



#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
