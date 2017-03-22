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


#Main method run script for analysis of SHEP21 CHIPRX data to make heatmaps, metas, and boxplots

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
    pipeline_dfci.summary(shep21_chiprx_dataFile)



    print('\n\n')
    print('#======================================================================')
    print('#========II. DEFINING ACTIVE GENES AND ENHANCERS IN SHEP21=============')
    print('#======================================================================')
    print('\n\n')

    #make_shep21_active()
    
    #bash_path = map_nb_enhancers(nb_all_chip_dataFile)
    #os.system('bash %s' % (bash_path))

    print('\n\n')
    print('#======================================================================')
    print('#==================III. MAKING +/- 5KB MYCN GFFs=======================')
    print('#======================================================================')
    print('\n\n')
    #make_shep21_mycn_landscape(nb_all_chip_dataFile)

    print('\n\n')
    print('#======================================================================')
    print('#===============IV. MAPPING MYCN GFFs FOR METAS AND HEATMAP============')
    print('#======================================================================')
    print('\n\n')
    
    #with and without spike in
    #map_shep21_for_heatmap(shep21_chiprx_dataFile,shep21_dataFile)

    print('\n\n')
    print('#======================================================================')
    print('#==================V. MAPPING MYCN GFFs FOR BOX PLOT===================')
    print('#======================================================================')
    print('\n\n')

    #mapping @ a 1 bin scale for the shep21 conserved mycn regions

    # gffList = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-0_+0.gff' % (gffFolder),
    #           '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-1kb_+1kb.gff' % (gffFolder),
    #            '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb.gff' % (gffFolder),
    #            '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb.gff' % (gffFolder),
    # ]

    # gffList = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0.gff' % (gffFolder),
    #           '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0.gff' % (gffFolder),
    #     ]

    # map_regions(shep21_chiprx_dataFile,gffList,names_list=[])
    # map_regions(shep21_dataFile,gffList,names_list=[])

    print('\n\n')
    print('#======================================================================')
    print('#===============VI. MAKING HEATMAPS AND METAS =========================')
    print('#======================================================================')
    print('\n\n')

    #set the output folder
    utils.formatFolder('%sfigures/5_chiprx_heatmaps/' % projectFolder,True)

    # #==========================================
    # #for shep21 mycn chiprx
    # plot_name = 'SHEP21_MYCN_RX'
    # names_list = ['SHEP21_0HR_MYCN_RX','SHEP21_2HR_MYCN_RX','SHEP21_24HR_MYCN_RX']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'red'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)

    # #for shep21 mycn regular chip
    # plot_name = 'SHEP21_MYCN_NOSPIKE'
    # names_list = ['SHEP21_0HR_MYCN_NOSPIKE','SHEP21_2HR_MYCN_NOSPIKE','SHEP21_24HR_MYCN_NOSPIKE']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'red'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)


    # #==========================================
    # #for shep21 h3k27ac chiprx
    # plot_name = 'SHEP21_H3K27AC_RX'
    # names_list = ['SHEP21_0HR_H3K27AC_RX','SHEP21_2HR_H3K27AC_RX','SHEP21_24HR_H3K27AC_RX']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'blue'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)

    # #for shep21 mycn regular chip
    # plot_name = 'SHEP21_H3K27AC_NOSPIKE'
    # names_list = ['SHEP21_0HR_H3K27AC_NOSPIKE','SHEP21_2HR_H3K27AC_NOSPIKE','SHEP21_24HR_H3K27AC_NOSPIKE']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'blue'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)


    # #==========================================
    # #for shep21 CTCF chiprx
    # plot_name = 'SHEP21_CTCF_RX'
    # names_list = ['SHEP21_0HR_CTCF_RX','SHEP21_2HR_CTCF_RX','SHEP21_24HR_CTCF_RX']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'black'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)


    # #==========================================
    # #for shep21 RNA Pol II chiprx
    # plot_name = 'SHEP21_POL2_RX'
    # names_list = ['SHEP21_0HR_POL2_RX','SHEP21_2HR_POL2_RX','SHEP21_24HR_POL2_RX']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'black'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)


    # #==========================================
    # #for shep21 RNA Pol II NOSPIKE
    # plot_name = 'SHEP21_POL2_NOSPIKE'
    # names_list = ['SHEP21_0HR_POL2_NOSPIKE','SHEP21_2HR_POL2_NOSPIKE','SHEP21_24HR_POL2_NOSPIKE']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'black'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)


    # #==========================================
    # #for shep21 H3K4ME3 chiprx
    # plot_name = 'SHEP21_H3K4ME3_RX'
    # names_list = ['SHEP21_0HR_H3K4ME3_RX','SHEP21_2HR_H3K4ME3_RX','SHEP21_24HR_H3K4ME3_RX']
    # gff_list = ['%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder),
    #             '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder),
    #             ]
    # plot_color = 'green'
    # makeHeatmap(names_list,gff_list,plot_name,plot_color)


    print('\n\n')
    print('#======================================================================')
    print('#=====================VII. MAKING BOXPLOTS ============================')
    print('#======================================================================')
    print('\n\n')
    utils.formatFolder('%sfigures/4_chiprx_plots/' % projectFolder,True)

    # #=============================================================================
    # #for nb mycn chiprx 
    # set_name = 'MYCN_CHIPRX'
    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_MYCN_RX','SHEP21_2HR_MYCN_RX','SHEP21_24HR_MYCN_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0'
    # names_list = ['SHEP21_0HR_MYCN_RX','SHEP21_2HR_MYCN_RX','SHEP21_24HR_MYCN_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_MYCN_RX','SHEP21_2HR_MYCN_RX','SHEP21_24HR_MYCN_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0'
    # names_list = ['SHEP21_0HR_MYCN_RX','SHEP21_2HR_MYCN_RX','SHEP21_24HR_MYCN_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # #=============================================================================
    # #for nb mycn chip no spike
    # set_name = 'MYCN_NOSPIKE'
    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_MYCN_NOSPIKE','SHEP21_2HR_MYCN_NOSPIKE','SHEP21_24HR_MYCN_NOSPIKE']
    # makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0'
    # names_list = ['SHEP21_0HR_MYCN_NOSPIKE','SHEP21_2HR_MYCN_NOSPIKE','SHEP21_24HR_MYCN_NOSPIKE']
    # makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_MYCN_NOSPIKE','SHEP21_2HR_MYCN_NOSPIKE','SHEP21_24HR_MYCN_NOSPIKE']
    # makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0'
    # names_list = ['SHEP21_0HR_MYCN_NOSPIKE','SHEP21_2HR_MYCN_NOSPIKE','SHEP21_24HR_MYCN_NOSPIKE']
    # makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)



    # #=============================================================================
    # #for nb h3k27ac chiprx 
    # set_name = 'H3K27AC_CHIPRX'
    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_H3K27AC_RX','SHEP21_2HR_H3K27AC_RX','SHEP21_24HR_H3K27AC_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0'
    # names_list = ['SHEP21_0HR_H3K27AC_RX','SHEP21_2HR_H3K27AC_RX','SHEP21_24HR_H3K27AC_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_H3K27AC_RX','SHEP21_2HR_H3K27AC_RX','SHEP21_24HR_H3K27AC_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0'
    # names_list = ['SHEP21_0HR_H3K27AC_RX','SHEP21_2HR_H3K27AC_RX','SHEP21_24HR_H3K27AC_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # #=============================================================================
    # #for nb H3K27ac chip no spike
    # set_name = 'H3K27AC_NOSPIKE'
    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_H3K27AC_NOSPIKE','SHEP21_2HR_H3K27AC_NOSPIKE','SHEP21_24HR_H3K27AC_NOSPIKE']
    # makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0'
    # names_list = ['SHEP21_0HR_H3K27AC_NOSPIKE','SHEP21_2HR_H3K27AC_NOSPIKE','SHEP21_24HR_H3K27AC_NOSPIKE']
    # makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_H3K27AC_NOSPIKE','SHEP21_2HR_H3K27AC_NOSPIKE','SHEP21_24HR_H3K27AC_NOSPIKE']
    # makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)

    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0'
    # names_list = ['SHEP21_0HR_H3K27AC_NOSPIKE','SHEP21_2HR_H3K27AC_NOSPIKE','SHEP21_24HR_H3K27AC_NOSPIKE']
    # makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)


    # #=============================================================================
    # #for nb ctcf chiprx 
    # set_name = 'CTCF_CHIPRX'
    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_CTCF_RX','SHEP21_2HR_CTCF_RX','SHEP21_24HR_CTCF_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0'
    # names_list = ['SHEP21_0HR_CTCF_RX','SHEP21_2HR_CTCF_RX','SHEP21_24HR_CTCF_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_CTCF_RX','SHEP21_2HR_CTCF_RX','SHEP21_24HR_CTCF_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0'
    # names_list = ['SHEP21_0HR_CTCF_RX','SHEP21_2HR_CTCF_RX','SHEP21_24HR_CTCF_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    #=============================================================================
    #for nb RNA Pol II chiprx 
    set_name = 'POL2_RX'
    gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb'
    names_list = ['SHEP21_0HR_POL2_RX','SHEP21_2HR_POL2_RX','SHEP21_24HR_POL2_RX']
    makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0'
    names_list = ['SHEP21_0HR_POL2_RX','SHEP21_2HR_POL2_RX','SHEP21_24HR_POL2_RX']
    makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb'
    names_list = ['SHEP21_0HR_POL2_RX','SHEP21_2HR_POL2_RX','SHEP21_24HR_POL2_RX']
    makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0'
    names_list = ['SHEP21_0HR_POL2_RX','SHEP21_2HR_POL2_RX','SHEP21_24HR_POL2_RX']
    makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    #=============================================================================
    #for nb RNA Pol II chiprx 
    set_name = 'POL2_NOSPIKE'
    gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb'
    names_list = ['SHEP21_0HR_POL2_NOSPIKE','SHEP21_2HR_POL2_NOSPIKE','SHEP21_24HR_POL2_NOSPIKE']
    makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)


    gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0'
    names_list = ['SHEP21_0HR_POL2_NOSPIKE','SHEP21_2HR_POL2_NOSPIKE','SHEP21_24HR_POL2_NOSPIKE']
    makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)


    gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb'
    names_list = ['SHEP21_0HR_POL2_NOSPIKE','SHEP21_2HR_POL2_NOSPIKE','SHEP21_24HR_POL2_NOSPIKE']
    makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)


    gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0'
    names_list = ['SHEP21_0HR_POL2_NOSPIKE','SHEP21_2HR_POL2_NOSPIKE','SHEP21_24HR_POL2_NOSPIKE']
    makeBoxPlot(shep21_dataFile,set_name,gff_name,names_list)



    # #=============================================================================
    # #for nb h3k4me3 chiprx 
    # set_name = 'H3K4ME3_CHIPRX'
    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_H3K4ME3_RX','SHEP21_2HR_H3K4ME3_RX','SHEP21_24HR_H3K4ME3_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0'
    # names_list = ['SHEP21_0HR_H3K4ME3_RX','SHEP21_2HR_H3K4ME3_RX','SHEP21_24HR_H3K4ME3_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb'
    # names_list = ['SHEP21_0HR_H3K4ME3_RX','SHEP21_2HR_H3K4ME3_RX','SHEP21_24HR_H3K4ME3_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)


    # gff_name = 'SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0'
    # names_list = ['SHEP21_0HR_H3K4ME3_RX','SHEP21_2HR_H3K4ME3_RX','SHEP21_24HR_H3K4ME3_RX']
    # makeBoxPlot(shep21_chiprx_dataFile,set_name,gff_name,names_list)






#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~DEFINING ACTIVE PROMOTERS IN SHEP21~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def make_shep21_active():

    mappedEnrichedFile = '%smappedEnriched/HG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_NB_TSS_H3K27AC.txt' % (projectFolder)
    setList = [['SHEP21_0HR_H3K27AC_NOSPIKE']]
    output = '%sHG19_SHEP21_H3K27AC_TRANSCRIBED.txt' % (geneListFolder)

    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~DEFINING ACTIVE ENHANCERS IN INDIVIDUAL NB SYSTEMS~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def map_nb_enhancers(nb_all_chip_dataFile):
    '''
    for enhancers in individual systems defined by k27ac
    '''
    dataDict = pipeline_dfci.loadDataTable(nb_all_chip_dataFile)
    namesList = dataDict.keys()

    print(namesList)

    parentFolder = '%senhancer_rose' % (projectFolder)
    parentFolder = utils.formatFolder(parentFolder,True)

    bashFileName = '%senhancer_rose/nb_enhancer_rose.sh' %(projectFolder)

    namesList = ['SHEP21_0HR_H3K27AC_NOSPIKE','BE2C_H3K27AC','KELLY_H3K27AC','NGP_H3K27AC']

    pipeline_dfci.callRose2(nb_all_chip_dataFile,macsEnrichedFolder,parentFolder,namesList,[],'',2500,'',bashFileName,maskFile)

    return bashFileName


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~DEFINING SHEP21 MYCN SITES THAT ARE CONSERVED IN NB~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_shep21_mycn_landscape(nb_all_chip_dataFile):

    '''
    finds mycn peaks in shep21 that are conserved in nb and segregates them into promoter or enhancer
    '''

    #first get the shep21 regions

    print('LOADING SHEP21 MYCN SITES')
    dataDict = pipeline_dfci.loadDataTable(nb_all_chip_dataFile)
    shep21_0hr_mycn_enriched_file = '%s%s' % (macsEnrichedFolder,dataDict['SHEP21_0HR_MYCN_NOSPIKE']['enrichedMacs'])
    shep21_0hr_mycn_bed = utils.parseTable(shep21_0hr_mycn_enriched_file,'\t')


    #now get the conserved NB MYCN regions
    nb_conserved_mycn_gff_file = '%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder)
    nb_conserved_mycn_collection = utils.gffToLocusCollection(nb_conserved_mycn_gff_file)

    print('LOADING SHEP21 ACTIVE ENHANCERS')
    #make a collection of enhancers
    shep21_enhancer_file = '%senhancer_rose/SHEP21_0HR_H3K27AC_NOSPIKE_ROSE/SHEP21_0HR_H3K27AC_NOSPIKE_peaks_AllEnhancers.table.txt' % (projectFolder)
    shep21_enhancer_collection = utils.makeSECollection(shep21_enhancer_file,'SHEP21_0HR_H3K27AC_NOSPIKE')

    #now get the active promoters
    print('LOADING SHEP21 ACTIVE PROMOTERS')
    startDict = utils.makeStartDict(annotFile)
    shep21_transcribed_file = '%sHG19_SHEP21_H3K27AC_TRANSCRIBED.txt' % (geneListFolder)
    shep21_transcribed_table = utils.parseTable(shep21_transcribed_file,'\t')
    transcribedList = [line[1] for line in shep21_transcribed_table]
    tssLoci = []
    for refID in transcribedList:
        tssLoci.append(utils.makeTSSLocus(refID,startDict,1000,1000))

    shep21_tss_collection = utils.LocusCollection(tssLoci,50)

    #now initialize the 6 gffs we will need
    shep21_mycn_conserved_gff = [] 
    shep21_mycn_conserved_gff_5kb = []
    shep21_mycn_conserved_gff_1kb = []

    shep21_mycn_conserved_promoter_gff = []
    shep21_mycn_conserved_promoter_gff_1kb = []
    shep21_mycn_conserved_promoter_gff_5kb = []

    shep21_mycn_conserved_enhancer_gff = []
    shep21_mycn_conserved_enhancer_gff_1kb = []
    shep21_mycn_conserved_enhancer_gff_5kb = []

    #and their respective file names
    shep21_mycn_conserved_gff_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-0_+0.gff' % (gffFolder)
    shep21_mycn_conserved_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-5kb_+5kb.gff' % (gffFolder)
    shep21_mycn_conserved_gff_1kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-1kb_+1kb.gff' % (gffFolder)

    shep21_mycn_conserved_promoter_gff_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0.gff' % (gffFolder)
    shep21_mycn_conserved_promoter_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder)
    shep21_mycn_conserved_promoter_gff_1kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-1kb_+1kb.gff' % (gffFolder)

    shep21_mycn_conserved_enhancer_gff_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0.gff' % (gffFolder)
    shep21_mycn_conserved_enhancer_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder)
    shep21_mycn_conserved_enhancer_gff_1kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-1kb_+1kb.gff' % (gffFolder)

    print('ITERATING THROUGH SHEP21 MYCN PEAKS')

    ticker = 0
    for line in shep21_0hr_mycn_bed:
        if ticker % 1000 == 0:
            print ticker
        ticker+=1
        peakID = '%s_%s' % ('SHEP21_0HR_MYCN_NOSPIKE',str(ticker))

        lineLocus = utils.Locus(line[0],line[1],line[2],'.',peakID)

        if nb_conserved_mycn_collection.getOverlap(lineLocus):

            gffLine = [line[0],peakID,peakID,line[1],line[2],'','.','',peakID]
            peakCenter = (int(line[1]) + int(line[2]))/2
            gffLine_5kb = [line[0],peakID,peakID,peakCenter - 5000,peakCenter + 5000,'','.','',peakID]
            #the 1kb is not a center +/- but a flank
            gffLine_1kb = [line[0],peakID,peakID,int(line[1]) - 1000,int(line[2]) + 1000,'','.','',peakID]

            shep21_mycn_conserved_gff.append(gffLine)
            shep21_mycn_conserved_gff_5kb.append(gffLine_5kb)
            shep21_mycn_conserved_gff_1kb.append(gffLine_1kb)

            #tss overlap should take precedence over enhancer overlap
            if shep21_tss_collection.getOverlap(lineLocus,'both'):
                shep21_mycn_conserved_promoter_gff.append(gffLine)
                shep21_mycn_conserved_promoter_gff_5kb.append(gffLine_5kb)
                shep21_mycn_conserved_promoter_gff_1kb.append(gffLine_1kb)
            #now check for enhancer overlap
            elif shep21_enhancer_collection.getOverlap(lineLocus,'both'):
                shep21_mycn_conserved_enhancer_gff.append(gffLine)
                shep21_mycn_conserved_enhancer_gff_5kb.append(gffLine_5kb)
                shep21_mycn_conserved_enhancer_gff_1kb.append(gffLine_1kb)

    #now write out the gffs
    utils.unParseTable(shep21_mycn_conserved_gff,shep21_mycn_conserved_gff_file,'\t')
    utils.unParseTable(shep21_mycn_conserved_gff_5kb,shep21_mycn_conserved_gff_5kb_file,'\t')
    utils.unParseTable(shep21_mycn_conserved_gff_1kb,shep21_mycn_conserved_gff_1kb_file,'\t')

    utils.unParseTable(shep21_mycn_conserved_promoter_gff,shep21_mycn_conserved_promoter_gff_file,'\t')
    utils.unParseTable(shep21_mycn_conserved_promoter_gff_5kb,shep21_mycn_conserved_promoter_gff_5kb_file,'\t')
    utils.unParseTable(shep21_mycn_conserved_promoter_gff_1kb,shep21_mycn_conserved_promoter_gff_1kb_file,'\t')

    utils.unParseTable(shep21_mycn_conserved_enhancer_gff,shep21_mycn_conserved_enhancer_gff_file,'\t')
    utils.unParseTable(shep21_mycn_conserved_enhancer_gff_5kb,shep21_mycn_conserved_enhancer_gff_5kb_file,'\t')
    utils.unParseTable(shep21_mycn_conserved_enhancer_gff_1kb,shep21_mycn_conserved_enhancer_gff_1kb_file,'\t')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~MAPPING SIGNAL AT SHEP21 MYCN SITES FOR METAS~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def map_shep21_for_heatmap(shep21_chiprx_dataFile,shep21_dataFile):

    '''
    map to both chiprx and regular chip bams
    to make quantification easier, all bams read lengths extended to 200
    chiprx is 150bp read length and chip is 75
    '''

    dataDict = pipeline_dfci.loadDataTable(shep21_chiprx_dataFile)

    #gff files
    shep21_mycn_conserved_gff_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-0_+0.gff' % (gffFolder)
    shep21_mycn_conserved_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_-5kb_+5kb.gff' % (gffFolder)

    shep21_mycn_conserved_promoter_gff_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-0_+0.gff' % (gffFolder)
    shep21_mycn_conserved_promoter_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_PROMOTER_-5kb_+5kb.gff' % (gffFolder)

    shep21_mycn_conserved_enhancer_gff_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-0_+0.gff' % (gffFolder)
    shep21_mycn_conserved_enhancer_gff_5kb_file = '%sHG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb.gff' % (gffFolder)

    #setting the list of gff's to plot
    gffList = [shep21_mycn_conserved_gff_5kb_file,shep21_mycn_conserved_promoter_gff_5kb_file,shep21_mycn_conserved_enhancer_gff_5kb_file]
    cellTypeList = ['SHEP21']
    mapList = [] # map everything

    #for the spike in without rpm
    mappedFolder_noRPM = utils.formatFolder('%smappedFolder_noRPM' % (projectFolder),True)
    pipeline_dfci.mapBams(shep21_chiprx_dataFile,cellTypeList,gffList,mappedFolder_noRPM,nBin = 200,overWrite =False,rpm=False,nameList = mapList,extension = 50)

    #for the spike in with rpm
    pipeline_dfci.mapBams(shep21_chiprx_dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,rpm=True,nameList = mapList,extension=50)

    #for the non spike in
    #note, this data is 75bp reads
    pipeline_dfci.mapBams(shep21_dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,rpm=True,nameList = mapList,extension=125)

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



#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
