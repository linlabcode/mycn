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


#Main method run script for generating heatmaps of enriched pathways from gsea analysis

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
genePlotFolder = '%sgene_plot/' % (projectFolder)
#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory
genomeDirectory = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder,maskFolder,genePlotFolder]

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
    for dataFile in chip_data_list:
        pipeline_dfci.summary(dataFile)


    print('\n\n')
    print('#======================================================================')
    print('#========================II. MAKING FIGURE GFF=========================')
    print('#======================================================================')
    print('\n\n')


    nb_figure_gff_path = make_nb_gff()
    
    #make the associated beds for plottings
    
    nb_mycn_conserved_gff = '%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder)
    canon_path,non_path = makeEboxBeds(nb_mycn_conserved_gff,name='')
    
    bed_string = ','.join([canon_path,non_path])
    print(bed_string)


    print('\n\n')
    print('#======================================================================')
    print('#=====================III. CALLING PLOTTING FUNCTIONS==================')
    print('#======================================================================')
    print('\n\n')


    #for the shep21 no spike system
    plot_shep21_genes(nb_figure_gff_path,bed_string)

    #for the shep21 chiprx system
    scale_path = '%sHG19_SHEP21_CHIPRX_SCALE_FACTORS.txt' % (tableFolder)
    plot_shep21_chiprx_genes(shep21_chiprx_dataFile,scale_path,nb_figure_gff_path,bed_string)

    #for the shep on system
    plot_shep_on_genes(shep_on_dataFile,nb_figure_gff_path,bed_string)

    #for the pan NB metas
    plot_nb_all_genes(nb_all_chip_dataFile,nb_figure_gff_path,bed_string)
    
    #for be2c only
    plot_be2c_genes(be2c_dataFile,nb_figure_gff_path,bed_string)

    #for atac
    pipeline_dfci.summary(atac_dataFile)
    plot_nb_atac_genes(atac_dataFile,nb_figure_gff_path,bed_string)

    #for p493-6
    pipeline_dfci.summary(p4936_young_dataFile)
    plot_p4936_genes(p4936_young_dataFile,nb_figure_gff_path,bed_string)


    #for mm1s
    pipeline_dfci.summary(mm1s_dataFile)
    plot_mm_genes(mm1s_dataFile,nb_figure_gff_path,bed_string)

#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAKING EBOX BEDS~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def makeEboxBeds(gff_path,name='',overwrite=False):

    '''
    makes an ebox bed for the corresponding gff
    '''

    if len(name) == 0:
        gff_name = gff_path.split('/')[-1].split('.')[0]
    else:
        gff_name = name

    #set output
    canon_path = '%s%s_CANON_EBOX.bed' % (bedFolder,gff_name)
    non_path = '%s%s_NONCANON_EBOX.bed' % (bedFolder,gff_name)

    #check to see if already done
    if not overwrite:
        if utils.checkOutput(canon_path,0.1,0.1) and utils.checkOutput(non_path,0.1,0.1):
            print('Found bed output at %s and %s' % (canon_path,non_path))
            return canon_path,non_path
        

    #for each region spit out canonical and non canonical eboxes
    canonBed = []
    nonBed = []

    #open up the repeat filtered gff
    region_gff = utils.parseTable(gff_path,'\t')

    ticker = 0
    for line in region_gff:

        if ticker %100 == 0:
            print(ticker)
        ticker+=1
        chrom = line[0]
        start = int(line[3])
        end = int(line[4])

        seq = string.upper(utils.fetchSeq(genomeDirectory,chrom,start,end,True))

        #get the canonical starts
        canonStarts = [start + match.start() for match in re.finditer('CACGTG',seq)]


        #get the non canonical starts
        nonStarts = [start + match.start() for match in re.finditer('CA..TG',seq)]

        #filter out the canonicals
        nonStarts = [x for x in nonStarts if canonStarts.count(x) == 0]

        #now fill out the bed
        if len(canonStarts) > 0:
            for ebox_start in canonStarts:
                canonBed.append([chrom,ebox_start,(ebox_start+6),'.'])

        if len(nonStarts) > 0:
            for ebox_start in nonStarts:
                nonBed.append([chrom,ebox_start,(ebox_start+6),'.'])

    print('FOUND %s CANONICAL EBOXES' % (len(canonBed)))
    print('FOUND %s NON CANONICAL EBOXES' % (len(nonBed)))
    utils.unParseTable(canonBed,canon_path,'\t')
    utils.unParseTable(nonBed,non_path,'\t')

    return canon_path,non_path



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~MAKING THE NB FIGURE GFF~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_nb_gff():

    '''
    making the gff of the nb figure plotting regions
    '''

    nb_figure_gff = [['chr1', 'RPL22', 'RPL22', '6257001', '6262000', '', '-', '', 'RPL22'],
                     ['chr1', 'RPL22_GENE', 'RPL22_GENE', '6241850', '6262300', '', '-', '', 'RPL22_GENE'],
                     ['chr5', 'NPM1_GENE', 'NPM1_GENE', '170810212', '170858084', '', '+', '', 'NPM1_GENE'],
                     ['chr12', 'CDK4_GENE', 'CDK4_GENE', '58137500', '58148000', '', '-', '', 'CDK4_GENE'],
                     ['chr2', 'ID2_ENHANCER', 'ID2_ENHANCER', '8817001', '8822000', '', '+', '', 'ID2_ENHANCER'],
                     ['chr4', 'HAND2_GENE', 'HAND2_GENE', '174436543', '174453771', '', '-', '', 'HAND2_GENE'],
                     ['chr4', 'HAND2_FULL', 'HAND2_FULL', '174426801', '174464794', '', '-', '', 'HAND2_FULL'],
                     ['chr7','TWIST1','TWIST1',19127919,19163227,'','+','','TWIST1'],
                     ['chr5','NSD1','NSD1',176558355,176718710,'','+','','NSD1'],
                     ['chr3','GATA2','GATA2',128218511,128198863,'','-','','GATA2'],
                     ['chr20','SRSF6','SRSF6',42085199,42090725,'','+','','SRSF6'],
                     ['chr9','BRD3','BRD3',136936744,136896333,'','-','','BRD3'],
                     ['chr9','PAX5','PAX5',36966922,37094885,'','-','','PAX5'],
                     ['chr19','TCF3','TCF3',1636079,1680263,'','-','','TCF3'],
                     ['chr11','TH','TH',2146785,2221909,'','-','','TH'],
                     ['chr6','IRF4','IRF4',314004,413545,'','+','','IRF4'],

                  ]
    
    nb_figure_gff_path = '%sHG19_NB_FIGURE_GENES.gff' % (gffFolder)

    utils.unParseTable(nb_figure_gff,nb_figure_gff_path,'\t')
 
    return nb_figure_gff_path

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING FOR SHEP21 NOSPIKE~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_shep21_genes(nb_figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep21 data
    '''

    #we will have a variety of different plot types
    #all nb_meta baseline
    #chiprx_scaled
    #chiprx w/o scaling
    #just shep21 nospike
    #shep on

    #first establish the plot folder
    plotFolder = utils.formatFolder('%sSHEP21_NOSPIKE/' % (genePlotFolder),True)
    plot_prefix = 'HG19_NB_FIGURE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    #for shep21_dataFile
    dataDict = pipeline_dfci.loadDataTable(shep21_dataFile)
    names_list = dataDict.keys()
    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (shep21_dataFile,bam_extension))

    #for shep21 we want meta of k27ac, pol2, mycn, and twist
    #individual of k27ac, pol2, mycn, and twist

    
    #first do individuals
    for plot_group in ['MYCN','TWIST','H3K27AC','POL2']:
        plotList = [name for name in dataDict.keys() if name.count(plot_group) > 0]
        plotName = '%s_SHEP21_%s_NOSPIKE' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(shep21_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    #now as metas


    plotList = ['SHEP21_0HR_MYCN_NOSPIKE',
                'SHEP21_2HR_MYCN_NOSPIKE',
                'SHEP21_24HR_MYCN_NOSPIKE',
                'SHEP21_0HR_H3K27AC_NOSPIKE',
                'SHEP21_2HR_H3K27AC_NOSPIKE',
                'SHEP21_24HR_H3K27AC_NOSPIKE',
                'SHEP21_0HR_TWIST',
                'SHEP21_2HR_TWIST',
                'SHEP21_24HR_B_TWIST',
                'SHEP21_0HR_POL2_NOSPIKE_R2',
                'SHEP21_2HR_POL2_NOSPIKE',
                'SHEP21_24HR_POL2_NOSPIKE',
                ]
    groupString = 'MYCN,MYCN,MYCN,H3K27AC,H3K27AC,H3K27AC,TWIST,TWIST,TWIST,POL2,POL2,POL2'

    plotName = '%s_SHEP21_NOSPIKE_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(shep21_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_SHEP21_NOSPIKE_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(shep21_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING FOR SHEP21 NOSPIKE~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_shep21_chiprx_genes(shep21_chiprx_dataFile,scale_path,nb_figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep21 chiprx data
    with both spikey normey and without spikey normey
    '''


    #we want a multiplicative scale factor for the data and to not have rpm on

    scale_table = utils.parseTable(scale_path,'\t')
    scale_dict = {}
    for line in scale_table[1:]:
        scale_dict[line[0]] = line[2]



    #first establish the plot folder
    plotFolder_scaled = utils.formatFolder('%sSHEP21_CHIPRX_SCALED/' % (genePlotFolder),True)
    plotFolder_rpm = utils.formatFolder('%sSHEP21_CHIPRX_RPM_NOSCALE/' % (genePlotFolder),True)
    plotFolder_raw = utils.formatFolder('%sSHEP21_CHIPRX_RAW_NOSCALE/' % (genePlotFolder),True)
    plot_prefix = 'HG19_NB_FIGURE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    #for shep21_dataFile
    dataDict = pipeline_dfci.loadDataTable(shep21_chiprx_dataFile)

    names_list = dataDict.keys()
    #initial check for consistency of read lengths
    # for name in names_list:
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,shep21_chiprx_dataFile,bam_extension))

    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (shep21_chiprx_dataFile,bam_extension))

    #for shep21 we want meta of k27ac, pol2, mycn, and twist
    #individual of k27ac, pol2, mycn, and twist

    
    #first do individuals rpm scaled
    for plot_group in ['MYCN','H3K4ME3','H3K27AC','POL2','CTCF']:
        plotList = [name for name in dataDict.keys() if name.count(plot_group) > 0]
        scaleList = [round(1/float(scale_dict[name]),4) for name in plotList]
        scaleList = [str(x) for x in scaleList]
        plot_scale_string =','.join(scaleList)

        #first raw no scaling
        plotName = '%s_SHEP21_%s_RX_RAW_NOSCALE' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(shep21_chiprx_dataFile,nb_figure_gff_path,plotName,plotFolder_raw,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=False,rxGenome = '')


        #first rpm no scaling
        plotName = '%s_SHEP21_%s_RX_RPM_NOSCALE' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(shep21_chiprx_dataFile,nb_figure_gff_path,plotName,plotFolder_rpm,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')


        #next w/ scaling
        plotName = '%s_SHEP21_%s_RX_SCALED' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(shep21_chiprx_dataFile,nb_figure_gff_path,plotName,plotFolder_scaled,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=False,rxGenome = '',scaleFactorString = plot_scale_string)
    #now as metas


    plotList = ['SHEP21_0HR_MYCN_NOSPIKE',
                'SHEP21_2HR_MYCN_NOSPIKE',
                'SHEP21_24HR_MYCN_NOSPIKE',
                'SHEP21_0HR_H3K27AC_NOSPIKE',
                'SHEP21_2HR_H3K27AC_NOSPIKE',
                'SHEP21_24HR_H3K27AC_NOSPIKE',
                'SHEP21_0HR_TWIST',
                'SHEP21_2HR_TWIST',
                'SHEP21_24HR_B_TWIST',
                'SHEP21_0HR_POL2_NOSPIKE_R2',
                'SHEP21_2HR_POL2_NOSPIKE',
                'SHEP21_24HR_POL2_NOSPIKE',
                ]
    groupString = 'MYCN,MYCN,MYCN,H3K27AC,H3K27AC,H3K27AC,TWIST,TWIST,TWIST,POL2,POL2,POL2'

    plotName = '%s_SHEP21_NOSPIKE_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(shep21_dataFile,nb_figure_gff_path,plotName,plotFolder_rpm,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_SHEP21_NOSPIKE_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(shep21_dataFile,nb_figure_gff_path,plotName,plotFolder_rpm,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING FOR SHEP ON SYSTEM~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_shep_on_genes(shep_on_dataFile,nb_figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sSHEP_ON/' % (genePlotFolder),True)
    plot_prefix = 'HG19_NB_FIGURE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(shep_on_dataFile)
    names_list = dataDict.keys()


    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (shep_on_dataFile,bam_extension))
    
    #first do individuals
    for plot_group in ['MYCN','H3K27AC']:
        plotList = [name for name in dataDict.keys() if name.count(plot_group) > 0]
        plotName = '%s_SHEP_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(shep_on_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    #now as metas
    plotList = ['SHEP_0HR_MYCN',
                'SHEP_2HR_MYCN',
                'SHEP_6HR_MYCN',
                'SHEP_0HR_H3K27AC',
                'SHEP_2HR_H3K27AC',
                'SHEP_6HR_H3K27AC',
                ]
    groupString = 'MYCN,MYCN,MYCN,H3K27AC,H3K27AC,H3K27AC'

    plotName = '%s_SHEP_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(shep_on_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_SHEP_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(shep_on_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING METAS ACROSS NB~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_nb_all_genes(nb_all_chip_dataFile,nb_figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sNB_ALL/' % (genePlotFolder),True)
    plot_prefix = 'HG19_NB_FIGURE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(nb_all_chip_dataFile)
    names_list = dataDict.keys()
    #initial check for consistency of read lengths
    # for name in names_list:
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,nb_all_chip_dataFile,bam_extension))

    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (nb_all_chip_dataFile,bam_extension))
    
    #first do individuals
    for plot_group in ['MYCN','H3K27AC']:
        plotList = [name for name in dataDict.keys() if name.count(plot_group) > 0]
        plotName = '%s_NB_ALL_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(nb_all_chip_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    #now as metas
    plotList = ['BE2C_MYCN',
                'KELLY_MYCN',
                'NGP_MYCN',
                'SHEP21_0HR_MYCN_NOSPIKE',
                'BE2C_H3K27AC',
                'KELLY_H3K27AC',
                'NGP_H3K27AC',
                'SHEP21_0HR_H3K27AC_NOSPIKE',
                ]
    groupString = 'MYCN,MYCN,MYCN,MYCN,H3K27AC,H3K27AC,H3K27AC,H3K27AC'

    plotName = '%s_NB_ALL_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(nb_all_chip_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_NB_ALL_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(nb_all_chip_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING BE2C ACROSS NB~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_be2c_genes(be2c_dataFile,nb_figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sBE2C/' % (genePlotFolder),True)
    plot_prefix = 'HG19_NB_FIGURE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(be2c_dataFile)
    names_list = dataDict.keys()
    print(names_list)

    # #initial check for consistency of read lengths
    # for name in names_list:
    #     print(name)
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,be2c_dataFile,bam_extension))


    print(dataDict[names_list[0]]['bam'])
    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (be2c_dataFile,bam_extension))
    
    #first do individuals except for twist using relative scaling

    plotList = [name for name in dataDict.keys() if name.count('TWIST') == 0 and name.count('INPUT') == 0]
    plotName = '%s_BE2C_RELATIVE' % (plot_prefix)
    print(plotName)
    pipeline_dfci.callBatchPlot(be2c_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    plotName = '%s_BE2C_UNIFORM' % (plot_prefix)
    print(plotName)
    pipeline_dfci.callBatchPlot(be2c_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    #now for twist
    plotList = ['BE2C_TWIST']
    twist_extension = 125
    plotName = '%s_BE2C_TWIST' % (plot_prefix)
    print(plotName)
    pipeline_dfci.callBatchPlot(be2c_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MULTIPLE',extension=twist_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING ATAC ACROSS NB~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_nb_atac_genes(atac_dataFile,nb_figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sNB_ATAC/' % (genePlotFolder),True)
    plot_prefix = 'HG19_NB_FIGURE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(atac_dataFile)
    names_list = dataDict.keys()
    print(names_list)

    #initial check for consistency of read lengths
    # for name in names_list:
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,atac_dataFile,bam_extension))



    print(dataDict[names_list[1]]['bam'])
    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]

    
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (atac_dataFile,bam_extension))
    
    #first do individuals
    for plot_group in ['ATAC']:
        plotList = [name for name in dataDict.keys() if name.count(plot_group) > 0 and name.count('MM1S') == 0]
        plotName = '%s_NB_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(atac_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    

    #now as metas
    plotList = ['BE2C_ATAC_rep1',
                'KELLY_ATAC',
                'NGP_ATAC',
                'SHEP21_ATAC',
                ]
    groupString = 'ATAC,ATAC,ATAC,ATAC'

    plotName = '%s_NB_ATAC_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(atac_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_NB_ATAC_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(atac_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING FOR P439-6 SYSTEM~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_p4936_genes(p4936_young_dataFile,nb_figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sP4936/' % (genePlotFolder),True)
    plot_prefix = 'HG19_NB_FIGURE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(p4936_young_dataFile)
    names_list = dataDict.keys()


    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (p4936_young_dataFile,bam_extension))
    
    #first do individuals
    for plot_group in ['MYC','H3K27AC','RNA_POL2','MAX']:
        plotList = [name for name in dataDict.keys() if name.count(plot_group) > 0]
        plotName = '%s_P4936_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(p4936_young_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING FOR P439-6 SYSTEM~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_mm_genes(mm1s_dataFile,nb_figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sMM1S/' % (genePlotFolder),True)
    plot_prefix = 'HG19_NB_FIGURE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(mm1s_dataFile)
    names_list = dataDict.keys()


    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (mm1s_dataFile,bam_extension))
    
    #first do individuals
    for plot_group in ['MYC','H3K27AC']:
        plotList = [name for name in dataDict.keys() if name.count(plot_group) > 0]
        plotName = '%s_MM1S_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(mm1s_dataFile,nb_figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')




    



#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
