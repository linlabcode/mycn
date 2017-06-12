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


#Main method run script for mycn code

#See README for additional information on downloading and installing dependencies

#Requires linlab pipeline set of utilities 
#

#Requires bamliquidator
#

#Requires 

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

#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory
genomeDirectory = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

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
u87_dataFile = '%sdata_tables/U87_TABLE.txt' % (projectFolder)
nb_all_chip_dataFile = '%sdata_tables/NB_ALL.txt' % (projectFolder)
p4936_young_dataFile = '%sdata_tables/P493-6_YOUNG_TABLE.txt' % (projectFolder)
sclc_dataFile = '%sdata_tables/SCLC_DATA_TABLE.txt' % (projectFolder)
shep21_dataFile = '%sdata_tables/SHEP21_TABLE.txt' % (projectFolder)
shep_on_dataFile = '%sdata_tables/SHEP_ON_TABLE.txt' % (projectFolder)

chip_data_list = [be2c_dataFile,mm1s_dataFile,nb_all_chip_dataFile,p4936_young_dataFile,sclc_dataFile,shep21_dataFile,shep_on_dataFile,u87_dataFile]


#note: all mouse analysis of THMYCN tumors are in a separate script

#CHIP-RX
shep21_chiprx_dataFile = '%sdata_tables/SHEP21_CHIPRX_TABLE.txt' % (projectFolder)

#RNA-Seq
be2c_rna_drug_dataFile = '%sdata_tables/BE2C_RNA_DRUG_TABLE.txt' % (projectFolder)
be2c_rna_twist_dataFile = '%sdata_tables/BE2C_RNA_TWIST_TABLE.txt' % (projectFolder)
shep21_rna_dataFile = '%sdata_tables/SHEP21_DOX_RNA_TABLE.txt' % (projectFolder)

all_data_list = [atac_dataFile,be2c_dataFile,mm1s_dataFile,nb_all_chip_dataFile,p4936_young_dataFile,sclc_dataFile,shep21_dataFile,shep_on_dataFile,u87_dataFile,shep21_chiprx_dataFile,be2c_rna_drug_dataFile,be2c_rna_twist_dataFile,shep21_rna_dataFile]

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
    #edit all of the data files to absolute path the
    for dataFile in chip_data_list:

        pipeline_dfci.summary(dataFile)


    print('\n\n')
    print('#======================================================================')
    print('#==========================II. CALLING MACS============================')
    print('#======================================================================')
    print('\n\n')

    #running peak finding using macs 1.4.2 on all chip datasets
    #this usually takes ~2-3 hours on a reasonably fast machine
    #a 3 hour time out on this entire operation is set
    #if peak calling takes longer than 3 hours, simply run the script again after completion


    # for dataFile in chip_data_list:

    #     run_macs(dataFile)



    print('\n\n')
    print('#======================================================================')
    print('#===================III. DEFINING ACTIVE GENES IN NB===================')
    print('#======================================================================')
    print('\n\n')

    
    # #here we will identify active promoters in various contexts as those with 
    # #an H3K27AC peak in the +/- 1kb tss region
    # #UCSC refseq annotations are used for all genes
    # #make_nb_active_gene_lists(nb_all_chip_dataFile)
    
    # make_active_gene_lists(mm1s_dataFile,p4936_young_dataFile,sclc_dataFile,shep_on_dataFile,u87_dataFile)

    print('\n\n')
    print('#======================================================================')
    print('#===============IV. DEFINING NB MYCN AND H3K27AC LANDSCAPE=============')
    print('#======================================================================')
    print('\n\n')

    # #for enhancers
    # enhancer_bashFileName,enhancer_region_map_path,namesList = define_enhancer_landscape(projectFolder,pipeline_dir,nb_all_chip_dataFile)

    # #runs only if no output detected
    # if not utils.checkOutput(enhancer_region_map_path,0,0):
    #     print(enhancer_bashFileName)
    #     os.system('bash %s' % (enhancer_bashFileName))


    # #for mycn
    # mycn_bashFileName,mycn_region_map_path,namesList = define_mycn_landscape(projectFolder,pipeline_dir,nb_all_chip_dataFile)

    # if not utils.checkOutput(mycn_region_map_path,0,0):
    #     print(mycn_bashFileName)
    #     os.system('bash %s' % (mycn_bashFileName))
    
    # #now we need to call the R script that creates the rank plots
    # if utils.checkOutput(mycn_region_map_path,1,30): #set a wait time for 30 minutes
    #     print('Found NB_MYCN meta_rose landscape and running rank plot R code')

    #     conserved_rank_path = '%smeta_rose/NB_MYCN/NB_MYCN_0KB_STITCHED_ENHANCER_REGION_RANK_CONSERVED.txt' % (projectFolder)
    #     if utils.checkOutput(conserved_rank_path,0,0):
    #         print('Identified NB rank conserved regions: %s' % (conserved_rank_path))
    #     else:
    #         print('Defining NB rank conserved regions')
    #         name_string = ','.join(namesList) #provides the dataset names used
    #         rank_script_path = '%sr_scripts/1_nb_mycn_rank.R' % (projectFolder)
    #         r_cmd = 'Rscript %s %s %s %s' % (rank_script_path,mycn_region_map_path,name_string,projectFolder)
    #         print(r_cmd)
    #         os.system(r_cmd)


    print('\n\n')
    print('#======================================================================')
    print('#==========V. MAPPING MYCN AND H3K27AC TO MYCN REGIONS=================')
    print('#======================================================================')
    print('\n\n')

    # #here we will first make a gff of conserved NB MYCN regions
    # #and then map MYCN and H3K27ac signal 

    # print('Making a gff and bed of conserved NB MYCN regions:')

    # mycn_gff_path,mycn_flank_gff_path = make_mycn_regions(conserved_rank_path) 
    
    # print('Mapping MYCN and H3K27AC signal')
    # gffList = [mycn_gff_path,mycn_flank_gff_path]
    #gffList = ['%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder),'%sHG19_NB_MYCN_CONSERVED_-500_+500.gff' % (gffFolder)]
    #pipeline_dfci.map_regions(nb_all_chip_dataFile,gffList,mappedFolder,signalFolder)

    print('\n\n')
    print('#======================================================================')
    print('#==================VI. CREATING NB MYCN STATS TABLE====================')
    print('#======================================================================')
    print('\n\n')
    

    # mycn_table_path = '%stables/HG19_NB_MYCN_CONSERVED_STATS_TABLE.txt' % (projectFolder)
    # if utils.checkOutput(mycn_table_path,0,0):
    #     print('Identified MYCN table %s' % (mycn_table_path))
    # else:
    #     print('Making MYCN stats table')        
    #     mycn_table_path = make_mycn_stats_table(nb_all_chip_dataFile,mycn_table_path)

    mycn_table_path = '%stables/HG19_NB_MYCN_CONSERVED_STATS_TABLE.txt' % (projectFolder)
    #mycn_table_path = make_mycn_stats_table(nb_all_chip_dataFile,mycn_table_path)
    print('\n\n')
    print('#======================================================================')
    print('#=================VII. MAKING VECTOR COMPARISON PLOTS==================')
    print('#======================================================================')
    print('\n\n')

    compare_script_path = '%sr_scripts/2_nb_mycn_vector_plots.R' % (projectFolder)
    r_cmd = 'Rscript %s %s %s' % (compare_script_path,mycn_table_path,projectFolder)
    print(r_cmd)
    os.system(r_cmd)


    print('\n\n')
    print('#======================================================================')
    print('#==================VIII. RANKING EBOXES IN MYCN PEAKS==================')
    print('#======================================================================')
    print('\n\n')

    # mycn_gff_path = '%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder)
    # ebox_rank_path = rank_eboxes(nb_all_chip_dataFile,mycn_gff_path,macsFolder,genomeDirectory,window = 100)
    
    # print(ebox_rank_path)

    # #now make the heatmap
    # ebox_heatmap_script_path = '%sr_scripts/3_nb_ebox_heatmap.R' % (projectFolder)
    # r_cmd = 'Rscript %s %s %s' % (ebox_heatmap_script_path,ebox_rank_path,projectFolder)
    # print(r_cmd)
    # os.system(r_cmd)

    print('\n\n')
    print('#======================================================================')
    print('#====================IX. MAPPING BE2C DATASETS TO TSS==================')
    print('#======================================================================')
    print('\n\n')

    # gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    # be2c_signal_path = pipeline_dfci.map_regions(be2c_dataFile,gffList,mappedFolder,signalFolder,[],False)
    
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
#~~~~~~~~~~~~~~~~~~~~~~~DEFINING NB ACTIVE GENES~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_nb_active_gene_lists(nb_all_chip_dataFile):

    pipeline_dfci.makeGeneGFFs(annotFile,gffFolder,species=genome.upper())


    dataDict = pipeline_dfci.loadDataTable(nb_all_chip_dataFile)
    setName = 'NB_TSS_H3K27AC'
    gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    cellTypeList = ['BE2C','KELLY','NGP','SHEP21']
    namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]

    pipeline_dfci.mapEnrichedToGFF(nb_all_chip_dataFile,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)

    #this is for the union
    mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_NB_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['BE2C_H3K27AC'],['KELLY_H3K27AC'],['NGP_H3K27AC'],['SHEP21_0HR_H3K27AC_NOSPIKE']]
    output = '%sgeneListFolder/HG19_NB_H3K27AC_ACTIVE_UNION.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)

    #this is for individual NB datasets
    namesList =['BE2C_H3K27AC','KELLY_H3K27AC','NGP_H3K27AC','SHEP21_0HR_H3K27AC_NOSPIKE']
    for name in namesList:
        mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_NB_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
        setList = [[name]]
        output = '%sgeneListFolder/HG19_%s_ACTIVE.txt' % (projectFolder,name)
        pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~DEFINING ACTIVE GENES IN OTHER LINES~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#now we need to make active gene lists in MM1S, P493-6, SCLC, and the SHEP ON system
def make_active_gene_lists(mm1s_dataFile,p4936_young_dataFile,sclc_dataFile,shep_on_dataFile,u87_dataFile):
    #first map to enriched for each dataset

    #for mm1s
    dataDict = pipeline_dfci.loadDataTable(mm1s_dataFile)
    setName = 'MM1S_TSS_H3K27AC'
    gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    cellTypeList = ['MM1S']
    namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]
    pipeline_dfci.mapEnrichedToGFF(mm1s_dataFile,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)

    #for p4936
    dataDict = pipeline_dfci.loadDataTable(p4936_young_dataFile)
    setName = 'P493-6_TSS_H3K27AC'
    gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    cellTypeList = ['P493-6']
    namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]
    pipeline_dfci.mapEnrichedToGFF(p4936_young_dataFile,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)

    #for sclc
    dataDict = pipeline_dfci.loadDataTable(sclc_dataFile)
    setName = 'SCLC_TSS_H3K27AC'
    gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    cellTypeList = ['H128','H2171']
    namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]
    pipeline_dfci.mapEnrichedToGFF(sclc_dataFile,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)

    #for shep on
    dataDict = pipeline_dfci.loadDataTable(shep_on_dataFile)
    setName = 'SHEP_ON_TSS_H3K27AC'
    gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    cellTypeList = ['SHEP']
    namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]
    pipeline_dfci.mapEnrichedToGFF(shep_on_dataFile,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)


    #for u87
    dataDict = pipeline_dfci.loadDataTable(u87_dataFile)
    setName = 'U87_TSS_H3K27AC'
    gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    cellTypeList = ['U87']
    namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]
    pipeline_dfci.mapEnrichedToGFF(u87_dataFile,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)


    #====================
    #now create the gene lists
    #this is for MM1S
    mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_MM1S_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['MM1S_H3K27AC_DMSO']]
    output = '%sgeneListFolder/HG19_MM1S_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)


    #this is for P493-6
    #here we will take the union of all datasets
    mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_P493-6_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['P493-6_T0_H3K27AC'],['P493-6_T1_H3K27AC'],['P493-6_T24_H3K27AC']]
    output = '%sgeneListFolder/HG19_P493-6_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)


    #this is for P493-6
    #here we will take the union of all datasets
    mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_P493-6_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['P493-6_T24_H3K27AC']]
    output = '%sgeneListFolder/HG19_P493-6_T24_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)



    #this is for SCLC
    #here we will take the union of all datasets
    mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_SCLC_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['H128_H3K27AC'],['H2171_H3K27AC']]
    output = '%sgeneListFolder/HG19_SCLC_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)


    #this is for SCLC
    #here we will take the union of all datasets
    mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_SCLC_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['H2171_H3K27AC']]
    output = '%sgeneListFolder/HG19_H2171_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)


    #this is for SHEP ON
    #here we will take the union of all datasets
    mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_SHEP_ON_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['SHEP_0HR_H3K27AC'],['SHEP_2HR_H3K27AC'],['SHEP_6HR_H3K27AC']]
    output = '%sgeneListFolder/HG19_SHEP_ON_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)

    #this is for U87
    #here we will take the union of all datasets
    mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_U87_TSS_H3K27AC.txt' % (mappedEnrichedFolder)
    #this setList variable defines overlap logic for promoters. In this case, it's asking for the union of all datasets
    setList = [['U87_H3K27AC']]
    output = '%sgeneListFolder/HG19_U87_H3K27AC_ACTIVE.txt' % (projectFolder)
    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~DEFINING NB H3K27AC ENHANCER LANDSCAPE~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_enhancer_landscape(projectFolder,pipeline_dir,nb_all_chip_dataFile):

    '''
    defines the NB enhancer baseline using H3K27ac chips from NGP, KELLY, BE2C, and SHEP21
    enhancers defined using auto optimized stitching of nearby regions
    w/ a 2.5kb tss exclusion
    uses the meta rose code and writes out a .sh file for reproducibility
    '''

    #For H3K27AC
    #with TSS exclusion and auto stitching

    dataDict = pipeline_dfci.loadDataTable(nb_all_chip_dataFile)
    analysisName = 'NB_H3K27AC'
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


    #the 4KB parameter is 
    region_map_path = '%s%s/%s_AllEnhancers.table.txt' % (roseFolder,analysisName,analysisName)
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
#~~~~~~~~~~~~~~~~MAKING NB CONSERVED MYCN REGION GFFS AND BEDS~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_mycn_regions(conserved_rank_path):

    '''
    takes conserved NB MYCN regions 
    then creates a bed and gff of regions
    '''

    conserved_rank_table = utils.parseTable(conserved_rank_path,'\t')
    mycn_gff = []
    mycn_flank_gff = []
    mycn_bed = []
    mycn_flank_bed = []

    for line in conserved_rank_table[1:]:
        locus_line = utils.Locus(line[1],line[2],line[3],'.')
        
        if int(line[3]) < int(line[2]):
            print('uh oh')
            print(line)
        gff_line = [line[1],line[0],'',line[2],line[3],'','.','',line[0]]
        bed_line = [line[1],line[2],line[3],line[0]]
        mycn_gff.append(gff_line)
        mycn_bed.append(bed_line)

        gff_flank_line = [line[1],line[0],'',int(line[2])-500,int(line[3])+500,'','.','',line[0]]
        bed_flank_line = [line[1],int(line[2])-500,int(line[3])+500,line[0]]
        mycn_flank_gff.append(gff_flank_line)
        mycn_flank_bed.append(bed_flank_line)
        
    mycn_gff_path = '%sHG19_NB_MYCN_CONSERVED_-0_+0.gff' % (gffFolder)
    mycn_flank_gff_path = '%sHG19_NB_MYCN_CONSERVED_-500_+500.gff' % (gffFolder)

    mycn_bed_path = '%sHG19_NB_MYCN_CONSERVED_-0_+0.bed' % (bedFolder)
    mycn_flank_bed_path = '%sHG19_NB_MYCN_CONSERVED_-500_+500.bed' % (bedFolder)

    #writing to disk
    utils.unParseTable(mycn_gff,mycn_gff_path,'\t')
    utils.unParseTable(mycn_flank_gff,mycn_flank_gff_path,'\t')

    utils.unParseTable(mycn_bed,mycn_bed_path,'\t')
    utils.unParseTable(mycn_flank_bed,mycn_flank_bed_path,'\t')

    print(mycn_gff_path)
    print(mycn_flank_gff_path)
    print(mycn_bed_path)
    print(mycn_flank_bed_path)
    return mycn_gff_path,mycn_flank_gff_path





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAKING NB MYCN STATS TABLE~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def nmers(n,baseList=[],nmerList = []):
    #print(baseList)
    if len(baseList) == 0:
        baseList =['A','T','G','C']
    if n == 0:
        return nmerList
    tempList = []
    if len(nmerList) == 0 and n == 1:
        return baseList
    if len(nmerList) == 0:
        nmerList = list(baseList)
        n = n -1
    for base in baseList:
        tempList = tempList + map(lambda x:x+base,nmerList)
    return nmers(n-1,baseList,tempList)


def make_mycn_stats_table(nb_all_chip_dataFile,outFile):

    '''
    making a table of conserved mycn peaks w/ some additional stats
    mycn and h3k27ac signal is avg. background normalized across 4 samples
    active tss defined as the union of all H3K27ac occupied promoters in NB
    active enhancers defined as the union of all H3K27ac sites outside of promoters
    '''
    dataDict = pipeline_dfci.loadDataTable(nb_all_chip_dataFile)

    print('SETTING UP OUTPUT TABLE')
    outTable = [['PEAK_ID','CHROM','START','STOP','LENGTH','ACTIVE_TSS_OVERLAP','ENHANCER_OVERLAP','CPG_ISLAND_OVERLAP','CPG_ISLAND_FRACTION','GC_FREQ','MYCN_RANK','AVG_MYCN_SIGNAL','AVG_H3K27AC_SIGNAL','CANON_EBOX_COUNT','NONCANON_EBOX_COUNT','TOTAL_EBOX_COUNT','CANON_EXP','NON_CANON_EXP','GABPA_COUNT','GABPA_EXP','GATA_COUNT','GATA_EXP']]

    dinuc = nmers(2,['A','T','G','C'])

    #input files
    mycnSignalFile = '%sHG19_NB_MYCN_CONSERVED_-0_+0_NB_ALL_SIGNAL.txt' % (signalFolder)
    h3k27acSignalFile = '%sHG19_NB_MYCN_CONSERVED_-500_+500_NB_ALL_SIGNAL.txt' % (signalFolder)
    mycnRankFile = '%smeta_rose/NB_MYCN/NB_MYCN_0KB_STITCHED_ENHANCER_REGION_RANK_CONSERVED.txt' % (projectFolder)
    activeGeneFile = '%sHG19_NB_H3K27AC_ACTIVE_UNION.txt' % (geneListFolder)
    #note, this is the ucsc hg19 cpg islands extended file
    #to download and format run ./beds/download_cpg.sh
    cpgFile = '%sbeds/hg19_cpg_islands.bed' % (projectFolder)
    enhancerFile = '%smeta_rose/NB_H3K27AC/NB_H3K27AC_AllEnhancers.table.txt' % (projectFolder)

    print('LOADING MYCN BINDING DATA')
    mycnSignalTable = utils.parseTable(mycnSignalFile,'\t')

    #making a signal dictionary for MYCN binding
    names_list = ['BE2C_MYCN','KELLY_MYCN','NGP_MYCN','SHEP21_0HR_MYCN_NOSPIKE']
    background_list = [dataDict[name]['background'] for name in names_list]
    header = mycnSignalTable[0]
    chip_columns = [header.index(name) for name in names_list]
    background_columns = [header.index(background_name) for background_name in background_list]
    
    mycn_sig_dict = {}
    #this only works if the first column are unique identifiers
    if len(mycnSignalTable) != len(utils.uniquify([line[0] for line in mycnSignalTable])):
        print('Error: Column 1 of must contain unique identifiers.' % (mycnSignalFile))
        sys.exit()
    for line in mycnSignalTable[1:]:
        line_sig = []
        for i in range(len(names_list)):
            line_sig.append(float(line[chip_columns[i]]) - float(line[background_columns[i]]))
        mycn_sig_dict[line[0]] = numpy.mean(line_sig)


    
    print('LOADING MYCN RANK DATA')
    mycnRankTable = utils.parseTable(mycnRankFile,'\t')

    print('LOADING H3K27AC BINDING DATA')
    h3k27acSignalTable = utils.parseTable(h3k27acSignalFile,'\t')
    #making a signal dictionary for background subtracted H3K27ac binding
    names_list = ['BE2C_H3K27AC','KELLY_H3K27AC','NGP_H3K27AC','SHEP21_0HR_H3K27AC_NOSPIKE']
    background_list = [dataDict[name]['background'] for name in names_list]
    header = h3k27acSignalTable[0]
    chip_columns = [header.index(name) for name in names_list]
    background_columns = [header.index(background_name) for background_name in background_list]
    
    h3k27ac_sig_dict = {}
    #this only works if the first column are unique identifiers
    if len(h3k27acSignalTable) != len(utils.uniquify([line[0] for line in h3k27acSignalTable])):
        print('Error: Column 1 of must contain unique identifiers.' % (h3k27acSignalFile))
        sys.exit()
    for line in h3k27acSignalTable[1:]:
        line_sig = []
        for i in range(len(names_list)):
            line_sig.append(float(line[chip_columns[i]]) - float(line[background_columns[i]]))
        h3k27ac_sig_dict[line[0]] = numpy.mean(line_sig)



    #making the cpg collection
    print('LOADING CPGS ISLANDS')
    cpgBed = utils.parseTable(cpgFile,'\t')
    cpgLoci = []
    for line in cpgBed:
        cpgLoci.append(utils.Locus(line[0],line[1],line[2],'.',line[-1]))
    cpgCollection = utils.LocusCollection(cpgLoci,50)
        
    #next make the tss collection of active promoters
    print('LOADING ACTIVE PROMOTERS')
    startDict = utils.makeStartDict(annotFile)
    activeTable = utils.parseTable(activeGeneFile,'\t')
    tss_1kb_loci = []
    for line in activeTable:
        tss_1kb_loci.append(utils.makeTSSLocus(line[1],startDict,1000,1000))
    tss_1kb_collection = utils.LocusCollection(tss_1kb_loci,50)


    #enhancer file
    print("LOADING ACTIVE ENHANCERS")
    enhancerTable = utils.parseTable(enhancerFile,'\t')
    print('STARTING WITH THE FOLLOWING NUMBER OF ENHANCERS IN NB')
    print(len(enhancerTable) - 6)
    enhancerLoci = []
    for line in enhancerTable:
        if line[0][0] != '#' and line[0][0] != 'R':
            try:
                lineLocus = utils.Locus(line[1],int(line[2]),int(line[3]),'.',line[0])
                enhancerLoci.append(lineLocus)
            except IndexError:
                print(line)
                sys.exit()
    enhancerCollection = utils.LocusCollection(enhancerLoci,50)

    print('CLASSIFYING MYCN PEAKS')
    ticker = 0
    for i in range(1,len(mycnSignalTable)):
        if ticker%100 == 0:
            print(ticker)
        ticker +=1

        line = mycnSignalTable[i]        

        mycn_signal = round(mycn_sig_dict[line[0]],4)
        h3k27ac_signal = round(h3k27ac_sig_dict[line[0]],4)
        
        peakID = line[0]
        locusString = line[1]
        chrom = locusString.split('(')[0]
        [start,stop] = [int(x) for x in line[1].split(':')[-1].split('-')]
        lineLocus = utils.Locus(chrom,start,stop,'.',peakID)
        
        tssOverlap = 0
        if tss_1kb_collection.getOverlap(lineLocus,'both'):
            tssOverlap = 1

        enhancerOverlap = 0
        if enhancerCollection.getOverlap(lineLocus,'both') and tssOverlap == 0:
            enhancerOverlap = 1

        cpgIslandOverlap = 0
        if cpgCollection.getOverlap(lineLocus,'both'):
            cpgIslandOverlap = 1

        #now do fractional cpgOverlap
        overlappingCpGLoci = cpgCollection.getOverlap(lineLocus,'both')
        overlappingBases = 0
        for locus in overlappingCpGLoci:
            cpgStart = max(locus.start(),lineLocus.start())
            cpgEnd = min(locus.end(),lineLocus.end())
            overlappingBases += (cpgEnd-cpgStart)
        overlapFraction = round(float(overlappingBases)/lineLocus.len(),2)
        
        #now get the seq
        lineSeq = string.upper(utils.fetchSeq(genomeDirectory,chrom,start,stop,True))
        gcFreq = round(float(lineSeq.count('GC') + lineSeq.count('CG'))/len(lineSeq),2)
            
        dinuc_dict = {}
        for nmer in dinuc:
            dinuc_dict[nmer] = float(lineSeq.count('GC'))/len(lineSeq)

        
        mycnRankLine = mycnRankTable[i]
        mycnRank = numpy.mean([float(x) for x in mycnRankLine[6:]])

        canonMatchList = re.findall('CACGTG',lineSeq)
        canon_count = len(canonMatchList)

        eboxMatchList = re.findall('CA..TG',lineSeq)
        ebox_count = len(eboxMatchList)

        non_canon_count = ebox_count-canon_count

        #get the expected values
        canon_exp = dinuc_dict['CA']*dinuc_dict['CG']*dinuc_dict['TG']*(len(lineSeq) - 5)
        canon_exp = round(canon_exp,2)
        notCG = 1- dinuc_dict['CG']
        non_exp = dinuc_dict['CA']*notCG*dinuc_dict['TG']*(len(lineSeq) - 5)
        non_exp = round(non_exp,2)



        #for gata and GABPA
        gabpaMatchList = re.findall('CGGAAG',lineSeq) + re.findall('CTTCCG',lineSeq)
        gabpa_count = len(gabpaMatchList)

        gabpa_exp_f = dinuc_dict['CG'] * dinuc_dict['GA'] * dinuc_dict['AG']*(len(lineSeq) - 5)
        gabpa_exp_r = dinuc_dict['CT'] * dinuc_dict['TC'] * dinuc_dict['CG']*(len(lineSeq) - 5)
        
        gabpa_exp = round(gabpa_exp_f,2) + round(gabpa_exp_r,2)

        gataMatchList = re.findall('GATAA',lineSeq) + re.findall('TTATC',lineSeq)
        gata_count = len(gataMatchList)

        an_freq = 1 - dinuc_dict['AA'] - dinuc_dict['AT'] - dinuc_dict['AG'] -dinuc_dict['AC']
        cn_freq = 1 - dinuc_dict['CA'] - dinuc_dict['CT'] - dinuc_dict['CG'] -dinuc_dict['CC']
        gata_exp_f = dinuc_dict['GA'] * dinuc_dict['TA'] * an_freq*(len(lineSeq) - 5)
        gata_exp_r = dinuc_dict['TT'] * dinuc_dict['AT'] * cn_freq*(len(lineSeq) - 5)
        gata_exp = round(gata_exp_f,2) + round(gata_exp_r,2)

        
        

        newLine = [peakID,chrom,start,stop,lineLocus.len(),tssOverlap,enhancerOverlap,cpgIslandOverlap,overlapFraction,gcFreq,mycnRank,mycn_signal,h3k27ac_signal,canon_count,non_canon_count,ebox_count,canon_exp,non_exp,gabpa_count,gabpa_exp,gata_count,gata_exp]
        outTable.append(newLine)

    utils.unParseTable(outTable,outFile,'\t')
    
    return outFile



#for the chiprx code
def filterPeaks(tabixFolder,mycTablePath,outputPath,repeatList = []):

    '''                                                                                             
    auto filters the 3 repeat classes LINE, LTR, Simple_repeat                                      
    outputs a bed in the format of                                                                  
    [PEAK_ID,CHROM, START,STOP,LENGTH, LINE, LTR, Simple_repeat]                                    
    '''

    if len(repeatList) == 0:
        repeatList = ['LINE','LTR','Simple_repeat']

    repeatTable = [['PEAK_ID','CHROM','START','STOP','LENGTH'] + repeatList]

    mycTable = utils.parseTable(mycTablePath,'\t')
    ticker =0
    for line in mycTable[1:]:
        if line[0][0] =='P':
            continue

        if ticker % 100 == 0:
            print ticker
        ticker +=1
        peak_ID = line[0]
        chrom = line[1]
        start = int(line[2])
        stop = int(line[3])
        length = line[4]
        locusString = '%s:%s-%s' % (chrom,start,stop)

        repeatFractions = []
        for repeatClass in repeatList:
            tabixGFF = '%shg19_%s_category_sorted.gff.gz' % (tabixFolder,repeatClass)

            tabixCmd = 'tabix %s %s' % (tabixGFF,locusString)

            tabix = subprocess.Popen(tabixCmd,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)

            tabixLines = tabix.stdout.readlines()

            tabixLines = [x.rstrip().split('\t') for x in tabixLines] #i think you get back essentially gff lines                                                                                      

            overlapFraction = 0.0
            for line in tabixLines:
                lineStart = int(line[3])
                lineStop = int(line[4])
                lineStart = max(start,lineStart)
                lineStop = min(stop,lineStop)
                overlapLength = lineStop - lineStart
                overlapFraction += float(overlapLength)/float(length)
            repeatFractions.append(round(overlapFraction,4))

        newLine = [peak_ID,chrom,start,stop,length] + repeatFractions
        repeatTable.append(newLine)

    utils.unParseTable(repeatTable,outputPath,'\t')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~RANKING MYCN PEAKS BY STRENGTH~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#this is adapted from the lin 2012 code, but uses actual binding strength as
#opposed to enrichment to quantify peaks

def rank_eboxes(nb_all_chip_dataFile,mycn_gff_path,macsFolder,genomeDirectory,window = 100):

    '''
    uses the  conserved MYCN sites and ranks eboxes within them
    by average background subtracted signal
    searches 100bp (window variable)  from mycn summits
    '''
    
    window = int(window)

    #bring in the conserved mycn region
    print('making gff of nb mycn summits')
    nb_mycn_gff = utils.parseTable(mycn_gff_path,'\t')

    nb_mycn_collection = utils.gffToLocusCollection(nb_mycn_gff,50)

    dataDict =pipeline_dfci.loadDataTable(nb_all_chip_dataFile)
    names_list = [name for name in dataDict.keys() if name.count('MYCN') == 1]
    names_list.sort()

    summit_loci = []
    #first makes a gff of all summits +/- 100bp for all nb mycn datasets
    for name in names_list:
        summit_bed_path = '%s%s/%s_summits.bed' % (macsFolder,name,name)
        summit_bed = utils.parseTable(summit_bed_path,'\t')
        for line in summit_bed:
            summit_locus = utils.Locus(line[0],int(line[1])-window,int(line[2])+window,'.',line[3])
            if len(nb_mycn_collection.getOverlap(summit_locus)) > 0:
                summit_loci.append(summit_locus)

    summit_collection =utils.LocusCollection(summit_loci,50)
    summit_merged_collection = summit_collection.stitchCollection()
    
    summit_gff = utils.locusCollectionToGFF(summit_merged_collection)
    summit_gff_path = '%sHG19_NB_MYCN_SUMMITS_-%s_+%s.gff' % (gffFolder,window,window)
    utils.unParseTable(summit_gff,summit_gff_path,'\t')

    #this is borrowed from above and maps chip-seq signal to the gff
    print('mapping to nb mycn summits and making signal dict')
    gffList = [summit_gff_path]
    summit_signal_path = pipeline_dfci.map_regions(nb_all_chip_dataFile,gffList)


    mycnSignalTable = utils.parseTable(summit_signal_path,'\t')

    #making a signal dictionary for MYCN binding
    names_list = ['BE2C_MYCN','KELLY_MYCN','NGP_MYCN','SHEP21_0HR_MYCN_NOSPIKE']
    background_list = [dataDict[name]['background'] for name in names_list]
    header = mycnSignalTable[0]
    chip_columns = [header.index(name) for name in names_list]
    background_columns = [header.index(background_name) for background_name in background_list]
    
    mycn_sig_dict = {}
    for line in mycnSignalTable[1:]:
        line_sig = []
        for i in range(len(names_list)):
            line_sig.append(float(line[chip_columns[i]]) - float(line[background_columns[i]]))
        region_id = line[1]
        coords = [int(x) for x in line[1].split(':')[-1].split('-')]
        line_length = coords[1]-coords[0]
        mycn_sig_dict[region_id] = numpy.mean(line_sig)*line_length

    #now for each region find the eboxes and then add up the signal
    print('making ebox ranking')
    ebox_list = ['CACGTG','CAGTTG','CAAGTG','CAGGTG','CAATTG','CAAATG','CATCTG','CAGCTG','CATGTG','CATATG']
    eboxDict = {}
    for ebox in ebox_list:
        eboxDict[ebox] = []
    ticker = 0
    for line in summit_gff:
        if ticker % 1000 == 0:
            print(ticker)
        ticker+=1

        chrom = line[0]
        sense = '.'

        start = int(line[3])
        end = int(line[4])
        region_id = '%s(%s):%s-%s' % (line[0],line[6],line[3],line[4])
        signal = mycn_sig_dict[region_id]

        sequenceLine = utils.fetchSeq(genomeDirectory,chrom,start,end,True)
        
        motifVector = []
        matches = re.finditer('CA..TG',str.upper(sequenceLine))
        if matches:
            for match in matches:
                motifVector.append(match.group())
        
        #count only 1 of each motif type per line
        #motifVector = utils.uniquify(motifVector)
        for motif in motifVector:
            if ebox_list.count(motif) > 0:
                eboxDict[motif].append(signal)
            else:
                eboxDict[utils.revComp(motif)].append(signal)


    eboxTable =[]
    eboxTableOrdered =[['EBOX','OCCURENCES','AVG_HEIGHT']]
    for ebox in eboxDict.keys():
        newLine = [ebox,len(eboxDict[ebox]),numpy.mean(eboxDict[ebox])]
        eboxTable.append(newLine)


    occurenceOrder = utils.order([line[2] for line in eboxTable],decreasing=True)
    
    for x in occurenceOrder:
        eboxTableOrdered.append(eboxTable[x])
    print(eboxTableOrdered)
    ebox_outfile = '%sHG19_NB_MYCN_CONSERVED_SUMMITS_-%s_+%s_EBOX_RANK.txt' % (tableFolder,window,window)
    utils.unParseTable(eboxTableOrdered,ebox_outfile,'\t')
    return ebox_outfile






#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
