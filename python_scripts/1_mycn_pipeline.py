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

#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory
genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder]

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
    for dataFile in chip_data_list:

        pipeline_dfci.summary(dataFile)

    print('\n\n')
    print('#======================================================================')
    print('#==========================II. CALLING MACS============================')
    print('#======================================================================')
    print('\n\n')


    for dataFile in chip_data_list:
        dataDict = pipeline_dfci.loadDataTable(dataFile)
        namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
        namesList.sort()
        print(namesList)
        #for name in namesList:
        #    print(dataDict[name]['bam'])
        pipeline_dfci.callMacs(dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9')
        os.chdir(projectFolder) # the silly call macs script has to change into the output dir
        #so this takes us back to the project folder


    #sym linking peak beds
    print('sym linking peak files to the macs enriched folder')
    file_string = '*_peaks.bed'
    source_dir = '%smacsFolder' % (projectFolder)
    dest_dir = '%smacsEnriched' % (projectFolder)
    utils.link_files(file_string,source_dir,dest_dir)


    print('\n\n')
    print('#======================================================================')
    print('#==============III. DEFINING NB MYCN AND H3K27AC LANDSCAPE=============')
    print('#======================================================================')
    print('\n\n')

    #for enhancers
    enhancer_bashFileName,enhancer_region_map_path,namesList = define_enhancer_landscape(projectFolder,pipeline_dir,nb_all_chip_dataFile)

    #runs only if no output detected
    if not utils.checkOutput(enhancer_region_map_path,0,0):
        print(enhancer_bashFileName)
        os.system('bash %s' % (enhancer_bashFileName))

    sys.exit()
    #for mycn
    mycn_bashFileName,mycn_region_map_path,namesList = define_mycn_landscape(projectFolder,pipeline_dir,nb_all_chip_dataFile)

    if not utils.checkOutput(mycn_region_map_path,0,0):
        print(bashFileName)
        os.system('bash %s' % (mycn_bashFileName))
    
    #now we need to call the R script that creates the rank plots
    if utils.checkOutput(region_map_path,1,30): #set a wait time for 30 minutes
        print('Found NB_MYCN meta_rose landscape and running rank plot R code')
        name_string = ','.join(namesList) #provides the dataset names used
        rank_script_path = '%sr_scripts/1_nb_mycn_rank.R' % (projectFolder)
        r_cmd = 'Rscript %s %s %s %s' % (rank_script_path,region_map_path,name_string,projectFolder)
        print(r_cmd)
        os.system(r_cmd)

    print('#======================================================================')
    print('#=================IV. CREATING NB MYCN STATS TABLE=====================')
    print('#======================================================================')

    
    

#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~DEFINING NB MYCN BINDING LANDSCAPE~~~~~~~~~~~~~~~~~~~~~
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

    #this is the expeceted region map output
    region_map_path = '%s%s/%s_0KB_STITCHED_ENHANCER_REGION_MAP.txt' % (roseFolder,analysisName,analysisName)
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
#~~~~~~~~~~~~~~~~~~~~~~~~RANKING MYCN PEAKS BY STRENGTH~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#this is adapted from the lin 2012 code, but uses actual binding strength as
#opposed to enrichment to quantify peaks

def rank_eboxes(data_file,signal_path,names_list,macsFolder,genomeDirectory):

    '''
    uses the repeat filtered conserved MYCN sites and ranks eboxes within them
    by average background subtracted signal
    '''

    print('ok')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAKING NB MYCN STATS TABLE~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


def makeMycnPeakTable():

    '''
    making a table of repeat filtered mycn peaks w/ some additional stats
    mycn and h3k27ac signal is avg. background normalized across 4 samples
    active tss defined as the union of all H3K27ac occupied promoters in NB
    active enhancers defined as the union of all H3K27ac sites outside of promoters
    repeat cutoffs defined as cumulative fraction > 0.2 of LINE,SIMPLE,LTR
    also produces a bed and gff
    '''

    print('SETTING UP OUTPUT TABLE')
    outTable = [['PEAK_ID','CHROM','START','STOP','LENGTH','ACTIVE_TSS_OVERLAP','ENHANCER_OVERLAP','CPG_ISLAND_OVERLAP','CPG_ISLAND_FRACTION','GC_FREQ','MYCN_RANK','AVG_MYCN_SIGNAL','AVG_H3K27AC_SIGNAL','CANON_EBOX_COUNT','NONCANON_EBOX_COUNT','TOTAL_EBOX_COUNT']]



    #setting up the output
    outFile = '%stables/HG19_NB_MYCN_CONSERVED_EBOX_TABLE.txt' % (projectFolder)

    #input files
    maskFile = '/raider/index/hg19/Masks/hg19_encode_blacklist.bed'
    mycnSignalFile = '%ssignalTables/HG19_NB_MYCN_CONSERVED_-0_+0_SIGNAL_TABLE.txt' % (projectFolder)
    h3k27acSignalFile = '%ssignalTables/HG19_NB_MYCN_CONSERVED_-0_+0_H3K27AC_SIGNAL_TABLE.txt' % (projectFolder)
    mycnRankFile = '%smeta_rose/NB_MYCN/NB_MYCN_0KB_STITCHED_ENHANCER_REGION_RANK_CONSERVED.txt' % (projectFolder)
    activeGeneFile = '/grail/projects/mycn_cyl/gffListFiles/HG19_NB_H3K27AC_TRANSCRIBED_UNION.txt'
    cpgFile = '%sbeds/cpg_islands.bed' % (projectFolder)
    enhancerFile = '%smeta_rose/NB_H3K27AC/NB_H3K27AC_AllEnhancers.table.txt' % (projectFolder)

    print('LOADING MASK')
    maskTable = utils.parseTable(maskFile,'\t')
    maskLoci = []
    for line in maskTable:
        lineLocus = utils.Locus(line[0],line[1],line[2],'.')
        maskLoci.append(lineLocus)
    maskCollection = utils.LocusCollection(maskLoci,50)

    print('LOADING MYCN BINDING DATA')
    mycnSignalTable = utils.parseTable(mycnSignalFile,'\t')
    print(len(mycnSignalTable))

    print('LOADING MYCN RANK DATA')
    mycnRankTable = utils.parseTable(mycnRankFile,'\t')

    print('LOADING H3K27AC BINDING DATA')
    h3k27acSignalTable = utils.parseTable(h3k27acSignalFile,'\t')


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
    tss_50kb_loci = []
    for line in activeTable:
        tss_1kb_loci.append(utils.makeTSSLocus(line[1],startDict,1000,1000))
        tss_50kb_loci.append(utils.makeTSSLocus(line[1],startDict,50000,50000))

    tss_1kb_collection = utils.LocusCollection(tss_1kb_loci,50)
    tss_50kb_collection = utils.LocusCollection(tss_50kb_loci,50)

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
            except IndexError:
                print(line)
                sys.exit()
            if not tss_1kb_collection.getOverlap(lineLocus):
                enhancerLoci.append(lineLocus)
    enhancerCollection = utils.LocusCollection(enhancerLoci,50)
    print("AFTER FILTERING THESE MANY ARE LEFT")
    print(len(enhancerCollection))


    print('CLASSIFYING MYCN PEAKS')
    ticker = 0
    for i in range(1,len(mycnSignalTable)):
        if ticker%100 == 0:
            print(ticker)
        ticker +=1

        line = mycnSignalTable[i]        

        signalVector = [float(x) for x in line[2:]]
        

        peakID = line[0]
        locusString = line[1]
        chrom = locusString.split('(')[0]
        [start,stop] = [int(x) for x in line[1].split(':')[-1].split('-')]
        lineLocus = utils.Locus(chrom,start,stop,'.',peakID)

        if maskCollection.getOverlap(lineLocus,'both'):
            continue
        newLine = [peakID,chrom,start,stop,lineLocus.len()]

        if tss_1kb_collection.getOverlap(lineLocus,'both'):
            newLine.append(1)
        else:
            newLine.append(0)

        #find distance to nearest gene with 1mb
        tssSearchLocus = utils.makeSearchLocus(lineLocus,999000,999000)
        overlappingPromoters = tss_1kb_collection.getOverlap(tssSearchLocus,'both')
        if len(overlappingPromoters) == 0:
            tssDistance = 1000000
        else:
            tssStartsList = [(locus.start() + locus.end())/2 for locus in overlappingPromoters]
            tssDistanceList = [abs((lineLocus.start()+lineLocus.end())/2 - tssStart) for tssStart in tssStartsList]
            tssDistance = min(tssDistanceList)
        newLine.append(tssDistance)


        if enhancerCollection.getOverlap(lineLocus,'both'):
            newLine.append(1)
        else:
            newLine.append(0)

        if twistCollection.getOverlap(lineLocus,'both'):
            newLine.append(1)
        else:
            newLine.append(0)
        
        if cpgCollection.getOverlap(lineLocus,'both'):
            newLine.append(1)
        else:
            newLine.append(0)

        #now do fractional cpgOverlap
        overlappingCpGLoci = cpgCollection.getOverlap(lineLocus,'both')
        overlappingBases = 0
        for locus in overlappingCpGLoci:
            cpgStart = max(locus.start(),lineLocus.start())
            cpgEnd = min(locus.end(),lineLocus.end())
            overlappingBases += (cpgEnd-cpgStart)
        overlapFraction = float(overlappingBases)/lineLocus.len()
        
        newLine.append(round(overlapFraction,2))

        #now get the seq
        lineSeq = string.upper(utils.fetchSeq(genomeDirectory,chrom,start,stop,True))
        
        gcFreq = float(lineSeq.count('GC') + lineSeq.count('CG'))/len(lineSeq)
        newLine.append(gcFreq)
            
        mycnRankLine = mycnRankTable[i]
        mycnRank = numpy.mean([float(x) for x in mycnRankLine[6:]])
        newLine.append(mycnRank)

        h3k27acLine = h3k27acSignalTable[i]
        h3k27acSignal = numpy.mean([float(x) for x in h3k27acLine[2:]])
        newLine.append(h3k27acSignal)


        #this is where we add the mycn signal
        newLine += signalVector


        eboxMatchList = re.findall('CA..TG',lineSeq)
        if len(eboxMatchList) == 0:
            newLine += [0]*12
        else:
            #first fix the eboxes
            for j in range(len(eboxMatchList)):
                ebox = eboxMatchList[j]
                if eboxList.count(ebox) == 0:
                    eboxMatchList[j] = utils.revComp(ebox)

            #now count for each ebox
            eboxCountVector = []
            for ebox in eboxList:
                eboxCountVector.append(eboxMatchList.count(ebox))
            
            newLine += eboxCountVector

            #now add the total number of eboxes
            newLine.append(len(eboxMatchList))

            #now get the cumulative score
            eboxScore = 0.0
            for ebox in eboxMatchList:
                eboxScore += eboxDict[ebox]

            newLine.append(eboxScore)

        #now find the overlapping and proximal genes
        overlappingGenes = [startDict[locus.ID()]['name'] for locus in tss_1kb_collection.getOverlap(lineLocus,'both')]
        overlappingGenes = utils.updfniquify(overlappingGenes)
 
        proximalGenes = [startDict[locus.ID()]['name'] for locus in tss_50kb_collection.getOverlap(lineLocus,'both')]
        proximalGenes = [gene for gene in proximalGenes if overlappingGenes.count(gene) == 0]
        proximalGenes = utils.uniquify(proximalGenes)


        overlappingString = string.join(overlappingGenes,',')
        proximalString = string.join(proximalGenes,',')
        
        newLine += [overlappingString,proximalString]

        outTable.append(newLine)

    utils.unParseTable(outTable,outFile,'\t')
    





#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
