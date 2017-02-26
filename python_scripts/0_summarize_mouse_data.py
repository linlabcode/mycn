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





#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
