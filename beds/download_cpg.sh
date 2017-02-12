#!/usr/bin/bash
#downloads and formats cpg island data into a bed

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz -O cpgIslandExt.txt.gz

gunzip cpgIslandExt.txt.gz

awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5}' ./cpgIslandExt.txt > hg19_cpg_islands.bed
