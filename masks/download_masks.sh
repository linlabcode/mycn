#!/usr/bin/bash

#links to download various masks used in the analysis

#encode blacklist
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz -O ./hg19_encode_blacklist.bed.gz

gunzip ./hg19_encode_blacklist.bed.gz
echo -e "chrM\t1\t16571\tMitochondria\t1000\t.\n" >> hg19_encode_blacklist.bed

#repeat masker
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz -O ./hg19_rmsk.txt.gz
gunzip ./hg19_rmsk.txt.gz

