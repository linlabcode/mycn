#!/usr/bin/bash

#links to download various masks used in the analysis

#encode blacklist
echo "downloading ENCODE blacklist"
see http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/hg19-blacklist-README.pdf
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz -O ./hg19_encode_blacklist.bed.gz
gunzip ./hg19_encode_blacklist.bed.gz
#adding the mitochondria genome to the blacklist
echo -e "chrM\t1\t16571\tMitochondria\t1000\t.\n" >> hg19_encode_blacklist.bed

#repeat masker
echo "downloading repeatmask"
http://www.repeatmasker.org/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz -O ./hg19_rmsk.txt.gz
gunzip ./hg19_rmsk.txt.gz

#now format into a bed
echo "pulling LINE, Simple_repeat, and LTR class into beds"
awk '{ if  ($12 == "LINE") print $6"\t"$7"\t"$8"\t"$12"\t"$2}' ./hg19_rmsk.txt > hg19_line_rmsk.bed
awk '{ if  ($12 == "Simple_repeat") print $6"\t"$7"\t"$8"\t"$12"\t"$2}' ./hg19_rmsk.txt > hg19_simple_repeat_rmsk.bed
awk '{ if  ($12 == "LTR") print $6"\t"$7"\t"$8"\t"$12"\t"$2}' ./hg19_rmsk.txt > hg19_ltr_rmsk.bed

#bgzip and tabix index
echo "bgzip and tabix indexing beds"
bgzip ./hg19_line_rmsk.bed
bgzip ./hg19_simple_repeat_rmsk.bed
bgzip ./hg19_ltr_rmsk.bed

tabix -p bed ./hg19_line_rmsk.bed.gz
tabix -p bed ./hg19_simple_repeat_rmsk.bed.gz
tabix -p bed ./hg19_ltr_rmsk.bed.gz


echo "all done setting up repeats"
