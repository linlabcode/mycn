#7_enhancer_invasion_plots.R


#produces line plots of mean changes in gene expression
#binned by overall myc/mycn load


#==================================================================
#===========================DEPENDENCIES===========================
#==================================================================

#==================================================================
#======================LOADING ARGUMENTS===========================
#==================================================================


args <- commandArgs()

#get a gene table from enhancer promoter analysis
peak_0_path = args[6]
peak_1_path = args[7]
peak_2_path = args[8]

print('Using the following peak paths from enhancer promoter analysis')
print(peak_0_path)
print(peak_1_path)
print(peak_2_path)

analysis_name = args[9]
print('analysis name:')
print(analysis_name)

sample_string = args[10]
sample_names = unlist(strsplit(sample_string,','))
print('setting sample names as:')
print(sample_names)

top = args[11]
top = as.numeric(top)
print('analyzing top genes:')
print(top)

#set the project folder
projectFolder = args[12]
print(projectFolder)

#see if there are spikeys involved
scale_path = args[13]
if(scale_path == 'NONE'){
	useScale = FALSE
	}else{
	useScale = TRUE
	print('Using scale table:')
	print(scale_path)
	}



#==================================================================
#========================HELPER FUNCTIONS==========================
#==================================================================







#==================================================================
#==================SAMPLE PARAMTERS FOR DEBUGGING==================
#==================================================================

#setwd('/Volumes/grail/projects/mycn_resub/mycn/')
# #get a set of gene tables from enhancer promoter analysis
# gene_0_path = './enhancerPromoter/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_0HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_0HR_MYCN_NOSPIKE_GENE_TABLE.txt'

# gene_1_path = './enhancerPromoter/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_2HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_2HR_MYCN_NOSPIKE_GENE_TABLE.txt'

# gene_2_path = './enhancerPromoter/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_24HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_24HR_MYCN_NOSPIKE_GENE_TABLE.txt'


# peak_0_path = './enhancerPromoter/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_0HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_0HR_MYCN_NOSPIKE_PEAK_TABLE.txt'

# peak_1_path = './enhancerPromoter/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_2HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_2HR_MYCN_NOSPIKE_PEAK_TABLE.txt'

# peak_2_path = './enhancerPromoter/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_24HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_24HR_MYCN_NOSPIKE_PEAK_TABLE.txt'

# #set the project folder
# projectFolder = './'
# print(projectFolder)

# sample_names = c('SHEP21_0HR_MYCN_NOSPIKE','SHEP21_2HR_MYCN_NOSPIKE','SHEP21_24HR_MYCN_NOSPIKE')

# analysis_name = 'SHEP21_0HR_MYCN_NOSPIKE_REGIONS_NO_WCE_SHEP21_MYCN_NOSPIKE'

# top=5000
#==================================================================
#========================DATA INPUT================================
#==================================================================




peak_table_0 = read.delim(peak_0_path,sep='\t')
peak_table_1 = read.delim(peak_1_path,sep='\t')
peak_table_2 = read.delim(peak_2_path,sep='\t')

if(useScale == TRUE){
	    scale_table = read.delim(scale_path,sep='\t')
	    }



#==================================================================
#===================CALCULATING AUC AT EACH TIME===================
#==================================================================

distal_rows = which(peak_table_0[,6]==0)
tss_rows = which(peak_table_0[,6]==1)

peak_0_auc = peak_table_0$LENGTH * peak_table_0$SIGNAL
peak_1_auc = peak_table_1$LENGTH * peak_table_1$SIGNAL
peak_2_auc = peak_table_2$LENGTH * peak_table_2$SIGNAL

peak_0_tss_auc = peak_0_auc[tss_rows] 
peak_1_tss_auc = peak_1_auc[tss_rows] 
peak_2_tss_auc = peak_2_auc[tss_rows]

tss_order =order(peak_0_tss_auc,decreasing=TRUE)
peak_0_tss_auc = peak_0_tss_auc[tss_order]
peak_1_tss_auc = peak_1_tss_auc[tss_order]
peak_2_tss_auc = peak_2_tss_auc[tss_order]

#for distal
peak_0_distal_auc = peak_0_auc[distal_rows]
peak_1_distal_auc = peak_1_auc[distal_rows]
peak_2_distal_auc = peak_2_auc[distal_rows]

distal_order =order(peak_0_distal_auc,decreasing=TRUE)

peak_0_distal_auc = peak_0_distal_auc[distal_order]
peak_1_distal_auc = peak_1_distal_auc[distal_order]
peak_2_distal_auc = peak_2_distal_auc[distal_order]


	

#==================================================================
#=======================SCALING SIGNAL=============================
#==================================================================



if(useScale == TRUE){
	print('scaling datasets with a multiplicative scale factor')
	#figure out the rpm scale factor #which I believe should be a multiplicative scale factor

	#for the peak 0 dataset
	peak_0_name = sample_names[1]
	scale_row = which(scale_table[,1] == peak_0_name)
	scale_factor_0 = scale_table[scale_row,4]
	print(peak_0_name)
	print(scale_factor_0)

	peak_0_tss_auc = peak_0_tss_auc *scale_factor_0
	peak_0_distal_auc = peak_0_distal_auc * scale_factor_0

	#for the peak 1 dataset
	peak_1_name = sample_names[2]
	scale_row = which(scale_table[,1] == peak_1_name)
	scale_factor_1 = scale_table[scale_row,4]
	print(peak_1_name)
	print(scale_factor_1)

	peak_1_tss_auc = peak_1_tss_auc *scale_factor_1
	peak_1_distal_auc = peak_1_distal_auc * scale_factor_1

	#for the peak 2 dataset
	peak_2_name = sample_names[3]
	scale_row = which(scale_table[,1] == peak_2_name)
	scale_factor_2 = scale_table[scale_row,4]
	print(peak_2_name)
	print(scale_factor_2)

	peak_2_tss_auc = peak_2_tss_auc *scale_factor_2
	peak_2_distal_auc = peak_2_distal_auc * scale_factor_2

	}


#==================================================================
#====================MAKING PEAK BOX PLOTS=========================
#==================================================================



#first vs. second time point fold change for top peaks
distal_1_fold = log2(peak_1_distal_auc/peak_0_distal_auc)
distal_1_fold = distal_1_fold[which(is.infinite(distal_1_fold)==FALSE)]
distal_1_fold = distal_1_fold[which(is.na(distal_1_fold)==FALSE)]

tss_1_fold = log2(peak_1_tss_auc/peak_0_tss_auc)
tss_1_fold = tss_1_fold[which(is.infinite(tss_1_fold)==FALSE)]
tss_1_fold = tss_1_fold[which(is.na(tss_1_fold)==FALSE)]

distal_2_fold = log2(peak_2_distal_auc/peak_0_distal_auc)
distal_2_fold = distal_2_fold[which(is.infinite(distal_2_fold)==FALSE)]
distal_2_fold = distal_2_fold[which(is.na(distal_2_fold)==FALSE)]

tss_2_fold = log2(peak_2_tss_auc/peak_0_tss_auc)
tss_2_fold = tss_2_fold[which(is.infinite(tss_2_fold)==FALSE)]
tss_2_fold = tss_2_fold[which(is.na(tss_2_fold)==FALSE)]

distal_2_1_fold = log2(peak_2_distal_auc/peak_1_distal_auc)
distal_2_1_fold = distal_2_fold[which(is.infinite(distal_2_fold)==FALSE)]
distal_2_1_fold = distal_2_fold[which(is.na(distal_2_fold)==FALSE)]

tss_2_1_fold = log2(peak_2_tss_auc/peak_1_tss_auc)
tss_2_1_fold = tss_2_fold[which(is.infinite(tss_2_fold)==FALSE)]
tss_2_1_fold = tss_2_fold[which(is.na(tss_2_fold)==FALSE)]


#set logical yMin and yMax
yMin = min(c(quantile(distal_1_fold,0.025),quantile(distal_2_fold,0.025),quantile(tss_1_fold,0.025),quantile(tss_2_fold,0.025),quantile(distal_2_1_fold,0.025),quantile(tss_2_1_fold,0.025)))
yMax = max(c(quantile(distal_1_fold,0.975),quantile(distal_2_fold,0.975),quantile(tss_1_fold,0.975),quantile(tss_2_fold,0.975),quantile(distal_2_1_fold,0.975),quantile(tss_2_1_fold,0.975)))

print('Using y limits:')
print(yMin)
print(yMax)

plot_title = paste(projectFolder,'figures/7_',analysis_name,'.pdf',sep='')

pdf(file=plot_title,width = 15,height = 5)
#for top peaks
par(mfrow=c(1,3))
plot_name_1 = paste(sample_names[2],' vs. ',sample_names[1],'\ntop ',top,' fold change',sep='')
boxplot(tss_1_fold[1:top],distal_1_fold[1:top],cex=0,names= c('TSS','DISTAL'),main=plot_name_1,ylim=c(yMin,yMax),ylab='log2 signal fold change (rpm)')
abline(h=0)
p1 = t.test(tss_1_fold[1:top],distal_1_fold[1:top])$p.value
text(0.5,1,p1)

plot_name_2 = paste(sample_names[3],' vs. ',sample_names[1],'\ntop ',top,' fold change',sep='')
boxplot(tss_2_fold[1:top],distal_2_fold[1:top],cex=0,names= c('TSS','DISTAL'),main=plot_name_2,ylim = c(yMin,yMax),ylab='log2 signal fold change (rpm)')
abline(h=0)
p2 = t.test(tss_2_fold[1:top],distal_2_fold[1:top])$p.value
text(0.5,1,p2)

plot_name_3 = paste(sample_names[3],' vs. ',sample_names[2],'\ntop ',top,' fold change',sep='')
boxplot(tss_2_fold[1:top],distal_2_fold[1:top],cex=0,names= c('TSS','DISTAL'),main=plot_name_3,ylim = c(yMin,yMax),ylab='log2 signal fold change (rpm)')
abline(h=0)
p3 = t.test(tss_2_1_fold[1:top],distal_2_1_fold[1:top])$p.value
text(0.5,1,p3)

#for all peaks
par(mfrow=c(1,3))
plot_name_1 = paste(sample_names[2],' vs. ',sample_names[1],'\nall fold change',sep='')
boxplot(tss_1_fold,distal_1_fold,cex=0,names= c('TSS','DISTAL'),main=plot_name_1,ylim=c(yMin,yMax),ylab='log2 signal fold change (rpm)')
abline(h=0)
p1 = t.test(tss_1_fold,distal_1_fold)$p.value
text(0.5,1,p1)

plot_name_2 = paste(sample_names[3],' vs. ',sample_names[1],'\nall fold change',sep='')
boxplot(tss_2_fold,distal_2_fold,cex=0,names= c('TSS','DISTAL'),main=plot_name_2,ylim = c(yMin,yMax),ylab='log2 signal fold change (rpm)')
abline(h=0)
p2 = t.test(tss_2_fold,distal_2_fold)$p.value
text(0.5,1,p2)

plot_name_3 = paste(sample_names[3],' vs. ',sample_names[2],'\nall fold change',sep='')
boxplot(tss_2_1_fold,distal_2_1_fold,cex=0,names= c('TSS','DISTAL'),main=plot_name_3,ylim = c(yMin,yMax),ylab='log2 signal fold change (rpm)')
abline(h=0)
p3 = t.test(tss_2_1_fold,distal_2_1_fold)$p.value
text(0.5,1,p3)

dev.off()


	

#==================================================================
#===========================OLD====================================
#==================================================================

# #loading the gene table
# gene_table_0 = read.delim(gene_0_path,sep='\t')
# gene_table_1 = read.delim(gene_1_path,sep='\t')
# gene_table_2 = read.delim(gene_2_path,sep='\t')

# gene_sig_0 = gene_table_0[,2]+ gene_table_0[,3]
# gene_0_order =order(gene_sig_0,decreasing=TRUE)

# plot(1:length(gene_sig_0),gene_sig_0[gene_0_order])

# #plot the avg fold change of tss and distal as a function of the order

# binSize = 100
# distal_1_fold = c()
# tss_1_fold = c()
# distal_2_fold = c()
# tss_2_fold = c()
# i = 1
# while((i+binSize) < length(gene_0_order)){
	
	# distal_fold = log2(mean(gene_table_1[gene_0_order[i:(i+binSize)],3])/mean(gene_table_0[gene_0_order[i:(i+binSize)],3]))
	# distal_1_fold = c(distal_1_fold, distal_fold)
	
	# distal_fold = log2(mean(gene_table_2[gene_0_order[i:(i+binSize)],3])/mean(gene_table_0[gene_0_order[i:(i+binSize)],3]))
	# distal_2_fold = c(distal_2_fold, distal_fold)
	

	# tss_fold = log2(mean(gene_table_1[gene_0_order[i:(i+binSize)],2])/mean(gene_table_0[gene_0_order[i:(i+binSize)],2]))
	# tss_1_fold = c(tss_1_fold, tss_fold)
	
	# tss_fold = log2(mean(gene_table_2[gene_0_order[i:(i+binSize)],2])/mean(gene_table_0[gene_0_order[i:(i+binSize)],2]))
	# tss_2_fold = c(tss_2_fold, tss_fold)	
	# i = i +binSize
# } 


# par(mfrow=c(1,2))
# plot(1:length(distal_1_fold),tss_1_fold,type='b',col='red',ylim = c(yMin,yMax))
# lines(1:length(distal_1_fold),distal_1_fold,type='b',col='blue')

# plot(1:length(distal_2_fold),tss_2_fold,type='b',col='red',ylim = c(yMin,yMax))
# lines(1:length(distal_2_fold),distal_2_fold,type='b',col='blue')







