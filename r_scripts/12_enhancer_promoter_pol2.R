#12_enhancer_promoter_pol2.R


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
gene_table_path = args[6]

#get the cell count norm exp path
exp_norm_path = args[7]

#get the raw exp path
exp_raw_path = args[8]

#setting group names
group_string = args[9]
group_names = unlist(strsplit(group_string,','))
print('group names:')
print(group_names)

#setting the top number of genes to run on
top_count = as.numeric(args[10])
print('running analysis on top genes:')
print(top_count)

#setting the top number of genes to run on
nBins = as.numeric(args[11])
print('Dividing genes into N bins:')
print(nBins)

#set the analysis name
analysis_name = args[12]
print('analysis name:')
print(analysis_name)

#set the x vector time points
x_string = args[13]
x_vector = unlist(strsplit(group_string,','))
print('x axis labels')
print(x_vector)

#set the color string
color_string = args[14]
color_vector = as.numeric(unlist(strsplit(color_string,',')))
color = rgb(color_vector[1],color_vector[2],color_vector[3],maxColorValue = 255)
print('color is:')
print(color_vector)


#set the project folder
projectFolder = args[15]
print(projectFolder)




#==================================================================
#========================HELPER FUNCTIONS==========================
#==================================================================


makeCorrectedMatrix <- function(signal_table,scale_table,sample_names){
	signal_matrix = as.matrix(signal_table[,3:ncol(signal_table)])
	rownames(signal_matrix) = signal_table[,1]
	colnames(signal_matrix) = colnames(signal_table)[3:ncol(signal_table)]
	#now correct for scale factor
	
	for(i in 1:ncol(signal_matrix)){	
		name = as.character(colnames(signal_matrix)[i])
		scale_row = which(scale_table[,1]==name)
		if(length(scale_row)>0){
			rpm_scale_factor = as.numeric(scale_table[scale_row,4])
		}else{
			rpm_scale_factor = 1
		}
		
		signal_matrix[,i] = signal_matrix[,i] * rpm_scale_factor
	}
	
	#can't assume that columns are ordered in any particular way
	
	#create a final matrix
	final_matrix = matrix(nrow = nrow(signal_matrix),ncol=length(sample_names))
	colnames(final_matrix) = sample_names
	rownames(final_matrix) = rownames(signal_matrix)
	
	for(i in 1:length(sample_names)){
		name = sample_names[i]
		name_col = which(colnames(signal_matrix)==name)
		final_matrix[,i] = signal_matrix[,name_col] #background subtraction doesn't make sense for chiprx
		final_matrix[which(final_matrix[,i]<0),i] <-0
	}
	return(final_matrix)
}


makeRawMatrix <- function(signal_table,sample_names){
	signal_matrix = as.matrix(signal_table[,3:ncol(signal_table)])
	rownames(signal_matrix) = signal_table[,1]
	colnames(signal_matrix) = colnames(signal_table)[3:ncol(signal_table)]

	#can't assume that columns are ordered in any particular way
	
	#create a final matrix
	final_matrix = matrix(nrow = nrow(signal_matrix),ncol=length(sample_names))
	colnames(final_matrix) = sample_names
	rownames(final_matrix) = rownames(signal_matrix)
	
	for(i in 1:length(sample_names)){
		name = sample_names[i]
		name_col = which(colnames(signal_matrix)==name)
		final_matrix[,i] = signal_matrix[,name_col]
		final_matrix[which(final_matrix[,i]<0),i] <-0
	}
	return(final_matrix)
}

formattedBoxPlot <- function(m,title){
	
	#makes nice formatted box plots w/ scaling and stats
	#get a logical y limit for a box plot
	yMax = max(apply(m,2,quantile,c(0.95)))
	boxplot(m,cex=0,ylim = c(0,yMax),ylab='reads',main=title)
	
	#add text for some pvalues want them to be sequential
	p_vector = c()
	i=1
	while(i < ncol(m)){
		p_vector = c(p_vector,signif(t.test(m[,i], m[,(i+1)])$p.value,2))
		i = i+1
	}
	
	text(1.5,yMax,paste(p_vector,'----',collapse='  '))	
	
}


makePol2Matrix <- function(pol2_tss_table,pol2_body_table,transcribed_table,geneList,cut=.1,columns = c(1,2,3)){
	
	subsetRows = c()
	nameList = c()
	for(name in geneList){
		#name is a common name, need to grab refseq IDs
		refseqIDs = as.character(transcribed_table[which(transcribed_table[,3] == name),2])
		if(length(refseqIDs) > 0){
			for(refID in refseqIDs){
				refRow = which(rownames(pol2_tss_table)==refID)
				subsetRows = c(subsetRows,refRow)
				nameList = c(nameList,name)
			}
			
		}
	}
	
	subPol2Matrix = as.matrix(cbind(pol2_tss_table[subsetRows,columns],pol2_body_table[subsetRows,columns]))
	rownames(subPol2Matrix) = nameList
	
	cutRows = which(apply(subPol2Matrix[,1:3],1,mean) > cut)
	return(subPol2Matrix[cutRows,])
}






#==================================================================
#==================SAMPLE PARAMTERS FOR DEBUGGING==================
#==================================================================
setwd('~/Dropbox/mycn_cyl')

signal_path_tss = './signalTables/HG19_TSS_ALL_-300_+300_SHEP21_CHIPRX_TABLE_SIGNAL.txt'
signal_path_body = './signalTables/HG19_BODY_ALL_+300_+3000_SHEP21_CHIPRX_TABLE_SIGNAL.txt'

#signal_path_tss = './signalTables/HG19_TSS_ALL_-300_+300_SHEP21_TABLE_SIGNAL.txt'
#signal_path_body = './signalTables/HG19_BODY_ALL_+300_+3000_SHEP21_TABLE_SIGNAL.txt'


scale_path = './tables/HG19_SHEP21_CHIPRX_SCALE_FACTORS.txt'

transcribed_path = './geneListFolder/HG19_SHEP21_0HR_H3K27AC_NOSPIKE_ACTIVE.txt'

gene_path = './enhancerPromoter_resub/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_SHEP21_0HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_SHEP21_0HR_MYCN_NOSPIKE_TOP_5000_ORDERED.txt'


analysis_name ='SHEP21_POL2_RX'

box_plot_path = paste('./figures/12_',analysis_name,'_BOX.pdf',sep='')

samples_string = 'SHEP21_0HR_POL2_RX,SHEP21_2HR_POL2_RX,SHEP21_24HR_POL2_RX'

#samples_string = 'SHEP21_0HR_POL2_NOSPIKE_R2,SHEP21_2HR_POL2_NOSPIKE,SHEP21_24HR_POL2_NOSPIKE'

sample_names= unlist(strsplit(samples_string,','))


#set the project folder
projectFolder = './'
print(projectFolder)

#==================================================================
#========================DATA INPUT================================
#==================================================================

pol2_tss_table = read.delim(signal_path_tss,sep='\t')
pol2_body_table = read.delim(signal_path_body,sep='\t')

scale_table = read.delim(scale_path,sep='\t')
transcribed_table = read.delim(transcribed_path,header=FALSE,sep='\t')

gene_table = read.delim(gene_path,sep='\t')

#==================================================================
#============MAKING SCALED AND UNSCALED SIGNAL MATRICIES===========
#==================================================================

pol2_tss_scaled = makeCorrectedMatrix(pol2_tss_table,scale_table,sample_names)
pol2_tss_raw = makeRawMatrix(pol2_tss_table,sample_names)

pol2_body_scaled = makeCorrectedMatrix(pol2_body_table,scale_table,sample_names)
pol2_body_raw = makeRawMatrix(pol2_body_table,sample_names)

#==================================================================
#========================GETTING GENE LISTS========================
#==================================================================

total_signal = gene_table[,2] + gene_table[,3]
total_order = order(total_signal,decreasing=TRUE)
gene_table_ordered = gene_table[total_order,]

top_genes = as.character(gene_table_ordered[1:1000,1])
mid_genes = as.character(gene_table_ordered[2001:3000,1])
bottom_genes = as.character(gene_table_ordered[4001:5000,1])
other_genes = setdiff(transcribed_table[,3],gene_table[,1])

#==================================================================
#===============MAKING FOLD MATRICIES FOR SCALED DATA==============
#==================================================================


pol2_table_top = makePol2Matrix(pol2_tss_scaled, pol2_body_scaled,transcribed_table, top_genes)
pol2_table_mid = makePol2Matrix(pol2_tss_scaled, pol2_body_scaled,transcribed_table, mid_genes)
pol2_table_bottom = makePol2Matrix(pol2_tss_scaled, pol2_body_scaled,transcribed_table, bottom_genes)
pol2_table_other = makePol2Matrix(pol2_tss_scaled, pol2_body_scaled,transcribed_table, other_genes)


top_tss_fold = log2(pol2_table_top[,c(2,3)]/pol2_table_top[,1])
top_body_fold = log2(pol2_table_top[,c(5,6)]/pol2_table_top[,4])

mid_tss_fold = log2(pol2_table_mid[,c(2,3)]/pol2_table_mid[,1])
mid_body_fold = log2(pol2_table_mid[,c(5,6)]/pol2_table_mid[,4])

bottom_tss_fold = log2(pol2_table_bottom[,c(2,3)]/pol2_table_bottom[,1])
bottom_body_fold = log2(pol2_table_bottom[,c(5,6)]/pol2_table_bottom[,4])

other_tss_fold = log2(pol2_table_other[,c(2,3)]/pol2_table_other[,1])
other_body_fold = log2(pol2_table_other[,c(5,6)]/pol2_table_other[,4])

par(mfrow=c(1,2))
boxplot(other_tss_fold[,1], top_tss_fold[,1],other_tss_fold[,2], top_tss_fold[,2], cex=0,main='Scaled TSS',col=c('grey','red'),ylim=c(-2.1,1.3),names=c('other 2hr','top 2hr','other 24hr','top 24hr'),las=2)

boxplot(other_body_fold[,1], top_body_fold[,1],other_body_fold[,2], top_body_fold[,2], cex=0,main='Scaled Body',col=c('grey','red'),ylim=c(-2.1,1.3),names=c('other 2hr','top 2hr','other 24hr','top 24hr'),las=2)





#==================================================================
#=================MAKING FOLD MATRICIES FOR RAW DATA===============
#==================================================================

pol2_table_top_raw = makePol2Matrix(pol2_tss_raw, pol2_body_raw,transcribed_table, top_genes)
pol2_table_mid_raw = makePol2Matrix(pol2_tss_raw, pol2_body_raw,transcribed_table, mid_genes)
pol2_table_bottom_raw = makePol2Matrix(pol2_tss_raw, pol2_body_raw,transcribed_table, bottom_genes)
pol2_table_other_raw = makePol2Matrix(pol2_tss_raw, pol2_body_raw,transcribed_table, other_genes)


top_tss_fold_raw = log2(pol2_table_top_raw[,c(2,3)]/pol2_table_top_raw[,1])
top_body_fold_raw = log2(pol2_table_top_raw[,c(5,6)]/pol2_table_top_raw[,4])

mid_tss_fold_raw = log2(pol2_table_mid_raw[,c(2,3)]/pol2_table_mid_raw[,1])
mid_body_fold_raw = log2(pol2_table_mid_raw[,c(5,6)]/pol2_table_mid_raw[,4])

bottom_tss_fold_raw = log2(pol2_table_bottom_raw[,c(2,3)]/pol2_table_bottom_raw[,1])
bottom_body_fold_raw = log2(pol2_table_bottom_raw[,c(5,6)]/pol2_table_bottom_raw[,4])

other_tss_fold_raw = log2(pol2_table_other_raw[,c(2,3)]/pol2_table_other_raw[,1])
other_body_fold_raw = log2(pol2_table_other_raw[,c(5,6)]/pol2_table_other_raw[,4])

par(mfrow=c(1,2))
boxplot(other_tss_fold_raw[,1], top_tss_fold_raw[,1],other_tss_fold_raw[,2], top_tss_fold_raw[,2], cex=0,main='Raw TSS',col=c('grey','red'),ylim=c(-1.5,1),names=c('other 2hr','top 2hr','other 24hr','top 24hr'),las=2)

boxplot(other_body_fold_raw[,1], top_body_fold_raw[,1],other_body_fold_raw[,2], top_body_fold_raw[,2], cex=0,main='Raw Body',col=c('grey','red'),ylim=c(-1.5,1),names=c('other 2hr','top 2hr','other 24hr','top 24hr'),las=2)






#==================================================================
#==========================MAKING BOX PLOTS========================
#==================================================================

pdf(file=box_plot_path,width = 6,height=6)

par(mfrow=c(1,2))

boxplot(other_tss_fold[,1], top_tss_fold[,1],other_tss_fold[,2], top_tss_fold[,2], cex=0,main='Scaled TSS',col=c('grey','red'),ylim=c(-2.1,1.3),names=c('other 2hr','top 2hr','other 24hr','top 24hr'),las=2)

boxplot(other_body_fold[,1], top_body_fold[,1],other_body_fold[,2], top_body_fold[,2], cex=0,main='Scaled Body',col=c('grey','red'),ylim=c(-2.1,1.3),names=c('other 2hr','top 2hr','other 24hr','top 24hr'),las=2)

par(mfrow=c(1,2))

boxplot(other_tss_fold_raw[,1], top_tss_fold_raw[,1],other_tss_fold_raw[,2], top_tss_fold_raw[,2], cex=0,main='Raw TSS',col=c('grey','red'),ylim=c(-1.5,1),names=c('other 2hr','top 2hr','other 24hr','top 24hr'),las=2)

boxplot(other_body_fold_raw[,1], top_body_fold_raw[,1],other_body_fold_raw[,2], top_body_fold_raw[,2], cex=0,main='Raw Body',col=c('grey','red'),ylim=c(-1.5,1),names=c('other 2hr','top 2hr','other 24hr','top 24hr'),las=2)

dev.off()


# pdf(file=box_plot_path,width = 6,height=6)

# boxplot(bottom_tss_fold[,2], mid_tss_fold[,2], top_tss_fold[,2],bottom_body_fold[,2], mid_body_fold[,2], top_body_fold[,2], cex=0,main='Scaled',col=c('grey','white','red'),ylim=c(-1.5,1.5),names=c('bottom 24hr','mid 24hr','top 24hr','bottom 24hr','mid 24hr','top 24hr'),las=2)

# boxplot(bottom_tss_fold_raw[,2], mid_tss_fold_raw[,2], top_tss_fold_raw[,2],bottom_body_fold_raw[,2], mid_body_fold_raw[,2], top_body_fold_raw[,2], cex=0,main='RPM',col=c('grey','white','red'),ylim=c(-1.5,1.5),names=c('bottom 24hr','mid 24hr','top 24hr','bottom 24hr','mid 24hr','top 24hr'),las=2)

# dev.off()















#NS,*, then 
boxplot(other_body_fold[,1],bottom_body_fold[,1], mid_body_fold[,1],top_body_fold[,1], other_body_fold[,2],bottom_body_fold[,2], mid_body_fold[,2], top_body_fold[,2], cex=0,main='BODY',col=c('grey','grey','white','red'),ylim=c(-2,0.8),names=c('other 2hr','bottom 2hr','mid 2hr','top 2hr','other 24hr','bottom 24hr','mid 2hr','top 24hr'),las=2)



par(mfrow=c(1,2))
boxplot(other_tss_fold[,1],bottom_tss_fold[,1], mid_tss_fold[,1], top_tss_fold[,1], other_tss_fold[,2],bottom_tss_fold[,2], mid_tss_fold[,2], top_tss_fold[,2], cex=0,main='TSS',col=c('grey','grey','white','red'),ylim=c(-2.3,1.2),names=c('other 2hr','bottom 2hr','mid 2hr','top 2hr','other 24hr','bottom 24hr','mid 2hr','top 24hr'),las=2)
#NS,*, then 
boxplot(other_body_fold[,1],bottom_body_fold[,1], mid_body_fold[,1],top_body_fold[,1], other_body_fold[,2],bottom_body_fold[,2], mid_body_fold[,2], top_body_fold[,2], cex=0,main='BODY',col=c('grey','grey','white','red'),ylim=c(-2,0.8),names=c('other 2hr','bottom 2hr','mid 2hr','top 2hr','other 24hr','bottom 24hr','mid 2hr','top 24hr'),las=2)