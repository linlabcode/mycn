#11_be2c_tss_clustering.R


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




makeRawMatrix <- function(signal_table,sample_names,useBackground=FALSE){
	signal_matrix = as.matrix(signal_table[,3:ncol(signal_table)])
	rownames(signal_matrix) = signal_table[,2]
	colnames(signal_matrix) = colnames(signal_table)[3:ncol(signal_table)]

	#can't assume that columns are ordered in any particular way
	
	#create a final matrix
	final_matrix = matrix(nrow = nrow(signal_matrix),ncol=length(sample_names))
	colnames(final_matrix) = sample_names
	rownames(final_matrix) = rownames(signal_matrix)
	
	for(i in 1:length(sample_names)){
		name = sample_names[i]
		name_col = intersect(grep(name,colnames(signal_matrix)),grep('MYC',colnames(signal_matrix)))
		background_col = setdiff(grep(name,colnames(signal_matrix)),grep('MYC',colnames(signal_matrix)))
		if(useBackground == TRUE){
			final_matrix[,i] = signal_matrix[,name_col]-signal_matrix[, background_col]

			}else{
			final_matrix[,i] = signal_matrix[,name_col]
			}
		final_matrix[which(final_matrix[,i]<0),i] <-0
	}
	return(final_matrix)
}



makeSampleHeatmap <- function(signal_matrix,color='red'){
	

	sampleDist = as.dist(1-cor(signal_matrix))
	sampleHC = hclust(sampleDist)
	sampleOrder = sampleHC$order
	
	#distance by samples
	sampleSimMatrix = cor(signal_matrix)	
	
	#Set the color spectrum
	colorSpectrum <- colorRampPalette(c("white",color))(100)	
	
	#setting a color data range
	minValue <- 0
	maxValue <- 1
	color_cuts <- seq(minValue,maxValue,length=100)
	color_cuts <- c(0, color_cuts,1)

	#add one extra min color to even out sampling
	colorSpectrum <- c(colorSpectrum[1],colorSpectrum)
	
	layout(matrix(data=c(1,2,2,2,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,4,4,4),ncol= 7))

	par(mar=c(8,9,5,9))
	plot(as.dendrogram(sampleHC),ylab='Distance')
	par(mar=c(4,2,2,2))

	plot(as.dendrogram(sampleHC),horiz=TRUE,xlab='Distance',leaflab='none')
	par(mar=c(6,4,4,2))

	image(1:ncol(sampleSimMatrix),1:nrow(sampleSimMatrix),t(sampleSimMatrix[sampleOrder, sampleOrder]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')

	par(mar=c(6,5,4,2))
	image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Similarity")
}








clusterSampleFile = paste('./figures/donowitz_differential_',sig_cut,'_pval_',fold_cut,'_fold_clusterSamples.pdf',sep='')

pdf(file = clusterSampleFile,width = 10,height =10)

layout(matrix(data=c(1,2,2,2,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,4,4,4),ncol= 7))

dev.off()





#==================================================================
#==================SAMPLE PARAMTERS FOR DEBUGGING==================
#==================================================================
setwd('~/Dropbox/mycn_cyl')

nb_signal_path_tss = './signalTables/HG19_MYC_TSS_REGIONS_-0_+0_NB_ALL_SIGNAL.txt'
nb_signal_path_distal = './signalTables/HG19_MYC_DISTAL_REGIONS_-0_+0_NB_ALL_SIGNAL.txt'

other_signal_path_tss = './signalTables/HG19_MYC_TSS_REGIONS_-0_+0_MYC_HIGH_TABLE_SIGNAL.txt'
other_signal_path_distal = './signalTables/HG19_MYC_DISTAL_REGIONS_-0_+0_MYC_HIGH_TABLE_SIGNAL.txt'



analysis_name ='MYC_HIGH'


nb_string = 'BE2C,KELLY,NGP,SHEP21_0HR'
other_string = 'H2171,MM1S,P493,U87'

nb_sample_names= unlist(strsplit(nb_string,','))
other_sample_names= unlist(strsplit(other_string,','))


#set the project folder
projectFolder = './'
print(projectFolder)

#==================================================================
#========================DATA INPUT================================
#==================================================================

nb_signal_table_tss = read.delim(nb_signal_path_tss,sep='\t')
nb_signal_table_distal = read.delim(nb_signal_path_distal,sep='\t')

other_signal_table_tss = read.delim(other_signal_path_tss,sep='\t')
other_signal_table_distal = read.delim(other_signal_path_distal,sep='\t')

#establish the output names
tss_output = paste(projectFolder,'figures/13_',analysis_name,'_TSS_SAMPLES.pdf',sep='')
distal_output = paste(projectFolder,'figures/13_',analysis_name,'_DISTAL_SAMPLES.pdf',sep='')


#==================================================================
#============MAKING RAW AND BACKGROUND CORRECTED MATRICIES=========
#==================================================================

nb_tss = makeRawMatrix(nb_signal_table_tss, nb_sample_names,FALSE)
nb_tss_background = makeRawMatrix(nb_signal_table_tss, nb_sample_names,TRUE)

nb_distal = makeRawMatrix(nb_signal_table_distal, nb_sample_names,FALSE)
nb_distal_background = makeRawMatrix(nb_signal_table_distal, nb_sample_names,TRUE)


other_tss = makeRawMatrix(other_signal_table_tss, other_sample_names,FALSE)
other_tss_background = makeRawMatrix(other_signal_table_tss, other_sample_names,TRUE)

other_distal = makeRawMatrix(other_signal_table_distal, other_sample_names,FALSE)
other_distal_background = makeRawMatrix(other_signal_table_distal, other_sample_names,TRUE)

#==================================================================
#=============AVG NB SIGNAL AND MAKING COMBINED MATRIX=============
#==================================================================

nb_tss_avg = apply(nb_tss,1,mean)
nb_distal_avg = apply(nb_distal,1,mean)

nb_tss_avg_background = apply(nb_tss_background,1,mean)
nb_distal_avg_background = apply(nb_distal_background,1,mean)

composite_tss = cbind(other_tss,nb_tss_avg)
composite_distal = cbind(other_distal,nb_distal_avg)

composite_tss_background = cbind(other_tss_background,nb_tss_avg_background)
composite_distal_background = cbind(other_distal_background,nb_distal_avg_background)



#==================================================================
#=======================MAKING SAMPLE HEATMAPS=====================
#==================================================================

pdf(file=tss_output,width = 10,height=10)
makeSampleHeatmap(composite_tss_background,'red')
dev.off()

pdf(file=distal_output,width = 10,height=10)

makeSampleHeatmap(composite_distal_background,'red')
dev.off()


