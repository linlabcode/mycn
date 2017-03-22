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

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
        }

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
                     function(x)
                       rgb(x[1], x[2], x[3], alpha=alpha))
}


compareVectors <- function(v1,v2,name1,name2,title,plot_color,nBins = 50,nIter = 100,yMin='',yMax=''){
	       	
	#first one as a line graph
	v1Order = order(v1)
	colorSpectrum <- colorRampPalette(c("green","black","black","red"))(100)

	#setting a color data range
	minValue <- -2
	maxValue <- 2
	color_cuts <- seq(minValue,maxValue,length=100)
	color_cuts <- c(min(v1,na.rm=TRUE), color_cuts,max(v1,na.rm=TRUE))


	#add one extra min color to even out sampling
	colorSpectrum <- c(colorSpectrum[1],colorSpectrum[1],colorSpectrum)

	colorVector = c()
	for(i in v1Order){
		delta = v1[i]
		color = colorSpectrum[max(which(color_cuts <= delta))]
		colorVector =c(colorVector,color)
	}
	
	#set the layout
	m =matrix(data = c(1,1,2,2,2,2),ncol=1)
	layout(m)
	
	
	plot(1:length(v1),v1[v1Order],ylim =c(-1.2*max(abs(v1)),1.2*max(abs(v1))),type='l',lwd=2,ylab=name1,xlab=paste(c('Regions ranked by',name1),sep=' '))
	lines(1:length(v1),v1[v1Order],type='h',lwd=3,col=colorVector)
	abline(h=0)
	
	#now for variable #2
	
	#matrix to bind the data
	binSize = length(v1)/nBins
	binMatrix = matrix(ncol = nBins,nrow=binSize)
	i = 1
	for(i in 1:nBins){
		start = 1 + (i-1)*binSize
		stop = i*binSize
		
		binMatrix[,i] = as.numeric(v2[v1Order[start:stop]])
		
		
	}
	
	meanMatrix = matrix(ncol = nBins,nrow=nIter)
	
	for(i in 1:nIter){
		for(j in 1:nBins){
			meanMatrix[i,j] = mean(sample(binMatrix[,j],binSize,replace=TRUE))
		}
	}
	
	meanVector = apply(meanMatrix,2,mean)
	upperError = apply(meanMatrix,2,quantile,probs=0.975) - apply(meanMatrix,2,mean)
	lowerError = apply(meanMatrix,2,mean) - apply(meanMatrix,2,quantile,probs=0.025)
	
	if(yMax == ''){
		yMax = 1.02* max(apply(meanMatrix,2,quantile,probs=0.975))
	}
	
	if(yMin == ''){
		yMin = 0.98 * min(apply(meanMatrix,2,quantile,probs=0.025))
	}
	plot(1:nBins,meanVector,ylim =c(yMin,yMax),ylab=name2,xlab='',cex=1,col=add.alpha(plot_color,.2),pch=16,xaxt='n',main=title)
	error.bar(1:nBins,meanVector,upperError,lowerError,length =0.05,col=add.alpha(plot_color,.2))	
	
	x=1:nBins
	lw1 = loess(meanVector~x)
	lines(x,lw1$fitted,col = plot_color,lwd=2)
	
	}









#==================================================================
#==================SAMPLE PARAMTERS FOR DEBUGGING==================
#==================================================================
setwd('~/Dropbox/mycn_cyl')

tss_signal_path = './signalTables/HG19_TSS_ALL_-1000_+1000_BE2C_TABLE_SIGNAL.txt'


cluster_heatmap_path = paste(outputFolder,genome,'_',analysisName,"_clusterSamples.pdf",sep='')
cluster_heatmap_path = './figures/11_be2c_tss_heatmap.pdf'


compare_plot_path = './figures/11_be2c_tss_compare.pdf'
#set the project folder
projectFolder = './'
print(projectFolder)

#==================================================================
#========================DATA INPUT================================
#==================================================================
tss_table = read.delim(tss_signal_path,sep='\t')

tss_matrix = as.matrix(tss_table[,c(3,4,5,6,9,10)])-tss_table[,7]
tss_matrix[which(tss_matrix < 0 )] <- 0 

tss_median_matrix = as.matrix(tss_matrix)
for(i in 1:ncol(tss_matrix)){
	
	tss_median_matrix[,i] = tss_matrix[,i]/median(tss_matrix[,i])
}


#==================================================================
#===================PERFORMING CLUSTERING==========================
#==================================================================

#just as in the clustering code
#distance by samples
sampleDist = as.dist(1-cor(tss_median_matrix,method='pearson'))

sampleHC = hclust(sampleDist)

#===================================================================
#=======================SAMPLE CLUSTER ORDERING=====================
#===================================================================

sampleOrder = sampleHC$order

#===================================================================
#=================MAKING SAMPLE PAIRWISE HEATMAP====================
#===================================================================
#distance by samples
sampleSimMatrix = cor(tss_median_matrix,method='pearson')


#Set the color spectrum
colorSpectrum <- colorRampPalette(c("white","white","red"))(100)

#setting a color data range
minValue <- .1
maxValue <- .9
color_cuts <- seq(minValue,maxValue,length=100)
color_cuts <- c(0, color_cuts,1)

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)

pdf(file = cluster_heatmap_path,width = 10,height =10)

layout(matrix(data=c(1,2,2,2,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,4,4,4),ncol= 7))
par(mar=c(8,9,5,9))
plot(as.dendrogram(sampleHC),ylab='Distance')
par(mar=c(4,2,2,2))

plot(as.dendrogram(sampleHC),horiz=TRUE,xlab='Distance',leaflab='none')
par(mar=c(6,4,4,2))

image(1:ncol(sampleSimMatrix),1:nrow(sampleSimMatrix),t(sampleSimMatrix[sampleOrder,sampleOrder]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')

par(mar=c(6,5,4,2))

image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Similarity")
dev.off()
#==================================================================
#===============COMPARING MYCN TO BRD4 AND OTHERS==================
#==================================================================
top = 20000

top_order = order(tss_matrix[,5],decreasing=TRUE)
top_rows = top_order[100:top]

pdf(file = compare_plot_path,width = 10,height =8)
compareVectors(tss_matrix[top_rows,5],tss_matrix[top_rows,1],'BE2C_MYCN','BE2C_BRD4','MYCN vs. BRD4 TOP20K TSS','purple',yMin=0,yMax=0.7)
compareVectors(tss_matrix[top_rows,5],tss_matrix[top_rows,6],'BE2C_MYCN','BE2C_POL2','MYCN vs. POL2 TOP20K TSS','black',yMin=0,yMax=0.7)

compareVectors(tss_matrix[top_rows,5],tss_matrix[top_rows,3],'BE2C_MYCN','BE2C_H3K27ME3','MYCN vs. H3K27ME3 TOP20K TSS','grey',yMin=0,yMax=0.7)
dev.off()
