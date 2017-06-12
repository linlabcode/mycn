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

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
        }


compareVectors <- function(v1,v2,name1,name2,title,nBins = 100,nIter = 1000,yMin='',yMax=''){
	
	#first one as a line graph
	
	v1Order = order(v1)
	colorSpectrum <- colorRampPalette(c("blue","black","red"))(100)

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
	plot(1:nBins,meanVector,ylim =c(yMin,yMax),ylab=name2,xlab='',cex=1,col=rgb(0.2,0.2,0.2,0.2),pch=16,xaxt='n',main=title)
	error.bar(1:nBins,meanVector,upperError,lowerError,length =0.01,col=rgb(0.2,0.2,0.2,0.2))	
	
	x=1:nBins
	lw1 = loess(meanVector~x)
	lines(x,lw1$fitted,col='orange',lwd=2)
	
	}








#==================================================================
#==================SAMPLE PARAMTERS FOR DEBUGGING==================
#==================================================================

setwd('~/Dropbox/mycn_cyl/')
twist_signal_path = './signalTables/HG19_SHEP21_0HR_TWIST_MYCN_INTERSECTION_-0_+0_SHEP21_TABLE_SIGNAL.txt'


#==================================================================
#========================DATA INPUT================================
#==================================================================


twist_table = read.delim(twist_signal_path,sep='\t')
twist_rows = which(apply(twist_table[,3:ncol(twist_table)],1,min) > 0)
twist_table = twist_table[twist_rows,]
twist_fold = log2(twist_table$SHEP21_2HR_TWIST/median(twist_table$SHEP21_2HR_TWIST)/twist_table$SHEP21_0HR_TWIST/median(twist_table$SHEP21_0HR_TWIST))
mycn_fold = log2(twist_table$SHEP21_2HR_MYCN_NOSPIKE/twist_table$SHEP21_0HR_MYCN_NOSPIKE)

mycn_2 = twist_table$SHEP21_2HR_MYCN_NOSPIKE/median(twist_table$SHEP21_2HR_MYCN_NOSPIKE)
mycn_0 = twist_table$SHEP21_0HR_MYCN_NOSPIKE/median(twist_table$SHEP21_0HR_MYCN_NOSPIKE)
mycn_fold = log2(mycn_2/mycn_0)


a = twist_table$SHEP21_2HR_TWIST/median(twist_table$SHEP21_2HR_TWIST)
b= twist_table$SHEP21_0HR_TWIST/median(twist_table$SHEP21_0HR_TWIST)

c = log2(a/b)


pdf(file='./figures/170321_twist_mycn_2hr.pdf',width = 8,height= 10)

compareVectors(mycn_fold,c,'twist 2 vs 0','mycn 2 vs 0','mycn change vs twist change @ 2 hours mycn shutdown',50)
dev.off()

twist_fold_24 = log2(twist_table$SHEP21_24HR_B_TWIST/twist_table$SHEP21_0HR_TWIST)
mycn_fold_24 = log2(twist_table$SHEP21_24HR_MYCN_NOSPIKE/twist_table$SHEP21_0HR_MYCN_NOSPIKE)

compareVectors(twist_fold_24,mycn_fold_24,'twist 2 vs 0','mycn 2 vs 0','mycn change vs twist change @ 2 hours mycn shutdown',50)

compareVectors(mycn_fold_24,twist_fold_24,'twist 2 vs 0','mycn 2 vs 0','mycn change vs twist change @ 2 hours mycn shutdown',50)


compareVectors(twist_fold,twist_table$SHEP21_2HR_MYCN_NOSPIKE,'twist 2 vs 0','mycn 2 vs 0','mycn change vs twist change @ 2 hours mycn shutdown',50)





