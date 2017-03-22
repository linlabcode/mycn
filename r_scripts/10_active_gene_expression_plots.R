
setwd('~/Dropbox/mycn_cyl')
library(ggplot2)
library(scatterplot3d)
library(rgl)
library(gplots)

#=========================================================
#========================HELPER FUNCTIONS=================
#=========================================================

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
        }





makeExpMatrix <- function(expTable,geneList,cut=10){
	
	subsetRows = c()
	for(name in geneList){

		expRow = which(rownames(expTable) == name)
		if(length(expRow) > 0){
			subsetRows = c(subsetRows,expRow)
		}
	}
	
	subExpMatrix = as.matrix(expTable[subsetRows,])
	
	cutRows = which(apply(subExpMatrix,1,max) > cut)
	return(subExpMatrix[cutRows,])
}

getExpRows <- function(expTable,geneList,cut=10){
	
	subsetRows = c()
	for(name in geneList){

		expRow = which(rownames(expTable) == name)
		if(length(expRow) > 0){
			subsetRows = c(subsetRows,expRow)
		}
	}
	
	subExpMatrix = as.matrix(expTable[subsetRows,])
	
	cutRows = which(apply(subExpMatrix,1,max) > cut)
	subsetRows[cutRows]
}




plotGene <- function(fullExpTable,gene,color){
	
	
	timePoints = c(0,2,4,6,8,16,24)
	gene_index= which(rownames(fullExpTable)==gene) 
	
	meanVector = c(mean(as.numeric(fullExpTable[gene_index,c(1,2)])))
	sdVector = c(sd(as.numeric(fullExpTable[gene_index,c(1,2)])))
	
	for(i in seq(3,18,3)){
		expMean = mean(as.numeric(fullExpTable[gene_index,i:(i+2)]))
		expSD = sd(as.numeric(fullExpTable[gene_index,i:(i+2)]))
		
		meanVector = c(meanVector,expMean)
		sdVector = c(sdVector,expSD)

	}
	maxIndex = which.max(meanVector)
	yMax = 1.1*(meanVector[maxIndex] + sdVector[maxIndex])
	
	minIndex = which.min(meanVector)
	yMin = .6*(meanVector[minIndex] - sdVector[minIndex])

	plot(timePoints,meanVector,type='b',ylim =c(yMin,yMax),ylab='Cell count normalized fpkm',main=gene,xlab='Hours after MYCN shutdown',cex=0,lwd=2,col=color)
	
	error.bar(timePoints,meanVector,sdVector,length =0.01,col=color,lwd=2)	

	
}




plotFoldGenes <- function(foldExpTable,gene_list,color_vector,yMin,yMax){
	
	
	timePoints = c(0,2,4,6,8,16,24)
	
	plot(0,0,xlim=c(0,max(timePoints)),ylim = c(yMin,yMax),ylab='Cell count normalized fpkm',main=gene,xlab='Hours after MYCN shutdown',cex=0,xaxt='n')
	axis(1,timePoints,timePoints)
	legend(0,0.5,gene_list,col=color_vector,lwd=2)
	for(i in 1:length(gene_list)){
		
		gene = gene_list[i]
		color = color_vector[i]
		
		gene_index= which(rownames(foldExpTable)==gene) 
	
		meanVector = c(mean(as.numeric(foldExpTable[gene_index,c(1,2)])))
		sdVector = c(sd(as.numeric(foldExpTable[gene_index,c(1,2)])))
	
		for(i in seq(3,18,3)){
			expMean = mean(as.numeric(foldExpTable[gene_index,i:(i+2)]))
			expSD = sd(as.numeric(foldExpTable[gene_index,i:(i+2)]))
		
			meanVector = c(meanVector,expMean)
			sdVector = c(sdVector,expSD)

		}
		lines(timePoints, meanVector,type='b',ylim =c(yMin,yMax),ylab='Cell count normalized fpkm',main=gene,xlab='Hours after MYCN shutdown',cex=0,lwd=2,col=color)
	
		error.bar(timePoints, meanVector, sdVector,length =0.01,col=color,lwd=2)	
	
	}
	
}
#========================================================
#========================DATA INPUT======================
#========================================================



transcribedTable = read.delim('gffListFiles/HG19_NB_H3K27AC_TRANSCRIBED_UNION.txt',sep='\t',header=FALSE)
transcribedGeneList = as.character(transcribedTable[,3])
#let's just look at the top N mycn peaks

#mean expression w/ replicate removed across all timepoints
expTable = read.delim('./expression/shep21_cufflinks_no_rep2/SHEP21_cuffnorm/output/SHEP21_all_fpkm_means.txt',header=TRUE)

#normExpTable = read.delim('./expression/SHEP21_all_fpkm_exprs_norm.txt',header=TRUE,sep='\t')
normExpTable = read.delim('./expression/shep21_cufflinks_no_rep2/SHEP21_cuffnorm/output/SHEP21_all_fpkm_exprs_norm.txt',sep='\t')
fullNormExpTable = rawExpTable[,seq(2,ncol(normExpTable))]

rawExpTable = read.delim('./expression/shep21_cufflinks_no_rep2/SHEP21_cuffnorm/output/SHEP21_all_fpkm_exprs_raw.txt',sep='\t')
fullRawExpTable = rawExpTable[,seq(1,ncol(rawExpTable))]

#converting this into a proper matrix
expMatrix = as.matrix(expTable)

exp_rows = getExpRows(expTable,transcribedGeneList,10)
transcribedExpMatrix = makeExpMatrix(expTable,transcribedGeneList,10)
#========================================================
#=========MAKING EXP MATRIX FOR NON SPIKEY DATA==========
#========================================================



rawExpMatrix = matrix(nrow= nrow(fullRawExpTable),ncol=7)
timePoints = c('0HR','2HR','4HR','6HR','8HR','16HR','24HR')
for(i in 1:length(timePoints)){
	
	cols = grep(as.character(timePoints[i]),as.character(colnames(fullRawExpTable)))
	rawExpMatrix[,i] = apply(fullRawExpTable[,cols],1,mean)
}
rownames(rawExpMatrix) = rownames(fullRawExpTable)
colnames(rawExpMatrix) = colnames(expMatrix)


transcribedExpMatrix = expMatrix[exp_rows,]
transcribedRawExpMatrix = rawExpMatrix[exp_rows,]

#========================================================
#==========================FIGURES=======================
#========================================================



pdf(file='./figures/170304_nb_mycn_expression_examples_norm.pdf',width = 5,height =4)
plotGene(fullExpTable,'HAND2','red')
plotGene(fullExpTable,'TWIST1','red')
plotGene(fullExpTable,'CDK4','red')
plotGene(fullExpTable,'ID1','red')
plotGene(fullExpTable,'ID2','red')
plotGene(fullExpTable,'NPM1','red')
plotGene(fullExpTable,'RPL22','red')

dev.off()


pdf(file='./figures/170304_nb_mycn_expression_examples_raw.pdf',width = 5,height =4)
plotGene(fullRawExpTable,'HAND2','red')
plotGene(fullRawExpTable,'TWIST1','red')
plotGene(fullRawExpTable,'CDK4','red')
plotGene(fullRawExpTable,'ID1','red')
plotGene(fullRawExpTable,'ID2','red')
plotGene(fullRawExpTable,'NPM1','red')
plotGene(fullRawExpTable,'RPL22','red')

dev.off()



foldExpTable = log2(fullExpTable/apply(fullExpTable[,c(1,2)],1,mean))
foldRawExpTable = log2(fullRawExpTable/apply(fullRawExpTable[,c(1,2)],1,mean))

pdf(file='./figures/170304_nb_mycn_expression_examples_fold_norm.pdf',width =5,height =4)
gene_list = c('NPM1','ID2')
color_vector =c('red','blue')
yMin= -1.2
yMax = 0.5
plotFoldGenes(foldExpTable,gene_list,color_vector,yMin,yMax)
#plotFoldGenes(foldRawExpTable,gene_list,color_vector,yMin,yMax)



gene_list = c('NPM1','HAND2')
color_vector =c('red','blue')
yMin= -1.2
yMax = 0.5
plotFoldGenes(foldExpTable,gene_list,color_vector,yMin,yMax)
#plotFoldGenes(foldRawExpTable,gene_list,color_vector,yMin,yMax)
dev.off()
#========================================================
#==================FOLD MATRIX BOX=======================
#========================================================

pdf(file='./figures/170304_shep21_expression_box_norm.pdf',width = 6,height = 5)
foldMatrix = log2(transcribedExpMatrix/transcribedExpMatrix[,1])
binaryMatrix = apply(foldMatrix,2,is.finite)
finite_rows = which(apply(binaryMatrix,1,min)>0)
boxplot(foldMatrix[finite_rows,2:7],cex=0,ylim=c(-1.8,1.3),ylab='Expression log2 fold change vs. 0hr',main ='active genes n=8,043')
abline(h=0)
dev.off()


pdf(file='./figures/170304_shep21_expression_box_raw.pdf',width = 6,height = 5)
foldMatrix = log2(transcribedRawExpMatrix/transcribedRawExpMatrix[,1])
binaryMatrix = apply(foldMatrix,2,is.finite)
finite_rows = which(apply(binaryMatrix,1,min)>0)
boxplot(foldMatrix[finite_rows,2:7],cex=0,ylim=c(-1.8,1.3),ylab='Expression log2 fold change vs. 0hr',main ='active genes n=8,043')
abline(h=0)
dev.off()
