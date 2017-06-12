#========================================================
#===========PLOTTING ENHANCER PROMOTER SIGNAL============
#========================================================

setwd('~/Dropbox/mycn_cyl/')
#for no spike


plotContrib <- function(geneTable){
	promoterSignal =geneTable[,2]
	enhancerSignal = geneTable[,3]
	totalSignal = promoterSignal + enhancerSignal
	
	totalOrder = order(totalSignal,decreasing=TRUE)[4:5002]
	
	enhancerContribution = enhancerSignal/(promoterSignal+enhancerSignal)
	
	
	distalLength = geneTable[,5]
	tssLength = geneTable[,4]
	
	
	enhancerVector = c()
	promoterVector = c()
	
	i = 1
	stepSize = 100
	
	while(i < (length(totalOrder) - stepSize/2)){
		
		enhancerVector = c(enhancerVector,mean(distalLength[totalOrder][i:(i+stepSize)]))
		promoterVector = c(promoterVector,mean(tssLength[totalOrder][i:(i+stepSize)]))
		i = i + stepSize/2
		
		
	}


	#enhancerContribVector = enhancerVector/(promoterVector+enhancerVector)
	yMax = max(0.8*max(enhancerVector[1:10]),0.8*max(promoterVector[1:10]))

	
	par(mfrow=c(4,1))
	plot(totalSignal[totalOrder],type='h',ylim =c(0,totalSignal[totalOrder][3]),col='red',xlab='Top 5,000 genes ranked by cumulative proximal MYCN signal',ylab='cumulative proximal MYCN signal/AUC (rpm)',main='TOP MYCN GENES',lwd=0.25)
	lines(enhancerSignal[totalOrder],type='h',col='blue',lwd=0.25)
	legend(3000,.6*max(totalSignal),c('promoter contribution','enhancer contribution'),fill=c('red','blue'))
	
	
	plot(1:length(enhancerVector), enhancerVector,type='p',col='blue',pch=16,ylab='contributing enhancer bases to total MYCN signal',xaxt='n',xlab='',ylim =c(0,7000),main='length')
	x=1:length(enhancerVector)
	lw1 = loess(enhancerVector ~x)
	lines(1:(length(enhancerVector)-1),lw1$fitted,col='blue',lwd=2)
	
	
	points(1:length(promoterVector), promoterVector,type='p',col='red',pch=16,ylab='contributing promoter bases to total MYCN signal',xaxt='n',xlab='')
	x=1:length(promoterVector)
	lw1 = loess(promoterVector ~x)
	lines(1:(length(promoterVector)-1),lw1$fitted,col='red',lwd=2)



	distalDensity = geneTable[,3]/geneTable[,5]
	distalDensity[is.na(distalDensity)] <- 0
	tssDensity = geneTable[,2]/geneTable[,4]
	tssDensity[is.na(tssDensity)] <- 0


	enhancerVector = c()
	promoterVector = c()
	
	i = 1
	stepSize = 100
	
	while(i < (length(totalOrder) - stepSize/2)){
		
		enhancerVector = c(enhancerVector,mean(distalDensity[totalOrder][i:(i+stepSize)]))
		promoterVector = c(promoterVector,mean(tssDensity[totalOrder][i:(i+stepSize)]))
		i = i + stepSize/2
		
		
	}
	yMax = max(0.8*max(enhancerVector[1:10]),0.8*max(promoterVector[1:10]))
	
	plot(1:length(enhancerVector), enhancerVector,type='p',col='blue',pch=16,ylab='MYCN density at contributing regions',xaxt='n',xlab='',ylim =c(0,yMax),main='density')
	x=1:length(enhancerVector)
	lw1 = loess(enhancerVector ~x)
	lines(1:(length(enhancerVector)-1),lw1$fitted,col='blue',lwd=2)
	
	
	points(1:length(promoterVector), promoterVector,type='p',col='red',pch=16,ylab='contributing promoter bases to total MYCN signal',xaxt='n',xlab='')
	x=1:length(promoterVector)
	lw1 = loess(promoterVector ~x)
	lines(1:(length(promoterVector)-1),lw1$fitted,col='red',lwd=2)
	
	distalAUC = geneTable[,3]

	tssAUC = geneTable[,2]
	
	
	
	enhancerVector = c()
	promoterVector = c()
	
	i = 1
	stepSize = 100
	
	while(i < (length(totalOrder) - stepSize/2)){
		
		enhancerVector = c(enhancerVector,mean(distalAUC[totalOrder][i:(i+stepSize)]))
		promoterVector = c(promoterVector,mean(tssAUC[totalOrder][i:(i+stepSize)]))
		i = i + stepSize/2
		
		
	}

	yMax = max(0.8*max(enhancerVector[1:10]),0.8*max(promoterVector[1:10]))

	
	
	plot(1:length(enhancerVector), enhancerVector,type='p',col='blue',pch=16,ylab='MYCN AUC at contributing regions (rpm)',xaxt='n',xlab='',ylim =c(0,yMax),main='AUC')
	x=1:length(enhancerVector)
	lw1 = loess(enhancerVector ~x)
	lines(1:(length(enhancerVector)-1),lw1$fitted,col='blue',lwd=2)
	
	
	points(1:length(promoterVector), promoterVector,type='p',col='red',pch=16,ylab='contributing promoter bases to total MYCN signal',xaxt='n',xlab='')
	x=1:length(promoterVector)
	lw1 = loess(promoterVector ~x)
	lines(1:(length(promoterVector)-1),lw1$fitted,col='red',lwd=2)

	

}




#========================================================
#========================DATA INPUT======================
#========================================================



geneTable = read.delim('./enhancerPromoter_resub/NB_MYCN_CONSERVED/NB_MYCN_CONSERVED_GENE_TABLE_LENGTH.txt',sep='\t')

pdf(file='./figures_resub/14_NB_MYCN_CONSERVED_ENHANCER_CONTRIBUTION.pdf',width =6,height = 12)

plotContrib(geneTable)

dev.off()





geneTable = read.delim('./enhancerPromoter_resub/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_SHEP21_0HR_MYCN_NOSPIKE/SHEP21_0HR_MYCN_NOSPIKE_REGIONS_SHEP21_0HR_MYCN_NOSPIKE_GENE_TABLE_LENGTH.txt',sep='\t')

pdf(file='./figures_resub/14_SHEP21_0HR_MYCN_NOSPIKE_ENHANCER_CONTRIBUTION.pdf',width =6,height = 12)

plotContrib(geneTable)

dev.off()
