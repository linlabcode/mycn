#2_nb_mycn_vector_plots.R
#comparing MYCN rank to various genomic features

#==================================================================
#======================LOADING ARGUMENTS===========================
#==================================================================


args <- commandArgs()
mycn_table_path = args[6]
print("Using mycn_table")
print(mycn_table_path)

projectFolder = args[7]





#=========================================================
#========================HELPER FUNCTIONS=================
#=========================================================

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
                     function(x)
                       rgb(x[1], x[2], x[3], alpha=alpha))
}


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
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


#========================================================
#========================DATA INPUT======================
#========================================================

mycn_table = read.delim(mycn_table_path,sep='\t',header=TRUE)
#print(mycn_table[1:5,])


#========================================================
#===================PLOTTING COMPARISONS=================
#========================================================

plot_path = paste(projectFolder,'figures/2_nb_mycn_comparisons.pdf',sep='')
print(plot_path)


pdf(file=plot_path,width=10,height =8)

#mycn and h3k27ac signal
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$AVG_MYCN_SIGNAL*mycn_table$LENGTH
name1 = 'MYCN RANK'
name2 = 'AVERAGE MYCN SIGNAL'
title = 'AVERAGE MYCN SIGNAL as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'red',nBins = 50,nIter = 100,yMin = 0,yMax=2500)


v1 = mycn_table$MYCN_RANK
v2 = mycn_table$AVG_H3K27AC_SIGNAL*mycn_table$LENGTH
name1 = 'MYCN RANK'
name2 = 'AVERAGE H3K27AC SIGNAL'
title = 'AVERAGE H3K27AC SIGNAL as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'blue',nBins = 50,nIter = 100,yMin = 0,yMax=2500)


#tss and enhancers
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$ACTIVE_TSS_OVERLAP*100
name1 = 'MYCN RANK'
name2 = 'ACTIVE TSS OVERLAP'
title = 'ACTIVE TSS overlap as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'red',nBins = 50,nIter = 100,yMin = 0,yMax=100)

v1 = mycn_table$MYCN_RANK
v2 = mycn_table$ENHANCER_OVERLAP*100
name1 = 'MYCN RANK'
name2 = 'ACTIVE ENHANCER OVERLAP'
title = 'ACTIVE ENHANCER overlap as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'blue',nBins = 50,nIter = 100,yMin = 0,yMax=100)


#canonical vs. non canonical eboxes
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$CANON_EBOX_COUNT/mycn_table$LENGTH * 1000
name1 = 'MYCN RANK'
name2 = 'CANON EBOX PER KB'
title = 'CANONICAL EBOX DENSITY as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'red',nBins = 50,nIter = 100,yMin = 0,yMax=1.2)


v1 = mycn_table$MYCN_RANK
v2 = mycn_table$NONCANON_EBOX_COUNT/mycn_table$LENGTH * 1000
name1 = 'MYCN RANK'
name2 = 'NON CANON EBOX PER KB'

title = 'NONCANONICAL EBOX DENSITY as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'black',nBins = 50,nIter = 100,yMin = 3,yMax=7)

#for GC
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$GC_FREQ
name1 = 'MYCN RANK'
name2 = 'GC fraction'

title = 'GC FRACTION as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'black',nBins = 50,nIter = 100,yMin = 0,yMax=0.3)


#for GABPA
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$GABPA_COUNT/mycn_table$LENGTH * 1000
name1 = 'MYCN RANK'
name2 = 'GABPA MOTIFS PER KB'

title = 'GABPA MOTIF DENSITY as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'red',nBins = 50,nIter = 100,yMin = 0,yMax=1.2)

#for GATA
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$GATA_COUNT/mycn_table$LENGTH * 1000
name1 = 'MYCN RANK'
name2 = 'GATA MOTIFS PER KB'

title = 'GATA MOTIF DENSITY as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'blue',nBins = 50,nIter = 100,yMin = 0,yMax=2)



#for CANON enrichment
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$CANON_EBOX_COUNT/(mycn_table$CANON_EXP+0.01)
name1 = 'MYCN RANK'
name2 = 'CANONICAL E-BOX FOLD ENRICHMENT'

title = 'CANONICAL E-Box enrichment as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'red',nBins = 50,nIter = 100,yMin = 0,yMax=20)


#for NONCANON enrichment
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$NONCANON_EBOX_COUNT/(mycn_table$NON_CANON_EXP+0.01)
name1 = 'MYCN RANK'
name2 = 'NON-CANONICAL E-BOX FOLD ENRICHMENT'

title = 'NON-CANONICAL E-Box enrichment as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'black',nBins = 50,nIter = 100,yMin = 0,yMax=5)


#for GABP enrichment
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$GABPA_COUNT/(mycn_table$GABPA_EXP+0.01)
name1 = 'MYCN RANK'
name2 = 'GABPA FOLD ENRICHMENT'

title = 'GABPA motif enrichment as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'red',nBins = 50,nIter = 100,yMin = 0,yMax=2)


#for GATA enrichment
v1 = mycn_table$MYCN_RANK
v2 = mycn_table$GATA_COUNT/(mycn_table$GATA_EXP+0.01)
name1 = 'MYCN RANK'
name2 = 'GATA FOLD ENRICHMENT'

title = 'GATA motif enrichment as a function of MYCN RANK at MYCN bound regions'
compareVectors(v1,v2,name1,name2,title,'blue',nBins = 50,nIter = 100,yMin = 0,yMax=2)



dev.off()
