#2_nb_mycn_vector_plots.R
#comparing MYCN rank to various genomic features

#==================================================================
#======================LOADING ARGUMENTS===========================
#==================================================================


args <- commandArgs()
ebox_rank_path = args[6]
print("Using EBOX RANK TABLE")
print(ebox_rank_path)

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


makeHeatmap <- function(m){


	colorSpectrum <- colorRampPalette(c("white","white","red"))(100)

	#setting a color data range
	minValue <- 0
	maxValue <- 1
	color_cuts <- seq(minValue,maxValue,length=100)

	#add one extra max color to even out sampling
	color_cuts <- c(color_cuts,1)

	layout(matrix(data=c(1,1,2),ncol=3))
	image(1:ncol(m),1:nrow(m),t(rev(m)),breaks=color_cuts,col=colorSpectrum,yaxt='n',xaxt='n',xlab='',ylab='')
	axis(2,1:nrow(m),rev(rownames(m)),las=2)

	image(1:2,color_cuts[1:100],t(matrix(data=color_cuts[1:100],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Relative EBOX Strength")	

	}





#========================================================
#========================DATA INPUT======================
#========================================================

ebox_rank_table = read.delim(ebox_rank_path,sep='\t',header=TRUE)
print(ebox_rank_table[1:5,])

#now turn into a matrix
m = as.matrix(ebox_rank_table[,3])
rownames(m) = as.character(ebox_rank_table[,1])
print(m)
m[,1] = m[,1]/max(m[,1])
print(m)


#========================================================
#===================PLOTTING COMPARISONS=================
#========================================================

plot_path = paste(projectFolder,'figures/3_ebox_heatmap.pdf',sep='')
print(plot_path)


pdf(file=plot_path,width = 3,height =5)
makeHeatmap(m)
dev.off()
