#151209_clusterHeatmap.R


#5_chiprx_heatmaps.R

#for a given set of 3 datsets e.g. mycn at enhancers 0,2,24 or ctcf at promoters 0,2,24
#makes box plots
#heatmaps
#metas

#==================================================================
#===========================DEPENDENCIES===========================
#==================================================================



#==================================================================
#======================LOADING ARGUMENTS===========================
#==================================================================


args <- commandArgs()

#getting mapped tables
mapped_string = args[6]
mapped_path_list = unlist(strsplit(mapped_string,','))
print("Using mapped files:")
print(mapped_path_list)

#getting scale table
scaleFile = args[7]
print("Using scale file")
print(scaleFile)

#setting sample names
name_string = args[8]
sample_names = unlist(strsplit(name_string,','))
print('sample names:')
print(sample_names)

#set the plot color
plot_color= args[9]
print(plot_color)

#set the region name
region_title = args[10]
print(region_title)

#set the plot name
plot_name= args[11]
print(plot_name)


#set rpm
rpm_string = args[12]
if(rpm_string == "TRUE"){
	rpm = TRUE
}else{
	rpm= FALSE
	}
print('RPM is:')
print(rpm_string)

#set the project folder
projectFolder = args[13]
print(projectFolder)


#==================================================================
#========================HELPER FUNCTIONS==========================
#==================================================================
getScaleVector <- function(names_list,scale_table,rpm=TRUE){
	#returns a multiplicative scale factor
	#if dataset not present in scale factor then returns 1
	scale_vector = c()
	for(name in names_list){
		row_index= which(scale_table[,1] == name)
		if(length(row_index) == 0){
			scale_vector = c(scale_vector,1)
		}else{
			if(rpm==TRUE){
				scale_vector = c(scale_vector,as.numeric(scale_table[row_index,4]))
			}else{
								scale_vector = c(scale_vector,1/as.numeric(scale_table[row_index,3]))
				
			}
		}
		
	}
	return(scale_vector)
}


plotMetas <- function(mapped_table_1,mapped_table_2,mapped_table_3,scale_table,names_list,plot_color,region_title,rpm=TRUE,scale_data = TRUE){
	
	if(scale_data == TRUE){
		scale_vector = getScaleVector(names_list,scale_table,rpm)
	}else{
		scale_vector = c(1,1,1)
	}
	meta_1 = apply(as.matrix(mapped_table_1[,3:202]),2,mean)*scale_vector[1]
	meta_2 = apply(as.matrix(mapped_table_2[,3:202]),2,mean)*scale_vector[2]
	meta_3 = apply(as.matrix(mapped_table_3[,3:202]),2,mean)*scale_vector[3]
	
	#find the yMax
	yMax = 1.05*max(c(meta_1, meta_2, meta_3))
	
	#figure out correct title
	
	if(scale_data == TRUE){
		y_title = 'scaled reads/bp'
	}else{
		if(rpm == TRUE){
			y_title = 'rpm/bp'
			}else{
				y_title = 'raw reads/bp'
			}
		}
			
	par(mfrow=c(1,3))
	plot(1:200, meta_1,type='l',ylim =c(0,yMax),col=plot_color,xlab='',ylab=y_title,main = names_list[1])
	plot(1:200, meta_2,type='l',ylim =c(0,yMax),col= plot_color,xlab='',ylab=y_title,main = names_list[2])
	plot(1:200, meta_3,type='l',ylim =c(0,yMax),col= plot_color,xlab='',ylab=y_title,main = names_list[3])
}



plotHeatmaps <- function(mapped_table_1,mapped_table_2,mapped_table_3,scale_table,plot_color,plot_path,names_list,rpm=TRUE,scale_data = TRUE){
	
	if(scale_data == TRUE){
		scale_vector = getScaleVector(names_list,scale_table,rpm)
	}else{
		scale_vector = c(1,1,1)
	}
	
	if(scale_data == TRUE){
		y_title = 'scaled reads/bp'
	}else{
		if(rpm == TRUE){
			y_title = 'rpm/bp'
			}else{
				y_title = 'raw reads/bp'
			}
	}
			
	matrix_1 = as.matrix(mapped_table_1[,3:202])*scale_vector[1]	
	matrix_2 = as.matrix(mapped_table_2[,3:202])*scale_vector[2]	
	matrix_3 = as.matrix(mapped_table_3[,3:202])*scale_vector[3]	
	
	#setting up color scaling
	colorSpectrum <- colorRampPalette(c("white",plot_color))(100)
	
	minValue = min(quantile(matrix_1,ra.rm=TRUE,prob=0.6),quantile(matrix_2,ra.rm=TRUE,prob=0.6),quantile(matrix_3,ra.rm=TRUE,prob=0.6))
	print('min value is')
	print(minValue)	
	maxValue = max(quantile(matrix_1,ra.rm=TRUE,prob=0.95),quantile(matrix_2,ra.rm=TRUE,prob=0.95),quantile(matrix_3,ra.rm=TRUE,prob=0.95))
	print('max value is')
	print(maxValue)
	color_cuts <- seq(minValue,maxValue,length=100)
	
	trueMin = min(min(matrix_1),min(matrix_2),min(matrix_3))
	trueMax = max(max(matrix_1), max(matrix_2), max(matrix_3))
	
	color_cuts <- c(trueMin, color_cuts,trueMax)
	colorSpectrum <- c(colorSpectrum[1],colorSpectrum)

	referenceOrder = order(apply(matrix_1,1,mean,na.rm=TRUE))	
	
	png(filename=plot_path,width = 1500,height = 1600)
	layout(matrix(data=c(1,1,1,1,1,2,3,3,3,3,3,4,5,5,5,5,5,6),nrow=1))
	image(1:ncol(matrix_1),1:nrow(matrix_1),t(matrix_1[referenceOrder,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",xlab="",ylab="",main=names_list[1])
	image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab=y_title)
	
	image(1:ncol(matrix_2),1:nrow(matrix_2),t(matrix_2[referenceOrder,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",xlab="",ylab="",main=names_list[2])
	image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab= y_title)
	
	image(1:ncol(matrix_3),1:nrow(matrix_3),t(matrix_3[referenceOrder,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",xlab="",ylab="",main=names_list[3])
	image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab= y_title)
	dev.off()	
}

#==================================================================
#================================DATA==============================
#==================================================================

#test data w/ hard paths
#setwd('/Volumes/grail/projects/mycn_resub/mycn/')

#getting the scale vectors
scale_table = read.delim(scaleFile,header=TRUE,sep='\t')

mapped_table_1 = read.delim(mapped_path_list[1])
mapped_table_2 = read.delim(mapped_path_list[2])
mapped_table_3 = read.delim(mapped_path_list[3])

#enhancers for mycn spike rpm
#enhancer_mycn_0hr_rpm = read.delim('mappedFolder/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb_SHEP21_0HR_MYCN_RX.gff',sep='\t',header=TRUE)
#enhancer_mycn_2hr_rpm = read.delim('mappedFolder/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb_SHEP21_2HR_MYCN_RX.gff',sep='\t',header=TRUE)
#enhancer_mycn_24hr_rpm = read.delim('mappedFolder/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb_SHEP21_24HR_MYCN_RX.gff',sep='\t',header=TRUE)

#enhancers for mycn spike raw (no rpm)
#enhancer_mycn_0hr_raw = read.delim('mappedFolder_noRPM/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb_SHEP21_0HR_MYCN_RX.gff',sep='\t',header=TRUE)
#enhancer_mycn_2hr_raw = read.delim('mappedFolder_noRPM/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb_SHEP21_2HR_MYCN_RX.gff',sep='\t',header=TRUE)
#enhancer_mycn_24hr_raw = read.delim('mappedFolder_noRPM/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb/HG19_SHEP21_0HR_MYCN_NOSPIKE_CONSERVED_ENHANCER_-5kb_+5kb_SHEP21_24HR_MYCN_RX.gff',sep='\t',header=TRUE)



#==================================================================
#======================PLOTTING MEAN METAS=========================
#==================================================================

#plot a meta for both scaled and non scaled data

plot_path = paste(projectFolder,'figures/5_',region_title,'_',plot_name,'_metas.pdf',sep='')
print(plot_path)
pdf(file=plot_path,width = 8,height = 6)

plotMetas(mapped_table_1,mapped_table_2,mapped_table_3,scale_table,sample_names,plot_color,rpm,scale_data=TRUE)

plotMetas(mapped_table_1,mapped_table_2,mapped_table_3,scale_table, sample_names,plot_color,rpm,scale_data = FALSE)
dev.off()
#==================================================================
#========================MAKING HEATMAPS===========================
#==================================================================

plot_path = paste(projectFolder,'figures/5_',region_title,'_',plot_name,'_scaled_heat.png',sep='')
print(plot_path)
plotHeatmaps(mapped_table_1,mapped_table_2,mapped_table_3,scale_table,plot_color,plot_path,sample_names,rpm,scale_data = TRUE)

plot_path = paste(projectFolder,'figures/5_',region_title,'_',plot_name,'_unscaled_heat.png',sep='')
print(plot_path)
plotHeatmaps(mapped_table_1,mapped_table_2,mapped_table_3,scale_table,plot_color,plot_path,sample_names,rpm,scale_data = FALSE)


