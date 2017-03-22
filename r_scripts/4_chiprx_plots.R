#4_chiprx_plots.R

#takes a signal table and a set of sample names
#and makes chiprx cell count normalized and background subtracted
#box plots w/ statistics

#==================================================================
#===========================DEPENDENCIES===========================
#==================================================================



#==================================================================
#======================LOADING ARGUMENTS===========================
#==================================================================


args <- commandArgs()

#getting signal table
signalFile = args[6]
print("Using signal file")
print(signalFile)

#getting scale table
scaleFile = args[7]
print("Using scale file")
print(scaleFile)

#setting sample names
name_string = args[8]
sample_names = unlist(strsplit(name_string,','))
print('sample names:')
print(sample_names)
nSamples = length(sample_names)

#setting background names
background_string = args[9]
background_names = unlist(strsplit(background_string,','))
print('background names:')
print(background_names)

#set the plot name
plotName = args[10]
print(plotName)
#set the project folder
projectFolder = args[11]
print(projectFolder)


#==================================================================
#========================HELPER FUNCTIONS==========================
#==================================================================

makeCorrectedMatrix <- function(signal_table,scale_table,sample_names,background_names){
	signal_matrix = as.matrix(signal_table[,3:ncol(signal_table)])
	rownames(signal_matrix) = signal_table[,2]
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
	
	#now background subtract
	#can't assume that columns are ordered in any particular way
	
	#create a final matrix
	final_matrix = matrix(nrow = nrow(signal_matrix),ncol=length(sample_names))
	colnames(final_matrix) = sample_names
	rownames(final_matrix) = rownames(signal_matrix)
	
	for(i in 1:length(sample_names)){
		name = sample_names[i]
		background = background_names[i]
		name_col = which(colnames(signal_matrix)==name)
		background_col = which(colnames(signal_matrix)==background)
		#final_matrix[,i] = signal_matrix[,name_col] - signal_matrix[,background_col]
		final_matrix[,i] = signal_matrix[,name_col] #background subtraction doesn't make sense for chiprx
		final_matrix[which(final_matrix[,i]<0),i] <-0
	}
	return(final_matrix)
}


makeRawMatrix <- function(signal_table,sample_names,background_names){
	signal_matrix = as.matrix(signal_table[,3:ncol(signal_table)])
	rownames(signal_matrix) = signal_table[,2]
	colnames(signal_matrix) = colnames(signal_table)[3:ncol(signal_table)]
	#now correct for scale factor
	

	#now background subtract
	#can't assume that columns are ordered in any particular way
	
	#create a final matrix
	final_matrix = matrix(nrow = nrow(signal_matrix),ncol=length(sample_names))
	colnames(final_matrix) = sample_names
	rownames(final_matrix) = rownames(signal_matrix)
	
	for(i in 1:length(sample_names)){
		name = sample_names[i]
		background = background_names[i]
		name_col = which(colnames(signal_matrix)==name)
		background_col = which(colnames(signal_matrix)==background)
		#final_matrix[,i] = signal_matrix[,name_col] - signal_matrix[,background_col]
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

#==================================================================
#================================DATA==============================
#==================================================================

#some hard coded paths during debugging.
#ctcf_intersect = read.delim('/Volumes/grail/projects/mycn_cyl/signalTables/HG19_H3K4ME3_ALL_-0_+0_signal.txt')
#scale_table = read.delim('/Volumes/grail/projects/mycn_cyl/tables/HG19_SHEP21_REVIEW_SCALE.txt')
#sample_names = c('SHEP21_0HR_H3K4ME3', 'SHEP21_2HR_H3K4ME3','SHEP21_24HR_H3K4ME3')
#background_names = c('SHEP21_0HR_INPUT', 'SHEP21_2HR_INPUT','SHEP21_24HR_INPUT')


#pull in the signal table and the scale table
signal_table = read.delim(signalFile)
scale_table = read.delim(scaleFile)





#==================================================================
#======================MAKING BOX PLOTS============================
#==================================================================

print('making uncorrected raw signal matrix')
raw_matrix = makeRawMatrix(signal_table,sample_names,background_names)


print('making corrected signal matrix')
corrected_matrix = makeCorrectedMatrix(signal_table,scale_table,sample_names,background_names)


plot_path = paste(projectFolder,'figures/4_chiprx_plots/4_',plotName,'.pdf',sep='')
print(plot_path)
pdf(file=plot_path,width = 10,height = 8)
par(mfrow=c(1,2))
plot_title = paste(plotName,'chiprx_scaled',sep=' ')
formattedBoxPlot(corrected_matrix,plot_title)
plot_title = paste(plotName,'before scaling',sep=' ')
formattedBoxPlot(raw_matrix,plot_title)
dev.off()



