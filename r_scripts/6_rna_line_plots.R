#6_rna_line_plots.R



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
x_vector = as.numeric(unlist(strsplit(x_string,',')))
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



makeExpMatrix <- function(gene_table,exp_table,group_names,cut=10){
	
	#first summarize the expression table into a proper mean matrix
	#with application of the cutoff
	m = matrix(nrow=nrow(exp_table),ncol = length(group_names))
	colnames(m) = group_names
	rownames(m)=rownames(exp_table)
	
	for(i in 1:length(group_names)){
		name = group_names[i]
		group_cols = grep(name,colnames(exp_table))
		m[,i] = apply(exp_table[,group_cols],1,mean)		
	}
	
	#now apply the cutoff 
	#max expression must be greater than or equal to cut
	cut_rows = intersect(which(apply(m,1,max)>=cut),which(apply(m,1,min) > 0))
	
	#now make sure each gene is actually present in the gene_table
	gene_list = as.character(gene_table[,1])
	
	final_rows = c()
	for(row in cut_rows){
		
		gene_name = as.character(rownames(m)[row])
		if(length(which(gene_list == gene_name)) > 0){
			final_rows = c(final_rows,row)
		}
		
	}
	
	final_m = m[final_rows,]
	
	return(final_m)
}


makeMeanCI <- function(m,nIter = 1000){
	#sample 1k times
	
	meanMatrix = matrix(ncol= ncol(m),nrow = nIter)
	for(i in 1:nIter){
		
		sampleRows = sample(1:nrow(m),nrow(m),replace=TRUE)
		meanMatrix[i,] = apply(m[sampleRows,],2,mean)
		
	}
	
	meanVector = apply(meanMatrix,2,mean)
	upperError = apply(meanMatrix,2,quantile,probs=0.975) - apply(meanMatrix,2,mean)
	lowerError = apply(meanMatrix,2,mean) - apply(meanMatrix,2,quantile,probs=0.025)
	
	return(list(mean=meanVector,ue = upperError,le =lowerError))
	
	
}

plotExpLines <- function(exp_matrix,gene_table,group_names,nBins,top_count,x_vector,analysis_name,plot_title,color,orderBy='total'){

	#set a few hard codies
	nIter = 1000 #number of iterations for resampling
		
		
	#first filter the gene_table to make an expressed_gene_table
	#w/ same ordering as exp_matrix
	expressed_genes = rownames(exp_matrix)
	expressed_rows = c()
	for(gene in expressed_genes){		
		gene_row = which(as.character(gene_table[,1])==gene)			
		expressed_rows = c(expressed_rows,gene_row)
	}
	
	#this should have same dimensionality and order as the exp_matrix
	expressed_gene_table = gene_table[expressed_rows,] 
	#next make a fold matrix
	fold_matrix = log2(exp_matrix/exp_matrix[,1])

	#next rank the gene_table by cumulative signal

	if(orderBy == 'promoter'){
		sig = as.numeric(expressed_gene_table[,2])		
	}else if(orderBy == 'enhancer'){
		sig = as.numeric(expressed_gene_table[,3])	
	}else{
		sig = as.numeric(expressed_gene_table[,2]) + as.numeric(expressed_gene_table[,3])
	}	
	
	
	sig_order = order(sig,decreasing=TRUE)
	expressed_gene_table_ordered = expressed_gene_table[sig_order,]	
	fold_matrix_ordered = fold_matrix[sig_order,]
	sig_ordered = sig[sig_order]
	
	#for each bin we need a mean vector w/ error bars
	#set the final matricies
	mean_matrix = matrix(nrow=nBins,ncol = length(group_names))
	ue_matrix = matrix(nrow=nBins,ncol = length(group_names))
	le_matrix = matrix(nrow=nBins,ncol = length(group_names))
	
	
	if(as.character(top_count) == 'all'){
		top_count = nrow(fold_matrix_ordered)
		
	}

	#important checkpoint
	if(nrow(fold_matrix_ordered) < top_count){
		print('ERROR not enough genes for analysis')
		print(nrow(fold_matrix_ordered))
		print(length(top_count))
	}	
	#set up binning
	binSize = top_count/nBins
	n= 1
	start = 1
	stop = binSize
	
	while(n <= nBins){
		
		bin_matrix = fold_matrix_ordered[start:stop,]
		bin_ci = makeMeanCI(bin_matrix,nIter)
		
		mean_matrix[n,] = bin_ci$mean
		ue_matrix[n,] = bin_ci$ue
		le_matrix[n,] = bin_ci$le
		
		start = start + binSize
		stop =  stop + binSize
		n = n+ 1
		
	}
	color_spectrum = colorRampPalette(c(color,rgb(0.75,0.75,0.75,maxColorValue=1)))(nBins)
	
	#now set appropriate y axis limits
	yMax = (max(mean_matrix) + max(ue_matrix)) + abs(mean(mean_matrix))
	yMin =  min(mean_matrix) - max(le_matrix)- abs(mean(mean_matrix))
	
	plot(x_vector,x_vector,cex=0,ylim = c(yMin,yMax),xlab='',xaxt='n',ylab='log2 fold change gene expression',main=plot_title)
	for(i in 1:nrow(mean_matrix)){
		lines(x_vector,mean_matrix[i,],type='b',col= color_spectrum[i],cex=0,lwd=2)
		error.bar(x_vector,mean_matrix[i,],ue_matrix[i,],le_matrix[i,],length=0.05,lwd=2,col= color_spectrum[i])
	}
	axis(1,x_vector,group_names,las=2)
	
	#legend(x_vector[2],mean(yMin+yMax),group_names,col=color_spectrum)
	
}






#==================================================================
#==================SAMPLE PARAMTERS FOR DEBUGGING==================
#==================================================================

# #get a gene table from enhancer promoter analysis
# gene_table_path = './enhancerPromoter/BE2C_MYCN/BE2C_MYCN_GENE_TABLE.txt'

# #get the cell count norm exp path
# exp_norm_path = './expression/BE2C_shTWIST_cuffnorm/output/BE2C_shTWIST_all_fpkm_exprs_norm.txt'

# #get the raw exp path
# exp_raw_path = './expression/BE2C_shTWIST_cuffnorm/output/BE2C_shTWIST_all_fpkm_exprs_raw.txt'

# #setting group names
# group_string = 'BE2C_shT_nodox,BE2C_shT_3HR,BE2C_shT_6HR,BE2C_shT_12HR,BE2C_shT_24HR,BE2C_shT_48HR'
# group_names = unlist(strsplit(group_string,','))
# print('group names:')
# print(group_names)

# #setting the top number of genes to run on
# top_count = 5000
# print('running analysis on top genes:')
# print(top_count)

# nBins = 3
# print('Dividing genes into N bins:')
# print(nBins)

# #set the analysis name
# analysis_name = 'BE2C_shTWIST'
# print('analysis name:')
# print(analysis_name)

# x_string = '0,3,6,12,24,48'
# x_vector = as.numeric(unlist(strsplit(x_string,',')))
# print('x axis labels')
# print(x_vector)

# #set the color string
# color_string = '175,0,0'
# color_vector = as.numeric(unlist(strsplit(color_string,',')))
# #color = rgb(color_vector[1],color_vector[2],color_vector[3],maxColorValue = 255)
# color=rgb(0,0,0.75,maxColorValue=1)
# print('color is:')
# print(color_vector)


# #set the project folder
# projectFolder = './'
# print(projectFolder)

#==================================================================
#========================DATA INPUT================================
#==================================================================

#loading the gene table
gene_table = read.delim(gene_table_path,sep='\t')

#loading the normalized expression table
exp_norm_table = read.delim(exp_norm_path,sep='\t')

#loading the raw expression table
exp_raw_table = read.delim(exp_raw_path,sep='\t')

#==================================================================
#===================MAKING FOLD EXP MATRIX=========================
#==================================================================

print('making norm exp matrix')
exp_norm_matrix = makeExpMatrix(gene_table,exp_norm_table,group_names,cut=10)

print('making raw exp matrix')
exp_raw_matrix = makeExpMatrix(gene_table,exp_raw_table,group_names,cut=10)



#==================================================================
#=================PLOTTING RAW AND NORM DATA=======================
#==================================================================


pdf_path = paste(projectFolder,'figures/6_mRNA_',analysis_name,'_n',as.character(nBins),'_top_',as.character(top_count),'.pdf',sep='')

pdf(file=pdf_path,width = 8,height =6)
#for the normalized by promoter,enhancer,total
plot_title = paste(analysis_name,'_norm_promoter',sep='')
print('plotting:')
print(plot_title)
plotExpLines(exp_norm_matrix,gene_table,group_names,nBins,top_count,x_vector,analysis_name,plot_title,color,orderBy ='promoter')

plot_title = paste(analysis_name,'_norm_enhancer',sep='')
print('plotting:')
print(plot_title)
plotExpLines(exp_norm_matrix,gene_table,group_names,nBins,top_count,x_vector,analysis_name,plot_title,color,orderBy ='enhancer')

plot_title = paste(analysis_name,'_norm_total',sep='')
print('plotting:')
print(plot_title)
plotExpLines(exp_norm_matrix,gene_table,group_names,nBins,top_count,x_vector,analysis_name,plot_title,color,orderBy ='total')

#for the raw by promoter,enhancer,total
plot_title = paste(analysis_name,'_raw_promoter',sep='')
print('plotting:')
print(plot_title)
plotExpLines(exp_raw_matrix,gene_table,group_names,nBins,top_count,x_vector,analysis_name,plot_title,color,orderBy ='promoter')

plot_title = paste(analysis_name,'_raw_enhancer',sep='')
print('plotting:')
print(plot_title)
plotExpLines(exp_raw_matrix,gene_table,group_names,nBins,top_count,x_vector,analysis_name,plot_title,color,orderBy ='enhancer')

plot_title = paste(analysis_name,'_raw_total',sep='')
print('plotting:')
print(plot_title)
plotExpLines(exp_raw_matrix,gene_table,group_names,nBins,top_count,x_vector,analysis_name,plot_title,color,orderBy ='total')

dev.off()

