#8_nes_heatmaps.R




#==================================================================
#===========================DEPENDENCIES===========================
#==================================================================

#==================================================================
#======================LOADING ARGUMENTS===========================
#==================================================================


args <- commandArgs()

# get the path to the table
nes_path = args[6]

#get the fdr cut
fdr_cut = as.numeric(args[7])

#get the enrichment score cutoff
enrich_cut = as.numeric(args[8])

#set the project folder
projectFolder = args[9]
print(projectFolder)


#==================================================================
#========================HELPER FUNCTIONS==========================
#==================================================================

plot_nes_heatmap <- function(nes_table,name,projectFolder,fdr_cut,enrich_cut){

	#first figure out the number of samples
	nSamples = (ncol(nes_table)-1)/2

	#next set up the path for the output figure
	figure_path = paste(projectFolder,'figures/8_',name,'_fdr_',fdr_cut,'_nes_',enrich_cut,'.pdf',sep='')
	print(figure_path)


	#next set up a matrix w/ just the fdr scores since we are doing binary encoding of the heatmap
	fdr_cols = seq(1,nSamples*2,2)+2
	nes_fdr = nes_table[,fdr_cols]
	rownames(nes_fdr)=as.character(nes_table[,1])
	pos_rows = c()
	neg_rows = c()
	for(i in 1:nSamples){
	
		cols = c(i*2,i*2+1)
		for(row in 1:nrow(nes_table)){
			if(nes_table[row,cols[1]] > enrich_cut && nes_table[row,cols[2]] < fdr_cut){
				pos_rows = c(pos_rows,row)
			}
		
		if(nes_table[row,cols[1]] < -1 * enrich_cut && nes_table[row,cols[2]] < fdr_cut){	
			neg_rows = c(neg_rows,row)
			}		
		}
	}

	neg_rows = unique(neg_rows)
	pos_rows = unique(pos_rows)

	#just do binary encoding
	nes_fdr_binary = nes_fdr
	for(i in 1:nrow(nes_fdr_binary)){
	
		for(j in 1:ncol(nes_fdr_binary)){
		      if(nes_fdr_binary[i,j] <= fdr_cut){
		      		nes_fdr_binary[i,j] = 0
			}else{
				nes_fdr_binary[i,j] = 1

			}
		}
	}

	#now for positive
	pos_binary= as.matrix(nes_fdr_binary[pos_rows,])
	pos_hclust = hclust(dist(pos_binary))
	pos_sample_hclust = hclust(dist(t(pos_binary)))
	pos_binary[pos_hclust$order, pos_sample_hclust$order]

	neg_binary= as.matrix(nes_fdr_binary[neg_rows,])
	neg_hclust = hclust(dist(neg_binary))
	neg_sample_hclust = hclust(dist(t(neg_binary)))
	neg_binary[neg_hclust$order, neg_sample_hclust$order]

	colorSpectrum <- c('red','white')
	color_cuts = c(0,0,1)

	pdf(file=figure_path,width =6,height = 10)

	colorSpectrum <- c('blue','white')

	image(1:ncol(neg_binary),1:nrow(neg_binary),t(neg_binary[neg_hclust$order,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='',main='enhancer enriched')
	axis(2,1:nrow(neg_binary),rownames(neg_binary)[neg_hclust$order],las=2)
	axis(1,1:ncol(neg_binary),colnames(neg_binary),las=2)

	colorSpectrum <- c('red','white')

	image(1:ncol(pos_binary),1:nrow(pos_binary),t(pos_binary[pos_hclust$order,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='',main='promoter enriched')
	axis(2,1:nrow(pos_binary),rownames(pos_binary)[pos_hclust$order],las=2)
	axis(1,1:ncol(pos_binary),colnames(pos_binary),las=2)
	dev.off()
}


#==================================================================
#==============================DATA_INPUT==========================
#==================================================================

#load in the nes table
nes_table = read.delim(nes_path,sep='\t')

#get the root name
name = unlist(strsplit(nes_path,'/'))
name = name[length(name)]

name = gsub('.txt','',name)

print(name)

plot_nes_heatmap(nes_table,name,projectFolder,fdr_cut,enrich_cut)


#==================================================================
#==================SAMPLE PARAMTERS FOR DEBUGGING==================
#==================================================================

#setwd('/Volumes/grail/projects/mycn_resub/mycn/')
#nes_path = './nes_tables/MYC_HIGH_NES.txt'
#nes_table = read.delim(nes_path,sep='\t')


#sigCut = 0.01
#enrichCut = 2


#go through pairs and
#nSamples = 5






