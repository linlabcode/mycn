#1_nb_mycn_rank.R

#takes the Meta rose output and identifies conserved regions using a rank cutoff

#==================================================================
#===========================DEPENDENCIES===========================
#==================================================================

library(ggplot2)


#==================================================================
#======================LOADING ARGUMENTS===========================
#==================================================================


args <- commandArgs()

signalFile = args[6]
print("Using signal file")
print(signalFile)

name_string = args[7]
print(name_string)
names_list = unlist(strsplit(name_string,','))
print(names_list)
nSamples = length(names_list)

projectFolder = args[8]

#==================================================================
#================================DATA==============================
#==================================================================

#input signal table and perform background subtraction
signalTable = read.table(signalFile,sep='\t',header=TRUE)

#formatting the signalMatrix
signalMatrix= as.matrix(signalTable[,7:10]) - as.matrix(signalTable[,11:14]) #hard coded to assume 4 samples
colnames(signalMatrix) = names_list
rownames(signalMatrix) = as.character(signalTable[,1])


rankTable = cbind(signalTable[,1:6],apply(signalMatrix,2,rank))

#input table of ranks
#high rank is good rank
rankFile = gsub('_MAP.txt','_RANK.txt',signalFile)
print('Writing rank table to:')
print(rankFile)
write.table(rankTable,file=rankFile,sep='\t',quote=FALSE,row.names=FALSE)

rankMatrix = as.matrix(rankTable[,7:ncol(rankTable)])
rownames(rankMatrix)=as.character(rankTable[,1])
colnames(rankMatrix) = names_list

mycnRankMatrix = rankMatrix

rankOrder = order(apply(rankMatrix,1,mean),decreasing=TRUE)

# #some simple diagnostic plotting that can be used to determine cutoff
# par(mfrow=c(1,4))
# plot(1:nrow(rankMatrix),rankMatrix[rankOrder,1],cex=0.5,pch=16,col=rgb(0.5,.5,.5,.5))
# abline(h=15000)
# abline(v=15000)
# plot(rankMatrix[rankOrder,2],cex=0.5,pch=16,col=rgb(0.5,.5,.5,.5))
# abline(h=15000)
# abline(v=15000)
# plot(rankMatrix[rankOrder,3],cex=0.5,pch=16,col=rgb(0.5,.5,.5,.5))
# abline(h=15000)
# abline(v=15000)
# plot(rankMatrix[rankOrder,4],cex=0.5,pch=16,col=rgb(0.5,.5,.5,.5))
# abline(h=15000)
# abline(v=15000)

#try a 15k cutoff

conservedRows = c()
cutOff = nrow(rankMatrix)-15000
for(i in 1:nrow(rankMatrix)){
	
	if(length(which(rankMatrix[i,] < cutOff)) <= 1){
		conservedRows = c(conservedRows,i)
		
	}
}

conservedRankMatrix = rankMatrix[conservedRows,]
conservedRankOrder = order(apply(conservedRankMatrix,1,mean),decreasing=TRUE)

#additional diagnostic plotting
# par(mfrow=c(1,4))
# plot(conservedRankMatrix[conservedRankOrder,1],cex=0.5,pch=16,col=rgb(0.5,.5,.5,.5))
# plot(conservedRankMatrix[conservedRankOrder,2],cex=0.5,pch=16,col=rgb(0.5,.5,.5,.5))
# plot(conservedRankMatrix[conservedRankOrder,3],cex=0.5,pch=16,col=rgb(0.5,.5,.5,.5))
# plot(conservedRankMatrix[conservedRankOrder,4],cex=0.5,pch=16,col=rgb(0.5,.5,.5,.5))

conservedSignalTable = signalTable[conservedRows,]
conservedRankTable = rankTable[conservedRows,]

conservedRankFile = gsub('_MAP.txt','_RANK_CONSERVED.txt',signalFile)
print('Writing conserved rank table to:')
print(conservedRankFile)
write.table(conservedRankTable,file=conservedRankFile,sep='\t',quote=FALSE,row.names=FALSE)

conservedSignalFile = gsub('_MAP.txt','_MAP_CONSERVED.txt',signalFile)
print('Writing conserved signal table to:')
print(conservedSignalFile)
write.table(conservedSignalTable,file=conservedSignalFile,sep='\t',quote=FALSE,row.names=FALSE)




#==================================================================
#=========================DENSITY PLOTS============================
#==================================================================

#for some reason ggplot doesn't like looping

#for shep21
print('Plotting rank plot for shep21')
pdf(file='./figures/1_rank_shep21_contour.pdf',width =6.5,height =5)
dataset = structure(list(x=1:length(rankOrder),y=rankMatrix[rankOrder,2]),class = "data.frame")

ggplot(dataset, aes(x, y)) + 
    geom_point(data = dataset,cex=0.5,pch=16,col=rgb(0.5,0.5,0.5,0.5)) +
    stat_density2d(aes(alpha=..level.., fill=..level..), size=2, bins=50, geom="polygon") +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_alpha(range = c(0.00, 0.5), guide = FALSE) +
    geom_density2d(colour=rgb(0.5,0.5,0.5,.8))

dev.off()    


#for be2c
print('Plotting rank plot for be2c')
pdf(file='./figures/1_rank_be2c_contour.pdf',width =6.5,height =5)
dataset = structure(list(x=1:length(rankOrder),y=rankMatrix[rankOrder,2]),class = "data.frame")

ggplot(dataset, aes(x, y)) + 
    geom_point(data = dataset,cex=0.5,pch=16,col=rgb(0.5,0.5,0.5,0.5)) +
    stat_density2d(aes(alpha=..level.., fill=..level..), size=2, bins=50, geom="polygon") +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_alpha(range = c(0.00, 0.5), guide = FALSE) +
    geom_density2d(colour=rgb(0.5,0.5,0.5,.8))

dev.off()    

#for kelly
print('Plotting rank plot for kelly')
pdf(file='./figures/1_rank_kelly_contour.pdf',width =6.5,height =5)
dataset = structure(list(x=1:length(rankOrder),y=rankMatrix[rankOrder,3]),class = "data.frame")

ggplot(dataset, aes(x, y)) + 
    geom_point(data = dataset,cex=0.5,pch=16,col=rgb(0.5,0.5,0.5,0.5)) +
    stat_density2d(aes(alpha=..level.., fill=..level..), size=2, bins=50, geom="polygon") +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_alpha(range = c(0.00, 0.5), guide = FALSE) +
    geom_density2d(colour=rgb(0.5,0.5,0.5,.8))

dev.off()    



#for ngp
print('Plotting rank plot for ngp')
pdf(file='./figures/1_rank_ngp_contour.pdf',width =6.5,height =5)
dataset = structure(list(x=1:length(rankOrder),y=rankMatrix[rankOrder,4]),class = "data.frame")

ggplot(dataset, aes(x, y)) + 
    geom_point(data = dataset,cex=0.5,pch=16,col=rgb(0.5,0.5,0.5,0.5)) +
    stat_density2d(aes(alpha=..level.., fill=..level..), size=2, bins=100, geom="polygon") +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_alpha(range = c(0.00, 0.5), guide = FALSE) +
    geom_density2d(colour=rgb(0.5,0.5,0.5,.8))

dev.off()    

print('all done w/ rank script')