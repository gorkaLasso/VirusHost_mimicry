library(gplots)
#library(RColorBrewer)
#############################
#    READING THE DATA
#############################

inputFile<-"matrix_molfun_0.01.txt"
heatmapPNG<-"matrix_molfun_0.01_heatmap.png"
orderedFile<-"matrix_molfun_0.01_orderedRows.txt"
clusteredFile<-"matrix_molfun_0.01_clustered.txt"

da<-read.table(inputFile,head=T,sep='\t')
rownames(da)<-da$term

#i = 0
#for (name in colnames(da)){
#  i=i+1
#  if (i>=2){
#    colnames(da)[i]=rownames(da)[i-1]
#  }
#}


#############################
#    CUSTOMIZING COLORS FOR HEATMAP
#############################

my_palette <- colorRampPalette(c("navyblue","white","red3"))(n=299)

col_breaks = c(seq(0,1.2,length=100),
  seq(1.4,2,length=100),
  seq(2.1,10,length=100))

#############################
#    PNG FOR HEATMAP
#############################

png(heatmapPNG,
  width = 25*500,
  height = 25*500,
  res = 500,
  pointsize = 10
)

#############################
#    HEATMAP
#############################

hm <- heatmap.2(data.matrix(da)[,-1],
  trace='none',
  revC='True',
  col=my_palette,
  symm=F,
  main="",
  breaks=col_breaks,
  margins=c(25,120),
  cexRow=1.0,
  cexCol=2.0,
  key='True',
  density.info="none",
 #lhei=c(2,4),lwid=c(2,3.5),
 keysize=0.5, key.par = list(cex=1),
  hclustfun=function(x) hclust(x,method="complete"),#this is the default hclustfun
  #distfun=function(x) as.dist(1-cor(data.matrix(da)[,-1])),#use pearson to compute similarity matrix
  distfun = function(x) dist(x,method = 'euclidean'),#this is the default distfun using euclidean
)

#save in a text file the row names in the same order as in the heatmap
rn <- rownames(da)[hm$rowInd]  
write.table(rn,file=orderedFile, sep="\t");

#save sorted matrix
write.table(data.matrix(hm$carpet),file= clusteredFile, row.names=TRUE, col.names=TRUE, sep='\t')


#############################
#    WRITING HEATMAP
#############################

dev.off()