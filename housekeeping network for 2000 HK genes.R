setwd("C:/Users/edoardo/Desktop/seattle experiments/joseph network")
#load the list of HK genes
l<-read.table("HK_genes.txt")
HKGenes<-as.character(l$V1)
#load the list of genes in the brca df
setwd("C:/Users/edoardo/Desktop/manuscript taramelli")
dataDir = "C:/Users/edoardo/Desktop/manuscript taramelli/Data/TCGA_Express_Data"
cancerID<-"brca"
inFileName = paste(cancerID,"_genes.tsv",sep="")
inFileName = paste(dataDir,inFileName,sep="/");
dfGenes    = read.csv(inFileName, sep="\t",header=F)

setwd("C:/Users/edoardo/Desktop/manuscript taramelli/KLD score")
#HK in the df
presenti<-HKGenes[!is.na(match(HKGenes,as.character(dfGenes$V1)))]
#HK not in the df
assenti<-HKGenes[is.na(match(HKGenes,as.character(dfGenes$V1)))]

set.seed(1)
sample<-sample(presenti,2000)
#print a list to give to string
write.csv(sample,"C:/Users/edoardo/Desktop/manuscript taramelli/sample 2000 HK genes.csv")
#remove the header and the first column
#give everithing to string
##http://string-db.org/cgi/input.pl?UserId=7RvKhUPTpGMf&sessionId=4FP98jX5X01K&input_page_active_form=multiple_identifiers

#filter only for experimental validated interaction
#save the file as a simple tabular text output
network<-read.csv("string_interactions sample 2000HK.csv")

col1<-as.character(network$X.node1)
col2<-as.character(network$node2)

#HK in col1
col1[!is.na(match(col1,HKGenes))]
#HK not in  col1
col1[is.na(match(col1,HKGenes))]

#HK in col2
col2[!is.na(match(col2,HKGenes))]
#HK not in  col1
col2[is.na(match(col2,HKGenes))]

#index of HK in col1
index1<-which(!is.na(match(col1,HKGenes)))
#index of HK in col2
index2<-which(!is.na(match(col2,HKGenes)))
#identify which pair is not HK-HK
r<-which(is.na(match(index1,index2)))
x<-which(c(T,T)==F)
if(!identical(r,x)){
  #the common index in both the column
  indexHKpair<-index1[-r]
  
  #who is missing
  index1[r]
  #chi sono i geni nell'index 1 (dovrebbero essere HK in quanto sono in index1 ma sono NA dal confronto con index2)
  col1[index1[r]]
  #prendo i nomi di chi sono i geni negli index2
  col2[index1[r]]
  #werifico se effettivamente sono coppie nel file network
  network[index1[r],]
  #confermo definitivamente che quelli in col2 non sono HK e quelli in index1 sono HK
  match(col2[index1[r]],HKGenes)
  match(col1[index1[r]],HKGenes)
} else {
  sum(index1!=index2)
  indexHKpair<-index1
}

#mi aspetto che la lunghezza di indexHKpair sia uguale a quella di index1 - la lunghezza di r
length(indexHKpair)
length(index1)-length(r)
#infine: tutti i geni in col1 non in indexHKpair non sono HK o hanno un partner non HK in col2
#devo lavorare con gli index specifici in quanto i geni sono ripetuti, solo le coppie sono sbagliate
index_r<-(1:length(col1))[-indexHKpair]
#determino in col1 gli HK rimossi
c1_HK_r<-(col1[index_r])[!is.na(match(col1[index_r],HKGenes))]
index_c1_HK_r<-index_r[!is.na(match(col1[index_r],HKGenes))]
#in col2 erano HK
match(col2[index_c1_HK_r],HKGenes)
#determino in col2 gli HK rimossi
c2_HK_r<-(col2[index_r])[!is.na(match(col2[index_r],HKGenes))]
index_c2_HK_r<-index_r[!is.na(match(col2[index_r],HKGenes))]
#in col1 erano HK
match(col1[index_c2_HK_r],HKGenes)
#######################################################################
#verifico che in col 1 indicizzata per indexHKpair siano tutti HK....
sum(is.na(match(col1[indexHKpair],HKGenes)))
#...ma soprattutto verifico che in col 2 indicizzata per indexHKpair siano tutti HK
sum(is.na(match(col2[indexHKpair],HKGenes)))

n<-network[indexHKpair,]
n[n$experimentally_determined_interaction==0.999,]

###############
###############
###############
dfGenes$V1[dfGenes$V1=="XIAP"]

nrow(network)
nrow(n)
write.csv(n,"C:/Users/edoardo/Desktop/manuscript taramelli/network filtrata per soli HK.csv")
