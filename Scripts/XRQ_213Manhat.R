#########################################################
### Burrito {Manhattanizer}
#########################################################
#
# Author: Sariel Hubner, sariel.hubner@botany.ubc.ca
# Ammended by: Rishi R. Masalia, rishimasalia@gmail.com, Jan 2016 - April 2017, Burke Lab, UGA
# Please don't distribute
#
# Summary: Burrito is a wrapper for the Manhattanizer (augmentation of Sariel's original Helimap pipeline)
# It treats traits and environments individually, and produces a Manhattan plot & list of significant SNPs.
# EMMAX is run with kinship and structure (PCA) to get associations
#
# Preferences are controlled in the configuration file ("Asada.config")
#
#########################################################
### Installing and loading required packages
#########################################################

if (!require("qqman")) {
  install.packages("qqman", dependencies = TRUE)
  library(qqman)
}
if (!require("CMplot")) {
  install.packages("CMplot", dependencies = TRUE)
  library(CMplot)
}

if (!require("SNPRelate")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("SNPRelate")
  library(SNPRelate)
}

if (!require("lme4")) {
  install.packages("lme4", dependencies = TRUE)
  library(lme4)
}

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

if (!require("metafor")) {
  install.packages("metafor", dependencies = TRUE)
  library(metafor)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

library(tools)

#########################################################
### Read and prepare files
#########################################################

setwd("/home/burke_lab/Desktop/XRQ_Nov2017_213/EnchiladaSuite/Burrito/") ##### change this to match your path to Enchilada

##########################################################33
config = read.table(file="Asada.config")

Variant_caller = config[1,1]
GROUP = as.character(config[2,1])
vcf = config[3,1]
trait = config[4,1]
time = config[5,1]
num_env = config[6,1]
num_rep = config[7,1]
prefix = config[8,1]
envNames = c(as.character(config[9,1]),as.character(config[11,1]),as.character(config[13,1]),as.character(config[15,1]),as.character(config[17,1]))
envFiles = c(as.character(config[10,1]),as.character(config[12,1]),as.character(config[14,1]),as.character(config[16,1]),as.character(config[18,1]))
kinship = config[19,1]
pca = config[20,1]
fdrcutoff = as.character(config[21,1])
color1 = as.character(config[22,1])
color2 = as.character(config[23,1])
comparison_name = envNames[1] 

#########read phenotypes and split to have data frame for each env and rep:

op <- options(stringsAsFactors=FALSE)
sam<-read.table(paste(prefix,".tfam",sep=""))
SAM<-sam[,1]
options(op)

#########################################################
### Association
#########################################################
num_env<-as.numeric(as.character(num_env))
num_rep<-as.numeric(as.character(num_rep))
print(num_env)
print(num_rep)

for(i in 1:num_env) {
cat ("\n..................Preparing the Pheno Files.............\n")
Pheno_in<-read.table(envFiles[i],header=TRUE)
Pheno_in[Pheno_in=="Na"]<-"NA"
PhenFile<-merge(SAM, Pheno_in, by.x=1,by.y=1,all.x=TRUE)
for(j in 1:num_rep) {
	options(op)
	write.table(data.frame(PhenFile[,1],PhenFile[,1], PhenFile[,j+1]),paste(trait,"_",envNames[i],".pheno",sep=""),row.names = FALSE,col.names = FALSE, quote = FALSE)
	EMMAX.env1<-paste("./bin/emmax-intel64 -v -d 10 -t ",prefix, " -p ", paste(trait,"_",envNames[i],".pheno",sep=""), " -c ", paste(prefix,".PCA_EV",sep=""), " -k ", paste(prefix,".aIBS.kinf",sep=""), " -o ", paste(trait,"_",envNames[i],sep=""), sep="")
	system (EMMAX.env1)
	}
}

cat ("Reading the map.......................")
snp.map = read.table(file="XRQ_213.map",header= T) 

numLoci = 640508 #New values with April2017
InferrMeff = 88009 #Generated from Lexuan using Beagle 
adj = (0.05/88009) #LOD6
cat("Total number of SNPs is: ", numLoci, "\n")
cat("Inferred Meff is: ", InferrMeff, "\n")
#adj<-0.05/sum(simpleMeff)
adjP<-format(adj, scientific = FALSE)

cat("Corrected threshold is: ", adjP, "\n")

#########################################################
### Summarize results
#########################################################
  cat("summarizing results.............\n")
  
  reps <- Sys.glob("*.ps")
  
  out.file<-data.frame()
  for (r in reps)
  {  
    Fname<- file_path_sans_ext(r)
    ps <- read.delim(r,header=FALSE)
    Sig<-data.frame()
    inFile<-read.delim(r,header=FALSE)
    print(r)
    Sig<-inFile[inFile[,4]<=adj,] #adj
    if (nrow(Sig)>0) {
      ifelse (nrow(Sig)>0, Sig[,5]<-Fname,NA)
      out.file <- rbind(out.file, Sig)
    }
  	
    ######## If there are NO significant Values in the current looped environment######################
     if (nrow(Sig)==0) { next }  
    ##################################################################

  }
if (nrow(Sig)==0){
  Wording = "There are no significant SNPs in this trait in this environment"
  print (Wording)
  cat(Wording, file=(paste(trait,"_",comparison_name,".sigsnps",sep="")), sep="\t", fill=FALSE, append=FALSE)
  reps <- Sys.glob("*.ps")
  ColSet<-c(color1,color2)
  
  pdf(paste(comparison_name,"_ManhattanPlot_",trait,".pdf",sep=""))
  if(num_env>=num_rep) {
    par(mfrow=(c(num_env,num_rep))) } else {par(mfrow=(c(num_rep,num_env)))}
  for(r in reps) {
    Fname<- file_path_sans_ext(r)
    ps <- read.delim(r,header=FALSE)
    pvalue = ps$V4
    ytop = -log10(min(pvalue))
    if (ytop > 10){ ytop = ytop+2} else { ytop = 12}
    ps1<- merge(snp.map,ps[,c(1,4)],by.x=1,by.y=1)
    colnames(ps1) = c("SNP_ID", "Chr", "Pos", "V4")
    print(r)
    NNN<-strsplit(Fname,"_")[[1]]
    envN<-NNN[2]
   #Print with qfdr:
    manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",main=paste(trait," ",comparison_name,sep=""),col=ColSet,suggestiveline = FALSE, genomewideline=-log10(adj), cex=0.6, cex.axis=0.8,ylim=c(2,ytop))
    }
    dev.off()
    
    cat("starting cleanup.............\n")
    system(paste("rm *.log",sep=""))
    system(paste("rm *.reml",sep=""))
    system(paste("rm *.pheno",sep=""))
    system(paste("rm *.txt",sep=""))
    cat("All Done!.............\n")
  quit()
} ##### This is a critical break; if there are 0 significant SNPs in this trait env , so it breaks

SigSNPs<-out.file
colnames(SigSNPs)<-c("SNP_ID","beta","SE_beta","Pvalue","DataSet")

phenotypes<-Sys.glob("*.pheno")
phenFile<-data.frame()
for (phen in phenotypes)
{
  Pname<- file_path_sans_ext(phen)
  print(phen)
  Var<-data.frame(Pname,"NA")
  trFile <- read.table(phen,header=FALSE)
  Var[,2]<-var(trFile[,3],na.rm=TRUE)
 # ifelse (is.na(Var), Var[,2]<-phen,NA)
  phenFile <- rbind(phenFile, Var)
}

Var.file<-phenFile
colnames(Var.file)<-c("DataSet","Phenotypic_Variance")
SumFile<-merge(SigSNPs,Var.file,by.x=5,by.y=1,all.x=TRUE)


SumFile$Phenotypic_Variance<-as.numeric(SumFile$Phenotypic_Variance)
SumFile$SE_beta<-as.numeric(SumFile$SE_beta)
PVE<-100*(SumFile$SE_beta^2/SumFile$Phenotypic_Variance)
Summary<-data.frame(SumFile[,1:6],PVE)

colnames(snp.map)<-c("SNP_ID","Chr","Pos")
Sum2<-merge(snp.map,Summary,by.x=1,by.y=2)
Sum2$Env = 0
for(val in 1:nrow(Sum2)){
 env_name = sub("Ha412HO_April2017:","",Sum2[val,4])
 env_name = sub(":.*","",env_name)
 Sum2$Env[val] <- env_name
 #Sum2$Chr = gsub("nXRQChr","",Sum2$Chr)
}
#Sum2$Chr = as.numeric(as.character(Sum2$Chr))

write.table(Sum2, paste(trait,"_",comparison_name,".sigsnps",sep=""),row.names = FALSE, quote = FALSE,na="NA",sep="\t")

#########################################################
### Plot results
#########################################################

if(pca=="yes"){
  #1) PCA
  zeva<-brewer.pal(12,"Set3")
  group<-read.table(file="../EnchiladaSuite/Burrito/XRQ_213.group",header=FALSE)
  evPCA<-read.table(paste("../EnchiladaSuite/Burrito/",prefix,".PCA_EV",sep=""),header=FALSE)
  OrPCA<-merge(evPCA,group,by.x=1,by.y=1)
  gr<-unique(OrPCA[,7])
  
  pdf(paste("../../Outputs/Plots/",prefix,"_PCA.pdf",sep=""))
  plot(OrPCA[,3],OrPCA[,4],main=paste(prefix," PCA",sep=""),xlab="PC1",ylab="PC2",col=zeva, pch=16)
  legend("bottomright",legend=gr,fill=zeva)
  
  plot(OrPCA[,4],OrPCA[,5],main=paste(prefix," PCA",sep=""),xlab="PC2",ylab="PC3",col=zeva, pch=16) 
  legend("bottomright",legend=gr,fill=zeva)
  
  plot(OrPCA[,5],OrPCA[,6],main=paste(prefix," PCA",sep=""),xlab="PC3",ylab="PC4",col=zeva, pch=16) 
  legend("bottomright",legend=gr,fill=zeva)
  
  dev.off()          
  
  
  #2) kinship
  aIBS <- read.delim(paste(prefix,".aIBS.kinf",sep=""), header=FALSE)
  row.names(aIBS)<-evPCA[,1]
  
  kinMat<-as.matrix(aIBS)
  
  pdf(paste("../../Outputs/Plots/",prefix," KinShip.pdf",sep=""))
  #heatmap.2(kinMat,main=paste(prefix," KinShip",sep=""),dendrogram='none',trace='none',labRow="",labCol="")
  heatmap.2(kinMat,main=paste(prefix," KinShip",sep=""),trace='none',labRow="",labCol="")
  
  dev.off()
  rm(kinMat)
  rm(aIBS)
} #Kinship and Structure plots

cat ("Running Manhattan Plots.........\n")
reps <- Sys.glob("*.ps")
ColSet<-c(color1,color2)

pdf(paste(comparison_name,"_ManhattanPlot_",trait,".pdf",sep=""))

for(r in reps) {
  Fname<- file_path_sans_ext(r)
  ps <- read.delim(r,header=FALSE)
  pvalue = ps$V4
  ytop = -log10(min(pvalue))
  if (ytop > 10){ ytop = ytop+2} else { ytop = 12}
  ps1<- merge(snp.map,ps[,c(1,4)],by.x=1,by.y=1)
  print(r)
  NNN<-strsplit(Fname,":")[[1]]
  envN<-NNN[2]
  repN<-NNN[3]
  
  #Print no qfdr:
  manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",main=paste(trait," ",comparison_name,sep=""),col=ColSet,suggestiveline = FALSE,genomewideline=-log10(adj), cex=0.6, cex.axis=0.8,ylim=c(2,ytop))
  
  #Print with qfdr:
  #manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",main=paste(trait," ",comparison_name,sep=""),col=ColSet,suggestiveline=-log10(qfdr),genomewideline=-log10(adj), cex=0.6, cex.axis=0.8,ylim=c(2,ytop))
  
  #Print w no qfdr:
  #manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",main=paste(trait," ",comparison_name,sep=""),col=ColSet,genomewideline=-log10(adj), cex=0.6, cex.axis=0.8,ylim=c(2,ytop))
  
}
dev.off()

#########################################################
### Cleanup 
#########################################################
cat("starting cleanup.............\n")
system(paste("rm *.log",sep=""))
system(paste("rm *.reml",sep=""))
system(paste("rm *.pheno",sep=""))
system(paste("rm *.txt",sep=""))
cat("All Done!.............\n")



