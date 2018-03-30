#########################################################
### Salsa
#########################################################
#
# Author: Rishi R. Masalia, rishimasalia@gmail.com
# 
# Summary: This script groups blocks of significant SNPs together and is run after Burrito. 
# 
# Preferences are controlled in the configuration file ("Tomato.config")
#
# Please set the current working directory, by replacing the "PATH" portion. 
#
# ******* Hard coded infomration!!! ****** If you are NOT using this for GWAS of 213 genotypes in the SAM population. You
# have to change some hard coded information. Harded coded infomration starts on lines: 
#
# 85 (Change the TPED file)
#
#########################################################
### Installing and loading required packages
#########################################################


setwd("/PATH/EnchiladaSuite/SALSA/")

if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("ggthemes")) {
  install.packages("ggthemes", dependencies = TRUE)
  library(ggthemes)
}
if (!require("IRanges")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("IRanges")
  library("IRanges")
}
library(lme4)
library(grid)
library(gridExtra)
library(gridBase)
library(coxme)

############ CONFIG set up
config = read.table(file="Tomato.config")
trait = config[1,1]
num_env = config[2,1]
prefix = config[3,1]
envFiles = c(as.character(config[5,1]),as.character(config[7,1]),as.character(config[9,1]),as.character(config[11,1]),as.character(config[13,1]))
envNames = c(as.character(config[4,1]),as.character(config[6,1]),as.character(config[8,1]),as.character(config[10,1]),as.character(config[12,1]))
thresh = config[14,1] #This is which thresholding value you want to use 'bon' for bonferonni SimpleM, or FDR for the fdr value

###################################
num_env<-as.numeric(as.character(num_env))
print(num_env)

picmat = NULL
##### Loop per Environment
for (t in 1:num_env){
  Environment = envNames[t]
  snp_file = read.table(file=envFiles[t],header=T) 
  if (nrow(snp_file)==0){
    print ("Skipping this file, no rows present.....................")
    break
  }
  chroms = as.factor(snp_file[,2])
  snp_file_good = snp_file[,c(1,2,3,5,6,7,4)]
  
  ############### Loop through chromosomes present
    for (c in levels(chroms)){############ Looks at a set of markers per chromosome and calculates pairwise R2 between them, output a table
    MarkMat = data.frame()
    rishi = snp_file[which(snp_file[,2]==c),]
    rishi = rishi[order(rishi$Pos),]
    rishi_markers = rishi[,1]
    markers = as.data.frame(rishi[,1])
    colnames(markers)= "markers"
    
    binmat=NULL
    ##### Run LDselect.pl ##############
    tped = fread(file="./XRQ_213.tped", header = F)
    colnames(tped)[2] = "markers"
    tmp.tped = merge(markers,tped,by = "markers")
    write.table(tmp.tped,file="tmp.tped",row.names = F,quote = F, col.names = F, sep="\t")
    
    marker_bin_mat = NULL
    tmp_bin_data = NULL
    system ("perl ConvertToPB.pl")
    system ("perl ldSelect.pl -pb Ready_tmp.tped -freq 0.05 -r2 0.49 >tmp.bins")
  
  
    ldselectthing = readLines("tmp.bins")
    #### Clean output and get number of bins #####
    NoBins = gsub("\tother_snps:.*","",ldselectthing[length(ldselectthing)-1])
    NoBins = gsub("Bin ","",NoBins)
    NoBins = as.numeric(as.character(gsub(" ","",NoBins)))
  
    ### create table of bins and markers ####
    sites = grep("total_sites:",readLines("tmp.bins"))
   
    for (i in 1:length(sites)){
      tmp_site = sites[i]
      tmpbin = gsub("\t.*","",ldselectthing[tmp_site])
      tmpbin = as.numeric(as.character(gsub("Bin ","",tmpbin)))
      tmp_nosites = gsub(".*\ttotal_sites: ","",ldselectthing[tmp_site])
      tmp_nosites = as.numeric(as.character(gsub("\t.*","",tmp_nosites)))
      tmp_row = cbind(tmpbin,tmp_nosites)
      binmat = rbind(binmat,tmp_row)
    
      tagline = tmp_site+1
      tmptag = gsub(".*: ","",ldselectthing[tagline])
      tmptag = gsub("\\s+$","",tmptag,perl=T)
      tmptag = strsplit(tmptag," ",perl=T) 
      tmptag = unlist(tmptag)
      if(length(tmptag)==0){next}
      else {
        for(b in 1:length(tmptag)){
        blah = cbind(as.numeric(as.character(tmptag[b])),tmpbin)
        marker_bin_mat = rbind(marker_bin_mat,blah)
        }
      }
      osline = tmp_site+2
      tmpos = gsub(".*: ","",ldselectthing[osline])
      tmpos = gsub("\\s+$","",tmpos,perl=T)
      tmpos = strsplit(tmpos," ",perl=T) 
      tmpos = unlist(tmpos)
      if(length(tmpos)==0){next}
      else {
        for(s in 1:length(tmpos)){
        blahos = cbind(as.numeric(as.character(tmpos[s])),tmpbin)
        marker_bin_mat = rbind(marker_bin_mat,blahos)
        }
      }
    } #End length sites
  
    colnames(marker_bin_mat) = c("Pos","Bin")
    ### Rematch with markers:
    tmp_bin_data = merge(snp_file_good,marker_bin_mat,by="Pos")
    tmp_bin_data = tmp_bin_data[order(tmp_bin_data$Pos),]
 
    ######################################################################################
    ########### Use PLINK to get R2 values for markers & Create graphs ###################
    
    if (length(rishi_markers)==1){ ##### Single Marker ##############
      for(v in 1:NoBins){
        subset = tmp_bin_data[which(tmp_bin_data$Bin==v),]
        start = as.character(subset[1,2])
        stop = as.character(subset[nrow(subset),2])
        ds = as.character(subset[1,7])
        chromosome = subset[1,3]
        binnybinbin = subset[1,8]
        picline = cbind(ds,chromosome,start,stop, binnybinbin)
        
        picmat = rbind(picmat,picline)
      }
    ##########End single marker loop
    } else { ######## More than one marker
      coffee = NULL
      bear = NULL
      mark_list = combn(rishi_markers,2) #makes pairwise combinations for markers extending by column
      for(r in 1:ncol(mark_list)){
        leftmark = as.character(mark_list[1,r])
        rightmark = as.character(mark_list[2,r])
        window<-paste("./bin/plink --file ",prefix," --ld ", leftmark, " ",rightmark, " --out tmp",sep="")
        system(window)   
        
        thing = readLines("tmp.log")
        the_line<- strsplit(thing[28],"", fixed = FALSE)
        r2value<-the_line[[1]][11:18]
        if (r2value[8]=="e"){
          r2value<-the_line[[1]][11:21]
          r2value<-paste(r2value[1:11],collapse="")
          r2value<- as.numeric(r2value)
          if (r2value>1){r2value=r2value*0.1}
          r2value<- format(r2value,scientific = FALSE)
          row = cbind(leftmark, rightmark, r2value, Environment)
          MarkMat = rbind(MarkMat, row)
        }else{
          r2value<-paste(r2value[1:8],collapse="")  
          r2value<- as.numeric(r2value)
          if (r2value>1){r2value=r2value*0.1}
          row = cbind(leftmark, rightmark, r2value, Environment)
          MarkMat = rbind(MarkMat, row)
        }
       } #End single marker in multiple marker set
      
      #### Create the diagonal ###
      digmat = as.data.frame(cbind(markers,markers))
      digmat$r = as.numeric(as.character(0))
      digmat$env = Environment
      colnames(digmat) = c("leftmark","rightmark","r2value","Environment")
      
      #### Create a LD heat map for all sig snps
      MarkMat$r2value=as.numeric(as.character(MarkMat$r2value))
      thing1 = as.data.frame(MarkMat)
      thing1 = rbind(thing1,digmat)
      thing1$leftvalue = as.factor(gsub(".*:","",thing1$leftmark))
      thing1$rightvalue = as.factor(gsub(".*:","",thing1$rightmark))
      bear = thing1[order(thing1$leftvalue, thing1$rightvalue),]
      
      ### Combine two dataframes: PLINK and LdSelect Bins ###
      test = tmp_bin_data[,c(1,8)]
      ###### see if bin is inbetween another bin #####
      test$Bin  = as.factor(test$Bin)
      binss = NULL
      for(o in 1:length(levels(test$Bin))){
        subbin = test[which(test$Bin==o),]
        binstart = subbin[1,1]
        binstop = subbin[nrow(subbin),1]
        bdiff = binstop - binstart
        btmp = cbind(o,binstart,binstop,bdiff)
        binss = rbind(binss,btmp)
      }
      binss = as.data.frame(binss)
      binss = binss[order(rev(binss$bdiff)),]
      sharpie=NULL
      for(y in 1:nrow(test)){
        checkmark = test[y,1]
        pen = NULL
        for (z in 1:nrow(binss)){
          bcstart = binss[z,2]
          bcstop =  binss[z,3]
          if (checkmark >= bcstart && checkmark <= bcstop){checkbin = binss[z,1]
          } else {checkbin = test[y,2]}
        checktmp = cbind(checkmark,checkbin)
        pen = rbind(pen,checktmp)
        pen = as.data.frame(pen)
        trubin = min(pen$checkbin)
        pentmp = cbind(checkmark,trubin)
        }
      sharpie = rbind(sharpie,pentmp)
      }
  
      test = sharpie
      colnames(test)[1] = "leftvalue"
      coffee = merge(bear,test,by="leftvalue")
      colnames(test)[1] = "rightvalue"
      coffee = merge(coffee,test,by="rightvalue")
      colnames(coffee)[7:8] = c("leftbin","rightbin")
      coffee = coffee[order(coffee$r2value),]
      coffee$rightvalue = as.character(coffee$rightvalue)
      coffee$leftvalue = as.character(coffee$leftvalue)
      
      #coffee$test<-factor(coffee$rightmark,levels=levels(coffee$leftmark))
      
      coffee$test<-factor(coffee$rightmark,levels=levels(coffee$leftmark)[nrow(snp_file):1])
      coffee$bla<-c(rep(NA,length(coffee$rightmark)))
      coffee$bla[coffee$rightmark==coffee$leftmark]<-coffee$leftbin[coffee$rightmark==coffee$leftmark]
      coffee$bla<-as.factor(coffee$bla)
      coffee$bla<-match(coffee$bla,levels(coffee$bla))
      
      ############## Create the heatplot ###############3
      #myCol<- c("aliceblue","pink2", "palevioletred1", "indianred2","red2","red4")
      #myColras<- c("aliceblue","blue1","olivedrab2","mediumorchid2", "turquoise1") #need unlimited colors
      rainbowcols = rainbow(10,s=0.5)
      myBreaks <- c(0, .1, .2, .4, .6, .8, 1)
      pdf(file=paste("../../Outputs/Plots/LD_Plots/",trait,"_",Environment,"_Chr_",c,"_LDplot.pdf",sep=""))
      print (
      ggplot(coffee,aes(x=leftmark, y=test))+
          geom_tile(aes(fill=factor(bla)))+
          geom_tile(aes(alpha=r2value))+
          #geom_tile(aes(fill=cut(r2value,breaks=myBreaks)),size=0.5)+
          #scale_fill_manual(values=myCol)+
          scale_fill_manual(drop=T,values=rainbowcols, na.translate=F,name="Bins")+
          coord_equal()+
          #geom_text(aes(label=round(r2value,3)), size = 4)+ #adds the R2 values into the box
          labs(x=NULL,y=NULL,title=paste("LD of chrom ",c," for ",trait," in ", Environment,sep = ""))+
          theme_classic()+
          theme(axis.ticks=element_blank())+
          theme(axis.text=element_text(size=7))+
          theme(axis.text.x=element_text(angle=45,hjust=1))+
          theme(legend.title=element_text(size=8))+
          theme(legend.text=element_text(size=6))
          #theme(legend.position="none")
     )
     dev.off()
      
    #### Create Pre_identified_Clusters file #############
      
     colnames(test)[1] = "Pos"
     tmp_bin = merge(tmp_bin_data,test,by="Pos")
     tmp_bin_data = tmp_bin[,c(1:7,9)]
     colnames(tmp_bin_data)[8] = "Bin"
     tmp_bin_data$Bin = as.factor(tmp_bin_data$Bin)
     tmp_bin_data$Bin<-as.factor(match(tmp_bin_data$Bin,levels(tmp_bin_data$Bin)))
     
      for(p in 1:length(levels(tmp_bin_data$Bin))){
        subset = tmp_bin_data[which(tmp_bin_data$Bin==p),]
        start = as.character(subset[1,2])
        stop = as.character(subset[nrow(subset),2])
        ds = as.character(subset[1,7])
        chromosome = subset[1,3]
        binnybinbin = subset[1,8]
        picline = cbind(ds,chromosome,start,stop, binnybinbin)
        
        picmat = rbind(picmat,picline)
      }
        
  } #End multi marker loop
  
  } #End Chromosome loop 
} #End environment loop
if(length(picmat)==0){ print ("No markers of interest")
} else {
  write.table(picmat,file=paste("./PICs/",trait,"_",Environment,".picfile",sep=""),row.names = F,quote = F, col.names = F, sep="\t")
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
