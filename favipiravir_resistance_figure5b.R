favipiravir_resistance_figure5b=function()
{
  #This is the code used to create figure 5b.  It requires msa which can be downloaded by using the following 2 commands "source("https://bioconductor.org/biocLite.R")" and then "biocLite("msa")"
  library(msa)
gb_file_list=list.files(pattern="*krpl.fasta")
transr=array(0,dim=c(4,4,length(gb_file_list)))
mutation_results=array(0,dim=c(length(gb_file_list)))
errgraph4=array(0,dim=c(length(gb_file_list),5))
cutoff=4
for (fn in 1:length(gb_file_list))
{
  dataf=readDNAStringSet(gb_file_list[fn])
  print(paste("Found",length(dataf),"sequences",sep=" "))
  
  oldl=length(dataf)
  
  barcodes=array(0,dim=c(length(dataf),1))
  barcodenums=barcodes
  hyphencontrol=2

  for(a in 2:length(dataf))
  {
    barcodes[a,1]=strsplit(attributes(attributes(dataf[a])$ranges)$NAMES,"-")[[1]][2+hyphencontrol]
    barcodenums[a,1]=as.integer(strsplit(attributes(attributes(dataf[a])$ranges)$NAMES,"-")[[1]][3+hyphencontrol])
  }
  
  
  print(paste("Removing bad barcodes from",gb_file_list[fn],sep=" "))
  goodbarcodes=which(substr(barcodes,5,5)=="A" & substr(barcodes,10,10)=="A" & substr(barcodes,15,15)=="A")
  barcodes=barcodes[c(1,goodbarcodes)] #Remove bad barcodes
  barcodenums=barcodenums[c(1,goodbarcodes)]
  dataf=dataf[c(1,goodbarcodes)]
  print(paste("Removed",oldl-length(dataf), "bad barcodes from",gb_file_list[fn],sep=" "))
  
  print(paste("Cutoff =", cutoff,sep=" "))
  b_cutoff=which(barcodenums>=cutoff)
  barcodes=barcodes[c(1,b_cutoff)] #Remove bad barcodes
  barcodenums=barcodenums[c(1,b_cutoff)]
  oldl=length(dataf)
  dataf=dataf[c(1,b_cutoff)]
  print(paste("Removed",oldl-length(dataf), "barcodes < cutoff from",gb_file_list[fn],sep=" "))
  
  datam=as.matrix(dataf)
  
  #Find the start and ends of the subamplicons
  errnums=0
  listerrs=array(0,dim=c(30*length(dataf),6))
  start1=which(datam[1,]!="-")[1] #A better way for finding the start and finish of the amplicon.
  end1=which(datam[1,]!="-")[245]
  for (a in 2:length(dataf))
  {
    temperrs=which(datam[1,start1:end1]!=datam[a,start1:end1])
    if (length(temperrs)>=1)
    {
      for(ab in 1:length(temperrs))
      {
        listerrs[errnums+ab,1]=datam[1,temperrs[ab]+start1-1] #Ref base
        listerrs[errnums+ab,2]=datam[a,temperrs[ab]+start1-1] #Mut base
        listerrs[errnums+ab,3]=temperrs[ab]+start1-1         #base location in alignment Nb! Insertions mean can't yet compare between alignments
        listerrs[errnums+ab,4]=a #Which sequence it came from in the alignment
        listerrs[errnums+ab,5]=barcodenums[a]
        listerrs[errnums+ab,6]=length(temperrs) #This is used to remove sequences with too many errors which are most likely due to misaligned sequences and not to mutations caused by favipiravir.
      }
      errnums=errnums+length(temperrs)
    }
  }
  listerrs=listerrs[1:errnums,]
  
  #This bit removes sequences with >error_removal number of errors. 
  #You can vary this by changing error_removal and get rid of it by setting error_removal to 0.
  error_removal=10 #Removes with this number of errors
  seqs_removed=0
  if (error_removal>1)
  {
    qa=as.integer(listerrs[,6])
  print(paste("Removing seqs with >",error_removal-1,"errors",sep=" "))
  freq1=data.frame(table(listerrs[,4]))
  freq2=data.frame(table(freq1$Freq))
  freq3=as.integer(as.matrix(freq2$Var1))
  freq4=which(freq3>=error_removal)
  seqs_removed=sum(freq2[freq4,]$Freq) #All this just calculates how many sequences were removed
  listerrs=listerrs[which(qa<error_removal),]
  }
  
  plot(listerrs[,4],listerrs[,3])
  points(listerrs[which(listerrs[,2]=="-"),4],listerrs[which(listerrs[,2]=="-"),3],col='red')
  points(listerrs[which(listerrs[,1]=="-"),4],listerrs[which(listerrs[,1]=="-"),3],col='green')
  print(paste("There were",dim(listerrs)[1],"errors",sep=" "))
  
  #Code to work out which types of errors e.g A->G
  rnames <- c("A", "C", "T", "G")
  cnames <- c("A", "C", "T", "G") 
  trans5=matrix(nrow=4,ncol=4,dimnames=list(rnames, cnames))
  trans5[1,1]=length(which(listerrs[,1]=="A"&listerrs[,2]=="A"))
  trans5[1,2]=length(which(listerrs[,1]=="A"&listerrs[,2]=="C"))
  trans5[1,3]=length(which(listerrs[,1]=="A"&listerrs[,2]=="T"))
  trans5[1,4]=length(which(listerrs[,1]=="A"&listerrs[,2]=="G"))
  trans5[2,1]=length(which(listerrs[,1]=="C"&listerrs[,2]=="A"))
  trans5[2,2]=length(which(listerrs[,1]=="C"&listerrs[,2]=="C"))
  trans5[2,3]=length(which(listerrs[,1]=="C"&listerrs[,2]=="T"))
  trans5[2,4]=length(which(listerrs[,1]=="C"&listerrs[,2]=="G"))
  trans5[3,1]=length(which(listerrs[,1]=="T"&listerrs[,2]=="A"))
  trans5[3,2]=length(which(listerrs[,1]=="T"&listerrs[,2]=="C"))
  trans5[3,3]=length(which(listerrs[,1]=="T"&listerrs[,2]=="T"))
  trans5[3,4]=length(which(listerrs[,1]=="T"&listerrs[,2]=="G"))
  trans5[4,1]=length(which(listerrs[,1]=="G"&listerrs[,2]=="A"))
  trans5[4,2]=length(which(listerrs[,1]=="G"&listerrs[,2]=="C"))
  trans5[4,3]=length(which(listerrs[,1]=="G"&listerrs[,2]=="T"))
  trans5[4,4]=length(which(listerrs[,1]=="G"&listerrs[,2]=="G"))
  print(trans5)
  print(paste("There were",length(which(listerrs[,1]=="-")),"insertions",sep=" "))
  print(paste("There were",length(which(listerrs[,2]=="-")),"deletions",sep=" "))
  numseqbases=sum(nchar(datam))-length(which(datam=="-"))-265-265*seqs_removed
  print(paste("The error rate was",signif(sum(trans5)/numseqbases,digits=3),"or 1 mutation per",round(numseqbases/sum(trans5)),"bases",sep=" "))
  print(paste(sum(trans5),"errors and",numseqbases,"bases",sep=" "))
  dat4=paste(strsplit(gb_file_list[fn],split="[.]")[[1]][1],"newerrlist",sep="-")
  #write.table(listerrs,dat4)
  transr[,,fn]=trans5
  mutation_results[fn]=sum(trans5)/numseqbases
  
   #For individual transition errors
  errgraph4[fn,1]=trans5[1,4]/numseqbases
  errgraph4[fn,2]=trans5[2,3]/numseqbases
  errgraph4[fn,3]=trans5[3,2]/numseqbases
  errgraph4[fn,4]=trans5[4,1]/numseqbases
  errgraph4[fn,5]=(sum(trans5)-trans5[1,4]-trans5[2,3]-trans5[3,2]-trans5[4,1])/numseqbases
}
print(paste("Overall Baseline was",(mutation_results[19]+mutation_results[20])/2*10000,sep=" "))

mutation_results=mutation_results-(mutation_results[19]+mutation_results[20])/2 #Normalizes


print(paste("A Baseline was",(errgraph4[19,1]+errgraph4[20,1])/2*10000*245/64,sep=" "))
print(paste("C Baseline was",(errgraph4[19,2]+errgraph4[20,2])/2*10000*245/53,sep=" "))
print(paste("U Baseline was",(errgraph4[19,3]+errgraph4[20,3])/2*10000*245/70,sep=" "))
print(paste("G Baseline was",(errgraph4[19,4]+errgraph4[20,4])/2*10000*245/58,sep=" "))

errgraph4[,1]=(errgraph4[,1]-(errgraph4[19,1]+errgraph4[20,1])/2)*10000*245/64 #The 245/64 adjusts for the number of As in the reference sequence.
errgraph4[,2]=(errgraph4[,2]-(errgraph4[19,2]+errgraph4[20,2])/2)*10000*245/53
errgraph4[,3]=(errgraph4[,3]-(errgraph4[19,3]+errgraph4[20,3])/2)*10000*245/70
errgraph4[,4]=(errgraph4[,4]-(errgraph4[19,4]+errgraph4[20,4])/2)*10000*245/58
errgraph4[,5]=(errgraph4[,5]-(errgraph4[19,5]+errgraph4[20,5])/2)*10000
print("The plot shows Figure 5B")
plot(c(0,0,10,10,100,100),mutation_results[13:18]*10000,col='blue')
points(c(0,0,10,10,100,100),mutation_results[1:6]*10000,col='red')
points(c(0,0,10,10,100,100),mutation_results[7:12]*10000,col='purple')
points(c(0,0,10,10,100,100),mutation_results[19:24]*10000,col='black')
print(mutation_results*10000)



#This gives the individual transitions in the order A C U G and then transversions found in supplemental figure 2
rnames <- c("KR-0A", "KR-0B", "KR-10A", "KR-10B","KR-100A", "KR-100B","KRPL-0A", "KRPL-0B", "KRPL-10A", "KRPL-10B","KRPL-100A", "KRPL-100B","PL-0A", "PL-0B", "PL-10A", "PL-10B","PL-100A", "PL-100B","WT-0A", "WT-0B", "WT-10A", "WT-10B","WT-100A", "WT-100B")
cnames <- c("A->G", "C->U", "U->C", "G->A","Transversions")
errgraphm=matrix(nrow=24,ncol=5,dimnames=list(rnames, cnames))
errgraphm[]=as.matrix(errgraph4,dimnames=list(rnames, cnames))
print(errgraphm)
return(mutation_results)
}