favipiravir_mutationbias_figure3=function()
{
  #This is the code used to create figure 3.  It is very similar to that used to make figure 2. The main difference is that some sequences with lots of errors (>10) are removed from the analysis but this doesn't qualitatively effect the analysis..
gb_file_list=list.files(pattern="*gb")[1:12]
transr=array(0,dim=c(4,4,length(gb_file_list)))
mutation_results=array(0,dim=c(length(gb_file_list)))
errgraph=array(0,dim=c(length(gb_file_list),2))
errgraph4=array(0,dim=c(length(gb_file_list),4))
cutoff=4
for (fn in 1:length(gb_file_list))
{
  dataf=readDNAStringSet(gb_file_list[fn])
  print(paste("Found",length(dataf),"sequences",sep=" "))
  
  oldl=length(dataf)
  
  barcodes=array(0,dim=c(length(dataf),1))
  barcodenums=barcodes
  hyphencontrol=2
  if(grepl("-",gb_file_list[fn],fixed=TRUE))
  {
    hyphencontrol=as.integer(length(gregexpr("-",gb_file_list[fn])[[1]]))
    print("Hyphen control on")
  }
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
  datam=datam[,-which(datam[1,]!="-")[c(150:153,237)]] #This is the first difference.  There are a lot of errors around the join of the forward and reverse reads which this removes and the polymorphism at position 237.
  
  #Find the start and ends of the subamplicons
  errnums=0
  listerrs=array(0,dim=c(30*length(dataf),6))
  start1=which(datam[1,]!="-")[1] #A better way for finding the start and finish the amplicon.
  end1=which(datam[1,]!="-")[265]
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
  
  #This bit removes sequences with >10 errors. This is the major difference between figure 2 and figure 3. 
  #You can vary this by changing error_removal andd get rid of it by setting error_removal to 0.
  error_removal=11 #Removes with this number of errors
  seqs_removed=0
  if (error_removal>1)
  {
    qa=as.integer(listerrs[,6])
  print(paste("Removing seqs with >",error_removal-1,"errors",sep=" "))
  freq1=data.frame(table(listerrs[,4]))
  freq2=data.frame(table(freq1$Freq))
  seqs_removed=sum(freq2[which(freq2$Var1==error_removal):dim(freq2)[1],]$Freq)
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
  
  errgraph[fn,1]=(sum(trans5)-trans5[1,4]-trans5[2,3]-trans5[3,2]-trans5[4,1])/numseqbases
  errgraph[fn,2]=(trans5[1,4]+trans5[2,3]+trans5[3,2]+trans5[4,1])/numseqbases
  #For individual transition errors
  errgraph4[fn,1]=trans5[1,4]/numseqbases
  errgraph4[fn,2]=trans5[2,3]/numseqbases
  errgraph4[fn,3]=trans5[3,2]/numseqbases
  errgraph4[fn,4]=trans5[4,1]/numseqbases
}
print(paste("Overall Baseline was",(mutation_results[1]+mutation_results[2]+mutation_results[3])/3*10000,sep=" "))

mutation_results=mutation_results-(mutation_results[1]+mutation_results[2]+mutation_results[3])/3 #Normalizes

print(paste("Transversion Baseline was",(errgraph[1,1]+errgraph[2,1]+errgraph[3,1])/3*10000,sep=" "))
print(paste("Transition Baseline was",(errgraph[1,2]+errgraph[2,2]+errgraph[3,2])/3*10000,sep=" "))

errgraph[,1]=errgraph[,1]-(errgraph[1,1]+errgraph[2,1]+errgraph[3,1])/3
errgraph[,2]=errgraph[,2]-(errgraph[1,2]+errgraph[2,2]+errgraph[3,2])/3

print(paste("A Baseline was",(errgraph4[1,3]+errgraph4[2,3]+errgraph4[3,3])/3*10000*265/100,sep=" "))
print(paste("C Baseline was",(errgraph4[1,4]+errgraph4[2,4]+errgraph4[3,4])/3*10000*265/41,sep=" "))
print(paste("U Baseline was",(errgraph4[1,1]+errgraph4[2,1]+errgraph4[3,1])/3*10000*265/51,sep=" "))
print(paste("G Baseline was",(errgraph4[1,2]+errgraph4[2,2]+errgraph4[3,2])/3*10000*265/73,sep=" "))

errgraph4[,1]=(errgraph4[,1]-(errgraph4[1,1]+errgraph4[2,1]+errgraph4[3,1])/3)*265/51 #The 265/51 adjusts for the number of As in the reference sequence.   This is where I messed up and flipped As and Ts etc.  All fixed now.
errgraph4[,2]=(errgraph4[,2]-(errgraph4[1,2]+errgraph4[2,2]+errgraph4[3,2])/3)*265/73
errgraph4[,3]=(errgraph4[,3]-(errgraph4[1,3]+errgraph4[2,3]+errgraph4[3,3])/3)*265/100
errgraph4[,4]=(errgraph4[,4]-(errgraph4[1,4]+errgraph4[2,4]+errgraph4[3,4])/3)*265/41
errgraph4[]=errgraph4[,c(3,4,1,2)]  #Puts it in +ve orientation
print("The plot shows Figure 3B")
plot(c(0,0,0,100,100,100,10,10,10,1,1,1),mutation_results*10000)
print(mutation_results*10000)

#This is the data for transversions vs transitions.  You can turn on the plots by removing the comments below.
#plot(c(0,0,0,100,100,100,10,10,10,1,1,1),errgraph*10000)
#points(c(0,0,0,100,100,100,10,10,10,1,1,1),errgraph*10000,col='red')

#This gives the individual transitions in the order A C U G 
rnames <- c("0A", "0B", "0C", "100A","100B", "100C","10A", "10B","10C", "1A","1B", "1C")
cnames <- c("A->G", "C->U", "U->C", "G->A")
errgraphm=matrix(nrow=12,ncol=4,dimnames=list(rnames, cnames))
errgraphm[]=as.matrix(errgraph4*10000,dimnames=list(rnames, cnames))
print(errgraphm)
print(errgraph*10000) #Note the order is the filenames so 0,100,10,1
}