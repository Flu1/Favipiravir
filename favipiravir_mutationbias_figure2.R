favipiravir_mutationbias_figure2=function()
{
  #This is the code that we used to make Figure 2 from the paper "Determining the Mutation Bias of Favipiravir in Influenza Using Next-generation Sequencing".  Please see the readme for more info.
source("https://bioconductor.org/biocLite.R") #This will make sure you have the tools you need to read in the files.
biocLite("msa")
library(msa)
  
ga_file_list=list.files(pattern="*ga")  #Finds all the files with ga in them.  Hopefullly 12.

transr=array(0,dim=c(4,4,length(ga_file_list))) #Sets up the arrays that will be populated with results
mutation_results=array(0,dim=c(length(ga_file_list)))
errgraph=array(0,dim=c(length(ga_file_list),4))
errgraph2=array(0,dim=c(length(ga_file_list),2))

cutoff=4 #Very important.  This is the cutoff.  You must have this many reads per barcode to be included.  In the paper we used 4.

for (fn in 1:length(ga_file_list))   #Sets up a loop to read each file and get the data
{
  dataf=readDNAStringSet(ga_file_list[fn])   #Reads in the fasta
  print(paste("Found",length(dataf),"sequences",sep=" "))
  
  oldl=length(dataf)
  
  barcodes=array(0,dim=c(length(dataf),1))
  barcodenums=barcodes
  hyphencontrol=0  #This is a bit that worries about if there were hyphens in the name which there shouldn't be for this.
  if(grepl("-",ga_file_list[fn],fixed=TRUE))
  {
    hyphencontrol=as.integer(length(gregexpr("-",ga_file_list[fn])[[1]]))
    print("Hyphen control on")
  }
  for(a in 2:length(dataf))
  {
    barcodes[a,1]=strsplit(attributes(attributes(dataf[a])$ranges)$NAMES,"-")[[1]][2+hyphencontrol]  #Pulls out the barcodes from the fasta
    barcodenums[a,1]=as.integer(strsplit(attributes(attributes(dataf[a])$ranges)$NAMES,"-")[[1]][3+hyphencontrol])  #Pulls out the number of reads for each barcode from the fasta
  }
 
  print(paste("Removing bad barcodes from",ga_file_list[fn],sep=" "))
  goodbarcodes=which(substr(barcodes,5,5)=="A" & substr(barcodes,10,10)=="A" & substr(barcodes,15,15)=="A") #If the barcodes don't have As in the right spot they get removed!
  barcodes=barcodes[c(1,goodbarcodes)] #Remove bad barcodes
  barcodenums=barcodenums[c(1,goodbarcodes)]
  dataf=dataf[c(1,goodbarcodes)]
  print(paste("Removed",oldl-length(dataf), "bad barcodes from",ga_file_list[fn],sep=" "))
  
  print(paste("Cutoff =", cutoff,sep=" "))
  b_cutoff=which(barcodenums>=cutoff)
  barcodes=barcodes[c(1,b_cutoff)] #Remove barcodes which don't have more reads than the cutoff
  barcodenums=barcodenums[c(1,b_cutoff)]
  oldl=length(dataf)
  dataf=dataf[c(1,b_cutoff)]
  print(paste("Removed",oldl-length(dataf), "barcodes < cutoff from",ga_file_list[fn],sep=" "))
  
  datam=as.matrix(dataf) #Turns the fasta into a matrix.  It's big but it lets you check errors.
  
  #Find the start and ends of the subamplicons.  This makes more sense when you have more than one subamplicon! There is an easier way of doing this.
  s_amp1="CGGGGAAAATATGCAACAAT" 
  end_amp1="AGGGTTTCACTTGGACTGGG"
  b=list(dim=c(2))
  s=array(dim=c(2))
  b[1]=aregexec(s_amp1,dataf[1],0.6)
  b[2]=aregexec(end_amp1,dataf[1],0.2)
  s[1]=b[1][[1]][1]
  s[2]=b[2][[1]][1]+attr(b[2][[1]],"match.length")-1
  print("Found these subamplicon splits")
  print(s)
  
  errnums=0
  listerrs=array(dim=c(30*length(dataf),5)) #Fails if more than 30 errors per sequence
  start1=s[1] #Finds the start of the sequence in the alignment
  end1=s[2] #Finds the end of the sequence in the alignment
  for (a in 2:length(dataf)) #Goes through all the sequences 1 by 1 looking for errors
  {
    temperrs=which(datam[1,start1:end1]!=datam[a,start1:end1]) #Finds the errors between the seqeunce and the reference
   
    if (length(temperrs)>=1) #If there were any errors then...
    {
      for(ab in 1:length(temperrs))  #For each error put in listerrs
      {
        listerrs[errnums+ab,1]=datam[1,temperrs[ab]+start1-1] #Reference base
        listerrs[errnums+ab,2]=datam[a,temperrs[ab]+start1-1] #Mutated base
        listerrs[errnums+ab,3]=temperrs[ab]+start1-1         #Base location in alignment Nb! Insertions mean can't yet compare between alignments
        listerrs[errnums+ab,4]=a #Which sequence it came from in the alignment
        listerrs[errnums+ab,5]=barcodenums[a]  #How many reads were from the particular barcode associated with this error
      }
      errnums=errnums+length(temperrs)
    }
  }
  listerrs=listerrs[1:errnums,]
  plot(listerrs[,4],listerrs[,3]) #Where were the errors in the alignment
  points(listerrs[which(listerrs[,2]=="-"),4],listerrs[which(listerrs[,2]=="-"),3],col='red')  #Deletions
  points(listerrs[which(listerrs[,1]=="-"),4],listerrs[which(listerrs[,1]=="-"),3],col='green')   #Insertions
  print(paste("There were",dim(listerrs)[1],"errors",sep=" "))
  
  #Code to work out which types of errors e.g A->G C->A etc.
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
  
  numseqbases=sum(nchar(datam))-length(which(datam=="-"))-245 #How many bases did we sequence? Minus the number of bases for the reference
  print(paste("The error rate was",signif(sum(trans5)/numseqbases,digits=3),"or 1 mutation per",round(numseqbases/sum(trans5)),"bases",sep=" ")) #Calculates the error rates.
  print(paste(sum(trans5),"errors and",numseqbases,"bases",sep=" "))
  dat4=paste(strsplit(ga_file_list[fn],split="[.]")[[1]][1],"newerrlist",sep="-")
  #write.table(listerrs,dat4)  #If you want to save the mutation data use this line.
  transr[,,fn]=trans5  
  mutation_results[fn]=sum(trans5)/numseqbases
  
 #For transition vs transversion
   errgraph2[fn,1]=(sum(trans5)-trans5[1,4]-trans5[2,3]-trans5[3,2]-trans5[4,1])/numseqbases
  errgraph2[fn,2]=(trans5[1,4]+trans5[2,3]+trans5[3,2]+trans5[4,1])/numseqbases
  #For individual transition errors
       errgraph[fn,1]=trans5[1,4]/numseqbases
  errgraph[fn,2]=trans5[2,3]/numseqbases
  errgraph[fn,3]=trans5[3,2]/numseqbases
  errgraph[fn,4]=trans5[4,1]/numseqbases
}
print(paste("Overall Baseline was",(mutation_results[1]+mutation_results[2])/2*10000,sep=" "))
mutation_results=mutation_results-(mutation_results[1]+mutation_results[2])/2 #Normalizes
print("The plot shows Figure 2C")
plot(c(0,0,100,100,10,10,1,1,20,20,50,50),mutation_results*10000)
print(mutation_results*10000)

#This is the data for transversions vs transitions.  You can turn on the plots by removing the comments below.
print(paste("Transversion Baseline was",(errgraph2[1,1]+errgraph2[2,1])/2*10000,sep=" "))
print(paste("Transition Baseline was",(errgraph2[1,2]+errgraph2[2,2])/2*10000,sep=" "))

errgraph2[,1]=errgraph2[,1]-(errgraph2[1,1]+errgraph2[2,1])/2
errgraph2[,2]=errgraph2[,2]-(errgraph2[1,2]+errgraph2[2,2])/2
#plot(c(0,0,100,100,10,10,1,1,20,20,50,50),errgraph2[,1]*10000)
#points(c(0,0,100,100,10,10,1,1,20,20,50,50),errgraph2[,2]*10000,col="red")
print(errgraph2*10000) #Note the order is the filenames so 0,100,10,1,20,50

#This gives the individual transitions in the order A C T G.  Adjusted for the number of bases in the referecnce
print(paste("A Baseline was",(errgraph[1,1]+errgraph[2,1])/2*10000*245/64,sep=" "))
print(paste("C Baseline was",(errgraph[1,2]+errgraph[2,2])/2*10000*245/53,sep=" "))
print(paste("U Baseline was",(errgraph[1,3]+errgraph[2,3])/2*10000*245/70,sep=" "))
print(paste("G Baseline was",(errgraph[1,4]+errgraph[2,4])/2*10000*245/58,sep=" "))

errgraph[,1]=(errgraph[,1]-(errgraph[1,1]+errgraph[2,1])/2)*245/64
errgraph[,2]=(errgraph[,2]-(errgraph[1,2]+errgraph[2,2])/2)*245/53
errgraph[,3]=(errgraph[,3]-(errgraph[1,3]+errgraph[2,3])/2)*245/70
errgraph[,4]=(errgraph[,4]-(errgraph[1,4]+errgraph[2,4])/2)*245/58
rnames <- c("0A", "0B", "100A", "100B","10A", "10B","1A", "1B","20A", "20B","50A", "50B")
cnames <- c("A->G", "C->U", "U->C", "G->A")
errgraphm=matrix(nrow=12,ncol=4,dimnames=list(rnames, cnames))
errgraphm[]=as.matrix(errgraph*10000,dimnames=list(rnames, cnames))
print(errgraphm)
}