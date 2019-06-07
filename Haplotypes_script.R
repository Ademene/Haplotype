### Set your work directory containing the files produced in steps 2) and 3).
setwd("~/Bureau/Work_Arthur/2_Article/Article1_Recomb_crypho/Article_Fungal_genetics_Biology/Resoumission/Haplotype_Test")
### Load the files from step 2)
files <- c("ABI_BAR", "ABI_BBE", "ABI_SAR", "BAR_SAR", "BBE_BAR", "BBE_SAR", "CER_ABI", "CER_BAR", "CER_BBE", "CER_CHA", "CER_CIR", "CER_DOI", "CER_SAR", "CHA_ABI", "CHA_BAR", "CHA_BBE", "CHA_CIR", "CHA_DOI", "CHA_SAR", "CIR_ABI", "CIR_BAR", "CIR_BBE", "CIR_SAR", "DOI_ABI", "DOI_BAR", "DOI_BBE", "DOI_CIR", "DOI_SAR", "STC_ABI", "STC_BAR", "STC_BBE", "STC_CER", "STC_CHA", "STC_CIR", "STC_DOI", "STC_SAR", "VG_ABI", "VG_BAR", "VG_BBE", "VG_CER", "VG_CHA", "VG_CIR", "VG_DOI", "VG_SAR", "VG_STC")
### Load the genome file from step 3)
Contig<-read.table("lengthCOntigPacbio.genome", header=TRUE, quote="\"")
table<-data.frame("Contig"=NA,"DEB"=NA,"FIN"=NA,"VAR"=NA)
### Create the haplotype file
for (j in 1:length(files)){
  var<-read.table(paste(files[j],'windowed.pi', sep="."), header = TRUE)
  table<-data.frame("Contig"=NA,"DEB"=NA,"FIN"=NA,"VAR"=NA)
  for (i in 1:dim(Contig)[1])
  {
    ContCible<-Contig[i,1]
    toto<-seq(1, Contig[i,2], by=10000)
    TabCom<-cbind("Contig"=as.character(ContCible),"DEB"=toto,"FIN"=toto+9999,"VAR"=0)
    CCompCi<-subset(var,var[,1]==paste("CrypaYVO003_",ContCible, sep=""))
    titi<-match(CCompCi[,2], TabCom, nomatch = 0)
    lign<-titi-dim(TabCom)[1]
    TabCom[lign,4]<-CCompCi[,4]
    table<-rbind(table,TabCom)
  }
  table <- table[-1,]
  varsubset<-table[complete.cases(table),]
  SNP<-c(1:(nrow(varsubset)))
  mydf<-data.frame(SNP,varsubset)
  mydf$Contig<-as.numeric(factor(mydf$Contig, levels=unique(as.factor(mydf$Contig))))
  mydf$DEB <- as.numeric(mydf$DEB)
  mydf$VAR <- as.numeric(mydf$VAR)
  assign(paste(files[j]),mydf)
}
Haplo <- matrix(nrow = nrow(mydf), ncol = 11)
colnames(Haplo)=(c("Position","VG","STC","CER","CHA","DOI","CIR","ABI","BBE","BAR","SAR"))
Haplo[,1]=mydf[,1]
Haplo[,2]="a"

for (i in 1:nrow(mydf)){
  if (VG_STC[i,5] < 2) Haplo[i,3]="a" else Haplo[i,3]="b"
  if (VG_CER[i,5] < 2) Haplo[i,4]="a"
  if (VG_CHA[i,5] < 2) Haplo[i,5]="a"
  if (VG_DOI[i,5] < 2) Haplo[i,6]="a"
  if (VG_CIR[i,5] < 2) Haplo[i,7]="a"
  if (VG_ABI[i,5] < 2) Haplo[i,8]="a"
  if (VG_BBE[i,5] < 2) Haplo[i,9]="a"
  if (VG_BAR[i,5] < 2) Haplo[i,10]="a"
  if (VG_SAR[i,5] < 2) Haplo[i,11]="a"
  if (STC_CER[i,5] < 2 & is.na(Haplo[i,4])==TRUE & (Haplo[i,3]=="b")==TRUE) Haplo[i,4]="b" 
  if (is.na(Haplo[i,4]==TRUE)) Haplo[i,4]="c"
  if (STC_CHA[i,5] < 2 & is.na(Haplo[i,5])==TRUE & (Haplo[i,3]=="b")==TRUE) Haplo[i,5]="b"
  if (STC_DOI[i,5] < 2 & is.na(Haplo[i,6])==TRUE & (Haplo[i,3]=="b")==TRUE) Haplo[i,6]="b"
  if (STC_CIR[i,5] < 2 & is.na(Haplo[i,7])==TRUE & (Haplo[i,3]=="b")==TRUE) Haplo[i,7]="b"
  if (STC_ABI[i,5] < 2 & is.na(Haplo[i,8])==TRUE & (Haplo[i,3]=="b")==TRUE) Haplo[i,8]="b"
  if (STC_BBE[i,5] < 2 & is.na(Haplo[i,9])==TRUE & (Haplo[i,3]=="b")==TRUE) Haplo[i,9]="b"
  if (STC_BAR[i,5] < 2 & is.na(Haplo[i,10])==TRUE & (Haplo[i,3]=="b")==TRUE) Haplo[i,10]="b"
  if (STC_SAR[i,5] < 2 & is.na(Haplo[i,11])==TRUE & (Haplo[i,3]=="b")==TRUE) Haplo[i,11]="b"
  if (CER_CHA[i,5] < 2 & is.na(Haplo[i,5])==TRUE & (Haplo[i,4]=="c")==TRUE) Haplo[i,5]="c"
  if (is.na(Haplo[i,5]==TRUE)) Haplo[i,5]="d"
  if (CER_DOI[i,5] < 2 & is.na(Haplo[i,6])==TRUE & (Haplo[i,4]=="c")==TRUE) Haplo[i,6]="c"
  if (CER_CIR[i,5] < 2 & is.na(Haplo[i,7])==TRUE & (Haplo[i,4]=="c")==TRUE) Haplo[i,7]="c"
  if (CER_ABI[i,5] < 2 & is.na(Haplo[i,8])==TRUE & (Haplo[i,4]=="c")==TRUE) Haplo[i,8]="c"
  if (CER_BBE[i,5] < 2 & is.na(Haplo[i,9])==TRUE & (Haplo[i,4]=="c")==TRUE) Haplo[i,9]="c"
  if (CER_BAR[i,5] < 2 & is.na(Haplo[i,10])==TRUE & (Haplo[i,4]=="c")==TRUE) Haplo[i,10]="c"
  if (CER_SAR[i,5] < 2 & is.na(Haplo[i,11])==TRUE & (Haplo[i,4]=="c")==TRUE) Haplo[i,11]="c"
  if (CHA_DOI[i,5] < 2 & is.na(Haplo[i,6])==TRUE & (Haplo[i,5]=="d")==TRUE) Haplo[i,6]="d" 
  if (is.na(Haplo[i,6]==TRUE)) Haplo[i,6]="e"
  if (CHA_CIR[i,5] < 2 & is.na(Haplo[i,7])==TRUE & (Haplo[i,5]=="d")==TRUE) Haplo[i,7]="d"
  if (CHA_ABI[i,5] < 2 & is.na(Haplo[i,8])==TRUE & (Haplo[i,5]=="d")==TRUE) Haplo[i,8]="d"
  if (CHA_BBE[i,5] < 2 & is.na(Haplo[i,9])==TRUE & (Haplo[i,5]=="d")==TRUE) Haplo[i,9]="d"
  if (CHA_BAR[i,5] < 2 & is.na(Haplo[i,10])==TRUE & (Haplo[i,5]=="d")==TRUE) Haplo[i,10]="d"
  if (CHA_SAR[i,5] < 2 & is.na(Haplo[i,11])==TRUE & (Haplo[i,5]=="d")==TRUE) Haplo[i,11]="d"
  if (DOI_CIR[i,5] < 2 & is.na(Haplo[i,7])==TRUE & (Haplo[i,6]=="e")==TRUE) Haplo[i,7]="e"
  if (is.na(Haplo[i,7]==TRUE)) Haplo[i,7]="f"
  if (DOI_ABI[i,5] < 2 & is.na(Haplo[i,8])==TRUE & (Haplo[i,6]=="e")==TRUE) Haplo[i,8]="e"
  if (DOI_BBE[i,5] < 2 & is.na(Haplo[i,9])==TRUE & (Haplo[i,6]=="e")==TRUE) Haplo[i,9]="e"
  if (DOI_BAR[i,5] < 2 & is.na(Haplo[i,10])==TRUE & (Haplo[i,6]=="e")==TRUE) Haplo[i,10]="e"
  if (DOI_SAR[i,5] < 2 & is.na(Haplo[i,11])==TRUE & (Haplo[i,6]=="e")==TRUE) Haplo[i,11]="e"
  if (CIR_ABI[i,5] < 2 & is.na(Haplo[i,8])==TRUE & (Haplo[i,7]=="f")==TRUE) Haplo[i,8]="f"
  if (is.na(Haplo[i,8]==TRUE)) Haplo[i,8]="g"
  if (CIR_BBE[i,5] < 2 & is.na(Haplo[i,9])==TRUE & (Haplo[i,7]=="f")==TRUE) Haplo[i,9]="f"
  if (CIR_BAR[i,5] < 2 & is.na(Haplo[i,10])==TRUE & (Haplo[i,7]=="f")==TRUE) Haplo[i,10]="f"
  if (CIR_SAR[i,5] < 2 & is.na(Haplo[i,11])==TRUE & (Haplo[i,7]=="f")==TRUE) Haplo[i,11]="f"
  if (ABI_BBE[i,5] < 2 & is.na(Haplo[i,9])==TRUE & (Haplo[i,8]=="g")==TRUE) Haplo[i,9]="g"
  if (is.na(Haplo[i,9]==TRUE)) Haplo[i,9]="h"
  if (ABI_BAR[i,5] < 2 & is.na(Haplo[i,10])==TRUE & (Haplo[i,8]=="g")==TRUE) Haplo[i,10]="g"
  if (ABI_SAR[i,5] < 2 & is.na(Haplo[i,11])==TRUE & (Haplo[i,8]=="g")==TRUE) Haplo[i,11]="g"
  if (BBE_BAR[i,5] < 2 & is.na(Haplo[i,10])==TRUE & (Haplo[i,9]=="h")==TRUE) Haplo[i,10]="h"
  if (is.na(Haplo[i,10]==TRUE)) Haplo[i,10]="i"
  if (BBE_SAR[i,5] < 2 & is.na(Haplo[i,11])==TRUE & (Haplo[i,9]=="h")==TRUE) Haplo[i,11]="h"
  if (BAR_SAR[i,5] < 2 & is.na(Haplo[i,11])==TRUE & (Haplo[i,10]=="i")==TRUE) Haplo[i,11]="i"
  if (is.na(Haplo[i,11]==TRUE)) Haplo[i,11]="j"
}

Haplo_sorted = cbind(Haplo[,1],Haplo[,5],Haplo[,3],Haplo[,2],Haplo[,4],Haplo[,6:11])
colnames(Haplo_sorted)=(c("Position","CHA","STC","VG","CER","DOI","CIR","ABI","BBE","BAR","SAR"))


#### Give the colors according to the haplo, from 1 to 4, with 1 most frequent and 4 the less frequent
### Here for the North American introduction
for (I in 1:nrow(Haplo_sorted)){
  Haplo<-table(Haplo_sorted[I,c(2:8)])
  Haplo<-sort(Haplo, decreasing = TRUE)
  if ((length(Haplo)==1)==TRUE) {
    Haplo1<-Haplo[1]
    for (J in 2:8){
      if ((labels(Haplo1)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"1"
    }
  }
  if ((length(Haplo)==2)==TRUE) {
    Haplo1<-Haplo[1]
    Haplo2<-Haplo[2]
    for (J in 2:8){
      if ((labels(Haplo1)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"1"
      if ((labels(Haplo2)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"2"
    }
  }
  if ((length(Haplo)==3)==TRUE) {
    Haplo1<-Haplo[1]
    Haplo2<-Haplo[2]
    Haplo3<-Haplo[3]
    for (J in 2:8){
      if ((labels(Haplo1)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"1"
      if ((labels(Haplo2)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"2"
      if ((labels(Haplo3)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"3"
    }
  }
  if ((length(Haplo)==4)==TRUE) {
    Haplo1<-Haplo[1]
    Haplo2<-Haplo[2]
    Haplo3<-Haplo[3]
    Haplo4<-Haplo[4]
    for (J in 2:8){
      if ((labels(Haplo1)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"1"
      if ((labels(Haplo2)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"2"
      if ((labels(Haplo3)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"3"
      if ((labels(Haplo4)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"4"
    }
  }
}

### Here for the Asian introduction

for (I in 1:nrow(Haplo_sorted)){
  Haplo<-table(Haplo_sorted[I,c(9:11)])
  Haplo<-sort(Haplo, decreasing = TRUE)
  if ((length(Haplo)==1)==TRUE) {
    Haplo1<-Haplo[1]
    for (J in 9:11){
      if ((labels(Haplo1)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"1"
    }
  }
  if ((length(Haplo)==2)==TRUE) {
    Haplo1<-Haplo[1]
    Haplo2<-Haplo[2]
    for (J in 9:11){
      if ((labels(Haplo1)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"1"
      if ((labels(Haplo2)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"2"
    }
  }
  if ((length(Haplo)==3)==TRUE) {
    Haplo1<-Haplo[1]
    Haplo2<-Haplo[2]
    Haplo3<-Haplo[3]
    for (J in 9:11){
      if ((labels(Haplo1)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"1"
      if ((labels(Haplo2)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"2"
      if ((labels(Haplo3)==Haplo_sorted[I,J])==TRUE) Haplo_sorted[I,J]<-"3"
    }
  }
}

### plot
plot.new(); plot.window(xlim = c(0,nrow(Haplo_sorted)), ylim = c(-2,21))
par(xpd=T)
for (j in 2:8){
  for (i in 1:nrow(Haplo_sorted)){
    if (Haplo_sorted[i,j]=="1") segments(x0 = i, y0 = 2*(j-2), x1 = i, y1 = 2*(j-1), col = "white")
    if (Haplo_sorted[i,j]=="2") segments(x0 = i, y0 = 2*(j-2), x1 = i, y1 = 2*(j-1), col = "black")
    if (Haplo_sorted[i,j]=="3") segments(x0 = i, y0 = 2*(j-2), x1 = i, y1 = 2*(j-1), col = "blue")
    if (Haplo_sorted[i,j]=="4") segments(x0 = i, y0 = 2*(j-2), x1 = i, y1 = 2*(j-1), col = "red")
    rect(xleft = 1, ybottom = 2*(j-2), xright = nrow(Haplo_sorted), ytop = 2*(j-1))
  }
}
for (j in 9:11){
  for (i in 1:nrow(Haplo_sorted)){
    if (Haplo_sorted[i,j]=="1") segments(x0 = i, y0 = 2*(j-2), x1 = i, y1 = 2*(j-1), col = "white")
    if (Haplo_sorted[i,j]=="2") segments(x0 = i, y0 = 2*(j-2), x1 = i, y1 = 2*(j-1), col = "black")
    if (Haplo_sorted[i,j]=="3") segments(x0 = i, y0 = 2*(j-2), x1 = i, y1 = 2*(j-1), col = "blue")
    rect(xleft = 1, ybottom = 2*(j-2), xright = nrow(Haplo_sorted), ytop = 2*(j-1))
  }
}

### Names of scaffolds ###
segments(x0 = 0,y0 = 20, x1 = 0, y1 = 21, col= "black")
segments(x0 = 3903,y0 = 20, x1 = 3903, y1 = 21, col= "black")
text(x = -225,y = 1, labels = "RE019")
text(x = -225,y = 3, labels = "RE092_1")
text(x = -225,y = 5, labels = "VG_1896")
text(x = -225,y = 7, labels = "RE092_2")
text(x = -225,y = 9, labels = "RE079")
text(x = -225,y = 11, labels = "RE103")
text(x = -225,y = 13, labels = "H13")
text(x = -225,y = 15, labels = "RE053")
text(x = -225,y = 17, labels = "RE043")
text(x = -225,y = 19, labels = "RE028")
text(x = -225,y = 20.5, labels = "scaffolds")
segments(x0 = -150,y0 = 14, x1 = 4029, y1 = 14, col= "red")
segments(x0 = 0,y0 = 21.5, x1 = 715, y1 = 21.5, col= "black")
segments(x0 = 716,y0 = 22.5, x1 = 778, y1 = 22.5, col= "black")
segments(x0 = 779,y0 = 21.5, x1 = 1287, y1 = 21.5, col= "black")
segments(x0 = 1288,y0 = 22.5, x1 = 1755, y1 = 22.5, col= "black")
segments(x0 = 1756,y0 = 21.5, x1 = 1882, y1 = 21.5, col= "black")
segments(x0 = 1883,y0 = 22.5, x1 = 2203, y1 = 22.5, col= "black")
segments(x0 = 2204,y0 = 21.5, x1 = 2237, y1 = 21.5, col= "black")
segments(x0 = 2238,y0 = 22.5, x1 = 2534, y1 = 22.5, col= "black")
segments(x0 = 2535,y0 = 21.5, x1 = 2817, y1 = 21.5, col= "black")
segments(x0 = 2818,y0 = 22.5, x1 = 2930, y1 = 22.5, col= "black")
segments(x0 = 2931,y0 = 21.5, x1 = 3044, y1 = 21.5, col= "black")
segments(x0 = 3045,y0 = 22.5, x1 = 3431, y1 = 22.5, col= "black")
segments(x0 = 3432,y0 = 21.5, x1 = 3708, y1 = 21.5, col= "black")
segments(x0 = 3709,y0 = 22.5, x1 = 3879, y1 = 22.5, col= "black")
text(x = 357,y = 22,cex = 0.6, labels = "MS1")
text(x = 747,y = 23,cex = 0.6, labels = "MS2")
text(x = 1033,y = 22,cex = 0.6, labels = "MS3")
text(x = 1521,y = 23,cex = 0.6, labels = "MS4")
text(x = 1819,y = 22,cex = 0.6, labels = "MS5")
text(x = 2043,y = 23,cex = 0.6, labels = "MS6")
text(x = 2220,y = 22,cex = 0.6, labels = "MS7")
text(x = 2386,y = 23,cex = 0.6, labels = "MS8")
text(x = 2676,y = 22,cex = 0.6, labels = "MS9")
text(x = 2874,y = 23,cex = 0.6, labels = "MS10")
text(x = 2987,y = 22,cex = 0.6, labels = "MS11")
text(x = 3238,y = 23,cex = 0.6, labels = "C02")
text(x = 3570,y = 22,cex = 0.6, labels = "RC05")
text(x = 3794,y = 23,cex = 0.6, labels = "C09")
segments(x0 = 100,y0 = 23.5, x1 = 200, y1 = 23.5, col= "black")
text(x = 150,y = 24,cex = 0.6, labels = "1 Mb")
for (i in 1:max(mydf[,2])){
  for (j in 1:nrow(mydf)){
    if (mydf[j,2] == i+1) segments(x0 = j,y0 = 20, x1 = j, y1 = 21, col= "black") & break
  }
}