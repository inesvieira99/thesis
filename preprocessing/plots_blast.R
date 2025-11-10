library(plyr)
library(ggplot2)
library(tidyr)

blastResFull<-read.table(header=TRUE,sep="\t","blastResultsFull.txt")
# exclude short hits
blastResFull<-blastResFull[-which(blastResFull$length<75),]

blastResFull<-blastResFull[which(duplicated(paste(blastResFull$qseqid,blastResFull$qaccver)) == FALSE),]

sapply(split(blastResFull$class,list(blastResFull$qseqid)),table)-> xx
lapply(xx,as.matrix)-> yy

combinedHits<-rbind.fill.matrix(
  t(yy$Lup2108/sum(yy$Lup2108)),
  t(yy$Lup2136/sum(yy$Lup2136)),
  t(yy$Lup2163/sum(yy$Lup2163)),
  t(yy$Lup2164/sum(yy$Lup2164)),
  t(yy$Lup2189/sum(yy$Lup2189)),
  t(yy$Lup2191/sum(yy$Lup2191)),
  t(yy$Lup2202/sum(yy$Lup2202)),
  t(yy$Lup2221/sum(yy$Lup2221)),
  t(yy$Lup2223/sum(yy$Lup2223)),
  t(yy$Lup2241/sum(yy$Lup2241)),
  t(yy$Lup2312/sum(yy$Lup2312)),
  t(yy$Lup2313/sum(yy$Lup2313)),
  t(yy$Lup2323/sum(yy$Lup2323)),
  t(yy$Lup2361/sum(yy$Lup2361)),
  t(yy$Lup2370/sum(yy$Lup2370)),
  t(yy$Lup2378/sum(yy$Lup2378)),
  t(yy$Lup2383/sum(yy$Lup2383)),
  t(yy$Lup2395/sum(yy$Lup2395)),
  t(yy$Lup2453/sum(yy$Lup2453)),
  t(yy$Lup2456/sum(yy$Lup2456)))

combinedHits[which(is.na(combinedHits))]<-0
# let s also hide anything with below 1% freq
combinedHits[which(combinedHits<0.001)]<-0
combinedHitsFLT<-combinedHits[,-which(colSums(combinedHits) == 0)]

namesSamples<-c("Lup2108","Lup2136","Lup2163","Lup2164","Lup2189","Lup2191","Lup2202","Lup2221", "Lup2223","Lup2241","Lup2312","Lup2313","Lup2323","Lup2361","Lup2370","Lup2378","Lup2383","Lup2395","Lup2453","Lup2456")

pdf("blastRes.pdf",w = 10, h = 10)
combinedHitsLong <- gather(data.frame(combinedHitsFLT, Sample = namesSamples), key = "Class", value = "Frequency", -Sample)

ggplot(combinedHitsLong, aes(x = Sample, y = Frequency, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rainbow(length(namesSamples))) +
  labs(title = "Classes, > 1%", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.margin = margin(1, 1, 1, 1, "cm"))


sapply(split(blastResFull$genus,list(blastResFull$qseqid)),table)-> xx
lapply(xx,as.matrix)-> yy

combinedHits<-rbind.fill.matrix(
  t(yy$Lup2108/sum(yy$Lup2108)),
  t(yy$Lup2136/sum(yy$Lup2136)),
  t(yy$Lup2163/sum(yy$Lup2163)),
  t(yy$Lup2164/sum(yy$Lup2164)),
  t(yy$Lup2189/sum(yy$Lup2189)),
  t(yy$Lup2191/sum(yy$Lup2191)),
  t(yy$Lup2202/sum(yy$Lup2202)),
  t(yy$Lup2221/sum(yy$Lup2221)),
  t(yy$Lup2223/sum(yy$Lup2223)),
  t(yy$Lup2241/sum(yy$Lup2241)),
  t(yy$Lup2312/sum(yy$Lup2312)),
  t(yy$Lup2313/sum(yy$Lup2313)),
  t(yy$Lup2323/sum(yy$Lup2323)),
  t(yy$Lup2361/sum(yy$Lup2361)),
  t(yy$Lup2370/sum(yy$Lup2370)),
  t(yy$Lup2378/sum(yy$Lup2378)),
  t(yy$Lup2383/sum(yy$Lup2383)),
  t(yy$Lup2395/sum(yy$Lup2395)),
  t(yy$Lup2453/sum(yy$Lup2453)),
  t(yy$Lup2456/sum(yy$Lup2456)))

combinedHits[which(is.na(combinedHits))]<-0
# let s also hide anything with below 1% freq
combinedHits[which(combinedHits<0.01)]<-0
combinedHitsFLT<-combinedHits[,-which(colSums(combinedHits) == 0)]

combinedHitsLong <- gather(data.frame(combinedHitsFLT, Sample = namesSamples), key = "Genus", value = "Frequency", -Sample)

ggplot(combinedHitsLong, aes(x = Sample, y = Frequency, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rainbow(length(namesSamples))) +
  labs(title = "Genus, > 1%", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.margin = margin(1, 1, 1, 1, "cm"))

dev.off()
---
#Grfico para blast unmapped 

library(plyr)
library(ggplot2)
library(tidyr)

blastResFull<-read.table(header=TRUE,sep="\t","blastUnResultsFull.txt")
# exclude short hits
blastResFull<-blastResFull[-which(blastResFull$length<75),]

blastResFull<-blastResFull[which(duplicated(paste(blastResFull$qseqid,blastResFull$qaccver)) == FALSE),]

sapply(split(blastResFull$class,list(blastResFull$qseqid)),table)-> xx
lapply(xx,as.matrix)-> yy

combinedHits<-rbind.fill.matrix(
  t(yy$Lup2108/sum(yy$Lup2108)),
  t(yy$Lup2136/sum(yy$Lup2136)),
  t(yy$Lup2163/sum(yy$Lup2163)),
  t(yy$Lup2164/sum(yy$Lup2164)),
  t(yy$Lup2189/sum(yy$Lup2189)),
  t(yy$Lup2191/sum(yy$Lup2191)),
  t(yy$Lup2202/sum(yy$Lup2202)),
  t(yy$Lup2221/sum(yy$Lup2221)),
  t(yy$Lup2223/sum(yy$Lup2223)),
  t(yy$Lup2241/sum(yy$Lup2241)),
  t(yy$Lup2312/sum(yy$Lup2312)),
  t(yy$Lup2313/sum(yy$Lup2313)),
  t(yy$Lup2323/sum(yy$Lup2323)),
  t(yy$Lup2361/sum(yy$Lup2361)),
  t(yy$Lup2370/sum(yy$Lup2370)),
  t(yy$Lup2378/sum(yy$Lup2378)),
  t(yy$Lup2383/sum(yy$Lup2383)),
  t(yy$Lup2395/sum(yy$Lup2395)),
  t(yy$Lup2453/sum(yy$Lup2453)),
  t(yy$Lup2456/sum(yy$Lup2456)))

combinedHits[which(is.na(combinedHits))]<-0
# let s also hide anything with below 1% freq
combinedHits[which(combinedHits<0.01)]<-0
combinedHitsFLT<-combinedHits[,-which(colSums(combinedHits) == 0)]

namesSamples<-c("Lup2108","Lup2136","Lup2163","Lup2164","Lup2189","Lup2191","Lup2202","Lup2221", "Lup2223","Lup2241","Lup2312","Lup2313","Lup2323","Lup2361","Lup2370","Lup2378","Lup2383","Lup2395","Lup2453","Lup2456")

pdf("blastUnRes.pdf",w = 10, h = 10)
combinedHitsLong <- gather(data.frame(combinedHitsFLT, Sample = namesSamples), key = "Class", value = "Frequency", -Sample)

ggplot(combinedHitsLong, aes(x = Sample, y = Frequency, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(title = "Classes, > 1%", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.margin = margin(1, 1, 1, 1, "cm"))


sapply(split(blastResFull$genus,list(blastResFull$qseqid)),table)-> xx
lapply(xx,as.matrix)-> yy

combinedHits<-rbind.fill.matrix(
  t(yy$Lup2108/sum(yy$Lup2108)),
  t(yy$Lup2136/sum(yy$Lup2136)),
  t(yy$Lup2163/sum(yy$Lup2163)),
  t(yy$Lup2164/sum(yy$Lup2164)),
  t(yy$Lup2189/sum(yy$Lup2189)),
  t(yy$Lup2191/sum(yy$Lup2191)),
  t(yy$Lup2202/sum(yy$Lup2202)),
  t(yy$Lup2221/sum(yy$Lup2221)),
  t(yy$Lup2223/sum(yy$Lup2223)),
  t(yy$Lup2241/sum(yy$Lup2241)),
  t(yy$Lup2312/sum(yy$Lup2312)),
  t(yy$Lup2313/sum(yy$Lup2313)),
  t(yy$Lup2323/sum(yy$Lup2323)),
  t(yy$Lup2361/sum(yy$Lup2361)),
  t(yy$Lup2370/sum(yy$Lup2370)),
  t(yy$Lup2378/sum(yy$Lup2378)),
  t(yy$Lup2383/sum(yy$Lup2383)),
  t(yy$Lup2395/sum(yy$Lup2395)),
  t(yy$Lup2453/sum(yy$Lup2453)),
  t(yy$Lup2456/sum(yy$Lup2456)))

combinedHits[which(is.na(combinedHits))]<-0
# let s also hide anything with below 1% freq
combinedHits[which(combinedHits<0.01)]<-0
combinedHitsFLT<-combinedHits[,-which(colSums(combinedHits) == 0)]

combinedHitsLong <- gather(data.frame(combinedHitsFLT, Sample = namesSamples), key = "Genus", value = "Frequency", -Sample)

ggplot(combinedHitsLong, aes(x = Sample, y = Frequency, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(title = "Genus, > 1%", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.margin = margin(1, 1, 1, 1, "cm"))

dev.off()

