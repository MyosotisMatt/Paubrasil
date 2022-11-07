#install.packages('devtools',dependencies=T)
#library(devtools)

#install pophelper package from GitHub
#remotes::install_github('royfrancis/pophelper')

library(pophelper)
library(ggplot2)
library(gridExtra)

setwd("C:/Users/edeli/OneDrive/Projets_Recherche/2018_Paubrasilia/Manuscripts/Scripts/07_admixture/data5mod2/")

#This file is important to sort the samples into the right order.
labels<-read.csv("../labels_group.txt",header=T, stringsAsFactors=F, sep = "\t")


#Put all the Q files in one folder and make it as working directory
#Read all the Q files
sfiles <- list.files(path=".", pattern='.Q',full.names=T)
slist <- sortQ(readQ(sfiles,indlabfromfile=T))

class(slist[5]$prunned_data.5.Q)

toto<-slist[5]$prunned_data.5.Q

rownames(toto)<-labels$Labels
rownames(toto)

labels2<-as.data.frame(labels$new.order)
x<-labels$new.order

toto2<-toto[match(labels$new.order,rownames(toto)),]

slist[5]$prunned_data.5.Q<-toto2



#Plot it

my_pal<-c("turquoise","orchid","cornflowerblue","yellowgreen","coral")

#pdf("Admixture.pdf",width = 10, height = 10, useDingbats = F)
p1<-plotQ(slist[c(2,3,4,5,6,7,8,9)],imgoutput="join",barsize = 1,panelspacer = 0.3,divtype = 1,returnplot=T,exportplot=F,showyaxis=T,basesize=10,
          grplab=labels,grplabangle=45,grplabheight = 0.2,grplabpos = 0.8,grplabsize=2,linesize=0.8,pointsize=3)#,clustercol = my_pal)

grid.arrange(p1$plot[[1]])


#pdf("Admixture_K5.pdf",width = 10, height = 10, useDingbats = F)

p2<-plotQ(slist[5],barsize = 1,panelspacer = 0.3,divtype = 1,returnplot=T,exportplot=F,showyaxis=T,basesize=10,
          grplab=labels2,grplabangle=90,grplabpos = 0.7,grplabsize=2.4,linesize=0.8,pointsize=3,clustercol = my_pal,sortind = "all")
grid.arrange(p2$plot[[1]])
#dev.off()

##Figure for cross validation
b<-read.table("CV.txt", header = T)
b$K <- factor(b$K , levels = b$K )

#pdf("Admixture_CV.pdf", height = 4, width = 6, useDingbats = F)
ggplot(data = b, aes(x = K, y = CV,group = 1)) + geom_point(size = 3) + geom_line(lty = 3) + 
  xlab("K") + ylab("Cross Validation error") + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                                                                    fill = NA,
                                                                                                                    size = 1))+xlab("K") + ylab("Cross validation error")
#dev.off()