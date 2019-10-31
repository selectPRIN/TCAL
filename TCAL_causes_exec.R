rm(list=ls())
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)
library(plyr)

source("TCALcauses_functions.R")
options(scipen=10)
   

#If you have multiple core desktop, you'd better parallelize, as bootstrap procedure might take time



no_cores <- 7
cl <- makeCluster(no_cores)
registerDoParallel(cl)
cty <- c("FRATNP","NOR","SWE","CHE",
         #"JPN",
         "AUS","ITA","FRATNP","NOR","SWE","CHE",
         #"JPN",
         "AUS","ITA")
SEx <- c("m","f")
clusterExport(cl,list("melt","CALDecompFunctionCauseboot"))
repli <- 20
res = foreach(i = 1:6, 
              .combine = "rbind") %dopar% {
                  j <- ifelse(i<=6,1,2)
                  x <- CALDecompFunctionCauseboot(cty[i],"JPN",Sex=SEx[j],R=repli)
                  x}

stopCluster(cl)

#Create 
FRAJAP <- res[1:repli,]
NORJAP <- res[(repli+1):(2*repli),]
SWEJAP <- res[(2*repli+1):(3*repli),]
CHEJAP <- res[(3*repli+1):(4*repli),]
AUSJAP <- res[(4*repli+1):(5*repli),]
ITAJAP <- res[(5*repli+1):(6*repli),]



swedat <- as.data.frame(SWEJAP)
nordat <- as.data.frame(NORJAP)
chedat <- as.data.frame(CHEJAP)
fradat <- as.data.frame(FRAJAP)
itadat <- as.data.frame(ITAJAP)
ausdat <- as.data.frame(AUSJAP)
colnames(swedat) <- colnames(nordat) <- colnames(chedat) <- c("INF","NEO","CVD",
          "SYM","MEN","NER", "DIA","DIG",
          "GEN","COS","RES",
          "EXT","LUNG","REST")
colnames(fradat) <- colnames(itadat) <- colnames(ausdat) <- c("INF","NEO","CVD",
          "SYM","MEN","NER", "DIA","DIG",
          "GEN","COS","RES",
          "EXT","LUNG","REST")
swedat2 <- melt(swedat)
nordat2 <- melt(nordat)
chedat2 <- melt(chedat)
fradat2 <- melt(fradat)
itadat2 <- melt(itadat)
ausdat2 <- melt(ausdat)



swe <- summarySE(swedat2, measurevar="value", groupvars="variable")
nor <- summarySE(nordat2, measurevar="value", groupvars="variable")
che <- summarySE(chedat2, measurevar="value", groupvars="variable")
fra <- summarySE(fradat2, measurevar="value", groupvars="variable")
ita <- summarySE(itadat2, measurevar="value", groupvars="variable")
aus <- summarySE(ausdat2, measurevar="value", groupvars="variable")


swe$Country <- "Sweden"
nor$Country <- "Norway"
che$Country <- "Switzerland"
fra$Country <- "France"
ita$Country <- "Italy"
aus$Country <- "Australia"

final <- rbind(swe,nor,che,fra,ita,aus)

final$CI <- 2.64*final$sd

final2 <- final[-which(final$variable=="MEN"|final$variable=="NER"|final$variable=="GEN"|final$variable=="COS"),]

 #Bubbles colors
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(final2, aes(x=variable, y=value, group=1, size = CI,colour=Country))+geom_point()+
ylab("Difference in TCALs: JPN vs other countries")+xlab("Causes of death")+ylim(-1,2)+scale_colour_manual(values=cbPalette)+ scale_size_continuous(range = c(5, 10))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+guides(size = guide_legend(order = 4))

