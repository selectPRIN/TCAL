#
#
#   This program includes the codes for distributing the causes of death aggregated in  
#   age-groups into single ages.             
#   
#   This is an adapted program from the original methodology for disagreagting data
#   by Silvia Rizzi. Adaptation by Stefano Mazzuco and Vladimir Canudas-Romo 
#
#

#Causes reported by WHO are the following


#Causes<-c("Infectious and parasitic diseases","Neoplasms","Diseases of the circulatory system",
#          "Symptoms not elsewhere classified","Mental and behavioural disorders","Diseases of the nervous system",
#          "Endocrine, nutritional and metabolic diseases","Diseases of the digestive system",
#          "Diseases of the genitourinary system","Congenital malformations","Diseases of the respiratory system",
#          "External causes")

#Here some specific causes are considered, the others are aggregated into REST category

# 2)  INF- Infectious and parasitic diseases,
# 3)  NEO-Neoplasms,
# 4)  CVD-Diseases of the circulatory system,
# 5)  SYM-Symptoms not elsewhere classified,
# 8)  DIA-Endocrine, nutritional and metabolic diseases,
# 9)  DIG-Diseases of the digestive system,
# 12) RES-Diseases of the respiratory system,
# 13) EXT-External causes
# REST calculated as a residual class of others causes not included in the above one


#Countries considered are the following, note that the code at their side is that used by WHO to identify them. So you need that code in order to refer to the right file


#France 4080,
#Netherlands 4210,
#Norway  4220,
#Sweden  4290,
#Switzerland 4300, 
#Japan	 3160,
#Australia   5020,
#New Zealand  5150,
#Italy 4180



rm(list=ls())

#ungroup package allows to estimate smoothed distributions from coarsely grouped data. See Pascariu, M. (Producer), Rizzi, S. (Producer), & Dańko, M. J. (Developer). (2018). ungroup: R package: Penalized Composite Link Model for Efficient Estimation of Smooth Distributions from Coarsely Binned Data. Computer programme, Retrieved from https://github.com/mpascariu/ungroup

library(ungroup)

# Cause specific data from WHO come with three different age groups. age0, age1, and age2 are the three groupings
age0<-c(0:4,seq(5,95,5))
age1<-c(0:4,seq(5,85,5))
age2<-c(0,1,seq(5,85,5))


cO<-c(4080,4210,4220,4290,4300,3160,5020,5150,4180)
D<-read.table("../ICD/FormatData.txt",header=TRUE,sep=",",fill=TRUE,skip=0)

for(i in 1:length(cO)){
    dat<-read.table(paste("ICD-",cO[i],".txt",sep=""),sep=",",header=TRUE)
    FD<-D[D$V1==cO[i],] #Data Format
    for(j in 1:c(1,2,3,4,5,8,9,12,13)){
        mal <- NULL
        fem <- NULL
        for(y in 1:length(table(dat$Year))){
            year <- as.numeric(names(table(dat$Year)))[y]
            VV<-FD[which(FD[,4]==year),5]
            if(VV==0&j>1){
                Males <- dat$Male[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1:24)]
                Females <- dat$Female[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1:24)]
                m <- pclm(age0[-c(1:5)],Males[-c(1:5)],nlast=26,control=list(lambda=10^7,deg=2))
                f <- pclm(age0[-c(1:5)],Females[-c(1:5)],nlast=26,control=list(lambda=10^7,deg=2))
                newMale <- c(Males[1:5],as.numeric(m$fitted))
                newFemale <- c(Females[1:5],as.numeric(f$fitted))
            }
            if(VV==0&j==1){
                Males <- dat$Male[dat$Year==year&dat$Cause==1&dat$Age!="Total"&dat$Age!="UNK"][c(1:24)]-tapply(dat[dat$Year==year&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Male,dat[dat$Year==year&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Age,sum)[c(1,2,5,8,11,14,3,4,6,7,9,10,12,13,15:24)]
                Females <- dat$Female[dat$Year==year&dat$Cause==1&dat$Age!="Total"&dat$Age!="UNK"][c(1:24)]-tapply(dat[dat$Year==1960&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Female,dat[dat$Year==1960&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Age,sum)[c(1,2,5,8,11,14,3,4,6,7,9,10,12,13,15:24)]
                m <- pclm(age0[-c(1:5)],Males[-c(1:5)],nlast=26,control=list(lambda=10^7,deg=2))
                f <- pclm(age0[-c(1:5)],Females[-c(1:5)],nlast=26,control=list(lambda=10^7,deg=2))
                newMale <- c(Males[1:5],as.numeric(m$fitted))
                newFemale <- c(Females[1:5],as.numeric(f$fitted))
            }
            if(VV==1&j>1){
                Males <- dat$Male[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1:22)]
                Females <- dat$Female[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1:22)]
                m <- pclm(age1[-c(1:5)],Males[-c(1:5)],nlast=26,control=list(lambda=10^7,deg=2))
                f <- pclm(age1[-c(1:5)],Females[-c(1:5)],nlast=26,control=list(lambda=10^7,deg=2))
                newMale <- c(Males[1:5],as.numeric(m$fitted))
                newFemale <- c(Females[1:5],as.numeric(f$fitted))
            }
            if(VV==1&j==1){
                Males <- dat$Male[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1:22)]-tapply(dat[dat$Year==year&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Male,dat[dat$Year==year&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Age,sum)[c(1,2,5,8,11,14,3,4,6,7,9,10,12,13,15:22)]
                Females <- dat$Female[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1:22)]-tapply(dat[dat$Year==1960&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Female,dat[dat$Year==1960&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Age,sum)[c(1,2,5,8,11,14,3,4,6,7,9,10,12,13,15:24)]
                m <- pclm(age1[-c(1:5)],Males[-c(1:5)],nlast=26,control=list(lambda=10^7,deg=2))
                f <- pclm(age1[-c(1:5)],Females[-c(1:5)],nlast=26,control=list(lambda=10^7,deg=2))
                newMale <- c(Males[1:5],as.numeric(m$fitted))
                newFemale <- c(Females[1:5],as.numeric(f$fitted))
            }
            if(VV==2&j>1){
                Males <- dat$Male[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1,2,6:22)]
                Females <- dat$Female[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1,2,6:22)]
                m <- pclm(age2[-1],Males[-1],nlast=26,control=list(lambda=10^7,deg=2))
                f <- pclm(age2[-1],Females[-1],nlast=26,control=list(lambda=10^7,deg=2))
                newMale <- c(Males[1],as.numeric(m$fitted))
                newFemale <- c(Females[1],as.numeric(f$fitted))
            }
            if(VV==2&j==1){
                Males <- dat$Male[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1,2,6:22)]-tapply(dat[dat$Year==year&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Male,dat[dat$Year==year&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Age,sum)[c(1,2,14,3,4,6,7,9,10,12,13,15:22)]
                Females <- dat$Female[dat$Year==year&dat$Cause==j&dat$Age!="Total"&dat$Age!="UNK"][c(1,2,6:22)]-tapply(dat[dat$Year==1960&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Female,dat[dat$Year==1960&dat$Cause%in%c(2,3,4,5,8,9,12,13)&dat$Age!="Total"&dat$Age!="UNK",]$Age,sum)[c(1,2,14,3,4,6,7,9,10,12,13,15:24)]
                m <- pclm(age2[-1],Males[-1],nlast=26,control=list(lambda=10^7,deg=2))
                f <- pclm(age2[-1],Females[-1],nlast=26,control=list(lambda=10^7,deg=2))
                newMale <- c(Males[1],as.numeric(m$fitted))
                newFemale <- c(Females[1],as.numeric(f$fitted))
            }
            fem <- cbind(fem,c(year,newFemale))
            mal <- cbind(mal,c(year,newMale))
        }
        
        nomem <- paste("../ICD_new2/Country",cO[i],"Cause",j,"m.txt",sep="")
        nomef <- paste("../ICD_new2/Country",cO[i],"Cause",j,"f.txt",sep="")
        write.csv(fem,file=nomef)
        write.csv(mal,file=nomem)
    }
}

            
            
                
        
    

