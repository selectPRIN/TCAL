library(reshape)
library(RColorBrewer)
#### here put the folder where you have your HMD data
CountryA<-"Country"
####  here put the folder where you want to include your results
FOLDER<-"Results"
#### here the ICD code
Data<-"ICD_new3"


#Country Codes
# France 4080, Netherlands 4210, Norway 4220, Sweden 4290, Switzerland 4300, 
# Japan	 3160, Australia 5020,New Zealand 5150, Italy 4180
###########################################################################
cO<-c(4080,4220,4290,4300,3160,5020,5150,4180)


#### here change the countries and  sex 
####you should get the comparison of the two
names<-c("FRATNP","NOR","SWE","CHE","JPN","AUS","NZL_NP","ITA")
names2<-c("France","Norway","Sweden","Switzerland","Japan","Australia","New Zealand","Italy")
Causes<-c("Infectious and parasitic diseases","Neoplasms","Diseases of the circulatory system",
          "Symptoms not elsewhere classified","Mental and behavioural disorders","Diseases of the nervous system",
          "Endocrine, nutritional and metabolic diseases","Diseases of the digestive system",
          "Diseases of the genitourinary system","Congenital malformations","Diseases of the respiratory system",
          "External causes")


################
## This function is needed to read Cause-specific n. deaths of the two countrys (cty1, cty2) you want to compare. The output is a dataset with Year, Age, n. deaths, Cause, Country
################

read.ICD <- function(cty1,cty2,sex="m"){
    cO<-c(4080,4220,4290,4300,3160,5020,5150,4180)
    names<-c("FRATNP","NOR","SWE","CHE","JPN","AUS","NZL_NP","ITA")
    cy1 <- cO[which(names==cty1)]
    cy2 <- cO[which(names==cty2)]

    B1<-t(read.csv(paste(Data,"/Country",cy1,"Cause1",sex,".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
    B2<-t(read.csv(paste(Data,"/Country",cy2,"Cause1",sex,".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])

    for (j in 2:14){
        B1C<-t(read.csv(paste(Data,"/Country",cy1,"Cause",j,sex,".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
        colnames(B1C)[2:112] <- as.character(c(0:110))
        colnames(B1C)[1] <- "Year"
        longB1C <- melt(as.data.frame(B1C),id.vars="Year")
        colnames(longB1C)[2] <- "Age"
        longB1C$Cause <- j
        longB1C$Country <- cty1
        B2C<-t(read.csv(paste(Data,"/Country",cy2,"Cause",j,sex,".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
        colnames(B2C)[2:112] <- as.character(c(0:110))
        colnames(B2C)[1] <- "Year"
        longB2C <- melt(as.data.frame(B2C),id.vars="Year")
        colnames(longB2C)[2] <- "Age"
        longB2C$Cause <- j
        longB2C$Country <- cty2
        if(j==2){doble <- rbind(longB1C,longB2C)}
        if(j>2){doble <- rbind(doble,rbind(longB1C,longB2C))}
    }
    return(doble)
}

        
        

    
bpy.colors <-
  function (n = 100, cutoff.tails = 0.1, alpha = 1)
  {
    n <- as.integer(n[1])
    if (n <= 0)
      return(character(0))
    
    if (cutoff.tails >= 1 || cutoff.tails < 0)
      stop("cutoff.tails should be in [0, 1]")
    i = seq(0.5 * cutoff.tails, 1 - 0.5 * cutoff.tails, length = n)
    r = ifelse(i < .25, 0, ifelse(i < .57, i / .32 - .78125, 1))
    g = ifelse(i < .42, 0, ifelse(i < .92, 2 * i - .84, 1))
    b = ifelse(i < .25, 4 * i, ifelse(i < .42, 1,
                                      ifelse(i < .92, -2 * i + 1.84, i / .08 - 11.5)))
    rgb(r, g, b, alpha)
  }



customAxis <- function() { 
n <- length(levels) 
y <- seq(min(levels), max(levels), length.out=n) 
rect(0, y[1:(n-1)], 1, y[2:n], col=WildColors) 
axis(4, at=y, labels=levels) 
} 



CALDecompFunction<-function(Mx1,Mx2,Y,Yb,Name1,Name2){
    CALlx<-c()
    CALlx1<-c()
    CALlx2<-c()
    PxCh<-c()

    YM<-Y-Yb

    for (x in 1:111){
        if (x <(YM+1)){
            px1<-c()
            px2<-c()
            for (z in 1:x){
                px1<-c(px1,Mx1[z,YM-x+z])
		px2<-c(px2,Mx2[z,YM-x+z])
            }
            pxCH<-c(log(px2/px1),rep(0,111-x))

            lx1<-prod(px1)
            lx2<-prod(px2)
        }
        if (x >(YM)){
            px1<-c()
            px2<-c()
            for (z in (x-YM+1):x){
                px1<-c(px1,Mx1[z,YM-x+z])
		px2<-c(px2,Mx2[z,YM-x+z])
            }
            px1<-c(rep(1,(x-YM)),px1)
            px2<-c(rep(1,(x-YM)),px2)
            pxCH<-c(log(px2/px1),rep(0,111-x))	 	

            lx1<-prod(px1)
            lx2<-prod(px2)
        }
        CALlx1<-c(CALlx1,lx1)
        CALlx2<-c(CALlx2,lx2)
        PxCh<-cbind(PxCh,pxCH)
    }
    CALlx<- t(matrix(rep((CALlx1+ CALlx2)/2,111),111))
    PxCh[is.na(PxCh)]<-0

    ## as Guillot calculates this plus a one for l(0)
    A1<-sum(c(1,CALlx1))+.5
    A2<-sum(c(1,CALlx2))+.5
    A3<-sum(CALlx2)-sum(CALlx1)
    A4<-sum(PxCh* CALlx)

    print(rbind(c(paste("CAL-",Name1),paste("CAL-",Name2),"Diff","est-Diff"),round(c(A1,A2,A3,A4),2)))
    return(A3)

}



CALDecompFunctionCause<-function(Mx1,Mx2,Y,Y1,Name1,Name2,CAU1,CAU2){
                                        # Mx1<-qx1
                                        # Mx2<-qx2
                                        # Y<-Y2
                                        # Name1<-Names[1]
                                        # Name2<-Names[2]
                                        # CAU1<-bb1c
                                        # CAU2<-bb2c
    CALlx<-c()
    CALlx1<-c()
    CALlx2<-c()
    PxCh<-c()
    PxChc<-c()

    YM<-Y-Y1
    for (x in 1:111){
        if (x <(YM+1)){
            px1<-c()
            px2<-c()

            px1c<-c()
            px2c<-c()
            for (z in 1:x){
                px1<-c(px1,Mx1[z,YM-x+z])
                px2<-c(px2,Mx2[z,YM-x+z])
                px1c<-c(px1c,Mx1[z,YM-x+z]^CAU1[z,YM-x+z])
                px2c<-c(px2c,Mx2[z,YM-x+z]^CAU2[z,YM-x+z])
            }
            pxCH<-c(log(px2/px1),rep(0,111-x))
            pxCHc<-c(log(px2c/px1c),rep(0,111-x))
            lx1<-prod(px1)
            lx2<-prod(px2)
        }
        if (x >(YM)){
            px1<-c()
            px2<-c()
            px1c<-c()
            px2c<-c()
            for (z in (x-YM+1):x){
                px1<-c(px1,Mx1[z,YM-x+z])
                px2<-c(px2,Mx2[z,YM-x+z])
                px1c<-c(px1c,Mx1[z,YM-x+z]^CAU1[z,YM-x+z])
                px2c<-c(px2c,Mx2[z,YM-x+z]^CAU2[z,YM-x+z])
            }

            px1<-c(rep(1,(x-YM)),px1)
            px2<-c(rep(1,(x-YM)),px2)

            px1c<-c(rep(1,(x-YM)),px1c)
            px2c<-c(rep(1,(x-YM)),px2c)

            pxCH<-c(log(px2/px1),rep(0,111-x))
            pxCHc<-c(log(px2c/px1c),rep(0,111-x))
            

            lx1<-prod(px1)
            lx2<-prod(px2)
        }
        CALlx1<-c(CALlx1,lx1)
        CALlx2<-c(CALlx2,lx2)

        PxCh<-cbind(PxCh,pxCH)
        PxChc<-cbind(PxChc,pxCHc)
    }

    CALlx<- t(matrix(rep((CALlx1+ CALlx2)/2,111),111))

    PxCh[is.na(PxCh)]<-0
    PxChc[is.na(PxChc)]<-0
    triangle <- PxChc* CALlx
    triangle[is.na(triangle)] <- 0


    ## as Guillot calculates this plus a one for l(0)
    A1<-sum(c(1,CALlx1))+.5
    A2<-sum(c(1,CALlx2))+.5
    A3<-sum(CALlx2)-sum(CALlx1)
    A4<-sum(PxCh* CALlx)
    A5<-sum(triangle)

    print(rbind(c(paste("CAL-",Name1),paste("CAL-",Name2),"Diff","est-Diff","Cause i"),round(c(A1,A2,A3,A4,A5),5)))
    return(A5)
}

count2rows2 <- function(cou){
    n <- sum(cou)
    p <- cou/n
    if(n>0){
        newdat <- cbind(c(2:14), Freq = rmultinom(1, n, p))
    }
    if(n==0){
        newdat <- cbind(c(2:14), Freq = rep(0, 12))
    }
    return(newdat[,2])
    }

CALDecompFunctionCauseboot<-function(cty1,cty2,Sex,R){

    data <- read.ICD(cty1,cty2,sex=Sex)#Reading data
                                        # First caclulating TCAL by cause

    if(Sex=="m"){NAME<-paste(CountryA,"/mltper_1x1_",cty1,".txt",sep="")}
    if(Sex=="f"){NAME<-paste(CountryA,"/fltper_1x1_",cty1,".txt",sep="")}

    A1<-read.table(NAME,header=TRUE,fill=TRUE,skip=1)
    if(Sex=="m"){NAME<-paste(CountryA,"/mltper_1x1_",cty2,".txt",sep="")}
    if(Sex=="f"){NAME<-paste(CountryA,"/fltper_1x1_",cty2,".txt",sep="")}
    A2<-read.table(NAME,header=TRUE,fill=TRUE,skip=1)
    B1<-t(read.csv(paste(Data,"/Country",cO[which(names==cty1)],"Cause1",Sex,".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
 

    B2<-t(read.csv(paste(Data,"/Country",cO[which(names==cty2)],"Cause1",Sex,".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
 
    Y1<-max(range(A1$Year)[1],range(A2$Year)[1])
    Y2<-min(range(A1$Year)[2],range(A2$Year)[2],2014)

    Y1<-max(range(A1$Year)[1],range(A2$Year)[1],range(B1[,1])[1],range(B2[,1])[1])
    Y2<-min(range(A1$Year)[2],range(A2$Year)[2],range(B1[,1])[2],range(B2[,1])[2],2014)
    commonY1 <- max(range(A1$Year)[1],range(A2$Year)[1],range(B1[,1])[1],range(B2[,1])[1])
    commonY2 <- min(range(A1$Year)[2],range(A2$Year)[2],range(B1[,1])[2],range(B2[,1])[2],2014)
    data <- data[data$Year %in% c(commonY1:commonY2),]#Ensuring we are working on the common range of years


    Year1<-Y1-1
    Year2<-Y2+1



    A2<-A2[(A2$Year>(Y1-1))&(A2$Year<(Y2+1)),]
    A1<-A1[(A1$Year>(Y1-1))&(A1$Year<(Y2+1)),]


    b1<-B1[(B1[,1]>(Y1-1))&(B1[,1]<(Y2+1)),-1]

    b2<-B2[(B2[,1]>(Y1-1))&(B2[,1]<(Y2+1)),-1]

    qx1<-matrix(1-A1$qx,111)
    qx2<-matrix(1-A2$qx,111)
    
 #   CALDecompFunctionCause(qx1,qx2,Y2,cty1,cty2,bb1c,bb2c)
    Repl <- matrix(NA,nrow=R,ncol=14)
                                            #Now the bootstrap Conf. Inf.
    for(i in 1:R){
        data$Strat <- paste(data$Age,data$Year,data$Country,sep=":")
        a <- tapply(data$value,data$Strat,count2rows2)#Resampling
        b <- as.data.frame(do.call(rbind, a))
        b$id <- row.names(b)
        c <- melt(b, id=c("id"))
        c$Cause <- as.numeric(gsub("V","",(c$variable)))+1
        newdata <- data
        newdata <- merge(data,c,by.x=c("Strat","Cause"),by.y=c("id","Cause"))
        newdata <- newdata[order(newdata$Country,newdata$Cause,newdata$Year,newdata$Age),]

 ##I have to create from data matrix bb1c and bb2c, which are matrices having (for country 1 and 2) the share of deaths, age-year-country-specific, with respect to total deaths. TOTAL1 and TOTAL2 have the total n. of deaths for each combination of age class (rows) and year(cols)

 ## Now the same matrices should be created having only the deaths (for country 1 and 2 ) due to specific cause of death

        for(j in range(data$Cause)[1]:range(data$Cause)[2]){
            SPEC1 <- t(matrix(newdata$value.y[newdata$Country==cty1&newdata$Cause==j],111))
            SPEC2 <- t(matrix(newdata$value.y[newdata$Country==cty2&newdata$Cause==j],111))

            ## Now CAU1 and CAU2 ive the share of deaths due to specific cause. Note that these are still in age classes, while in the original code data had been split into one year age classes (and qx are in one year age classes)
            CAU1 <- t(SPEC1/b1)
            CAU2 <- t(SPEC2/b2)
            print(c(i,j))

            Repl[i,(j-1)] <- CALDecompFunctionCause(qx1,qx2,Y2,Y1,cty1,cty2,CAU1,CAU2)
            Repl[i,14] <- CALDecompFunction(qx1,qx2,Y2,Y1,cty1,cty2)-sum(Repl[i,c(1,2,3,4,7,8,11,12)])
            }
    print(paste("Iteration n",i,sep=":"))
    }
    return(Repl)
}


## This function is used to plot (with ggplot) results.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
 
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

        # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )
 
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
 
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
 
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
 
    return(datac)
}
