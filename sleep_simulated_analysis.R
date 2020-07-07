####
library(survival)
library(mstate)
library(frailtyEM)
source("cox_markov_test.R") #Rename the function
source("simulation_code.R")
source("plotMarkovTest.R")
  

set.seed(2072020)

#Simulate the dataset
N<-70
v <- c(0.5,1.5,0.4,0.1,0.3,0.1)
simdata <-  sim_data_sleep(N,v=v,frailties=TRUE,a=c(42.39,  0.08 , 4.79,  0.78 , 1.92,  3.92),b=c(0.9, 1.5,1,1.5,1.5,0.8))
simdata <- collapse_data(simdata,times=seq(0,10,by=30/3600))

#Convert the data to mstate format
  simdata1 <- simdata
  simdata1$status <- 1*(simdata1$from==1 & simdata1$to==2) + 1*(simdata1$from==2 & simdata1$to==3) + 1*(simdata1$from==3 & simdata1$to==1)
  simdata1$to[simdata1$from==1]<-2
  simdata1$to[simdata1$from==2]<-3
  simdata1$to[simdata1$from==3]<-1
  simdata2 <- simdata
  simdata2$status <- 1*(simdata2$from==1 & simdata2$to==3) + 1*(simdata2$from==2 & simdata2$to==1) + 1*(simdata2$from==3 & simdata2$to==2)
  simdata2$to[simdata2$from==1]<-3
  simdata2$to[simdata2$from==2]<-1
  simdata2$to[simdata2$from==3]<-2
  msdata <- rbind(simdata1,simdata2)
  msdata <- msdata[order(msdata$id,msdata$entry),]
  msdata$trans <- as.numeric(factor(paste(msdata$from,msdata$to,sep=" ")))

  
  #Fit each of the intensities
  fit_trans <- coxph(Surv(entry,exit,status)~strata(trans),control=coxph.control(timefix=TRUE),data=msdata,x=TRUE,model=TRUE)
  #Commenges-Andersen test across all transitions
  caALL <- ca_test(fit_trans,id=msdata$id)
  caALL
  
  #Fit a frailty model to the Awake -> Non-REM transition.
  fit1 <- coxph(Surv(entry,exit,status)~frailty(id),control=coxph.control(timefix=TRUE),data=msdata,subset=(trans==1))
  #W/o frailty
  fit1_0 <- coxph(Surv(entry,exit,status)~1,control=coxph.control(timefix=TRUE),data=msdata,subset=(trans==1),x=TRUE)
  
  simdata$offs1 <- fit1$frail[simdata$id] #Frailty for the 1->2 transition
  
  #Fit a frailty model to the Non-REM -> Awake transition.
  fit3 <- coxph(Surv(entry,exit,status)~frailty(id),control=coxph.control(timefix=TRUE),data=msdata,subset=(trans==3))
  #w/0 frailty
  fit3_0 <- coxph(Surv(entry,exit,status)~1,control=coxph.control(timefix=TRUE),data=msdata,subset=(trans==3),x=TRUE)
  
  simdata$offs3 <- fit3$frail[simdata$id] #Frailty for the 2->1 transition

  #Obtain the weights functions for the 1->2 transition
  owm1 <- weights_multiple(simdata,grid=tseq,from=1,to=2,min_time=0)
  opw_ind1 <- weights_matrix(simdata,grid=tseq,from=1,to=2,min_time=0,other_weights=list(function(x) mean(abs(x),na.rm=TRUE),function(x) max(abs(x),na.rm=TRUE)))
 
  tmat <- transMat(x = list( c(2, 3), c(1,3), c(1,2) ),
                   names = c("Awake", "REM", "Non-REM"))
  B<-1000
  #Version without correcting for frailty
  ct1_0 <- cox_markov_test(simdata,formula=NULL,tfrom=1,tto=2,trans=tmat,grid = tseq,B = B,fn=opw_ind1,fn2=list(function(x) weighted.mean(abs(x),w=owm1,na.rm=TRUE), function(x) mean(abs(x),na.rm=TRUE),function(x) max(abs(x),na.rm=TRUE)))
  #Version treating estimated frailties as offsets.
  ct1 <- cox_markov_test(simdata,formula="offset(offs1)",tfrom=1,tto=2,trans=tmat,grid = tseq,B = B,fn=opw_ind1,fn2=list(function(x) weighted.mean(abs(x),w=owm1,na.rm=TRUE), function(x) mean(abs(x),na.rm=TRUE),function(x) max(abs(x),na.rm=TRUE)))
  
  #Figure S16 in the Supplementary Materials
   plot.MarkovTest(ct1, tseq, what="states", idx=1:50, qsup=3, states=c("Awake","Non-REM","REM"),
                   xlab="Hours since sleep onset", ylab="Log-rank test statistic",
                   main="Awake -> non-REM")
   
   
   