rweib <- function(n, a, b)
{
  # Generates from Weibull distribution with rate a (lambda) and shape b (alpha);
  # rate a and shape b correspond to lambda and alpha in Klein & Moeschberger, p. 38, resp.
  return((-log(runif(n))/a)^(1/b))
}

rcweib <- function(n, a, b, s)
{
  # Generates from conditional Weibull distribution with rate a (lambda) and shape b (alpha),
  # conditional on having survived until time s
  return((s^b - log(runif(n))/a)^(1/b))
}

#Convert data from truly continuous to actually panel observed at a series of times
collapse_data <- function(dataset,times) {
  #dataset: dataframe with entries entry, exit, from, to, id
  #times: vector of times at which the processes will be observed
  
  #Matrix of state occupied by each patient at each time
  occmat <- array(NA,c(length(unique(dataset$id)),length(times)))
  for (i in 1:length(times)) {
    tvals <- which(dataset$entry <= times[i] & dataset$exit > times[i])
    if (length(tvals)>0) {
      occmat[cbind(dataset$id[tvals],i)] <- dataset$from[tvals]
    }
  }
  
  #Now create etm format data from this assuming it is continuously observed
  simdat_red <- data.frame()
  for (j in 1:dim(occmat)[1]) {
    ttimes <- (1+which(diff(occmat[j,])!=0))
    ltime <- max(which(!is.na(occmat[j,])))
    if (ltime > max(ttimes)) {
      entries <- c(0, times[ttimes])
      exits <- times[c(ttimes,ltime)]
      froms <- occmat[j,c(ttimes-1,ltime-1)]
      tos <- c(occmat[j,c(ttimes)],99) #99 = censoring indicator
    }else{
      entries <- c(0, times[ttimes[1:(length(ttimes)-1)]])
      exits <- times[ttimes]
      froms <- occmat[j,c(ttimes-1)]
      tos <- occmat[j,c(ttimes)]
    }
    simdat_red <- rbind(simdat_red,data.frame(entry=entries,exit=exits,from=froms,to=tos,id=rep(j,length(tos))))
  }
  return(simdat_red)
}


#Function to generate a dataset similar to the sleeping behaviour data
#Assume Markov with Weibull intensities conditional on independent Gamma frailties
#except after t=6 tendency to stay awake if become awake
sim_data_sleep <- function(N,v,frailties=TRUE,a,b,twaking=6,pwake=0.0125) {
  #v = Gamma frailty variances
  #a = Weibull rates
  #b = Weibull shapes
  #twaking = time at which start to have possibility of staying away
  #pwake = probability of staying awake at any sojourn into awake
  simdata <- NULL
  C <- runif(N,7.5,8.5)
  
  if (frailties) {
    Z1 <- rgamma(N,1/v[1],1/v[1])
    Z2 <- rgamma(N,1/v[2],1/v[2])
    Z3 <- rgamma(N,1/v[3],1/v[3])
    Z4 <- rgamma(N,1/v[4],1/v[4])
    Z5 <- rgamma(N,1/v[5],1/v[5])
    Z6 <- rgamma(N,1/v[6],1/v[6])
  }else{
    Z1<-Z2<-Z3<-Z4<-Z5<-Z6<-rep(1,N)
  }
  
  for (j in 1:N) {
    awake <- 0
    state <- 2 #Everyone starts in Non-REM...
    time <- 0
    cond <- 0
    frs <- tos <- en <- ex <- NULL
    while (time < C[j]) { #NB: No absorbing state
      if (state==1) {
        if (time < twaking) {
          T2 <- rcweib(1,Z1[j]*a[1],b[1] ,time)
          T3 <- rcweib(1,Z2[j]*a[2],b[2] ,time)
          Tp <- min(T2,T3)
        }
        if (time > twaking) {
          if (awake==0) awake <- 1*(runif(1)< pwake) #1/80 chance of staying awake
          if (awake==0) {
            T2 <- rcweib(1,Z1[j]*a[1],b[1] ,time)
            T3 <- rcweib(1,Z2[j]*a[2],b[2] ,time)
          }else{
            T2 <- 10 #Set time beyond the censoring point
            T3 <- 12
          }
          Tp <- min(T2,T3)
        }
        nstate <- 3 - 1*(T2<T3)
      }
      if (state==2) {
        T1 <- rcweib(1,Z3[j]*a[3],b[3] ,time)
        T3 <- rcweib(1,Z4[j]*a[4],b[4] ,time)
        Tp <- min(T1,T3)
        nstate <- 3 - 2*(T1<T3)    
      }
      if (state==3) {
        T1 <- rcweib(1,Z5[j]*a[5],b[5] ,time)
        T2 <- rcweib(1,Z6[j]*a[6],b[6] ,time)
        Tp <- min(T1,T2)
        nstate <- 1+1*(T2<T1)    
      }
      
      if (Tp < C[j]) {
        frs <- c(frs,state)
        tos <- c(tos,nstate)
        en <- c(en,time)
        ex <- c(ex,Tp)
      }else{
        frs <- c(frs,state)
        tos <- c(tos,99)
        en <- c(en,time)
        ex <- c(ex,C[j])  
      }
      time <- Tp
      state <- nstate
    }
    simdata <- rbind(simdata,data.frame(entry=en,exit=ex,from=frs,to=tos,id=rep(j,length(tos))))
  }
  return(simdata)
}
