#' @importFrom Matrix bdiag
#' @export

generate.data = function(n,Mu0.list,beta,beta0,sigma,rate=.8,num.differ=7){

  K0 = length(Mu0.list)
  N <- rep(n,K0)
  p = length(Mu0.list[[1]])
  if(K0==3){
  A.list <-Power.law.network3(p,s=10,umin1=0.8,umax1=1,umin2=1.8,umax2=2,umin3=2.8,umax3=3,num.differ=num.differ,Sparse=TRUE,m=1,Dnum=3,rate=1)

  Theta01 <- A.list$A1
  Theta02 <- A.list$A2
  Theta03 <- A.list$A3
  sigma01 <- solve(Theta01)
  sigma02 <- solve(Theta02)
  sigma03 <- solve(Theta03)
  Sigma0.list <- list(sigma01,sigma02,sigma03)
  Theta0.list <- list(Theta01,Theta02,Theta03)
  }else if(K0==2){
    A.list <-Power.law.network1(p,s=10,umin1=0.5,umax1=1,umin2=.5,umax2=1,num.differ=num.differ,Sparse=TRUE,m=1,Dnum=3,rate=1)

    Theta01 <- A.list$A1
    Theta02 <- A.list$A2
    sigma01 <- solve(Theta01)
    sigma02 <- solve(Theta02)
    Sigma0.list <- list(sigma01,sigma02)
    Theta0.list <- list(Theta01,Theta02)

}
  # set.seed(3)
  Mu0=matrix(0,K0,p);L0=NULL;Theta0=array(0, dim = c(p, p, K0));data=NULL
  Y=NULL
  for (k in 1:K0) {
    Mu0[k,] <- Mu0.list[[k]]
    L0 <- c(L0,rep(k,N[k]))
  }
  for (k in 1:K0) {
    Theta0[,,k] <- as.matrix(Theta0.list[[k]])
  }
  for (k in 1:K0) {
    bb=mvrnorm(N[k],Mu0[k,],Sigma0.list[[k]])
    data <- rbind(data,bb)
    Y=c(Y,bb%*%as.matrix(beta[k,],ncol=1)+beta0[k]+rnorm(N[k],0,sigma[k]))
  }
  n_all = dim(data)[1]


  ti=exp(Y)
  ns=dim(data)[1]
  delta=rep(0,ns)





  i=1
  niter=1
  rrrr=NULL
  time <- ti


    delta_all=time_all=NULL
    ntmp=1
    for(kk in 1:K0){
      i=1
      rrrr=NULL


      niter=0
      delta=rep(0,N[kk])
      while((abs(sum(delta)/N[kk]-rate)>0.02)&niter<=50){
        niter=niter+1
        censoring.time<-rgamma(round(N[kk]*0.5), shape =10, scale =10)
        cc<-sample(c(censoring.time,rep(exp(9999),N[kk]-round(N[kk]*0.5))),N[kk])
        time<- pmin(ti[(ntmp):(ntmp+N[kk]-1)], cc)
        delta = time == ti[(ntmp):(ntmp+N[kk]-1)]

      }
      if(niter==51){

        ss=c(0.05,seq(0.1,0.2,0.01),seq(0.3,1,0.1),seq(1,1000,1),1200,2000)

        timm=ccmm=deltam=NULL
        for(i in ss){
          censoring.time<-rgamma(round(N[kk]*0.5), shape =i, scale =i)
          cc<-sample(c(censoring.time,rep(exp(9999),N[kk]-round(N[kk]*0.5))),N[kk])
          time<- pmin(ti[(ntmp):(ntmp+N[kk]-1)], cc)


          ccmm=rbind(ccmm,cc)
          timm=rbind(timm,time)


          delta = time == ti[(ntmp):(ntmp+N[kk]-1)]
          deltam=rbind(deltam,delta)
          rrrr=append(rrrr,sum(delta)/N[kk])
          if(abs(sum(delta)/N[kk]-rate)<=0.02){
            break;
          }
        }
      }
      if(i==2000){
        tmr=rrrr-rate
        tmr[tmr<0]=5
        temp=which.min((tmr))
        censoring.time=ccmm[temp,]
        time=timm[temp,]
        delta=deltam[temp,]
      }
      delta_all=c(delta_all,delta)
      time_all=c(time_all,time)
      ntmp=ntmp+N[kk]

    }

  delta=delta_all
  time=time_all
  censoring=delta
  ct=cbind(time,censoring)







  whole.data=list()
  whole.data$L0=L0
  whole.data$Mu0=Mu0
  whole.data$Theta0=Theta0
  whole.data$data=data
  whole.data$n_all=n_all
  whole.data$K0=K0
  whole.data$Y=Y
  whole.data$ct=ct
  return(whole.data)
}

sim_subi<-function(subi,pp,hh=0.1){
  for (i in 1:pp) {
    subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+hh
  }
  if(det(subi)<=0) print('wrong with subi')


  return(subi)

}
