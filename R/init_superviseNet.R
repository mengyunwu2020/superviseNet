#' @export
init_superviseNet<-function(
    ct,
    xx,
    K,
    lambda_mu,
    lambda_b,
    tau_0=0.01,
    v_0,
    v_1=1,
    p_2,
    lambda_s,
    l.m,
    member_input,
    eps=1e-2,
    maxiter=50,
    threshold=1e-3,
    eps_z=1e-5,
    l.update=TRUE,tau1=0.001,init_time=40,thres_gamma=0.5){

n=nrow(ct)
p=ncol(xx)
init_set=list()
for (m in 1:init_time){

  if (m==1){
    class_old<-sample(1:K,n,replace=T)
    init_set[[1]]<-class_old
  } else {
    aa=0
    while (aa<0.5*n){
      for (ii in 1:(m-1)){
        class_old<-sample(1:K,n,replace=T)
        aa=sum(abs(class_old-init_set[[ii]]))
        if (aa>0.5*n){
          break
        }
      }
    }
  }

  init_set[[m]]<-class_old
}

  bic_init<-matrix(1e+10,init_time,1)
  lap.list=list()
  for (m in c(1:init_time)){
    class_old=init_set[[m]]
    Thetab <- array(0, dim = c(p, p, K))
    Sb <- array(0, dim = c(p, p, K))
    for(k in 1:K)
    {
      Sb[,,k]  <- cov(xx[class_old == k, , drop = FALSE])

      if(det(Sb[,,k])<1e-4) Sb[,,k]=Sb[,,k]+diag(0.01,p,p)
      Thetab[,,k] <- solve(Sb[,,k])
    }



    l.m=NULL
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        tmm=1/sqrt(Thetab[i,j,]^2+tau1^2)
        tmp.m=-matrix(tmm,ncol=1)%*%matrix(tmm,nrow=1)
        diag(tmp.m)=(K-1)*(tmm^2)
        tmp=tmp.m
        l.m=c(l.m,c(tmp))
      }
    }

    lap.list[[m]]=l.m


    PP <- tryCatch({
      superviseNet(ct,xx,K,lambda_mu=lambda_mu,
                   lambda_b=lambda_b,tau_0=tau_0,
                   v_0=v_0,v_1=v_1,p_2=p_2,
                   lambda_s=lambda_s,l.m=l.m,
                   member_input=class_old,eps=eps,maxiter=maxiter,
                   threshold=threshold, eps_z=eps_z)
    }, error = function(e) {
      message(paste("Error in m =", m, ":", e$message))
      message("Skipping to next iteration...")
      return(NULL)
    })

    if (is.null(PP)) {
      next
    }


    mu_hat=PP$mu
    mmem=PP$member
    gamma_w=PP$gamma
    Theta_hat=PP$Omega
    A.hat2=array(NA,dim=c(p,p,K))
    for(k in 1:K){
      A.hat2[,,k]=abs(gamma_w[[k]])>thres_gamma
    }

    L.mat=matrix(0,dim(xx)[1],K)
    for(il in 1:dim(xx)[1]){
      L.mat[il, mmem[il]]=1
    }

    prob = PP$prob;
    residual2=PP$residual2
    beta_hat=PP$beta
    beta0h=PP$beta0
    sigma_hat=PP$sigma

    tmp= AdapBIC(ct, xx, residual2, mu_hat, (Theta_hat*A.hat2),beta0h, beta_hat, sigma_hat, L.mat=L.mat,pi_vec=prob)


    bic_init[m]<- tmp$bic




  }


  id=which.min(bic_init)[1]

  class_old=init_set[[id]]

  laplace.m= lap.list[[id]]


  return(list(class_old=class_old,laplace.m=laplace.m))
}
