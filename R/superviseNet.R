#' Supervised Bayesian joint graphical model for simultaneous network estimation and subgroup identification (superviseNet).
#'
#' @param ct  A two-column matrix with the first column being the survival time and the second column being the censoring indicator. The indicator is a binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param xx Input matrix of p genetic measurements consisting of n rows. Each row is an observation vector.
#' @param K Number of subgroups.
#' @param v_0 The spike prior parameter.
#' @param v_1 The slab prior parameter.
#' @param maxiter Maximum of number of iterations.
#' @param p_2 The probability parameter of Bernoulli prior on binary latent indicators \eqn{\gamma_{k,jl}}'s.
#' @param lambda_b A non-negative tuning parameter of laplace prior on subgroup regression parameters.
#' @param lambda_mu A non-negative tuning parameter of laplace prior on subgroup mean parameters.
#' @param tau_0 A non-negative tuning parameter of exponential prior on diagonal elements of precision matrices.
#' @param eps Tolerance for the EM algorithm.
#' @param threshold A small constant that thresholds the final precision matrix estimator.
#' @param member_input  A vector indicating initialized subgroup memberships of each subjects.
#' @param eps_z A small regularized constant used in the risk ratio computation.
#' @param lambda_s A non-negative tuning parameter controlling the similarity across different networks.
#' @param l.m A vector composed of similarity-based matrix.
#' @param l.update Whether update the similarity matrix.
#' @param tau1 A small constant used in the computation of similarity matrix.
#' @return A list
#' @import MASS
#' @importFrom igraph ba.game
#' @importFrom igraph get.edgelist
#' @export
#'
#' @examples
#' n <- 150
#' p <- 100
#' K <- 3
#' set.seed(1)
#' mu01 = mu02 = mu03 = rep(0,p)
#' beta=matrix(0,K,p)
#' beta[1,1:5]=2
#' beta[2,1:5]=-2
#' beta[3,3:7]=1
#' beta0=c(0,0,0)
#' sigma=c(.01,.01,0.01)
#' rate=.8
#' Mu0.list <- list(mu01,mu02,mu03)
#' whole.data <- generate.data(n,Mu0.list,beta,beta0,sigma)
#' whole.data$beta0=beta
#' xx=whole.data$data
#' ct=whole.data$ct
#' set.seed(1)
#' class_old<-sample(1:K,nrow(xx),replace=TRUE)
#' Thetab <- array(0, dim = c(p, p, K))
#' for(k in 1:K)
#' {
#' Thetab[,,k] <- diag(1,p,p)
#' }
#' tau1=0.001
#' l.m=NULL
#' for(i in 1:(p-1)){
#' for(j in (i+1):p){
#' tmm=1/sqrt(Thetab[i,j,]^2+tau1^2)
#' tmp.m=-matrix(tmm,ncol=1)%*%matrix(tmm,nrow=1)
#' diag(tmp.m)=(K-1)*(tmm^2)
#' tmp=tmp.m
#' l.m=c(l.m,c(tmp))
#' }
#' }
#' res= superviseNet(ct,xx,K,lambda_mu=sqrt(dim(ct)[1]*log(p))/2,lambda_b=sqrt(dim(ct)[1]*log(p))/2,v_0=0.057,v_1=1,p_2=0.85,lambda_s=.01,l.m=l.m,member_input=class_old,eps=1e-2,maxiter=50,threshold=1e-3,eps_z=1e-5,l.update=TRUE,tau1=tau1,tau_0=0.01)




superviseNet<-function(
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
    l.update=TRUE,tau1=0.001)
{

  l.m=lambda_s*l.m

  if(length(lambda_b)==1){
    lambda_b=matrix(lambda_b,K,p)
  }


  set.seed(1)
  n <- as.integer(dim(xx)[1])
  p <- as.integer(dim(xx)[2])

  logcensortime=z_esti = log(ct[,1])

  censoring <- ct[,2]

  memb = member_input

  if(length(unique(memb))!=K){
    print('Wrong with the member_input! Random initialization will be used.')
    print(n)
    memb=sample(1:K,n,replace = TRUE)
  }
  Mu = matrix(0,nrow =K, ncol=p)
  prob =rep(0,K)
  for(l in 1:K){
    Mu[l,] = apply(xx[memb==l,], 2, mean)
    prob[l]=sum(memb==l)/n
  }


  Omega = list()
  sigma=rep(0,K)
  beta0=rep(0,K)
  beta=matrix(0,K,p)
  S=list()
  set.seed(1)


  for(k in 1:K){
    S[[k]]  <- diag(1,p,p)

    if(det(S[[k]])<1e-4) S[[k]]=S[[k]]+diag(0.01,p,p)
    Omega[[k]] <- solve(S[[k]])
  }


  for(k in 1:K){
    beta[k,]=0
    beta0[k]=0
    sigma[k]=as.numeric(sqrt(sum((z_esti[memb==k])^2)/nrow(xx[memb == k,])))
  }


  sigma <-1/sigma
  for (k in 1:K) {
    beta0[k]<-beta0[k]*sigma[k]
    beta[k,] <- beta[k,]*sigma[k]
  }
  G_mat = matrix(0,n,K)
  for(jj in 1:n) G_mat[jj, memb[jj]]=1
  W_l= S
  mu = Mu


  # EM algorithm
  tn = 0
  diff_mu = 10
  diff_omega = 10

  nK= apply(G_mat,2,sum)


  eta<-matrix(0,n,K)
  for (k in 1:K){
    eta[,k]=xx%*%beta[k,]+beta0[k]
  }
  #updata gamma_w
  gamma_w=Gamma_int(v_0,v_1, p_2,p,k=K,Theta_l=Omega)

  z_esti=matrix(rep(z_esti,K),ncol=K,byrow = F)
  while(tn < maxiter)
  {

    prob.old = prob
    mu.old = mu
    Theta_l= Omega
    Omega.old = Omega
    G_mat.old = G_mat
    nK.old = nK
    beta.old=beta
    beta0.old=beta0
    sigma.old=sigma
    z_esti.old=z_esti


    # E-step
    for(i in 1:n){
      if(censoring[i]==0){
        for(k in 1:K){
          z_esti[i,k]=0
          mu_tmp=(sum(xx[i,]*beta.old[k,])+beta0.old[k])/sigma.old[k]
          sd_tmp=1/sigma.old[k]
          fc=dnorm(logcensortime[i],mu_tmp,sd_tmp)
          sc=1-pnorm(logcensortime[i],mu_tmp,sd_tmp)
          if(sc==0) sc=eps_z
          tmp1=(mu_tmp+sd_tmp^2*fc/sc)
          if(is.infinite(tmp1)){
            z_esti[i,k]=2^1023*sign(tmp1)
          }else{
            z_esti[i,k]=tmp1
          }

        }
      }
    }

    det_vec=rep(0,K)
    for(kk in 1:K){
      det_vec[kk]=det(Omega.old[[kk]])
      while(det_vec[kk]<1e-4){
        Omega.old[[kk]]=Omega.old[[kk]]+diag(0.01,p,p)
        det_vec[kk]=det(Omega.old[[kk]])
      }
    }

    tmpp=apply(mu.old, 1,function(x) (t(xx) -x))
    Omega.old2=indtemp=list()
    for(kind in 1:K){
      temp=abs(Omega.old[[kind]])>0
      diag(temp)=0
      indtemp[[kind]]=which(colSums(temp)!=0)
      Omega.old2[[kind]]=Omega.old[[kind]][indtemp[[kind]],indtemp[[kind]]]
    }



    ptmp = c()
    for(j in 1:K){
      probtmp=rep(0,n)
      for(i in 1:n){
        probtmp[i]=0
        for(l in 1:K){
          tm2=tmpp[((i-1)*p+1):(i*p),j]
          tm1=tmpp[((i-1)*p+1):(i*p),l]
          tmpind1=setdiff(1:p,indtemp[[j]])
          tmpind2=setdiff(1:p,indtemp[[l]])

          if(j==l){
            con2=con=0
          }else{
            con=matrix(tm2[indtemp[[j]]],nrow=1)%*%Omega.old2[[j]]%*%matrix(tm2[indtemp[[j]]],ncol=1)+sum(diag(Omega.old[[j]])[tmpind1]*(tm2[tmpind1])^2)-matrix(tm1[indtemp[[l]]],nrow=1)%*%Omega.old2[[l]]%*%matrix(tm1[indtemp[[l]]],ncol=1)-sum(diag(Omega.old[[l]])[tmpind2]*(tm1[tmpind2])^2)
            con=con[1,1]

          }

          tmp=exp((con)*0.5)
          if(is.infinite(tmp)) tmp=1e100
          if(j==l){
            mml=1}else{mml=det_vec[l]/det_vec[j]}
          if(is.infinite(mml)) mml=1e100
          if(is.na(mml)) mml=1


          if(censoring[i]==0){
            con3=(eps_z+1-pnorm(logcensortime[i],eta[i,l]/sigma.old[l],1/sigma.old[l]))/(eps_z+1-pnorm(logcensortime[i],eta[i,j]/sigma.old[j],1/sigma.old[j]))
          }else{
            con2=((sigma.old[j]*z_esti[i,1]-eta[i,j])^2)-((sigma.old[l]*z_esti[i,1]-eta[i,l])^2)
            con3=(sigma.old[l]/sigma.old[j])*exp((con2)*0.5)
          }
          probtmp[i]=probtmp[i]+prob.old[l]*(mml)^0.5*tmp*con3
        }
        probtmp[i]=prob.old[j]/probtmp[i]
      }
      ptmp = cbind(ptmp,probtmp)
    }




    ptmp=ptmp/rowSums(ptmp)
    G_mat = ptmp



    #M-step: prob
    prob= colMeans(G_mat)
    if(sum(is.na(prob))>0){
      warning('something wrong with probability calculation!')
      prob[which(is.na(prob),arr.ind = T)]=1/K
    }

    if(sum(prob)!=1) prob=prob/sum(prob)
    nK= apply(G_mat,2,sum)




    #M-setp: Omega
    S_l=list()
    for (k.ind in 1:K) {
      L_ikx = sqrt(G_mat[,k.ind])*t(t(xx) - mu.old[k.ind,])
      S_l[[k.ind]]= t(L_ikx) %*% L_ikx / nK[k.ind]
    }

    lam.m=array(0,c(p,p,K))
    for(k in 1:K){
      lam.m[,,k]=gamma_w[[k]]/v_1 + (1-gamma_w[[k]])/v_0
    }


    tmp=W.ADMM(S_l, lam.m, nK, l.m, epsilon = 1e-5, maxiter = 100, rho = 1,rho.incr = 1.2,rho.max = 1e10,K=K,tau_0=tau_0)

    res1=tmp$W
    for(k in 1:K){

      Omega[[k]]=res1[,,k]
    }
    Theta_l=Omega
    gamma_w=Gamma_int(v_0,v_1,p_2,p,k=K,Theta_l=Omega)





    #update mu
    for(k in 1:K){
      for(j in 1:p){
        g2=sum(G_mat[,k]*(xx%*%matrix(Omega[[k]][,j],ncol=1)))-prob[k]*n*sum(mu.old[k,-j]*Omega[[k]][-j,j])

        if ((abs(g2)>lambda_mu)){
          mu[k,j]=(g2-sign(g2)*lambda_mu)/((prob[k]*n*Omega[[k]][j,j]))
        } else {
          mu[k,j]=0
        }
      }
    }



    #M-step
    #1. update sigma  2. update beta 3. update beta0
    a<-rep(1,K)
    b<-rep(1,K)
    c2<-rep(1,K)
    for (k in 1:K) {
      tmp=sum(G_mat[,k]*(z_esti[,k]^2))
      if(is.infinite(tmp)) tmp=2^1023
      a[k]=tmp
      tmp1=sigma[k]*logcensortime-eta[,k]
      tmp2=1-pnorm(tmp1,0,1)
      tmp2[tmp2==0]=eps_z
      tmp2=tmp1-dnorm(tmp1,0,1)/(tmp2)

      b[k]=-sum(G_mat[,k]*(1-censoring)*(tmp2)*z_esti[,k])+sum(G_mat[,k]*z_esti[,k]*eta[,k])
      c2[k]= sum(G_mat[,k]*(1-censoring)*(tmp2)*eta[,k])+sum(censoring*G_mat[,k])
      if(a[k]==0) a[k]=eps_z
      sigma[k]= (b[k])/(2*a[k])+(sqrt(b[k]^2+4*a[k]*c2[k]))/(2*a[k])
    }


    zeta <- matrix(0,p,K)

    for (k in 1:K) {
      for (j in 1:p) {
        zeta[j,k]= sum(G_mat[,k]*xx[,j]^2)
        sm=sum(G_mat[,k]*(sigma[k]*z_esti[,k]-eta[,k])*xx[,j])+zeta[j,k]*beta[k,j]



        if (abs(sm)>lambda_b[k,j]){
          beta[k,j]=(sm-sign(sm)*lambda_b[k,j])/(zeta[j,k])
        } else {
          beta[k,j]=0
        }

        beta[abs(beta)<threshold]=0
        eta[,k]=eta[,k]+xx[,j]*beta[k,j]-xx[,j]*beta.old[k,j]
      }
      beta0[k]=sum(G_mat[,k]*(sigma[k]*z_esti[,k]-xx%*%beta[k,]))/nK[k]
    }





    tn = tn + 1

    diff_mu = norm(mu.old-mu,type="2")/(norm(mu,type="2")+0.001)
    diff_omega=ss=0
    diff_diagnal= 0
    for(k in 1:K){
      diff_omega=diff_omega+sqrt(sum((Omega.old[[k]]-Omega[[k]])^2))

      diff_diagnal=0
      ss=ss+sqrt(sum(Omega[[k]]^2))
    }
    diff_omega =diff_omega/ss
    diff_beta = norm(beta-beta.old,type="2")/(norm(beta,type="2")+0.001)

    diff_z= norm(z_esti-z_esti.old,type="2")/(norm(z_esti,type="2")+0.001)

    diff_sigma = sqrt(sum((sigma-sigma.old)^2))/(sqrt(sum(sigma^2))+0.001)



    if(max(diff_mu,diff_beta)<eps){
      break;
    }

    if(l.update){
    Thetab <- array(0, dim = c(p, p, K))
    for(k in 1:K)
    {
      Thetab[,,k] <- Omega[[k]]
    }


    l.m=NULL
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        tmm=1/sqrt(Thetab[i,j,]^2+tau1^2)
        tmp.m=-matrix(tmm,ncol=1)%*%matrix(tmm,nrow=1)
        diag(tmp.m)=(K-1)*(tmm^2)

        if(det(tmp.m)<0){ tmp.m=tmp.m+diag(.001,K,K)}
        l.m=c(l.m,c(tmp.m))
      }
    }

    l.m=lambda_s*l.m
    }

  }

  for (k in 1:K) {
    beta[k,]<-beta[k,]/sigma[k]
    beta0[k]<- beta0[k]/sigma[k]
  }

  sigma=1/sigma

  residual2=matrix(0,K,n)
  for(k in 1:K){
    residual2[k,]=(z_esti[,k]-xx%*%beta[k,]-beta0[k])^2
  }

  tmp=array(NA,dim=c(p,p,K))
  for(kk in 1:K){
    tmp[,,kk]=Omega[[kk]]*(abs(Omega[[kk]])>threshold)
    diag(tmp[,,kk])=diag(Omega[[kk]])
  }

  member = apply(G_mat,1,which.max)
  if(class(member)=='list'){
    member=as.numeric(member)
  }
  return(list(mu=mu, Omega=tmp,prob=prob, G_mat= G_mat,member=member,beta=beta,residual2=residual2,beta0=beta0,sigma=sigma,z=z_esti,gamma=gamma_w))





}







Gamma_int<-function(v_0, v_1, p_2, p, k,Theta_l){

  P_l=list();
  for(i in 1:k) {
    P_l[[i]] = (1.0/(1.0 + v_1/v_0*exp(-abs((Theta_l[[i]]))/v_0 + abs((Theta_l[[i]]))/v_1)*(1-p_2)/p_2));
  }


  return(P_l)

}

