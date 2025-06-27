#' Supervised heterogeneous network estimation via survival-based Bayesian graphical models (superviseNet) along path.
#'
#' @param ct  A two-column matrix with the first column being the survival time and the second column being the censoring indicator. The indicator is a binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param xx Input matrix of p genetic measurements consisting of n rows. Each row is an observation vector.
#' @param Kseq A user supplied vector of number of subgroups.
#' @param lambda_mu A user supplied non-negative tuning parameters of laplace prior on subgroup mean parameters.
#' @param v1 A user supplied slab prior parameters.
#' @param v0 A user supplied spike prior parameters.
#' @param lambda_beta A user supplied non-negative tuning parameters of laplace prior on subgroup regression parameters.
#' @param lambda_sim A user supplied non-negative tuning parameters controlling the similarity across different networks.
#' @param tau_0 A user supplied non-negative tuning parameter of exponential prior on diagonal elements of precision matrices.
#' @param p2 A user supplied probability parameters of Bernoulli prior on binary latent indicator $gamma_{k,jl}$.
#' @param l.m A vector composed of similarity-based matrix.
#' @param member_input A vector indicating initialized subgroup memberships of each subjects.
#' @param eps Tolerance for the EM algorithm.
#' @param maxiter Maximum of number of iterations.
#' @param threshold A small constant that thresholds the final precision matrix estimator.
#' @param eps_z A small regularized constant used in the risk ratio computation.
#' @param l.update Whether update the similarity matrix.
#' @param tau1 A small constant used in the computation of similarity matrix.
#' @return A list
#' @import MASS
#' @importFrom igraph ba.game
#' @importFrom igraph get.edgelist
#' @export
#'
#' @examples
#' n <- 150              # The sample size of each subgroup
#' p <- 100             # The dimension of the precision matrix
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
#' res=superviseNetpath(ct,xx,K,lambda_mu=sqrt(dim(ct)[1]*log(p))/2,lambda_beta=sqrt(dim(ct)[1]*log(p))/2,tau_0=0.01,v0=seq(0.057,.06,length.out=2),v1=1,p2=0.85,lambda_sim=c(0.1,0.05,0.01),l.m=l.m,member_input=class_old,eps=1e-2,maxiter=50,threshold=1e-3,eps_z=1e-5)
superviseNetpath= function(ct, xx, Kseq,lambda_mu=0,v1=1,v0,lambda_beta,lambda_sim,tau_0=0.01,p2,l.m,member_input,eps =1e-2,maxiter=50,threshold=1e-3,eps_z=1e-5,l.update=TRUE,tau1=0.001){


  L1 = length(v1)
  L2 = length(v0)
  L3 = length(lambda_beta)
  L4 = length(lambda_sim)
  L5 = length(p2)
  L6 = length(tau_0)
  p=ncol(xx)
  LL = L1*L2*L3*L4*L5*L6

  bic=rep(1e100,length(Kseq))
  lt=0
  resultall=vector('list',length(Kseq))
  oobic=oobic2=oobic3=rep(10e10,length(Kseq))


  for(K in Kseq){
    cat('the ',K,'th kseq is running','\n')
    lt=lt+1

    Mu_hat.list = list()
    Theta_hat.list = list()
    prob.list = list()
    L.mat.list = list()
    member.list = list()
    beta.list=sigma.list=list()
    aBIC=aicp=ADBIC=rep(10^10, LL)
    lam=matrix(0,LL,6)
    ltl=0

    for(jj in 1:L1){
      for(ll in 1:L2){
        for(mm in 1:L3){
          for(ss in 1:L4){
            for(qq in 1:L5){
              for(yy in 1:L6){
          ltl=ltl+1

          v_1=v1[jj]
          v_0=v0[ll]
          lambda_b=lambda_beta[mm]
          lambda_s=lambda_sim[ss]
          p_2=p2[qq]
          tau0=tau_0[yy]
          lam[ltl,]=c(v_1,v_0,lambda_b,lambda_s,p_2,tau0)

          cat('the ',ltl,'th lams is running','\n')
          print(lam[ltl,])

          PP=superviseNet(ct,xx, K=K,
                                v_0=v_0,
                                v_1=v_1,
                                maxiter=maxiter,
                                p_2=p_2,
                                lambda_b=lambda_b,lambda_mu=lambda_mu,tau_0=tau0,
                                threshold=threshold,eps=eps,member_input=member_input,eps_z=eps_z,lambda_s = lambda_s,l.m =l.m,l.update=l.update
          )

            mu_hat=PP$mu
            mmem=PP$member
            Theta_hat=PP$Omega

            L.mat=matrix(0,dim(xx)[1],K)
            for(il in 1:dim(xx)[1]){
              L.mat[il, mmem[il]]=1
            }

            prob = PP$prob;
            residual2=PP$residual2
            beta_hat=PP$beta
            beta0=PP$beta0
            sigma_hat=(PP$sigma)^2

          prob.list[[ltl]]=prob;
          Mu_hat.list[[ltl]]=mu_hat;
          Theta_hat.list[[ltl]]=Theta_hat;
          L.mat.list[[ltl]] = L.mat
          beta.list[[ltl]]=c(beta0,beta_hat)
          sigma.list[[ltl]]=sigma_hat


          tmp=AdapBIC(ct, xx, residual2, mu_hat, Theta_hat, beta0,beta_hat, sigma_hat, L.mat=L.mat,pi_vec=prob)

          aicp[ltl]=tmp$aic

}
          }
          }
        }

      }
    }

    indtmp=which.min(aicp)


    if(length(indtmp)!=1){
      indtmp=ifelse(length(indtmp)==0,1,indtmp[1])
    }

    Opt_aAIC = aicp[indtmp]
    
    colnames(lam)=c('v_1','v_0','lambda_b','lambda_s','p_2','tau0')
    
    Opt_lambda_aic = lam[indtmp,]



    result_adap.bic = list(K=K, Opt_lambda_aic =Opt_lambda_aic,Mu_hat.list=Mu_hat.list,Theta_hat.list=Theta_hat.list,
                           L.mat.list=L.mat.list, Opt_aAIC=Opt_aAIC, Opt_num_aic=indtmp, aicp=aicp,
                           prob.list=prob.list,
                           beta.list= beta.list, sigma.list=sigma.list,
                           gamma=gamma,lam= lam,l.m=l.m,member_input=member_input)
    resultall[[lt]]=result_adap.bic

  }

  return(resultall)
}

