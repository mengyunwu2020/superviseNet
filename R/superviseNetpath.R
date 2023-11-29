#' Supervised heterogeneous network estimation via survival-based Bayesian graphical models (superviseNet) along path.
#'
#' @param ct  A two-column matrix with the first column being the survival time and the second column being the censoring indicator. The indicator is a binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param xx Input matrix of p genetic measurements consisting of n rows. Each row is an observation vector.
#' @param Kseq A user supplied vector of number of subgroups.
#' @param lambda_mu A user supplied non-negative tuning parameters of laplace prior on subgroup mean parameters
#' @param v1 A user supplied slab prior parameters.
#' @param v0 A user supplied spike prior parameters.
#' @param lambda_beta A user supplied non-negative tuning parameters of laplace prior on subgroup regression parameters
#' @param lambda_sim A user supplied non-negative tuning parameters controlling the similarity across different networks.
#' @param p2 	A user supplied probability parameters of Bernoulli prior on binary latent indicator $gamma_{k,jl}$
#' @param laplace.m Laplace matrix
#' @param member_input  A vector indicating initialized subgroup memberships of each subjects.
#' @param eps Tolerance for the EM algorithm. The default value is 1e-3.
#' @param maxiter Maximum of number of iterations.
#' @param threshold A small constant that thresholds the final precision matrix estimator.
#' @param traces Whether to trace intermediate results.
#' @param eps_z A small regularized constant used in the risk ratio computation.
#'
#' @return A list
#' @import MASS
#' @importFrom igraph ba.game
#' @importFrom igraph get.edgelist
#' @export
#'
#' @examples
#' n <- 150
#' p <- 100
#' K0 <- 3
#' repli_result=list()
#' num.method=10
#' set.seed(1)
#' mue <-0
#' nonnum <- 2
#' mu03 <- c(rep(mue,nonnum),rep(-mue,nonnum),rep(-0,p-2*nonnum))#rep(0,p)
#' mu01 <- c(rep(-mue,2*nonnum),rep(-0,p-2*nonnum))
#' mu02 <- c(rep(mue,2*nonnum),rep(0,p-2*nonnum))
#' beta=matrix(0,K0,p)
#' beta[1,1:5]=2
#' beta[2,1:5]=-2
#' beta[3,3:7]=1
#' beta0=c(0,0,0)
#' sigma=c(.01,.01,0.01)
#' Mu0.list <- list(mu01,mu02,mu03)
#' set.seed(1)
#' whole.data <- generate.data(n,Mu0.list,beta,beta0,sigma)
#' whole.data$beta0=beta
#' xx=whole.data$data
#' ct=whole.data$ct
#' set.seed(4)
#' class_old<-sample(1:K0,nrow(xx),replace=TRUE)
#' Thetab <- array(0, dim = c(p, p, K0))
#' Sb <- array(0, dim = c(p, p, K0))
#' for(k in 1:K0)
#' {
#' Sb[,,k]  <- cov(xx[class_old == k, , drop = FALSE])
#' if(det(Sb[,,k])<1e-4) Sb[,,k]=Sb[,,k]+diag(0.01,p,p)
#' Thetab[,,k] <- solve(Sb[,,k])
#' }
#' tau1=0.001
#' laplace.m=NULL
#' for(i in 1:(p-1)){
#' for(j in (i+1):p){
#' tmm=1/sqrt(Thetab[i,j,]^2+tau1^2)
#' tmp.m=-matrix(tmm,ncol=1)%*%matrix(tmm,nrow=1)
#' diag(tmp.m)=(K0-1)*(tmm^2)
#' tmp=tmp.m+diag(tau1,K0,K0)
#' laplace.m=c(laplace.m,c(tmp))
#' }
#' }
#' res=superviseNetpath(ct,xx,K0,lambda_mu=.01,lambda_beta=0.06,v0=c(0.01,0.017),v1=1,p2=0.9,lambda_sim=0.005,laplace.m=laplace.m,member_input=class_old,eps=1e-2,maxiter=50,threshold=1e-3, traces=TRUE, eps_z=1e-5)
superviseNetpath= function(ct, xx, Kseq,lambda_mu=0,v1,v0,lambda_beta,lambda_sim,p2,laplace.m,
                                 member_input,eps =1e-2,
                                 maxiter=50,threshold=1e-3,
                                 traces=TRUE,eps_z=1e-5){


  L1 = length(v1)
  L2 = length(v0)
  L3 = length(lambda_beta)
  L4 = length(lambda_sim)
  L5 = length(p2)
  p=ncol(xx)

  LL = L1*L2*L3*L4*L5

  bic=rep(1e100,length(Kseq))
  lt=0
  resultall=vector('list',length(Kseq))
  oobic=oobic2=oobic3=rep(10e10,length(Kseq))
  for(K in Kseq){
    cat('the ',K,'th kseq is running','\n')
    lt=lt+1

    # gamma=list()
    Mu_hat.list = list()
    Theta_hat.list = list()
    prob.list = list()
    L.mat.list = list()
    member.list = list()
    beta.list=sigma.list=list()
    aBIC=aicp=ADBIC=rep(10^10, LL)
    lam=matrix(0,LL,5)
    ltl=0

    for(jj in 1:L1){
      for(ll in 1:L2){
        for(mm in 1:L3){
          for(ss in 1:L4){
            for(qq in 1:L5){

          ltl=ltl+1

          v_1=v1[jj]
          v_0=v0[ll]
          lambda_b=lambda_beta[mm]
          lambda_s=lambda_sim[ss]
          p_2=p2[qq]
          print(c(v_1,v_0,lambda_b,lambda_s,p_2))
          lam[ltl,]=c(v_1,v_0,lambda_b,lambda_s,p_2)

          cat('the ',ltl,'th lams is running','\n')
          print(lam[ltl,])

          PP=superviseNet(ct,xx, K=K,
                                v_0=v_0,
                                v_1=v_1,
                                maxiter=maxiter,
                                p_2=p_2,
                                lambda_b=lambda_b,lambda_mu=lambda_mu,
                                threshold=threshold,
                                traces=traces,eps=eps,member_input=member_input,eps_z=eps_z,lambda_s = lambda_s,laplace.m =laplace.m
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

          # gamma[[ltl]]=PP$gamma
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

    indtmp=which.min(aicp)


    if(length(indtmp)!=1){
      indtmp=ifelse(length(indtmp)==0,1,indtmp[1])
    }

    Opt_aAIC = aicp[indtmp]
    Opt_lambda_aic = lam[indtmp,]


    result_adap.bic = list(K=K, Opt_lambda_aic =Opt_lambda_aic,Mu_hat.list=Mu_hat.list,Theta_hat.list=Theta_hat.list,
                           L.mat.list=L.mat.list, Opt_aAIC=Opt_aAIC, Opt_num_aic=indtmp, aicp=aicp,
                           prob.list=prob.list,
                           beta.list= beta.list, sigma.list=sigma.list,
                           gamma=gamma,lam= lam,laplace.m=laplace.m,member_input=member_input)
    resultall[[lt]]=result_adap.bic

  }

  return(resultall)
}

