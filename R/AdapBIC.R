AdapBIC = function(ct, data, residual2, mu_hat, Theta_hat,beta0, beta_hat, sigma, L.mat=NULL,pi_vec=NULL){

  n = nrow(data)
  K = nrow(mu_hat)

  logcensortime = log(ct[,1])
  if(!(is.null(L.mat))){
    pi_vec = apply(L.mat, 2, sum)/n
  }

  fit.error_mat = matrix(0, n, K)

  f.mat= fy.mat=matrix(0,n,K)
  for(k in 1:K) {
    f.mat[,k]=f.den.vec(data, as.numeric(mu_hat[k,]),Theta_hat[,,k])
    fy.mat[,k]=(1/(2*pi*sigma[k]))^(0.5)*exp(-residual2[k,]/(2*sigma[k]))
    for(ii in 1:n){
      if(ct[ii,2]==0){
      mu_tmp=sum(data[ii,]*beta_hat[k,])+beta0[k]
      sd_tmp=sqrt(sigma[k])
      fy.mat[ii,k]=1-pnorm(logcensortime[ii],mu_tmp,sd_tmp)
      }
    }
    fit.error_mat[,k] = pi_vec[k] * f.mat[,k] * fy.mat[,k]
  }
  fit0 = apply(fit.error_mat, 1, sum)
  if(length(which(fit0==0))>0){
  fit.error = sum(log( fit0 + min(fit0[fit0>0]) ))
  }else{
    fit.error = sum(log( fit0))
  }
  fit.error = - 2*fit.error

  # degrees of freedom
  for(i in 1:K){
    Theta_hat[upper.tri(Theta_hat[, , i], diag = T)] = 0
  }
  p=dim(data)[2]

  dfomega =  length(which(Theta_hat != 0))
  dfmu = length(which(mu_hat!=0))
  dfbeta = length(which(beta_hat!=0))


  bic = fit.error + log(n) * (dfomega+dfmu+dfbeta)


  P = list()
  P$fit.error = fit.error

  P$bic=bic
  P$f.mat=f.mat
  return(P)
}
