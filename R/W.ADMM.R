soft <-function(a,lam11){
  out <- sign(a)*max(0, abs(a)-lam11)
  return(out)
}




softfun<-function(a,lam,rho,laplace.m,zz,max.iter=100,eps=1e-2,tau_0){
  p=dim(a)[1]
  K=dim(a)[3]

  out=zz
  s=0
  for(i in 1:(p-1)){
    for(j in (i+1):p){

      laptmp=matrix(laplace.m[(s+1):(s+K^2)],K,K)

      lap=diag(rho,K,K)+laptmp

      con=diag(lap)


      iter=0
      while(iter<max.iter){
        zs=out[i,j,]


        for(k in 1:K){

          tmp=(-sum(lap[k,-k]*out[i,j,-k])+rho*a[i,j,k])[1]
          if(is.infinite(tmp)){
            tmp=sign(sum(lap[k,-k]*out[i,j,-k]))*2^1023
          }
          out[i,j,k]=soft(tmp[1],lam[i,j,k])/con[k]
        }

        diff=sum((zs-out[i,j,])^2)
        if(diff<eps) break;
        iter=iter+1
      }
      for(k in 1:K)
        out[j,i,k]=out[i,j,k]
      s=s+K^2
    }
  }

  for(k in 1:K){
    diag(out[,,k]) <-diag(a[,,k])+tau_0/rho
  }

  return(out)

}


expand = function(A, rho, n){
  A[which(is.infinite(A))]=2^1023*sign(A[which(is.infinite(A))])
  edecomp <- eigen(A)
  D <- edecomp$values
  U <- edecomp$vectors
  D2 <- 0.5*(D + sqrt(D^2 + 4*n/rho ))

  D2[which(is.infinite(D2))]=2^1023*sign(D2[which(is.infinite(D2))])
  Omega <- U %*% diag(D2) %*% t(U)
  Omega
}


W.ADMM<-function(S, lam,  num, laplace.m,
                 epsilon = 1e-5, maxiter = 500, rho = 0.1,
                 rho.incr = 1.2,rho.max = 1e10,K,tau_0){


  p = dim(S[[1]])[1]
  W1 = Z1 = diag(1,p,p)
  A1 = matrix(0, p, p)

  W = Z = A = array(NA,c(p,p,K))
  for(k in 1:K){
    W[,,k]=Z[,,k]=W1
    A[,,k]=A1
  }

  diff_vec=NULL
  for (iii in 1:maxiter){

    W_prev = W
    Z_prev = Z

    temp=array(0,c(p,p,K))
    for(k in 1:K){

      W[,,k] = expand((-num[k]*S[[k]]/rho-A[,,k]/rho+Z[,,k]), rho, num[k])

      temp[,,k] = 1/rho*A[,,k]+W[,,k]
    }

    Z = softfun(temp,lam,rho,laplace.m,Z_prev,max.iter=100,eps=1e-2,tau_0)



    A  =  A + rho*(W - Z)


    diff_value_1 =  sum(abs(Z - Z_prev))
    diff_value_2 =  sum(abs(W - W_prev))
    diff_value_3 = sum(abs(W - Z))
    norm_value  =  sum(abs(Z))
    if (max(diff_value_1, diff_value_2, diff_value_3) <= norm_value*epsilon){
      break
    }
    diff_vec=c(diff_vec,max(diff_value_1, diff_value_2, diff_value_3))

    rho = min(rho*rho.incr,rho.max)


  }


  out = list(W= Z, W.1 = W, diff_value_1=diff_value_1,iii=iii)



}

