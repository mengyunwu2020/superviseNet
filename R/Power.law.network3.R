
Power.law.network3 = function(p,s=10,umin=0.4,umax=0.7,num.differ,Sparse=TRUE,m,Dnum=3,rate,option='ne',N,K0){



  pp=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  submatrix=submatrix2=submatrix3=list()


  for (ss in 1:num.differ) {
    g = ba.game(pp, m=m,directed = FALSE)
    Eg= as.data.frame(get.edgelist(g))
    subi = diag(1,pp,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }


    submatrix[[ss]]=subi


    g = ba.game(pp, m=m,directed = FALSE)
    Eg = as.data.frame(get.edgelist(g))
    subi = diag(1,pp,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
    submatrix2[[ss]] = subi

    g = ba.game(pp, m=m,directed = FALSE)
    Eg = as.data.frame(get.edgelist(g))
    subi = diag(1,pp,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
    submatrix3[[ss]] = subi
  }
  if(num.differ!=(s*rate)){
    if(!Sparse){


      for (ss in (num.differ+1):(s*rate)) {
        g = ba.game(pp, m=m,directed = FALSE)
        Eg= as.data.frame(get.edgelist(g))
        subi = diag(1,pp,pp)
        for (q in 1:dim(Eg)[1]) {
          i=Eg[q,1];j=Eg[q,2]
          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi[i,j]=ij;subi[j,i]=ij
        }

        submatrix[[ss]]=submatrix2[[ss]]=subi
      }

    }else{
      for (ss in (num.differ+1):(s*rate)) {
        g = ba.game(pp, m=m,directed = FALSE)
        Eg= as.data.frame(get.edgelist(g))
        subi = subi2 = subi3=diag(1,pp,pp)
        for (q in 1:dim(Eg)[1]) {
          i=Eg[q,1];j=Eg[q,2]
          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi[i,j]=ij;subi[j,i]=ij

          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi2[i,j]=ij;subi2[j,i]=ij
          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi3[i,j]=ij;subi3[j,i]=ij
        }


        submatrix[[ss]]=subi

        subi=subi2


        submatrix2[[ss]]=subi

        submatrix3[[ss]]=subi3

      }


    }
  }

  if(ss==1){
    A=submatrix[[1]]
  }else{
    A=submatrix[[1]]
    for (ss in 2:(s*rate)) {
      A=bdiag(A,submatrix[[ss]])
    }
    A=bdiag(A,diag(1,(p*(1-rate)),(p*(1-rate))))
  }
  A = as.matrix(A)

  if(ss==1){
    A2=submatrix2[[1]]
  }else{
    A2=submatrix2[[1]]
    for (ss in 2:(s*rate)) {
      A2=bdiag(A2,submatrix2[[ss]])
    }
    A2=bdiag(A2,diag(1,(p*(1-rate)),(p*(1-rate))))
  }
  A2 = as.matrix(A2)


  if(ss==1){
    A3=submatrix3[[1]]
  }else{
    A3=submatrix3[[1]]
    for (ss in 2:(s*rate)) {
      A3=bdiag(A3,submatrix3[[ss]])
    }
    A3=bdiag(A3,diag(1,(p*(1-rate)),(p*(1-rate))))
  }
  A3 = as.matrix(A3)


  if(option=='ne'){
    ds <- rowSums(abs(A)) * 0.8
    diag(A) <- ds
    for (i in 1:p) {
      for (j in 1:p){
        A[i, j] <- A[i, j] / sqrt(ds[i]) / sqrt(ds[j])
      }
    }
  }
  if(option=='ne'){
    ds <- rowSums(abs(A2)) * 0.8
    diag(A2) <- ds
    for (i in 1:p) {
      for (j in 1:p){
        A2[i, j] <- A2[i, j] / sqrt(ds[i]) / sqrt(ds[j])
      }
    }
  }

  if(option=='ne'){
    ds <- rowSums(abs(A3)) * 0.8
    diag(A3) <- ds
    for (i in 1:p) {
      for (j in 1:p){
        A3[i, j] <- A3[i, j] / sqrt(ds[i]) / sqrt(ds[j])
      }
    }
  }






  return(list(A1=A,A2=A2,A3=A3))
}
