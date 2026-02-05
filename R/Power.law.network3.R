




Power.law.network1 = function(p,s=10,umin1=0.4,umax1=0.7,umin2=0.4,umax2=0.7,num.differ,Sparse=TRUE,m,Dnum=3,rate){

  pp=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  submatrix=submatrix2=list()


  for (ss in 1:num.differ) {
    g = ba.game(pp, m=m,directed = FALSE)
    Eg= as.data.frame(get.edgelist(g))
    subi = diag(1,pp,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin1,umax1),runif(1,-umax1,-umin1)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }

    submatrix[[ss]]=subi


    g = ba.game(pp, m=m,directed = FALSE)
    Eg = as.data.frame(get.edgelist(g))
    subi = diag(1,pp,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin1,umax1),runif(1,-umax1,-umin1)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }


    submatrix2[[ss]] = subi


  }



  if(num.differ!=(s*rate)){
    if(!Sparse){


      for (ss in (num.differ+1):(s*rate)) {
        g = ba.game(pp, m=m,directed = FALSE)
        Eg= as.data.frame(get.edgelist(g))
        subi = diag(1,pp,pp)
        for (q in 1:dim(Eg)[1]) {
          i=Eg[q,1];j=Eg[q,2]
          ij = sample(c(runif(1,umin1,umax1),runif(1,-umax1,-umin1)))[1]
          subi[i,j]=ij;subi[j,i]=ij
        }
        submatrix[[ss]]=submatrix2[[ss]]=subi
      }

    }else{
      for (ss in (num.differ+1):(s*rate)) {
        g = ba.game(pp, m=m,directed = FALSE)
        Eg= as.data.frame(get.edgelist(g))
        subi = subi2 = diag(1,pp,pp)
        for (q in 1:dim(Eg)[1]) {
          i=Eg[q,1];j=Eg[q,2]
          ij = sample(c(runif(1,umin1,umax1),runif(1,-umax1,-umin1)))[1]
          subi[i,j]=ij;subi[j,i]=ij

          ij = runif(1,umin2,umax2)
          subi2[i,j]=abs(ij)*sign(subi[j,i]);subi2[j,i]=abs(ij)*sign(subi[j,i])

        }

        submatrix[[ss]]=subi

        subi=subi2




        submatrix2[[ss]]=subi

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

    diag(A)=0
    tmp=0.01+abs(min(eigen(A)$values))
    A=A+diag(tmp,p,p)


    diag(A2)=0
    tmp=0.01+abs(min(eigen(A2)$values))
    A2=A2+diag(tmp,p,p)

  return(list(A1=A,A2=A2))
}

Power.law.network3 = function(p,s=10,umin1=0.4,umax1=0.7,umin2=0.4,umax2=0.7,umin3=0.4,umax3=0.7,num.differ,Sparse=TRUE,m,Dnum=3,rate){



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
      ij = runif(1,umin1,umax1)
      subi[i,j]=ij;subi[j,i]=ij
    }


    submatrix[[ss]]=subi


    g = ba.game(pp, m=m,directed = FALSE)
    Eg = as.data.frame(get.edgelist(g))
    subi = diag(1,pp,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = runif(1,umin2,umax2)
      subi[i,j]=ij;subi[j,i]=ij
    }


    submatrix2[[ss]] =subi






    g = ba.game(pp, m=m,directed = FALSE)
    Eg = as.data.frame(get.edgelist(g))
    subi = diag(1,pp,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = runif(1,umin3,umax3)
      subi[i,j]=ij;subi[j,i]=ij
    }

    submatrix3[[ss]] = subi

    tmp_sign=matrix(0,pp,pp)
    tmp=sample(c(1,-1),pp*(pp-1)/2,replace = T)
    tmp_sign[upper.tri(tmp_sign)]=tmp
    tmp_sign=(t(tmp_sign)+tmp_sign)


    submatrix[[ss]]=submatrix[[ss]]*tmp_sign

    submatrix2[[ss]]=submatrix2[[ss]]*tmp_sign
    submatrix3[[ss]]=submatrix3[[ss]]*tmp_sign
  }
  if(num.differ!=(s*rate)){
    if(!Sparse){


      for (ss in (num.differ+1):(s*rate)) {
        g = ba.game(pp, m=m,directed = FALSE)
        Eg= as.data.frame(get.edgelist(g))
        subi = diag(1,pp,pp)
        for (q in 1:dim(Eg)[1]) {
          i=Eg[q,1];j=Eg[q,2]
          ij = sample(c(runif(1,umin1,umax1),runif(1,-umax1,-umin1)))[1]
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
          ij = sample(c(runif(1,umin1,umax1),runif(1,-umax1,-umin1)))[1]
          subi[i,j]=ij;subi[j,i]=ij

          ij = (runif(1,umin2,umax2))
          subi2[i,j]=abs(ij)*sign(subi[j,i]);subi2[j,i]=abs(ij)*sign(subi[j,i])
          ij = (runif(1,umin3,umax3))
          subi3[i,j]=abs(ij)*sign(subi[j,i]);subi3[j,i]=abs(ij)*sign(subi[j,i])
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
  diag(A)=1
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
  diag(A2)=1

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
  diag(A3)=1


    diag(A)=0
    tmp=0.01+abs(min(eigen(A)$values))
    A=A+diag(tmp,p,p)


    diag(A2)=0
    tmp=0.01+abs(min(eigen(A2)$values))
    A2=A2+diag(tmp,p,p)


    diag(A3)=0
    tmp=0.01+abs(min(eigen(A3)$values))
    A3=A3+diag(tmp,p,p)

  return(list(A1=A,A2=A2,A3=A3))
}
