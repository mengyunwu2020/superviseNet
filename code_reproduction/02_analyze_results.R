# ==============================================================================
# Script: 02_analyze_results.R
# Purpose: Analyze simulation results


Esti.error = function(mu_hat, Theta_hat, beta_hat, Mu0, Theta0, beta0, K0, L0, prob,Linput,beta0hat){

  FailureK=FALSE
  p = dim(mu_hat)[2]
  K_hat = dim(mu_hat)[1]
  n_all = n=dim(Linput)[1]

  num = rep(0,K0)
  for (k in 1:K0) {
    error1 = apply((t(beta_hat) - beta0[k,])^2,2,sum)
    nums= which(error1 == min(error1))[1]
    num[k] = nums
  }
  if(length(unique(num))!=K_hat){
    errortmpp=matrix(0,K0,K0)
    for (k in 1:K0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      errork2=apply((t(beta_hat) - beta0[k,])^2,2,sum)
      errork =errork2+ errork
      errortmpp[k,]=errork
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }


    if(length(unique(num))!=K_hat){
      bb1=apply(errortmpp,2,which.min)
      sb=bb1
      bm=as.data.frame(table(bb1))
      tmp=bm[bm['Freq']==1,'bb1']
      bm=as.numeric(as.vector(tmp))
      bm=which(bb1==bm)

      bb1=setdiff(c(1:K0),bb1)
      for(kk in bb1){
        if(length(bm)==0){
          im=which(errortmpp[kk,]==min(errortmpp[kk,]))
          sb[kk]=im
        }else{
          im=which(errortmpp[kk,]==min(errortmpp[kk,-bm]))
          sb[kk]=im
        }
        bm=c(bm,im)
      }
      num=order(sb)

      FailureK=TRUE
    }
  }
  beta_hat=beta_hat[num,]
  mu_hat = mu_hat[num,]
  Theta_hat = Theta_hat[,,num]
  prob=prob[num]


  temp1=(Theta_hat - Theta0)^2
  temp1=as.numeric(apply(temp1, 3, sum))
  PME = sum(sqrt(temp1))/K0

  TPk = (apply((Theta_hat!=0) + (Theta0!=0) == 2,3,sum) - p) / (apply(Theta0!=0,3,sum) - p)
  TPk[(is.na(TPk))] = 1
  TPR_omega = sum(TPk) / K0
  FPk = apply((Theta_hat!=0) + (Theta0==0) == 2,3,sum) / apply(Theta0==0,3,sum)
  FPk[(is.na(FPk))] = 0
  FPR_omega = sum(FPk[(!is.na(FPk))]) / K0



  L.mat=Linput
  member = apply(L.mat,1,which.max)

  if(class(member)=='list'){
    member=as.numeric(member)
  }

  aa = L0
  cap_matrix0 = matrix(0,n_all,n_all)
  for(i in 1:(n_all-1)){
    for (j in (i+1):(n_all)) {
      cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
    }
  }
  aa = member
  cap_matrix1 = matrix(0,n_all,n_all)
  for(i in 1:(n_all-1)){
    for (j in (i+1):(n_all)) {
      cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
    }
  }
  CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)

  if(length(member)==0) member=rep(0,n_all)
  if(is.na(sum(member))) member[is.na(member)]=1

  index =c(K_hat,CE,PME,TPR_omega,FPR_omega)
  names(index) = c("K","CE","PME","TPR_omega","FPR_omega")

  return(index)
}





result_file_demo <- "output/sim_results_DEMO.RData"
result_file_full <- "output/sim_results_FULL.RData"

if(file.exists(result_file_demo)) {
  load(result_file_demo)
  message("Loaded results from: ", result_file_demo)
} else if(file.exists(result_file_full)) {
  load(result_file_full)
  message("Loaded results from: ", result_file_full)
} else {
  stop("No result file found in ../output/. Please run 01_simulation_main.R first.")
}

runtimes <- sapply(sim_results, function(x) {
  if(is.list(x) && !is.null(x$runtime)) return(x$runtime)
  else return(NA)
})

cat("\n--- Simulation Summary ---\n")
cat("Number of Replications:", length(sim_results), "\n")
cat("Average Runtime per Rep:", mean(runtimes, na.rm=TRUE), "seconds\n")

valid_indices <- which(!sapply(sim_results, function(x) all(is.na(x))))

if(length(valid_indices) > 0) {
  for(ii in 1:length(valid_indices)){
  idx <- valid_indices[ii]
  cat(paste0("\n--- Structure of the Result (",idx, "th Replication) ---\n"))
  ind=sim_results[[idx]]$fit[[1]]$Opt_num_bic




  PP=sim_results[[idx]]$fit[[1]]
  Theta_hat=PP$Theta_hat.list[[ind]]
  p=dim(Theta_hat)[1]
  K0=dim(Theta_hat)[3]

  mu_hat=PP$Mu_hat.list[[ind]]
  mmem=apply(PP$L.mat.list[[ind]],1,which.max)



  gamma_w=PP$gamma.list[[ind]]
  A.hat2=array(NA,dim=c(p,p,K0))
  for(k in 1:K0){
    A.hat2[,,k]=abs(gamma_w[[k]])>0.5
  }

  L.mat=matrix(0,dim(xx)[1],K0)
  for(il in 1:dim(xx)[1]){
    L.mat[il, mmem[il]]=1
  }

  prob = PP$prob.list[[ind]];
  residual2=PP$residual.list[[ind]];
  beta_hat=PP$beta.list[[ind]][,-1]
  beta0h=PP$beta.list[[ind]][,1]
  sigma_hat=PP$sigma.list[[ind]]


  res=Esti.error(mu_hat, Theta_hat*A.hat2, beta_hat,sim_results[[idx]]$whole.data$Mu0,sim_results[[idx]]$whole.data$Theta0,sim_results[[idx]]$whole.data$beta, K0, L0=whole.data$L0, prob,Linput=L.mat,beta0h)
  print(res)
}
} else {
  cat("\nWarning: No valid results found (all runs failed?).\n")
}
