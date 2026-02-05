# ==============================================================================
# Script: 01_simulation_main.R
# Purpose: Run simulation studies for superviseNet
# ==============================================================================

library(superviseNet)

# ------------------------------------------------------------------------------
# 1. Setting
# ------------------------------------------------------------------------------
RUN_MODE <- "DEMO"
# RUN_MODE <- "FULL"

n <- 150
p <- 100
K <- 2
maxiter <- 50

if (RUN_MODE == "DEMO") {
  message(sprintf(">>> Running in DEMO mode (n=%d, p=%d, maxiter=%d) <<<", n, p, maxiter))
  message("Note: Running only 1 replications for quick verification.")
  n_rep <- 1
} else {
  message(sprintf(">>> Running in FULL mode (n=%d, p=%d, maxiter=%d) <<<", n, p, maxiter))
  n_rep <- 100
}

# ------------------------------------------------------------------------------
# 2. Results
# ------------------------------------------------------------------------------
if(!dir.exists("output")) dir.create("output")
sim_results <- list()

# ------------------------------------------------------------------------------
# 3. Simultion
# ------------------------------------------------------------------------------
for (r in 1:n_rep) {
  set.seed(r)
  cat(sprintf("\n[Repetition %d / %d] Initializing...\n", r, n_rep))

  mu01 = mu02 = rep(0, p)
  Mu0.list <- list(mu01, mu02)

  beta = matrix(0, K, p)
  beta[1, 1:5] = 2
  beta[2, 1:5] = -2

  beta0 = c(0, 0)
  sigma = c(0.01, 0.01)

  whole.data <- superviseNet::generate.data(n, Mu0.list, beta, beta0, sigma)

  xx = whole.data$data
  ct = whole.data$ct
  whole.data$beta=beta
  whole.data$beta0=beta0


  lambda_pro=list()
  lambda_pro$v_0=.045#seq(0.04,.05,length.out=3)
  lambda_pro$lambda_sim=c(.01)
  init=init_superviseNet(
    ct,
    xx,
    K,
    lambda_mu=sqrt(dim(ct)[1]*log(p))/2,
    lambda_b=sqrt(dim(ct)[1]*log(p))/2,
    tau_0=0.01,
    v_0=lambda_pro$v_0[1],
    v_1=1,
    p_2=0.85,
    lambda_s=lambda_pro$lambda_sim[1],
    l.m,
    member_input,
    eps=1e-2,
    maxiter=20,
    threshold=1e-3,
    eps_z=1e-5,
    l.update=TRUE,tau1=0.001,init_time=2,thres_gamma=0.5)


  class_old=init$class_old
  l.m=init$laplace.m


  cat(sprintf("[Repetition %d / %d] Running superviseNet algorithm...\n", r, n_rep))
  start_time <- Sys.time()

  tryCatch({
    fit<-superviseNet::superviseNetpath(ct,xx,K,lambda_mu=sqrt(dim(ct)[1]*log(p))/2,
                                        lambda_beta=sqrt(dim(ct)[1]*log(p))/2,
                                        tau_0=0.01,v0=lambda_pro$v_0,
                                        v1=1,p2=0.85,lambda_sim=lambda_pro$lambda_sim,
                                        l.m=l.m,member_input=class_old,eps=1e-2,
                                        maxiter=maxiter,threshold=1e-3,eps_z=1e-5)

    end_time <- Sys.time()
    fit$runtime <- as.numeric(end_time - start_time, units = "secs")

    sim_results[[r]] <-list(fit=fit,whole.data=whole.data)
    cat(sprintf("   Done. Runtime: %.2f secs.\n", fit$runtime))

  }, error = function(e) {
    message("   Error in this run: ", e$message)
    sim_results[[r]] <- NA
  })
}

# ------------------------------------------------------------------------------
# 4. Save results
# ------------------------------------------------------------------------------
output_filename <- paste0("output/sim_results_", RUN_MODE, ".RData")
save(sim_results, file = output_filename)

message("\n=======================================================")
message(sprintf("Simulation (%s mode) All Done!", RUN_MODE))
message("Results saved to: ", output_filename)
message("=======================================================")
