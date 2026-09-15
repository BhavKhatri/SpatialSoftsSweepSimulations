#install.packages("./abcWasserstein", repos = NULL, type="source")

library(R.matlab)
library(T4transport)
library(transport)
library(ptw)
library(abcWasserstein)
library(abc)
library(future.apply)
plan(multisession)  # or multiprocess on older R

# Simulated parameters and summary statistics.
input_10_222 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=10_nsample=222.mat")
input_10_2218 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=10_nsample=2218.mat")
input_10_22180 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=10_nsample=22180.mat")
input_50_222 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=50_nsample=222.mat")
input_50_2218 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=50_nsample=2218.mat")
input_50_22180 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=50_nsample=22180.mat")
input_100_222 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=100_nsample=222.mat")
input_100_2218 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=100_nsample=2218.mat")
input_100_22180 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=100_nsample=22180.mat")
input_500_222 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=500_nsample=222.mat")
input_500_2218 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=500_nsample=2218.mat")
input_500_22180 <- readMat("./GambiaeRange_RandomSampling_N=1.0e7_5.0e9_D=5.0_1000.0_s=0.01_1.0_Tmin=1940.0_Tmax=2000.0_nspatialsites=500_nsample=22180.mat")

datasets <- list(
  input_10_222 = input_10_222,
  input_10_2218 = input_10_2218,
  input_10_22180 = input_10_22180,
  input_50_222 = input_50_222,
  input_50_2218 = input_50_2218,
  input_50_22180 = input_50_22180,
  input_100_222 = input_100_222,
  input_100_2218 = input_100_2218, 
  input_100_22180 = input_100_22180,
  input_500_222 = input_500_222,
  input_500_2218 = input_500_2218,
  input_500_22180 = input_500_22180
)

# Loop through the different datasets.
counter=0
for (d in datasets) {
  
  counter=counter+1
  
  # Subset only the summary statistics.
  sust <- lapply(d[["n"]], function(x) x[[1]])
  length(sust)
  
  # Subset only the parameters.
  parameters <- as.data.frame(rbind(as.matrix(d$N),as.matrix(d$s),as.matrix(d$D),as.matrix(d$t)))
  parameters<- data.frame(t(parameters))
  colnames(parameters)<-c("N","s","D","t")

  # Create a storage data frame for results
  results <- data.frame(
    iteration = 1:500,
    orig_param1 = NA_real_,
    median_param1 = NA_real_,
    ci_low_param1 = NA_real_,
    ci_high_param1 = NA_real_,
    orig_param2 = NA_real_,
    median_param2 = NA_real_,
    ci_low_param2 = NA_real_,
    ci_high_param2 = NA_real_,
    orig_param3 = NA_real_,
    median_param3 = NA_real_,
    ci_low_param3 = NA_real_,
    ci_high_param3 = NA_real_,
    orig_param4 = NA_real_,
    median_param4 = NA_real_,
    ci_low_param4 = NA_real_,
    ci_high_param4 = NA_real_,
    orig_param5 = NA_real_,
    median_param5 = NA_real_,
    ci_low_param5 = NA_real_,
    ci_high_param5 = NA_real_
  )

  # Define a function to pad smaller matrices on the reight side, so that the Wasserstein distance can be computed.
  pad_right <- function(x, target_cols) {
    x <- as.matrix(x)
    nc <- NCOL(x)
    if (nc >= target_cols) return(x)
    cbind(x, matrix(0, nrow = NROW(x), ncol = target_cols - nc))
  }
  
  # Define a function to calculate the Wasserstein distance.
  wasserstein_pad <- function(x, y, p = 2) {
    tgt <- max(NCOL(x), NCOL(y))
    x2 <- pad_right(x, tgt)
    y2 <- pad_right(y, tgt)
    transport::wasserstein(transport::pp(x2), transport::pp(y2), p = 2)
  }

  # Loop through the datasets calculating the distances and using them to perform ABC.
  i <- 1
  for (p in 1:500) {

  dist <- future_sapply(sust, function(su, tgt) 
    wasserstein_pad(su, tgt, p = 2), tgt = sust[[p]])
  
    # Run abc for the i-th line of sust and parameters
    post.parest <- abc::abc(
      target = 0, 
      param = parameters[-p, ], 
      sumstat = dist[-p], 
      tol = 0.01,
      transf=c("none","none","none","none","none"), 
      prior.range=rbind(c(10000000,5000000000),c(1e-2,1.0),c(5,1000.0),c(150,900),c(8,900)),
      method = "rejection"
    )

    # Extract original values from parameters
    orig_vals <- parameters[p, ]
    
    # Compute confidence intervals (e.g., 95% CI)
    ci_vals <- apply(post.parest$unadj.values, 2, quantile, probs = c(0.025, 0.975))
    
    # Store results
    results$orig_param1[i] <- orig_vals[[1]]
    results$median_param1[i] <- median(post.parest$unadj.values[,1])
    results$ci_low_param1[i] <- ci_vals[1, 1]
    results$ci_high_param1[i] <- ci_vals[2, 1]
    
    results$orig_param2[i] <- orig_vals[[2]]
    results$median_param2[i] <- median(post.parest$unadj.values[,2])
    results$ci_low_param2[i] <- ci_vals[1, 2]
    results$ci_high_param2[i] <- ci_vals[2, 2]
    
    results$orig_param3[i] <- orig_vals[[3]]
    results$median_param3[i] <- median(post.parest$unadj.values[,3])
    results$ci_low_param3[i] <- ci_vals[1, 3]
    results$ci_high_param3[i] <- ci_vals[2, 3]
    
    results$orig_param4[i] <- orig_vals[[4]]
    results$median_param4[i] <- median(post.parest$unadj.values[,4])
    results$ci_low_param4[i] <- ci_vals[1, 4]
    results$ci_high_param4[i] <- ci_vals[2, 4]
    
    results$orig_param5[i] <- orig_vals[[5]]
    results$median_param5[i] <- median(post.parest$unadj.values[,5])
    results$ci_low_param5[i] <- ci_vals[1, 5]
    results$ci_high_param5[i] <- ci_vals[2, 5]
    
    i <- i+1
    
  }

# Save to output file.
write.csv(results,paste("SimpleRejection_CI_results_",names(datasets[counter]),".csv",sep=""))
