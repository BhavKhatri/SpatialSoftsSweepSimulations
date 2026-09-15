library(abc)

# Simulated parameters and summary statistics.
input_10_2218 <- read.csv("../../../nsites=10_nsample=2218.csv")
input_10_22180 <- read.csv("../../../nsites=10_nsample=22180.csv")
input_50_2218 <- read.csv("../../../nsites=50_nsample=2218.csv")
input_50_22180 <- read.csv("../../../nsites=50_nsample=22180.csv")
input_100_2218 <- read.csv("../../../nsites=100_nsample=2218.csv")
input_100_22180 <- read.csv("../../../nsites=100_nsample=22180.csv")
input_500_2218 <- read.csv("../../../nsites=500_nsample=2218.csv")
input_500_22180 <- read.csv("../../../nsites=500_nsample=22180.csv")

datasets <- list(
  input_10_2218 = input_10_2218,
  input_10_22180 = input_10_22180,
  input_50_2218 = input_50_2218,
  input_50_22180 = input_50_22180,
  input_100_2218 = input_100_2218, 
  input_100_22180 = input_100_22180,
  input_500_2218 = input_500_2218,
  input_500_22180 = input_500_22180
)

# Loop through the different datasets.
counter=0
for (d in datasets) {
  
  # Subset only the parameters.
  parameters<-d[c(1:4)]
  counter=counter+1
  parameters["s_t"]=parameters$s*parameters$t
  
  # Subset only the summary statistics.
  sust<-d[c(6,9:12)]

  # Identify rows with NaN values
  parameters <- parameters[is.finite(rowSums(sust)),]
  sust <- sust[is.finite(rowSums(sust)),]
  

  # Perform 100 iterations of cross-validation with different thresholds to select the value with the highest accuracy.
  cv.parest <- cv4abc(parameters, sust, nval=100, tol =c(.05,.01,.005), transf=c("none","none","none","none","none"), prior.range=rbind(c(10000000,5000000000),c(1e-2,1.0),c(5,1000.0),c(150,900),c(8,900)), method = "rejection")
  print(names(datasets[counter]))
  print(summary(cv.parest))
  

  # Create a storage data frame for results
  results <- data.frame(
    iteration = 1:nrow(sust),
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
  
  for (i in 1:nrow(sust)) {
    # Run abc for the i-th line of sust and parameters
    post.parest <- abc(
      target = sust[i, ], 
      param = parameters[-i, ], 
      sumstat = sust[-i, ], 
      tol = 0.005,
      transf=c("none","none","none","none","none"), 
      prior.range=rbind(c(10000000,5000000000),c(1e-2,1.0),c(5,1000.0),c(150,900),c(8,900)),
      method = "rejection"
    )
    
    # Extract original values from parameters
    orig_vals <- parameters[i, ]
    
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
    
  }
  # Save to output file.
  write.csv(results,paste("SimpleRejection_CI_results_",names(datasets[counter]),".csv",sep=""))
}