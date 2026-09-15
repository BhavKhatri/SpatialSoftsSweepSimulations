library(abc)

# Simulated parameters and summary statistics.
input <- read.csv("../../../SampleSummaryStats_Area=GambiaeRange_mu=2e-09_M=9600.csv")

datasets <- list(
  input = input
)

# Subset only the parameters.
parameters<-input[c(1:4)]

# Subset only the summary statistics.
sust<-input[c(6,9:12)]
  
# Identify rows with NaN values
parameters <- parameters[is.finite(rowSums(sust)),]
sust <- sust[is.finite(rowSums(sust)),]
  
# Vector of empirical summary statistics.
emp<-c(10,259,0.295,0.9780,0.629)


# Run abc with a neural network approach.
post.parest <- abc(
      target = emp, 
      param = parameters[1:5000,], 
      sumstat = sust[1:5000,], 
      tol = .05,
      transf=c("log","log","log","log","log"), 
      prior.range = rbind(c(10000000,5000000000),c(1e-2,1.0),c(5,1000.0),c(150,900),c(8,900)), 
      method = "neuralnet"
    )

# Print parameter estimates.
summary(post.parest)

# Save the parameter estimates.
write.csv(post.parest$adj.values,"NN_posterior_5K.csv")





