library(ggplot2)
library(gridExtra)
library(coda)

###### Importing Data ####
data <- airquality
# Only retain the variables 
ozone<- cbind(data[,1:4])


# delete casewise missing values

ozone <- na.omit(ozone)

#Center Variables
Ozone <- ozone$Ozone
Solar.R <- ozone$Solar.R-mean(ozone$Solar.R)
Wind <- ozone$Wind - mean(ozone$Wind)
Temp <- ozone$Temp - mean(ozone$Temp)

### Define Variables
teller <- 0
N <- nrow(ozone)
# no. of iterations
XG <-10000
#no. of chains
nchains <- 2
# list to store coefficient estimates, residuals and predicted values
mcmc <- list(matrix(ncol = 4, nrow = XG, dimnames = list(c(),c("beta0","beta1", "beta2", "sigma"))), matrix(ncol = 4, nrow= XG, dimnames = list(c(),c("beta0","beta1", "beta2", "sigma"))))
predicted <- list(matrix(ncol = 111, nrow=XG), matrix(ncol=111, nrow=XG))
residuals <- list(matrix(ncol = 111, nrow=XG), matrix(ncol=111, nrow=XG))

#### MCMC Algorithm ####

set.seed(1234)
for (c in 1:nchains){
  #set intial values
  mcmc[[c]][1,1] <- rnorm(1,0,1)
  mcmc[[c]][1,2] <- rnorm(1,0,1)
  mcmc[[c]][1,3] <- rnorm(1,0,1)
  mcmc[[c]][1,4] <- runif(1,1,10)
  
  #Compute intial residuals
  predicted[[c]][1,] <- mcmc[[c]][1,1]+mcmc[[c]][1,2]*Solar.R+mcmc[[c]][1,3]*Wind
  residuals[[c]][1,] <- Ozone - predicted[[c]][1,]
  
 
 ### Gibbs for beta0 
  
for (i in 2:XG){  
  ## calculate mean 
  muBeta0 <-  (sum(Ozone-(mcmc[[c]][i-1,2]*Solar.R)-(mcmc[[c]][i-1,3]*Wind))/mcmc[[c]][i-1,4]^2 + 0/1000)/(N/(mcmc[[c]][i-1,4]^2) + 1/1000)
  ##calculate variance
  varBeta0 <- 1/((N/mcmc[[c]][i-1,4]^2) + 1/1000)
  mcmc[[c]][i,1]  <- rnorm(1, mean = muBeta0, sd = sqrt(varBeta0))
  
  ## Gibbs sampler for Beta 1
  muBeta1 <- (sum(Solar.R*(Ozone-mcmc[[c]][i,1]-mcmc[[c]][i-1,3]*Wind)) /mcmc[[c]][i-1,4]^2+0/1000)/(sum(Solar.R^2)/mcmc[[c]][i-1, 4]^2 + 0.13/5) ### last part indicates the prior
  varBeta1 <- 1/(sum(Solar.R^2)/mcmc[[c]][i-1, 4]^2 + 1/5)
  mcmc[[c]][i,2] <- rnorm(1, muBeta1, sqrt(varBeta1))

  
  ### Metropolis Hastings sampler
  propfunc <- rnorm(1, mcmc[[c]][i-1,3],1.35)
  # draw a u
  u <- runif(1,0,1)
  #predict scores with current beta1 coefficient
  predCurrent<- mcmc[[c]][i,1]+ mcmc[[c]][i,2]*Solar.R+propfunc*Wind
  #predict scores with  beta1 from the previous iteration coefficient
  predPrevious <- mcmc[[c]][i,1]+mcmc[[c]][i,2]*Solar.R+mcmc[[c]][i-1,3]*Wind
  # calculate Acceptance probability
  r <- exp((sum(dnorm(Ozone, mean = predCurrent, sd = mcmc[[c]][i-1,4], log = T)) +
              dunif(x = propfunc, min = -100, max = 100, log = T)) -
             (sum(dnorm(Ozone, mean = predPrevious, sd = mcmc[[c]][i-1,4], log = T)) +
                dunif(x = mcmc[[c]][i-1,3], min = -100, max = 100, log = T)))
  # decide which coefficient to keep
  if (r > u) {
    mcmc[[c]][i,3] <- propfunc
    ## teller used to calculate accpetacne rate
    teller <- teller+1
  } else {
    mcmc[[c]][i,3] <- mcmc[[c]][i-1,3]
  }
  # Gibbs sampler for beta2
  # derive mu
  # muBeta2 <- (sum(Wind*(Ozone-mcmc[[c]][i,1]-mcmc[[c]][i,2]*Solar.R))/mcmc[[c]][i-1,4]+0/1000)/(sum(Wind^2)/mcmc[[c]][i-1,4] + 1/1000)
  # derive the variance
  # varBeta2 <- 1/(sum(Wind^2)/mcmc[[c]][i-1, 4] + 1/1000)
  # draw a sample 
  # mcmc[[c]][i,3] <- rnorm(1, muBeta2, varBeta2)
  
  
  # Gibbs sampler for Sigma
  alpha <- N/2 + 0.001
  beta <- sum((Ozone-(mcmc[[c]][i,1]+mcmc[[c]][i,2]*Solar.R+mcmc[[c]][i,3]*Wind))^2)/2 + 0.001
  
  # Square root is taken to get sigma instead of variance
  mcmc[[c]][i,4] <- sqrt(1/rgamma(1, shape = alpha, rate = beta))
  
  #### Residuals & Predicted Values ###
  predicted[[c]][i,] <- mcmc[[c]][i,1]+mcmc[[c]][i,2]*Solar.R+mcmc[[c]][i,3]*Wind
  residuals[[c]][i,] <- Ozone - predicted[[c]][i,] 
  }
}



### Burin in period ### Remove first 1000 iterations
mcmc.results <- lapply(mcmc, function(x){
  x[-(1:1000),]
})


#### Traceplot ####

chain1<- as.data.frame(mcmc.results[[1]])
chain2<- as.data.frame(mcmc.results[[2]])
library(ggplot2)

plot.0 <-  ggplot(data =chain1, aes(x=1001:10000, y=beta0))+
  geom_line(color="red")+
  geom_line(data = chain2,aes(x=1001:10000, y=beta0),color= "blue") + xlab("Iteration")+ ylab("Estimated Beta 0")

plot.1 <-  ggplot(data =chain1, aes(x=1001:10000, y=beta1))+
  geom_line(color="red")+
  geom_line(data = chain2,aes(x=1001:10000, y=beta1),color= "blue")+ xlab("Iteration")+ ylab("Estimated Beta 1")

plot.2 <-  ggplot(data =chain1, aes(x=1001:10000, y=beta2))+
  geom_line(color="red")+
  geom_line(data = chain2,aes(x=1001:10000, y=beta2),color= "blue")+ xlab("Iteration")+ ylab("Estimated Beta 2")

plot.sigma <-  ggplot(data =chain1, aes(x=1001:10000, y=sigma))+
  geom_line(color="red")+
  geom_line(data = chain2,aes(x=1001:10000, y=sigma),color= "blue")+ xlab("Iteration")+ ylab("Estimated Sigma")

grid.arrange(plot.0, plot.1, plot.2, plot.sigma, name="Traceplots")

### Autocorraltion ####

### Calculate Autocorrelation
autocor <- matrix(nrow =31, ncol=4)
for (l in 1:4){
  for (i in 1:31){
  
    autocor[i,l]<-mean(cor(chain1[i:9000,l], chain1[1:(9000-i+1),l]), cor(chain2[i:9000,l], chain2[1:(9000-i+1),l]))
  } 
}

acceptrate<-teller/(2*XG)
autocor <- as.data.frame(autocor)

### Plot Autocorrelation
a.plot.0 <- ggplot(autocor, aes(x = 0:30, y=autocor[,1])) + geom_point(size=1)+geom_segment(aes(x = 0:30, y = 0,xend= 0:30, yend= autocor[,1]))+ xlab("Lag - Beta0")+ ylab("Autocorrelation")
a.plot.1 <- ggplot(autocor, aes(x = 0:30, y=autocor[,2])) + geom_point(size=1)+geom_segment(aes(x = 0:30, y = 0,xend= 0:30, yend= autocor[,2]))+ xlab("Lag - Beta1")+ ylab("Autocorrelation")
a.plot.2 <- ggplot(autocor, aes(x = 0:30, y=autocor[,3])) + geom_point(size=1)+geom_segment(aes(x = 0:30, y = 0,xend= 0:30, yend= autocor[,3]))+ xlab("Lag - Beta2")+ ylab("Autocorrelation")
a.plot.sigma <- ggplot(autocor, aes(x = 0:30, y=autocor[,4])) + geom_point(size=1)+geom_segment(aes(x = 0:30, y = 0,xend= 0:30, yend= autocor[,4]))+ xlab("Lag - Sigma")+ ylab("Autocorrelation")

## Assesing Convergence

grid.arrange(a.plot.0, a.plot.1, a.plot.2, a.plot.sigma)




### Gelman Rubin Statistic ####
#### function to plot gelman rubin plot

gel.plot <- function(mcmc){
  
  gel.plot.stat <- matrix(ncol = ncol(mcmc[[1]]), nrow = (nrow(mcmc[[1]])-1000)/10)
  for (i in 1:(((nrow(mcmc[[1]])-1000)/10))){
    mcmc.short <- lapply(mcmc, function(x){x[1:(1000+i*10),]})
    gel.plot.stat[i,] <-gelman(mcmc.short)
  }
  
  par(mfrow=c(2,2))
  labels <- c("beta0","beta1", "beta2", "sigma")
  gel.plot.stat <- as.data.frame(gel.plot.stat)
  A <- ggplot(gel.plot.stat, aes(x=seq(1001, 10000, by=10 ), y=gel.plot.stat[,1])) + geom_line(color="black")+ xlab("Iteration")+ ylab("GR-Statistic") + ggtitle("Beta0") + geom_hline(yintercept = 1, color="red") 
  B <- ggplot(gel.plot.stat, aes(x=seq(1001, 10000, by=10 ), y=gel.plot.stat[,2])) + geom_line(color="black")+ xlab("Iteration")+ ylab("GR-Statistic") + ggtitle("Beta1") + geom_hline(yintercept = 1, color="red")
  C <- ggplot(gel.plot.stat, aes(x=seq(1001, 10000, by=10 ), y=gel.plot.stat[,3])) + geom_line(color="black")+ xlab("Iteration")+ ylab("GR-Statistic") + ggtitle("Beta2") + geom_hline(yintercept = 1, color="red")
  D <- ggplot(gel.plot.stat, aes(x=seq(1001, 10000, by=10 ), y=gel.plot.stat[,4])) + geom_line(color="black")+ xlab("Iteration")+ ylab("GR-Statistic") + ggtitle("Sigma") + geom_hline(yintercept = 1, color="red")
  grid.arrange(A,B,C,D)
}

#### function to calculate GR statistic
gelman <- function(mcmclist){
  gel.stat <- c()
  for(i in 1:(ncol(mcmclist[[1]]))){
  gel.dat<-rbind(mcmclist[[1]], mcmclist[[2]])
  GM <- mean(gel.dat[,i])

  theta1 <- mean(mcmclist[[1]][,i])
  theta2 <- mean(mcmclist[[2]][,i])

  M <- 2

  N <- nrow(mcmclist[[1]])
  B <- (N/M-1)*((theta1-GM)^2+ (theta2-GM)^2)
  W <-   sum((mcmclist[[1]][,i]-theta1)^2, (mcmclist[[2]][,i]-theta2)^2 )/(2*(N-1))

   gel.stat[i]<-((((XG-1)/XG) * W) + (B/XG)) /W
  }
  
  return(gel.stat)
}
### call the functions
gel.plot(mcmc)



#### Sumary of Statistics ####
sum.res <- rbind(chain1, chain2)
sum.m <- matrix(ncol= 4, nrow = ncol(chain1), dimnames = list(c("beta0", "beta1", "beta2", "sigma"), c("Mean", "SD", "2.5%", "97.5%")))

for (i in 1:ncol(chain1)){
  sum.m[i,1:2] <- c(mean(sum.res[,i]),sd(sum.res[,i]))
  sum.m[i,3:4]<- quantile(sum.res[,i], c(0.025, 0.975))
} 
sum.m

###################### Posterior Predictive Check #######################################


# Function to test for heteroscedasicity 
hsc <- function(x){
  # first the residuals are split into 11 seperate chunks, the input vector should already be sorted according to the predicted value
  split.x <- split((x), ceiling(seq_along(x)/10))
  # for this dataset the last chunk would only contain one datapoint, so it is just included in the last chunk
  
  
  split.x[[11]] <- x[101:111]
  split.x[[12]] <- NULL
  
  ## get the variance of each chunk
  var.x <- sapply(split.x, var)

  # the discrepancy measure is returned
  # it is the ratio of max and min variance of the chunks
  return( max(var.x)/min(var.x))
  
}



# extract predicted values 

prediction.chain1 <-as.matrix(predicted[[1]])
prediction.chain2 <- as.matrix(predicted[[2]])

#extract  residuals for the observed data 
residuals.chain1 <- as.matrix(residuals[[1]])
residuals.chain2 <- as.matrix(residuals[[2]])


### the residuals for the observed values need to be sorted according to the predicted values ####
#create empty matrixes to store in orderd residuals and normal residuals
order1 <- matrix(nrow = nrow(residuals.chain1), ncol = 111)
order2 <- matrix(nrow = nrow(residuals.chain1), ncol = 111)
normal.residuals1 <- matrix(nrow= (nrow(residuals.chain1)), ncol=111)
normal.residuals2 <- matrix(nrow= (nrow(residuals.chain1)), ncol=111)



# residuals are orderd and
# normal distributed residuals are being simulated as well 

set.seed(1234)
for (i in 1:nrow(residuals.chain1)){
  
  # Sorting residuals accroding to predictions 
  order.prediction.1 <- order(prediction.chain1[i,])
  res.1 <- residuals.chain1[i,]
  order1[i,] <- res.1[order.prediction.1]

  order.prediction.2 <- order(prediction.chain2[i,])
  res.2 <- residuals.chain2[i,]
  order2[i,] <- res.2[order.prediction.2]
  ## Generating residuals from a normal distribution with sigma = e.sigma from the simulated data
  normal.residuals1[i,] <- rnorm(111, 0, mcmc[[1]][i,4])
  normal.residuals2[i,] <- rnorm(111, 0, mcmc[[2]][i,4])
}



# calculate the disprepancy measure over all residuals of the observed data
discrepancy.obs1 <- apply(order1,1, hsc)

#calculate the dispreancy measure over all simulated residuals 
discrepancy.sim1 <- apply(normal.residuals1, 1,hsc)


# p-value chain 1
sum(discrepancy.sim1 > discrepancy.obs1)/XG

# for second chain
discrepancy.obs2 <- apply(order2,1, hsc)
discrepancy.sim2 <- apply(normal.residuals2, 1,hsc)

#p.value chain 2
sum(discrepancy.sim2 > discrepancy.obs2)/XG

par(mfrow=c(1,1))
plot(density(discrepancy.obs1), lwd = 2, col="red", main="Distribution of Discrepancy Meassures", xlab = "Discrepancy Meassure")
lines(density(discrepancy.sim1), lwd =2, col ="blue")
legend(20,0.2, legend=c("Observerd residuals", "Simulated residuals"),
       col= c("red", "blue"), lty = 1)




################### DIC #############

logLikelihood = function(theta) {
  #Get the individual parameters out of theta.
  
  mu = (theta[1] + theta[2]*Solar.R + theta[3]*Wind)

  sigma2 = theta[4]
  
  
  #sum of log likelihoods = log of product of likelihoods
  sum( dnorm(Ozone, mu,(sigma2), log=TRUE) )
}

calculateDIC = function(y, theta_post) {
  #Calculate L
  theta_hat = apply(theta_post, 2, mean)
  #Dhat = -2*logLikelihood(theta_hat)
  L = logLikelihood(theta_hat)
  #Calculate P
  S = nrow(theta_post) #S = number of iterations
  #Add up the log likelihoods of each iteration
  
 llSum <- sum(apply(theta_post,1, logLikelihood))
  P = 2 * (L - (1 / S * llSum))
  
  #Calculate DIC
  DIC = -2 * (L - P)
  
  #Return the results
  list(DIC=DIC, P=P, L = L)
}

out <-calculateDIC(Ozone, as.matrix(sum.res))








### Frequentist approach ####


fit <- lm(Ozone~Solar.R + Wind, baindat)
summary(fit)
sd(fit$residuals)
confint(fit)
sum.m
