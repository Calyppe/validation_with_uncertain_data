# MSc Thesis
# DSM with uncertain data
# RQ 1: Validation metrics
# Author: Chloe Girka
# Date: 2020-11-03          Increase epsilon and measurement errors


#############################
#### Indicator functions ####
#############################

# Mean Error (ME)
EstimateME <- function(prediction, observation) {
  mean(observation-prediction)
} 

# Mean Squared Error (MSE)
EstimateMSE <- function(prediction, observation, error_variance) {
  mean((prediction-observation)^2)-mean(error_variance)
}

# Model Efficiency Coefficient (MEC)
EstimateMEC <- function(prediction, observation, error_variance) {
  (1 - (sum((prediction-observation)^2) / sum((observation-mean(observation))^2))) * 
    (1+ ((mean(error_variance)/(1/length(prediction)*sum((observation-mean(observation))^2)+length(prediction)*mean(error_variance)))))
}

# Lin's Concordance Correlation Coefficient (CCC)
EstimateCCC <- function(prediction, observation, error_variance) {
  (2* (((1/length(observation)*sum(observation*prediction))-(mean(observation)*(mean(prediction))))/
         (var(observation)+var(prediction)+((mean(observation)-mean(prediction))^2)))) * 
    ((var(observation)+var(prediction)+((mean(observation)-mean(prediction))^2))
     / (var(observation)-mean(error_variance)+var(prediction)+((mean(observation)-mean(prediction))^2)))
}

# Prediction Interval Coverage Probability (PICP)
EstimatePICP <-function(interval, observation){
  q <- 0
  for (i in 1:length(observation)){
    if (observation[i] >= min(interval[i,])) {
      if (observation[i] <= max(interval[i,])){
        q <- q+1
      }
    }
  }
  q/length(observation)
}


#################################
#### Creating synthetic data ####
#################################

sample_size <- 1000
seed <- 1

set.seed(seed)

# Creating covariates
x1 <- sample(1:100, sample_size, replace=TRUE) + rnorm(mean=0, sd=1, n=sample_size) # picking integers between 1 and 100 and adding small variation
x2 <- sample(1:100, sample_size, replace=TRUE) + rnorm(mean=0, sd=1, n=sample_size)

# Creating the relationship between X1, X2 and Y
variation <- rnorm(mean=0, sd=6, n=sample_size)
y_true <- 2*x1 - 0.01*x1*x1 + 0.5 + x2 + variation

indices_method <- sample(1:sample_size, floor(sample_size/4), replace=FALSE)    # Drawing a portion of the observations
var_e1 <- 5
var_e2 <- 16
var_e <- rep(0, sample_size)
var_e[indices_method] <- var_e1                                                 # Giving the subsample measurement error variance 1
var_e[-indices_method] <- var_e2                                                # Giving the rest of the observations measurement error variance 2

measurement_error <- rnorm(mean=0, sd=sqrt(var_e), n=sample_size)               # Computing individual errors
y_obs <- y_true + measurement_error                                             # Computing observed y
udf <- data.frame(x1, x2, var_e, y_true, y_obs)                                 # Creating "uncertain" dataframe containing the covariates,
                                                                                # the measurement error variances, true y and observed y

###################
#### Modelling ####
###################

nb_runs <- 50

ME <- rep(0, nb_runs)
ME_obs <- rep(0, nb_runs)
MSE <- rep(0, nb_runs)
MSE_obs <- rep(0, nb_runs)
MSE_adjusted <- rep(0, nb_runs)
MEC <- rep(0, nb_runs)
MEC_obs <- rep(0, nb_runs)
MEC_adjusted <- rep(0, nb_runs)
CCC <- rep(0, nb_runs)
CCC_obs <- rep(0, nb_runs)
CCC_adjusted <- rep(0, nb_runs)
PICP <- rep(0, nb_runs)
PICP_obs <- rep(0, nb_runs)
PICP_adjusted <- rep(0, nb_runs)

for (i in 1:nb_runs) {


# Splitting data into training and validation set
train_ind <- sample(1:nrow(udf),0.8*nrow(udf))
train_data <- udf[train_ind,]

# Fitting
lm.fit <- lm(y_obs~x1+I(x1^2)+x2, data=train_data)

# Prediction
udf$y_sim <- predict(lm.fit, udf)
test_data <- udf[-train_ind,]


#####################################
#### Assessing model performance ####
#####################################

# Calculating ME
ME[i] <- EstimateME(test_data$y_sim, test_data$y_true)      # Real ME
ME_obs[i] <- EstimateME(test_data$y_sim, test_data$y_obs)   # Observed ME


# Calculating MSE

MSE[i] <- (1/length(test_data$y_true))*sum((test_data$y_sim-test_data$y_true)^2)                       # Real MSE, calculated from the true values of y

MSE_obs[i] <- (1/length(test_data$y_true))*sum((test_data$y_sim-test_data$y_obs)^2)                    # Observed MSE, using the same formula but on the observed data

MSE_adjusted[i] <- EstimateMSE(test_data$y_sim, test_data$y_obs, test_data$var_e)                           # Adjusted MSE, accounting for measurement errors


# Calculating MECs

MEC[i] <- 1 - (sum((test_data$y_sim-test_data$y_true)^2) / sum((test_data$y_true-mean(test_data$y_true))^2))   # Real MEC, calculated from the true values of y

MEC_obs[i] <- 1 - (sum((test_data$y_sim-test_data$y_obs)^2) / sum((test_data$y_obs-mean(test_data$y_obs))^2))  # Observed MEC, using the same formula but on the observed data

MEC_adjusted[i] <- EstimateMEC(test_data$y_sim, test_data$y_obs, test_data$var_e)                                   # Adjusted MEC, accounting for measurement errors


# Calculating CCC

CCC[i] <- 2* ((( (1/length(test_data$y_true)) * sum(test_data$y_true*test_data$y_sim)) - (mean(test_data$y_true) * (mean(test_data$y_sim))))/
         (var(test_data$y_true) + var(test_data$y_sim) + ((mean(test_data$y_true)-mean(test_data$y_sim))^2)))               # Real CCC, calculated from the true values of y

CCC_obs[i] <- 2* (((1/length(test_data$y_obs)*sum(test_data$y_obs*test_data$y_sim))-(mean(test_data$y_obs)*(mean(test_data$y_sim))))/
         (var(test_data$y_obs)+var(test_data$y_sim)+((mean(test_data$y_obs)-mean(test_data$y_sim))^2)))                     # Observed CCC, using the same formula but on the observed data

CCC_adjusted[i] <- EstimateCCC(test_data$y_sim, test_data$y_obs, test_data$var_e)                                                       # Adjusted MEC, accounting for measurement errors


# Calculating PICP
pred_interval <- predict(lm.fit, test_data, interval='prediction', level=0.95)

PICP[i] <- EstimatePICP(pred_interval, test_data$y_true)            # PICP using true values of Y
PICP_obs[i] <- EstimatePICP(pred_interval, test_data$y_obs)         # PICP calculated using observed values, ignoring uncertainties


error_sd <- rep(0, nrow(test_data))

for (j in 1:nrow(test_data)) {
  error_sd[j] <- (max(pred_interval[j,])-min(pred_interval[j,]))/(2*1.96)
}
mini <- test_data$y_sim - 1.96*sqrt(error_sd^2 + test_data$var_e)
maxi <- test_data$y_sim + 1.96*sqrt(error_sd^2 + test_data$var_e)
adjusted_interval <- data.frame(mini, maxi)
PICP_adjusted[i] <- EstimatePICP(adjusted_interval, test_data$y_obs)

}

mean(ME)
mean(ME_obs)
mean(MSE)
mean(MSE_obs)
mean(MSE_adjusted)
mean(MEC)
mean(MEC_obs)
mean(MEC_adjusted)
mean(CCC)
mean(CCC_obs)
mean(CCC_adjusted)
mean(PICP)
mean(PICP_obs)
mean(PICP_adjusted)
