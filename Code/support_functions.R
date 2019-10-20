# -- dwj : October 29th 2018

# -----------------------------------------------------------------
# Functions for calibration curves

Calibrate.Bins <- function( obs, pred, ref.level, break.levels = 11 ) {
  # Computes a calibration curve by equal-width bins over the predicted probability range
  #
  # Arguments:
  #   obs   : vector of N observed events (factor)
  #   pred  : vector of N predicted probabilities (e.g. from a predictive model)
  #   ref.level : a string indicating which value in obs should be treated as the reference class (i.e. the "positive" class)
  #   break.levels : the number of bins to use
  #
  # Returns:
  #   calib.df : a data.frame of the bin middle points (for plotting), the factor of bin intervals, actual probability (for each bin) from obs, 
  #               the mean predicted probability for each bin, and a measure of spread for the bin (standard deviation)
  
  bin.obs       <- cut( pred, breaks = break.levels, ordered = TRUE )
  bin.levels    <- levels( bin.obs )
  bin.midpoints <- midpoints( bin.levels )
  
  calib.df <- data.frame( Prob.bins.mids = bin.midpoints,
                          Prob.bins      = bin.levels,
                          P.Y            = rep(NA, length(bin.levels)),
                          Pred.Y         = rep(NA, length(bin.levels)),
                          Spread.Pred.Y  = rep(NA, length(bin.levels))
  )
  
  # -- for each interval in bin.levels
  for ( i in 1:length( bin.levels ) ) {
    
    # -- for the ith bin, compute the observed probability of events
    observed   <- obs[ which( bin.obs == bin.levels[i] ) ]
    P.observed <- length( which( observed == ref.level ) ) / length( observed )
    
    # -- for the same bin, compute the mean predicted probability
    predicted  <- pred[ which(  bin.obs == bin.levels[i] ) ]
    
    # -- store away in calib.df
    calib.df$P.Y[i]           <- P.observed
    calib.df$Pred.Y[i]        <- mean( predicted )
    calib.df$Spread.Pred.Y[i] <- sd( predicted )
    
  }
  
  return( calib.df )
}

Calibrate.Lowess <- function( obs, pred, ref.level, break.levels ) {
  # Computes calibration curve, but using lowess scatterplot smoother
  # break.levels, here, is used as the fraction of points used to 
  # to compute the smoother "span"
  #
  # Arguments: as for Calibrate.Bins
  # 
  # Returns:
  #   as for Calibrate.Bins, but returns one value for Pred.Y (the x-axis of the calibration curve) 
  #   corresponding to one point in the P.Y (observed probability)
  
  f.val <- 1/break.levels
  
  bin.obs <- ifelse( obs == ref.level, 1, 0 )
  calib.low <- as.data.frame( lowess( pred, bin.obs, iter = 0, f = f.val ) )

  calib.df <- data.frame( P.Y     = calib.low$y,
                          Pred.Y  = calib.low$x
  )  
  
  return( calib.df )
}

# -----------------------------------------------------------------
# Functions to generate simulated data with correlation structure
GenSimulatedData <- function( N.sim.samples, new.cov ) {
  # -- create some simulated data
  x.1 <- rnorm(N.sim.samples, 15, 5)
  x.234 <- scale(matrix( rnorm(3 * N.sim.samples), ncol=3 ))
  x.1234 <- cbind(scale(x.1),x.234)
  
  # -- covariance
  cov.1 <- var(x.1234)
  LL.1  <- solve(chol(cov.1))
  new.x <-  x.1234 %*% LL.1
  
  # - desired correlations
  LL.2 <- chol(new.cov)
  cor.data <- new.x %*% LL.2 * sd(x.1) + mean(x.1)
  return( cor.data )
}

GenSimDichot <- function( sim.data ) {
  # -- takes a data set simulated by GenSimulatedData
  #    and turns it into data scenario we might see in a real life study 
  #    So, we preserve the structure of the data, but we do some transforms to 
  #    simulate the kinds of problems seen in the literature
  #    The idea is to preserve, but "hide", the underlying statistical structure of sim.data
  #
  # -- In this scenrio, what we produce is a data set
  #    where we have a biomarker B, and Age, DAP (duration attenuated psychotic symptoms) and a dichotomised outcome
  #    variable which might be e.g. transition to full blown psychosis.
  
  # apply some transformations to make the data look real/plausible
  df <- data.frame( B   = sim.data[,1],
                    Age = sim.data[,2] * 1.2,
                    DAP = sim.data[,3] / 10,
                    Y   = sim.data[,4] * 5 )
 
  # Now, assume a clinician has labelled the classes based on a threshold for Y
  df$Y <- ifelse( df$Y > 90, 1, 0 )  # 1 is the positive class
  
  return( df )
}


# ---------------------------------------
# Stan / model building code

TrainModel <- function( training.data, testing.data ) {
  
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
  
  # build the model : we'll seperate out all variables to make the model explicit (rather than use matrices)
  train_B   <- training.data$B
  train_Age <- training.data$Age
  train_DAP <- training.data$DAP
  train_Y   <- training.data$Y
  
  test_B   <- testing.data$B
  test_Age <- testing.data$Age
  test_DAP <- testing.data$DAP
  test_Y   <- testing.data$Y
  
  train_N  <- length( train_B )
  test_N   <- length( test_B )
  
  fit <- stan(file = "logreg_model.stan",
              data = list(train_B, train_Age, train_DAP, train_Y, train_N,
                          test_B, test_Age, test_DAP, test_N),
              chains = 4, iter = 4000)
  
  # # create a convenient dataframe for storage / reloading that doesn't depend on Rstan
  mcmc.output <- list(
    # posterior samples for parameters : e.g. P( theta | y )
    alpha    = extract(fit, 'alpha')$alpha,
    beta_B   = extract(fit, 'beta_B')$beta_B,
    beta_Age = extract(fit, 'beta_Age')$beta_Age,
    beta_DAP = extract(fit, 'beta_DAP')$beta_DAP,

    # samples from the PPD rows = samples, cols = patients
    PPD = extract(fit, 'test_Y_prob')$test_Y_prob,
    
    # retrieve the training performance 
    y_pred.train = extract(fit, 'train_Y_prob')$train_Y_prob
    
  )
  
  # save the complete fit object -- so we don't have to use Stan to run the example as a markdown notebook
  saveRDS( mcmc.output, file = "../Data/MCMC-output.rds" )
}



