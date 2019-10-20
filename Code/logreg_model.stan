data {
  // Define training data
  int<lower = 1> train_N;
  vector[train_N] train_B;
  vector[train_N] train_Age;
  vector[train_N] train_DAP;
  int<lower = 0, upper = 1> train_Y[train_N];
  
  // Define test data (we want predictions for)
  int<lower = 1> test_N;
  vector[test_N] test_B;
  vector[test_N] test_Age;
  vector[test_N] test_DAP;
}
parameters {
  real alpha;
  real beta_B;
  real beta_Age;
  real beta_DAP;
}
model {
  train_Y ~ bernoulli_logit(alpha + beta_Age*train_Age + beta_B*train_B + beta_DAP*train_DAP );
  // priors on parameters
  alpha ~ normal(0, 1);
  beta_B ~ normal(0, 1);   
  beta_Age ~ normal(0, 1);
  beta_DAP ~ normal(0, 1);
}
generated quantities {
  vector[test_N] test_Y;
  vector[test_N] test_odds;
  vector[test_N] test_Y_prob;
  
  vector[train_N] train_odds;
  vector[train_N] train_Y_prob;
  
  for(i in 1:test_N) {
    // using samples from the posterior distribution p( theta | y ) 
    test_Y[i] = bernoulli_rng(inv_logit( alpha + beta_B*test_B[i] + beta_Age*test_Age[i] + beta_DAP*test_DAP[i] ));
    // compute odds, so we can record prob( Y = 1 | x )
    test_odds[i] = exp( alpha + beta_B*test_B[i] + beta_Age*test_Age[i] + beta_DAP*test_DAP[i] );
    // probability of Y = 1 (positive case)
    test_Y_prob[i] = test_odds[i] / (test_odds[i] + 1); 
  }
  
  // we may also want the predictions for the training set so we can 
  // derive an operating threshold to compare with validation set performance
  for(i in 1:train_N) {
    train_odds[i] = exp( alpha + beta_B*train_B[i] + beta_Age*train_Age[i] + beta_DAP*train_DAP[i] );
    train_Y_prob[i] = train_odds[i] / (train_odds[i] + 1); 
  }
  
  
}

