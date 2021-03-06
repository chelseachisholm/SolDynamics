data {
  int<lower=0> nsite;           // Number of populations
  int<lower=0> nyear;           // Number of years
  real<lower=0> C[nyear, nsite]; // Abundance
  vector[nyear] year;           // Year covariate
}

transformed data {
  vector[nyear] year_squared;

  year_squared = year .* year;
}

parameters {
  real mu;
  vector[nsite] alpha;
  real eps[nyear];
  real beta[2];
  real<lower=0,upper=2> sd_alpha;
  real<lower=0,upper=1> sd_year;
  real<lower=0> sigma;
}

transformed parameters {
  vector[nsite] log_lambda[nyear];

  // Linear predictor including random site and random year effects
  for (i in 1:nyear)
    log_lambda[i] = alpha
                  + beta[1] * year[i]
                  + beta[2] * year_squared[i]
                  + eps[i];
}

model {
  // Priors

  // Random site effects
  alpha ~ normal(mu, sd_alpha);

  // Hyperparameter 1
  mu ~ normal(0, 10);

  // Hyperparameter 2
  //  sd_alpha ~ uniform(0, 2); // Implicitly defined

  beta ~ normal(0, 10);

  // Hyperparameter 3
  //  sd_year ~ uniform(0, 1); // Implicitly defined

  // Random year effects
  eps ~ normal(0, sd_year);

  //did not specify priors for sigma, as uniform flat is better according to Gelman

  // Likelihood
  for (i in 1:nyear) {
    // Distribution for random part
    // Link function
    C[i] ~ normal(log_lambda[i], sigma);
  }
}
