
data {
  int<lower=1> nsite;                   // Number of sites
  int<lower=1> nyear;                   // Number of years
  int<lower=0,upper=2> y[nsite, nyear]; // Number of detections
}

transformed data {
  int ny_minus_1 = nyear - 1;
}

parameters {
  real<lower=0,upper=1> psi1; 
	vector<lower=0,upper=1>[ny_minus_1] phi;
	vector<lower=0,upper=1>[ny_minus_1] gamma;


}

model {
  // Priors
  phi ~ normal(0.5, 0.5);
  gamma ~ normal(0.5, 0.5);
  // Flat priors Uniform(0, 1) are implicitly used on psi1

  // Likelihood
  for (i in 1:nsite) {
    //vector[2] z[nyear];
    // First year
    //y[i, 1] = psi1;
    //y[i, 2] = (1 - psi1);

    for (k in 2:nyear) {
    	y[i,k] ~ bernoulli(y[i,k-1] * phi[k-1] + (1 - y[i,k-1]) * gamma[k-1])
    ;
    } //k
  } //i

    
} //model

generated quantities {
  // Sample and population occupancy, growth rate and turnover
  vector<lower=0,upper=1>[nyear] psi;        // Occupancy probability
  //vector<lower=0,upper=1>[ny_minus_1] phi;   // Survival probability
  //vector<lower=0,upper=1>[ny_minus_1] gamma; // Colonization probability
  // int<lower=0,upper=nsite> n_occ[nyear];     // Number of occupancy

psi[1] = psi1;
for (k in 2:nyear) 
	psi[k] = psi[k - 1] * phi[k - 1] + (1 - psi[k - 1]) * gamma[k - 1];
	
}
