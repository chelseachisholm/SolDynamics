
data {
  int<lower=1> nsite;                   // Number of sites
  int<lower=1> nyear;                   // Number of years
  real<lower=0,upper=2> y[nsite, nyear]; // Number of detections
  matrix[nsite, nsite] D; // Distance matrix
}

transformed data {
  int ny_minus_1 = nyear - 1;
}

parameters {
  	real<lower=0,upper=1> psi1; 
	vector<lower=0,upper=1>[ny_minus_1] phi;
	matrix[nsite,ny_minus_1] gamma;
	real a;
	real y2;
	
}

model {
  // Priors
  phi ~ normal(0.5, 0.5);
  a ~ uniform(0, 10);
  y2 ~ uniform(0, 100);
  // Flat priors Uniform(0, 1) are implicitly used on psi1

  // Likelihood
  for (i in 1:nsite) {
    for (k in 2:nyear) {
    	S[i,k-1] = sum((y'[k-1]) * exp((-1) * a * D[i,]))
    	gamma[i,k-1] = pow(S[i,k-1], 2)/(pow(S[i,k-1], 2) + y2)
    	y[i,k] ~ bernoulli(y[i,k-1] * phi[k-1] + (1 - y[i,k-1]) * gamma[i,k-1])
    ;
    } //k
  } //i

    
} //model

generated quantities {
  // Sample and population occupancy, growth rate and turnover
  vector<lower=0,upper=1>[nyear] psi;        // Occupancy probability

psi[1] = psi1;
for (k in 2:nyear) 
	psi[k] = psi[k - 1] * phi[k - 1] + (1 - psi[k - 1]) * gamma[k - 1];
	
}
