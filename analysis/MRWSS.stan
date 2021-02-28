data {
  int<lower=1> N;            // total # observations
  int<lower=1> year[N];      // year index of each obs, starting at 1
  int<lower=1> turtle[N];    // turtle ID, starting at 1
  vector[N] y;               // observations
}

transformed data {
  int<lower=1> N_year = max(year);      // # years
  int<lower=1> N_turtle = max(turtle);  // # turtles
}

parameters {
  real theta1;               // initial value of shared random walk
  real mu;                   // mean annual drift 
  vector[N_year] omega_std;  // annual shared process errors (Z-scored)
  real<lower=0> sigma_omega; // SD of shared process error
  vector[N] nu_std;          // unique process errors (Z-scored)
  real<lower=0> sigma_nu;    // SD of unique process error
  real<lower=0> tau;         // observation error SD
}

transformed parameters {
  vector[N_year] theta;      // shared random walk
  vector[N] x;               // states
  
  // States include shared RW component plus independent process errors
  // Implies a multivariate RW with compound symmetric process covariance matrix
  // (Q "equalvarcov" in MARSS terminology)
  
  // Shared random walk process
  theta[1] = theta1;
  for(t in 2:N_year)
    theta[t] = theta[t-1] + mu + omega_std[t]*sigma_omega;
    
  // Independent process errors
  x = theta[year] + nu_std*sigma_nu;
}

model {
  // Priors
  theta1 ~ normal(0,100);
  mu ~ normal(0,10);
  sigma_omega ~ student_t(3,0,10);
  omega_std ~ std_normal();       // implies omega ~ N(0, sigma_omega) 
  sigma_nu ~ student_t(3,0,10);
  nu_std ~ std_normal();          // implies nu ~ N(0, sigma_nu)
  tau ~ student_t(3,0,10);
  
  // Likelihood
  y ~ normal(x, tau);
}

generated quantities {
  vector[N] LL;  // pointwise log-likelihood
  
  for(i in 1:N) LL[i] = normal_lpdf(y[i] | x, tau);
}
