data {
  int<lower=1> N;             // total # observations
  int<lower=1> year[N];       // year index of each obs, starting at 1
  int<lower=1> turtle[N];     // turtle ID, starting at 1
  vector[N] y;                // observations
}

transformed data {
  int<lower=1> N_year = max(year);      // # years
  int<lower=1> N_turtle = max(turtle);  // # turtles
  real y_bar = mean(y);                 // sample grand mean
  real<lower=0> s = sd(y);              // sample SD
  vector[N] z = (y - y_bar)/s;          // Z-scored data
}

parameters {
  // "0" denotes parameters defined for Z-scored data
  real alpha01;                // initial value of random walk
  real mu0;                    // drift 
  vector[N_year] omega0_std;   // process errors (Z-scored)
  real<lower=0> sigma_alpha0;  // process error SD
  vector[N_turtle] theta0_std; // turtle effects (Z-scored)
  real<lower=0> sigma_theta0;  // among-turtle SD
  real<lower=0> sigma0;        // observation error SD
}

transformed parameters {
  vector[N_year] alpha0;       // random walk
  vector[N_turtle] theta0;     // turtle effects
  vector[N] z_hat;             // predicted values of Z-scored data
  
  // Random walk process
  alpha0[1] = alpha01;
  for(t in 2:N_year)
    alpha0[t] = alpha0[t-1] + mu0 + omega0_std[t]*sigma_alpha0;

  // Turtle effects
  theta0 = sigma_theta0*theta0_std;
  
  // Predicted values
  z_hat = alpha0[year] + theta0[turtle];
}

model {
  // Priors
  alpha01 ~ normal(0,1);
  mu0 ~ normal(0,1);
  sigma_alpha0 ~ student_t(3,0,2);
  omega0_std ~ std_normal();        // implies omega ~ N(0, sigma_omega) 
  sigma_theta0 ~ student_t(3,0,2);
  theta0_std ~ std_normal();        // implies theta ~ N(0, sigma_theta)
  sigma0 ~ student_t(3,0,2);
  
  // Likelihood
  z ~ normal(z_hat, sigma0);
}

generated quantities {
  // shift and rescale parameters to apply to raw data
  vector[N_year] alpha = s*alpha0 + y_bar;
  real mu = s*mu0;
  real sigma_alpha = s*sigma_alpha0;
  vector[N_turtle] theta = s*theta0;
  real sigma_theta = s*sigma_theta0;
  real sigma = s*sigma0;
  vector[N] y_hat = s*z_hat + y_bar;
  vector[N] LL;  // pointwise log-likelihood
  
  for(i in 1:N) LL[i] = normal_lpdf(z[i] | z_hat[i], sigma0);
}
