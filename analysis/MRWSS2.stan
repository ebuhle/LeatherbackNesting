data {
  int<lower=1> N;              // total # observations
  int<lower=1> year[N];        // year index of each obs, starting at 1
  int<lower=1> turtle[N];      // turtle ID, starting at 1
  vector[N] y;                 // observations
}

transformed data {
  int<lower=1> N_year = max(year);      // # years
  int<lower=1> N_turtle = max(turtle);  // # turtles
}

parameters {
  vector[N_turtle] x1;         // initial states
  real mu;                     // mean annual drift 
  real<lower=0> sigma;         // process error SD
  real<lower=-1,upper=1> rho;  // between-turtle process error correlation
  vector[N_turtle] w[N_year];  // annual process error vectors
  real<lower=0> tau;           // observation error SD
}

transformed parameters {
  corr_matrix[N_turtle] R;     // process error correlation matrix 
  // matrix[N_turtle,N_turtle] L; // Cholesky factor of covariance matrix
  vector[N_turtle] x[N_year];  // states
  
  // Compound symmetric process error correlation matrix 
  // ("equalvarcov" in MARSS terminology)
  R = diag_matrix(rep_vector(1,N_turtle));
  for(i in 2:N_turtle)
    for(j in 1:(i-1))
    {
      R[i,j] = rho;
      R[j,i] = rho;
    }
    
  // L = diag_pre_multiply(rep_vector(sigma,N_turtle), cholesky_decompose(R));
  
  // Construct states using multivariate Matt trick
  x[1] = x1;
  for(t in 2:N_year)
    x[t] = x[t-1] + mu + w[t];
}

model {
  vector[N] yhat;  // elements of x matched with elements of y
  
  // Priors
  x1 ~ normal(0,100);
  mu ~ normal(0,10);
  sigma ~ student_t(3,0,10);
  // for(t in 1:N_year) w_std[t] ~ normal(0,1); // implies w[t] ~ MVN(0, L*L')
  w ~ multi_normal(rep_vector(0, N_turtle), quad_form_diag(R, rep_vector(sigma,N_turtle)));
  tau ~ student_t(3,0,10);
  
  // Likelihood
  for(i in 1:N) yhat[i] = x[year[i]][turtle[i]];
  y ~ normal(yhat, tau);
}

generated quantities {
  vector[N] LL;  // pointwise log-likelihood
  
  for(i in 1:N) LL[i] = normal_lpdf(y[i] | x[year[i]][turtle[i]], tau);
}
