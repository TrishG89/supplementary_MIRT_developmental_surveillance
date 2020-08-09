data {
  int<lower=1> J; // Number of items
  int<lower=1> I; // Number of subjects
  int<lower=1> N; // Number of observations
  int<lower=1,upper=6> K; //Number of latent dimensions
  int<lower=1, upper=J> jj[N]; // Item for n
  int<lower=1, upper=I> ii[N]; // Examinee for n
  int<lower=1,upper=K> kk[N]; //Dimension for observation n
  int<lower=0, upper=1> y[N]; // Binary response for n
  matrix[6,6]W;
}
parameters {
  vector<lower=0>[J] alpha; // Discrimination for item J
  vector[J] beta; // Difficulty for item J
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[K] mu_theta;
  cov_matrix[K] Sigma_theta;
  matrix[K,I] delta;
}
transformed parameters {
  cholesky_factor_cov[K] L;
  matrix[K,I] theta;
  matrix[K,I] mat_mu_theta;
  for (i in 1:I){
    mat_mu_theta[,i]= mu_theta;
  }
  L = cholesky_decompose(Sigma_theta);
  theta = mat_mu_theta + L * delta;
  
}
model {
  vector[N] eta;
  alpha ~ normal(mu_alpha, sigma_alpha); // Prior for discrimination parameter
  beta ~ normal(0,1); // Prior for difficulty parameter
  mu_alpha ~ normal(1,1)T[0,];
  sigma_alpha ~ inv_gamma(1,1);
  mu_theta ~ std_normal();
  Sigma_theta ~ inv_wishart(7,W[,]); 
  for (k in 1:K){
  delta[k,] ~ std_normal();
  }
   for (n in 1:N)
    eta[n] = alpha[jj[n]] * (theta[kk[n],ii[n]] - beta[jj[n]]); // 2PL Model
    
   
  y ~ bernoulli_logit(eta);
}
generated quantities {
   matrix[K,K] cor_sigma;
   for(i in 1:K){
     for(j in 1:K){
         cor_sigma[i,j] = Sigma_theta[i,j]/sqrt(Sigma_theta[i,i]*Sigma_theta[j,j]);
         }
      }
}
