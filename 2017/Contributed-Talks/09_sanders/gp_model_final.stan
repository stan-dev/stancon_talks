
data {
	int<lower=1> N1;
	vector[N1] x1;
	int z1[N1];
	int<lower=1> N2;
	vector[N2] x2;
	real<lower=0> alpha_rho;
	real<lower=0> beta_rho;
}
transformed data {
	int<lower=1> N;
	vector[N1+N2] x;
	// cov_exp_quad wants real valued inputs
	real rx[N1+N2]; 
	real rx1[N1];
	real rx2[N2];
	
	N = N1 + N2;
	x = append_row(x1, x2);
	
	rx = to_array_1d(x);
	rx1 = to_array_1d(x1);
	rx2 = to_array_1d(x2);
}
parameters {
	vector[N1] y_tilde1;
	real<lower=0> eta_sq;
	real<lower=1> inv_rho;
	real<lower=0> sigma_sq;
	real mu_0;
	real mu_b;
	real<lower=0> NB_phi_inv;
}
model {
	vector[N1] mu1;
	vector[N1] y1;
	matrix[N1,N1] Sigma1;
	matrix[N1,N1] L1;
	
	// Calculate mean function
	mu1 = mu_0 + mu_b * x1;

	// GP hyperpriors
	eta_sq ~ cauchy(0, 1);
	sigma_sq ~ cauchy(0, 1);
	inv_rho ~ gamma(alpha_rho, beta_rho); // Gamma prior with mean of 4 and std of 2
	
	// Calculate covariance matrix using new optimized function
	Sigma1 = cov_exp_quad(rx1, sqrt(eta_sq), sqrt(0.5) * inv_rho);
	for (n in 1:N1) Sigma1[n,n] = Sigma1[n,n] + sigma_sq;
	
	// Decompose
	L1 = cholesky_decompose(Sigma1);
	// We're using a the non-centered parameterization, so rescale y_tilde
	y1 = mu1 + L1 * y_tilde1;

	// Mean model priors
	mu_0 ~ normal(0, 2);
	mu_b ~ normal(0, 0.2);
	
	// Negative-binomial prior
	// For neg_binomial_2, phi^-1 controls the overdispersion.  
	// phi^-1 ~ 0 reduces to the poisson.  phi^-1 = 1 represents variance = mu+mu^2
	NB_phi_inv ~ cauchy(0, 5);
	
	// Generate non-centered parameterization
	y_tilde1 ~ normal(0, 1);
	
	// Likelihood
	z1 ~ neg_binomial_2_log(y1, inv(NB_phi_inv));
}

generated quantities {
	vector[N1] y1;
	vector[N2] y2;
	vector[N] y;
	int z_rep[N];
	
	{
		// Don't save these parameters
		matrix[N,N] Sigma;
		matrix[N,N] L;
		vector[N] y_tilde;
		
		Sigma = cov_exp_quad(rx, sqrt(eta_sq), sqrt(0.5) * inv_rho);
		for (n in 1:N) Sigma[n,n] = Sigma[n,n] + sigma_sq;
		
		for (n in 1:N1) y_tilde[n] = y_tilde1[n];
		for (n in (N1 + 1):N) y_tilde[n] = normal_rng(0,1);
		
		// Decompose
		L = cholesky_decompose(Sigma);
		y = mu_0 + mu_b * x + L * y_tilde;

		for (n in 1:N1) y1[n] = y[n];
		for (n in 1:N2) y2[n] = y[N1+n];
				
		for (n in 1:N) z_rep[n] = neg_binomial_2_log_rng(y[n], inv(NB_phi_inv));
	}
}
