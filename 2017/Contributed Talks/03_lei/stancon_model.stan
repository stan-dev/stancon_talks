data {
	int<lower=0> N ; # Total number of data points
	int<lower=0> K_m; # Number of films
	int<lower=0> K_mg; # Number of film groups
	int<lower=0> K_ads; # Number of ad platforms

	int<lower=0,upper=1> y[N] ; # Have seen target film

	int<lower=0,upper=1> x_parent[N] ; # Parent code
	real x_age[N] ; # Age
	real<lower=0> x_gender[N] ; # Gender
	int<lower=0> x_film[N] ; # Film
	matrix[N, K_ads] x_ads ; # Ad impressions per platform

	int<lower=1, upper=K_mg> x_mg[K_m] ; # Film groups
}

transformed data {
	real realN;
	real logit_mean_y;
	realN = N;
	
	logit_mean_y = logit(sum(y) / realN);
}

parameters {
	real b;
	
	real<lower=0> v_age_sigma;
	vector[K_mg] v_age_g;
	vector<lower=0>[K_mg] v_age_sigma_g;
	vector[K_m] v_age;
		
	real<lower=0> v_gender_sigma;
	vector[K_mg] v_gender_g;
	vector<lower=0>[K_mg] v_gender_sigma_g;
	vector[K_m] v_gender;
	
	real<lower=0> v_parent_sigma;
	vector[K_mg] v_parent_g;
	vector<lower=0>[K_mg] v_parent_sigma_g;
	vector[K_m] v_parent;
	
	real<lower=0> v_ad_sigma;
	vector[K_ads] v_ad_platform_mean;
	vector<lower=0>[K_ads] v_ad_platform_sigma;
	matrix[K_ads, K_mg] v_ad_platform_g;
	matrix<lower=0>[K_ads, K_mg] v_ad_platform_sigma_g;
	matrix[K_m,K_ads] v_ad_platform_film;
	
	real<lower=0> v_film_sigma;
	vector<lower=0>[K_mg] v_film_g_sigma;
	vector[K_mg] v_film_g;
	vector[K_m] v_film_s;
}

model {
	real y_pred[N];
	
	b ~ normal(0, 0.5);
	
	v_age_sigma ~ cauchy(0, 0.5);
	v_age_g ~ normal(0, 1);
	v_age_sigma_g ~ cauchy(0, 1);
	v_age ~ normal(0, 1);
	
	v_gender_sigma ~ cauchy(0, 0.5);
	v_gender_g ~ normal(0, 1);
	v_gender_sigma_g ~ cauchy(0, 1);
	v_gender ~ normal(0, 1);
	
	v_parent_sigma ~ cauchy(0, 0.5);
	v_parent_g ~ normal(0, 1);
	v_parent_sigma_g ~ cauchy(0, 1);
	v_parent ~ normal(0, 1);
	
	v_ad_sigma ~ cauchy(0, 0.5);
	v_ad_platform_mean ~ normal(0, 1);
	v_ad_platform_sigma ~ cauchy(0, 1);
	
  to_vector(v_ad_platform_g) ~ normal(0, 1);
  to_vector(v_ad_platform_sigma_g) ~ cauchy(0, 1);
  to_vector(v_ad_platform_film) ~ normal(0, 1);
	
	v_film_sigma ~ cauchy(0,1);
	v_film_g_sigma ~ cauchy(0,1);
	v_film_g ~ normal(0,1);
	v_film_s ~ normal(0,1);

	for (n in 1:N) {
		int g;
		int m;
		real comb_ad_platform;
		
		m = x_film[n];
		g = x_mg[m];
		
		comb_ad_platform = 0;
		for (k in 1:K_ads) {
			comb_ad_platform = comb_ad_platform +
				v_ad_sigma * (v_ad_platform_mean[k] + v_ad_platform_sigma[k] * (v_ad_platform_g[k,g]  + v_ad_platform_sigma_g[k,g] * v_ad_platform_film[m,k])) * x_ads[n,k];
		}
		
		y_pred[n] = logit_mean_y + b + 
			(v_age_sigma * (v_age_sigma_g[g] * v_age[m] + v_age_g[g])) * x_age[n] +		
			(v_gender_sigma * (v_gender_sigma_g[g] * v_gender[m] + v_gender_g[g])) * x_gender[n] + 
			(v_parent_sigma * (v_parent_sigma_g[g] * v_parent[m] + v_parent_g[g])) * x_parent[n] + 
			comb_ad_platform +
			v_film_sigma * (v_film_g_sigma[g] * v_film_s[m] + v_film_g[g]);
	}
	
	y ~ bernoulli_logit(y_pred);
}
	
generated quantities {
  real log_lik[N]; # Log-likelihood of each data point given a posterior sample
  real theta[N]; # The probabilities of p(y=1|x) for each data point and MCMC sample

	for (n in 1:N) {
		int g;
		int m;
		real t_i;
		real comb_ad_platform;

		m = x_film[n];
		g = x_mg[m];

		comb_ad_platform = 0;
		for (k in 1:K_ads) {
			comb_ad_platform = comb_ad_platform +
				v_ad_sigma * (v_ad_platform_mean[k] + v_ad_platform_sigma[k] * (v_ad_platform_g[k,g]  + v_ad_platform_sigma_g[k,g] * v_ad_platform_film[m,k])) * x_ads[n,k];
		}

		t_i = logit_mean_y + b + 
			(v_age_sigma * (v_age_sigma_g[g] * v_age[m] + v_age_g[g])) * x_age[n] +		
			(v_gender_sigma * (v_gender_sigma_g[g] * v_gender[m] + v_gender_g[g])) * x_gender[n] + 
			(v_parent_sigma * (v_parent_sigma_g[g] * v_parent[m] + v_parent_g[g])) * x_parent[n] + 
			comb_ad_platform +
			v_film_sigma * (v_film_g_sigma[g] * v_film_s[m] + v_film_g[g]);

		log_lik[n] = bernoulli_logit_lpmf( y[n] | t_i );
		theta[n] = inv_logit( t_i );
	}
}
