data {
  int<lower=1> N;                   // number of observations
  int<lower=1> T;                   // number of tracts
  int<lower=1> G;                   // number of grades
  int<lower=1> Y;                   // number of years
  int<lower=1> C;                   // number of cohorts
  int num_students[N];              // response var
  int<lower=1,upper=T> tract[N];    // tract id's
  int<lower=1,upper=G> grade[N];    // grade id's
  int<lower=1,upper=Y> year[N];     // year id's
  int<lower=1,upper=C> cohort[N];   // cohort id's
}

parameters {
  real mu[T];
  real alpha_grade[G, T];
  real alpha_year[Y, T];
  real alpha_cohort[C, T];
  real<lower=0> tau_year[Y];
  real<lower=0> tau_grade[G];
  real<lower=0> tau_cohort[C];
  real<lower=0> sigma;
}

transformed parameters {
  real theta[N];
  for(n in 1:N)
    theta[n] = sigma * mu[tract[n]] + 
               tau_year[  year[n]]   * alpha_year[  year[n],   tract[n]] + 
               tau_grade[ grade[n]]  * alpha_grade[ grade[n],  tract[n]] + 
               tau_cohort[cohort[n]] * alpha_cohort[cohort[n], tract[n]];
}

model {
    mu ~ normal(0, 1);
    
    for(g in 1:G) alpha_grade[g,]  ~  normal(0, 1);
    for(y in 1:Y) alpha_year[y,]   ~  normal(0, 1);
    for(c in 1:C) alpha_cohort[c,] ~  normal(0, 1);
    
    num_students ~ poisson_log(theta);
    
}
