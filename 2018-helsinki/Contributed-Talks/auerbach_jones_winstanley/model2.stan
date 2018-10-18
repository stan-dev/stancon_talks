data {
  int<lower=1> N;                   // number of observations
  int<lower=1> T;                   // number of tracts
  int<lower=1> Z;                   // number of zones
  int<lower=1> H;                   // number of neighorhoods
  int<lower=1> L;                   // number of land
  int<lower=1> G;                   // number of grades
  int<lower=1> Y;                   // number of years
  int<lower=1> C;                   // number of cohorts
  int num_students[N];              // response var
  int<lower=1,upper=T> tract[N];    // tract id's
  int<lower=1,upper=Z> zone[T];     // zone id's
  int<lower=1,upper=H> nbhd[T];     // neighborhood id's
  int<lower=1,upper=L> land[T];     // land id's
  int<lower=1,upper=G> grade[N];    // grade id's
  int<lower=1,upper=Y> year[N];     // year id's
  int<lower=1,upper=C> cohort[N];   // cohort id's
}

parameters {
  real mu[T];
  real<lower=0> sigma;
  real<lower=0> tau_grade_zone[G];
  real<lower=0> tau_grade_nbhd[G];
  real<lower=0> tau_grade_land[G];
  real<lower=0> tau_year_zone[Y];
  real<lower=0> tau_year_nbhd[Y];
  real<lower=0> tau_year_land[Y];
  real<lower=0> tau_cohort_zone[C];
  real<lower=0> tau_cohort_nbhd[C];
  real<lower=0> tau_cohort_land[C];
  real beta_grade_zone[G, Z];
  real beta_grade_nbhd[G, H];
  real beta_grade_land[G, L];
  real beta_year_zone[Y, Z];
  real beta_year_nbhd[Y, H];
  real beta_year_land[Y, L];  
  real beta_cohort_zone[C, Z];
  real beta_cohort_nbhd[C, H];
  real beta_cohort_land[C, L];
}

transformed parameters {
  real theta[N];
  for(n in 1:N)
    theta[n] = sigma             * mu[tract[n]]                                + 
      tau_grade_zone[grade[n]]   * beta_grade_zone[grade[n],   zone[tract[n]]] +
      tau_grade_nbhd[grade[n]]   * beta_grade_nbhd[grade[n],   nbhd[tract[n]]] +
      tau_grade_land[grade[n]]   * beta_grade_land[grade[n],   land[tract[n]]] + 
      tau_year_zone[year[n]]     * beta_year_zone[year[n],     zone[tract[n]]] +
      tau_year_nbhd[year[n]]     * beta_year_nbhd[year[n],     nbhd[tract[n]]] +
      tau_year_land[year[n]]     * beta_year_land[year[n],     land[tract[n]]] +
      tau_cohort_zone[cohort[n]] * beta_cohort_zone[cohort[n], zone[tract[n]]] +
      tau_cohort_nbhd[cohort[n]] * beta_cohort_nbhd[cohort[n], nbhd[tract[n]]] +
      tau_cohort_land[cohort[n]] * beta_cohort_land[cohort[n], land[tract[n]]];
}

model {

 mu ~ normal(0, 1);

  for(g in 1:G) {
    beta_grade_zone[g, ]  ~ normal(0, 1);
    beta_grade_nbhd[g, ]  ~ normal(0, 1);
    beta_grade_land[g, ]  ~ normal(0, 1);
  }
  
  for(y in 1:Y) {
    beta_year_zone[y, ]   ~ normal(0, 1);
    beta_year_nbhd[y, ]   ~ normal(0, 1);
    beta_year_land[y, ]   ~ normal(0, 1);
  }
  
  for(c in 1:C) {
    beta_cohort_zone[c, ] ~ normal(0, 1);
    beta_cohort_nbhd[c, ] ~ normal(0, 1);
    beta_cohort_land[c, ] ~ normal(0, 1);
  }
  
tau_grade_zone          ~ gamma(2, 1/G);
tau_grade_nbhd          ~ gamma(2, 1/G);
tau_grade_land          ~ gamma(2, 1/G);
tau_year_zone           ~ gamma(2, 1/Y);
tau_year_nbhd           ~ gamma(2, 1/Y);
tau_year_land           ~ gamma(2, 1/Y);
tau_cohort_zone         ~ gamma(2, 1/Z);
tau_cohort_nbhd         ~ gamma(2, 1/Z);
tau_cohort_land         ~ gamma(2, 1/Z);
  
  num_students            ~ poisson_log(theta);
}
  