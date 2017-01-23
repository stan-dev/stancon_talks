data {
  int<lower=1> N_train;// number of training data points
  int<lower=1> N;      // number of total data points
  int<lower=1> G;      // number of groupings
  int<lower=3> J[G];   // group sizes
  int COND[N];         // index for surface condition/inclement weather
  int CITY[N];         // index for city
  int SLIM[N];         // index for posted speed limit
  int SIGN[N];         // index for signs and signals i.e. school zone/work zone
  int LGHT[N];         // index for light and time
  int BLTE[N];         // index for built environment
  int CITYxSLIM[N];
  int LGHTxSLIM[N];
  int BLTExSLIM[N];
  int CONDxSLIM[N];
  int LGHTxCOND[N];
  int LGHTxCONDxSLIM[N];
  vector[N] EXPR;      // population exposed
  int count[N];        // number of pedestrian deaths
}
transformed data {
  vector[N] offset;
  offset = log(EXPR);
}
parameters {
  real offset_e;
  vector[J[1]] COND_eta;
  vector[J[2]] CITY_eta;
  vector[J[3]] SLIM_eta;
  vector[J[4]] SIGN_eta;
  vector[J[5]] LGHT_eta;
  vector[J[6]] BLTE_eta;
  vector[J[7]] CITYxSLIM_eta;
  vector[J[8]] LGHTxSLIM_eta;
  vector[J[9]] BLTExSLIM_eta;
  vector[J[10]] CONDxSLIM_eta;
  vector[J[11]] LGHTxCOND_eta;
  vector[J[12]] LGHTxCONDxSLIM_eta;
  vector<lower=0>[G + 1] sds;
  vector[N_train] cell_eta;
  real mu;
}
transformed parameters {
  vector[J[1]]  COND_e;
  vector[J[2]]  CITY_e;
  vector[J[3]]  SLIM_e;
  vector[J[4]]  SIGN_e;
  vector[J[5]]  LGHT_e;
  vector[J[6]]  BLTE_e;
  vector[J[7]]  CITYxSLIM_e;
  vector[J[8]]  LGHTxSLIM_e;
  vector[J[9]]  BLTExSLIM_e;
  vector[J[10]] CONDxSLIM_e;
  vector[J[11]] LGHTxCOND_e;
  vector[J[12]] LGHTxCONDxSLIM_e;
  vector[N_train]     cell_e;
  vector[N_train]     mu_indiv;

  COND_e           = sds[1]   * COND_eta;
  CITY_e           = sds[2]   * CITY_eta;
  SLIM_e           = sds[3]   * SLIM_eta;
  SIGN_e           = sds[4]   * SIGN_eta;
  LGHT_e           = sds[5]   * LGHT_eta;
  BLTE_e           = sds[6]   * BLTE_eta; 
  CITYxSLIM_e      = sds[7]   * CITYxSLIM_eta;
  LGHTxSLIM_e      = sds[8]   * LGHTxSLIM_eta;
  BLTExSLIM_e      = sds[9]   * BLTExSLIM_eta;
  CONDxSLIM_e      = sds[10]  * CONDxSLIM_eta;
  LGHTxCOND_e      = sds[11]  * LGHTxCOND_eta;
  LGHTxCONDxSLIM_e = sds[12]  * LGHTxCONDxSLIM_eta;
  cell_e           = sds[G+1] * cell_eta;

  for(n in 1:N_train)
  mu_indiv[n] = mu + offset_e * offset[n]
                   + COND_e[COND[n]]
                   + CITY_e[CITY[n]]
                   + SLIM_e[SLIM[n]]
                   + SIGN_e[SIGN[n]]
                   + LGHT_e[LGHT[n]]
                   + BLTE_e[BLTE[n]]
                   + CITYxSLIM_e[CITYxSLIM[n]]
                   + LGHTxSLIM_e[LGHTxSLIM[n]]
                   + BLTExSLIM_e[BLTExSLIM[n]]
                   + CONDxSLIM_e[CONDxSLIM[n]]
                   + LGHTxCOND_e[LGHTxCOND[n]]
                   + LGHTxCONDxSLIM_e[LGHTxCONDxSLIM[n]]
                   + cell_e[n];
}
model {
  COND_eta           ~ normal(0,1);
  CITY_eta           ~ normal(0,1);
  SLIM_eta           ~ normal(0,1);
  SIGN_eta           ~ normal(0,1);
  LGHT_eta           ~ normal(0,1);
  BLTE_eta           ~ normal(0,1);
  CITYxSLIM_eta      ~ normal(0,1);
  LGHTxSLIM_eta      ~ normal(0,1);
  BLTExSLIM_eta      ~ normal(0,1);
  CONDxSLIM_eta      ~ normal(0,1);
  LGHTxCOND_eta      ~ normal(0,1);
  LGHTxCONDxSLIM_eta ~ normal(0,1);
  cell_eta           ~ normal(0,1);
  offset_e           ~ normal(0,1);
  sds                ~ normal(0,1);
  mu                 ~ normal(0,10);
  
 for (n in 1:N_train){
    target += poisson_log_lpmf(count[n] | mu_indiv[n]);
    target += -log1m_exp(-exp(mu_indiv[n]));
 }
}
generated quantities {
  real COND_sd;
  real CITY_sd;
  real SLIM_sd;
  real SIGN_sd;
  real LGHT_sd;
  real BLTE_sd;
  real CITYxSLIM_sd;
  real LGHTxSLIM_sd;
  real BLTExSLIM_sd;
  real CONDxSLIM_sd;
  real LGHTxCOND_sd;
  real LGHTxCONDxSLIM_sd;
  real cell_sd;
  vector[N - N_train] mu_indiv_pred;
  vector[N - N_train] cell_e_pred;

  COND_sd           = sd(COND_e);
  CITY_sd           = sd(CITY_e);
  SLIM_sd           = sd(SLIM_e);
  SIGN_sd           = sd(SIGN_e);
  LGHT_sd           = sd(LGHT_e);
  BLTE_sd           = sd(BLTE_e);
  CITYxSLIM_sd      = sd(CITYxSLIM_e);
  LGHTxSLIM_sd      = sd(LGHTxSLIM_e);
  BLTExSLIM_sd      = sd(BLTExSLIM_e);
  CONDxSLIM_sd      = sd(CONDxSLIM_e);
  LGHTxCOND_sd      = sd(LGHTxCOND_e);
  LGHTxCONDxSLIM_sd = sd(LGHTxCONDxSLIM_e);
  cell_sd           = sd(cell_e);
  
  for (n in 1:(N - N_train)){
    cell_e_pred[n]     = normal_rng(0, sds[G+1]);
    mu_indiv_pred[n] = mu + offset_e * offset[N_train + n] 
                          + COND_e[COND[N_train + n]]
                          + CITY_e[CITY[N_train + n]]
                          + SLIM_e[SLIM[N_train + n]]
                          + SIGN_e[SIGN[N_train + n]]
                          + LGHT_e[LGHT[N_train + n]]
                          + BLTE_e[BLTE[N_train + n]]
                          + CITYxSLIM_e[CITYxSLIM[N_train + n]]
                          + LGHTxSLIM_e[LGHTxSLIM[N_train + n]]
                          + BLTExSLIM_e[BLTExSLIM[N_train + n]]
                          + CONDxSLIM_e[CONDxSLIM[N_train + n]]
                          + LGHTxCOND_e[LGHTxCOND[N_train + n]]
                          + LGHTxCONDxSLIM_e[LGHTxCONDxSLIM[N_train + n]]
                          + cell_e_pred[n];
 }
}
