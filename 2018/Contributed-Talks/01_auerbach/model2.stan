functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) 
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]); 
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) / 
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) + 
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int N;                    // number of data points
  int num_knots;            // number of knots
  vector[num_knots] knots;  // location of knots
  int spline_degree;        // the degree of spline
  int<lower=0> D;           // number of days
  int<lower=0> O;           // number of officers
  int<lower=0> W;           // number of day of week
  int<lower=0> M;           // number of months
  int<lower=0> H;           // number of holidays
  int day[N];
  int officer[N]; 
  int week[N]; 
  int month[N];
  int holiday[N];
  int y[N];                 // number of tickets written per officer per day
  real X[N];                // relative productivity of officer in prev. period
  int z[D];                 // number of crashes per day
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; 
  matrix[num_basis, N] B;            
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; 
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
    B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  B[num_knots + spline_degree - 1, N] = 1; 
}

parameters {
  row_vector[num_basis] a_raw; 
  real a0;
  real<lower=0> tau;
  
  real beta_0;
  real beta_1[O];
  real beta_2[W];
  real beta_3[M];
  real beta_4[H];
  vector[D] epsilon;
  
  real<lower=0> sigma;
  real<lower=0> sigma_0;
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;
  real<lower=0> sigma_3;
  real<lower=0> sigma_4;
}

transformed parameters {
  row_vector[num_basis] a;
  vector[N] lambda;
  a[num_basis] = 0;
  for (i in 1:(num_basis - 1))
      a[num_basis - i] = a[num_basis - i + 1] + a_raw[num_basis - i] * tau;   
  lambda = a0 * to_vector(X) + to_vector(a * B);
}

model {
  // Prior
  a_raw ~ normal(0, 1);
  a0    ~ normal(0, 1);
  tau   ~ normal(0, 1);
  
  beta_0  ~ normal(0, 1);
  beta_1  ~ normal(0, 1);
  beta_2  ~ normal(0, 1);
  beta_3  ~ normal(0, 1);
  beta_4  ~ normal(0, 1);
  epsilon ~ normal(0, 1);
  
  sigma   ~ chi_square(1);
  sigma_0 ~ chi_square(1);
  sigma_1 ~ chi_square(1);
  sigma_2 ~ chi_square(1);
  sigma_3 ~ chi_square(1);
  sigma_4 ~ chi_square(1);
  
  //Likelihood
  for(n in 1:N) 
    y[n] ~ poisson_log(lambda[n] + 
                        sigma   * epsilon[day[n]] +
                        sigma_1 * beta_1[officer[n]] +
                        sigma_2 * beta_2[week[n]] + 
                        sigma_3 * beta_3[month[n]] +
                        sigma_4 * beta_4[holiday[n]]);
  z ~ poisson_log(sigma_0 * epsilon + beta_0);
}
