data {
  // numbers of things
  int<lower=1> I;  // items
  int<lower=1> S;  // subjects
  int<lower=1> R;  // responses
  int<lower=1> G;  // grades
  // training data
  int<lower=1,upper=I> item[R];
  int<lower=1,upper=S> subject[R];
  int<lower=1,upper=G> grade[R];
}
parameters {
  // parameters
  ordered[G-1] difficulty[I];
  vector[S] ability;
  // hyperparameters
  real mu_first_difficulty;
  real<lower=0> sigma_ability;
  real<lower=0> sigma_first_difficulty;
  real<lower=0> sigma_step;
}
model {
  // priors
  ability ~ normal(0, sigma_ability);
  difficulty[1:I, 1] ~ normal(mu_first_difficulty,
                              sigma_first_difficulty); // prior for easiest grades
  for (i in 1:I){
    difficulty[i, 2:G-1]
    - difficulty[i, 1:G-2] ~ normal(0, sigma_step);  // priors for steps between grades
  }
  // data model
  for (response in 1:R){
    grade[response] ~ ordered_logistic(ability[subject[response]],
                                       difficulty[item[response]]);}
}
