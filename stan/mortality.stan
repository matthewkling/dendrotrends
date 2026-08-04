
data {
      int<lower=0> N;
      int<lower=0> P; // n preds, excl intercept
      array[N] int y;
      matrix[N, P] x; // predictors
      vector[N] t;
}

transformed data {
      matrix[N, P + 1] x1 = append_col(rep_vector(1, N), x);
}

parameters {
      vector[P+1] beta;
}

model {
      y ~ bernoulli(1 - (1 - 1/(1 + exp(- x1 * beta))) .^ t);

      beta[1] ~ normal(0, 10);
      for(i in 2:(P+1)) beta[i] ~ normal(0, 1);
}
