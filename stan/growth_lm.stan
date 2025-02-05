
data {
      int<lower=0> N;
      int<lower=0> P; // n preds, excl intercept
      vector[N] y;
      matrix[N, P] x; // predictors
}

transformed data {
      matrix[N, P + 1] x1 = append_col(rep_vector(1, N), x);
}

parameters {
      vector[P+1] beta;
      real<lower=0> sigma;
}

model {
      y ~ normal(x1 * beta, sigma);
}
