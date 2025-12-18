
functions {
      // https://discourse.mc-stan.org/t/custom-likelihood-for-zero-inflated-lognormal-model/15247

      // Arguments:
      // y: the response value
      // mu: mean parameter of the lognormal distribution
      // sigma: sd parameter of the lognormal distribution
      // zeta: zero-inflation probability
      // Returns:
      // a scalar to be added to the log posterior

      // zero-inflated lognormal log-PDF of a single response
      real hurdle_lognormal_lpdf(real y, real mu, real sigma, real zeta) {
            if (y == 0) {
                  return bernoulli_lpmf(0 | zeta);
            } else {
                  return bernoulli_lpmf(1 | zeta) + lognormal_lpdf(y | mu, sigma);
            }
      }

}

data {
      int<lower=0> N; // n records
      int<lower=0> P; // n predictor vars
      vector[N] y; // recruits
      vector[N] t; // years
      matrix[N, P] x; // predictors
}

transformed data {
      matrix[N, P + 1] x1 = append_col(rep_vector(1, N), x);
}

parameters {
      row_vector[P+1] zeta;
      row_vector[P+1] mu;
      real<lower=0> sigma;
}

model {
      for(i in 1:N) {
            real zeta_i = inv_logit(dot_product(zeta, x1[i,]));
            zeta_i = 1 - ((1 - zeta_i) ^ t[i]); // convert from annual prob to multiyear

            real mu_i = dot_product(mu, x1[i,]);
            target += hurdle_lognormal_lpdf(y[i] | mu_i, sigma, zeta_i);
      }

      zeta[1] ~ normal(0, 10);
      for(i in 2:(P+1)) zeta[i] ~ normal(0, 1);

      mu[1] ~ normal(0, 10);
      for(i in 2:(P+1)) mu[i] ~ normal(0, 1);
      sigma ~ normal(0, 1);
}
