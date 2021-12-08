// generated with brms 2.10.0
functions {

  //prior for a unit_vector estimate of mean angle
   real von_mises_unitvector_lpdf(vector in_vec, real mu, real kappa) {
        real angle = atan2(in_vec[2], in_vec[1]);
        //extreme kappa correction for numeric stability
     if (kappa < 100) {
       return von_mises_lpdf(angle | mu, kappa);
     } else {
       real sd = sqrt(1 / kappa) +1e-16; // avoid sd ~= 0
       return normal_lpdf(angle | mu, sd);
     }
   }
  /* von Mises log-PDF of a single response
   * for kappa > 100 the normal approximation is used
   * for reasons of numerial stability
   * Args:
   *   y: the response vector between -pi and pi
   *   mu: location parameter vector
   *   kappa: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real von_mises_real_lpdf(real y, real mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(y | mu, kappa);
     } else {
       real sd = sqrt(1 / kappa) +1e-16; // avoid sd ~= 0
       return normal_lpdf(y | mu, sd);
     }
   }
   // calculate the mean angle of a circular distribution (in radians)
   real mean_circular(vector y, int N){
      real sumsin = 0;
      real sumcos = 0;
     for (n in 1:N){
       sumsin += sin(y[n]);
       sumcos += cos(y[n]);
     }
     sumsin = sumsin/N;
     sumcos = sumcos/N;
     return(atan2(  sumsin, sumcos) );
   }
    // calculate the mean vector length for a circular distribution
   real rho_circular(vector y, int N){
      real sumsin = 0;
      real sumcos = 0;
     for (n in 1:N){
       sumsin += sin(y[n]);
       sumcos += cos(y[n]);
     }
     sumsin = sumsin/N;
     sumcos = sumcos/N;
     return(sqrt(  sumsin^2 + sumcos^2)/N );
   }
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_kappa;  // number of population-level effects
  matrix[N, K_kappa] X_kappa;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_kappa_1;
  vector[N] Z_2_kappa_2;
  int<lower=1> NC_2;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
      // int Kc = K - 1;
      // matrix[N, Kc] Xc;  // centered version of X without an intercept
      // vector[Kc] means_X;  // column means of X before centering
      // int Kc_kappa = K_kappa - 1;
      // matrix[N, Kc_kappa] Xc_kappa;  // centered version of X_kappa without an intercept
      // vector[Kc_kappa] means_X_kappa;  // column means of X_kappa before centering
      // for (i in 2:K) {
      //   means_X[i - 1] = mean(X[, i]);
      //   Xc[, i - 1] = X[, i] - means_X[i - 1];
      // }
      // for (i in 2:K_kappa) {
      //   means_X_kappa[i - 1] = mean(X_kappa[, i]);
      //   Xc_kappa[, i - 1] = X_kappa[, i] - means_X_kappa[i - 1];
      // }
  //X is just two columns,
  //one of which is Intecept = 1,
  //the other of which is light_typeul = 0 || 1
  // It might be easier to just use the 2nd
  vector[N] Xc; // column of just the light_typeul
  int Kc_kappa = K_kappa - 1;  //kappa can be centred
  matrix[N, Kc_kappa] Xc_kappa;  // centered version of X_kappa without an intercept
  vector[Kc_kappa] means_X_kappa;  // column means of X_kappa before centering
  for (i in 2:K_kappa) {
    means_X_kappa[i - 1] = mean(X_kappa[, i]);
    Xc_kappa[, i - 1] = X_kappa[, i] - means_X_kappa[i - 1];
  }
  Xc = X[,2]; //I think this should work for mu
}
parameters {
  //all effects on mean vector need to be between -pi and pi!
  // since this is an angle, it may need a mu_vec transform
  // vector<lower=-pi(), upper=pi()>[Kc] b;  // population-level effects
  unit_vector [2] b_vec;  //  I'll use b_vec as a temporary contrast
      // temporary intercept for centered predictors
      // real Intercept;
  unit_vector [2] mu_vec;  //  I'll use mu_vec as a temporary intercept
      // vector[Kc_kappa] b_kappa;  // population-level effects
  vector[Kc_kappa] b_kappa;  // population-level effects
  // temporary intercept for centered predictors
  real Intercept_kappa;
    // vector<lower=0>[M_1] sd_1;  // group-level standard deviations
    // // standardized group-level effects
    // vector[N_1] z_1[M_1];
    // vector<lower=0>[M_1] sd_1;  // group-level standard deviations
    vector<lower=0>[M_1] sd_1;  // logkappa for group-level means
  // there is only one random effect (animal),
  //starting values may fall outside these of (-pi,pi)
  // vector<lower=-1, upper=1>[N_1] z_1[M_1];  // unscaled group-level effects
  array[M_1, N_1] unit_vector[2] z_1;  // standardized group-level effects
  // cholesky factor of correlation matrix; Unused?
  // cholesky_factor_corr[M_1] L_1;
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  matrix[M_2, N_2] z_2;  // standardized group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_2] L_2;
}
transformed parameters {
  // some kind of angle difference
  real b = atan2(b_vec[2], b_vec[1]);
  // some kind of angle
  real temp_angle = atan2(mu_vec[2], mu_vec[1]);
  // actual group-level effects
  // matrix[N_1, M_1] r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)';
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  // actual group-level effects
  matrix[N_2, M_2] r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)';
  // using vectors speeds up indexing in loops
  vector[N_2] r_2_kappa_1 = r_2[, 1];
  vector[N_2] r_2_kappa_2 = r_2[, 2];
  for (n in 1:N_1) {
  r_1_1[n] = atan2(z_1[n,1][2], z_1[n,1][1]);  //  Convert unit vector to arc
  r_1_2[n] = atan2(z_1[n,2][2], z_1[n,2][1]);  //  Convert unit vector to arc
  }
}
model {
  //population level circular mean estimate for mu_vec prior
  real cmean = mean_circular(Y,N);
  // initialize linear predictor term mu
  vector[N] mu = temp_angle + Xc * b;
  // initialize linear predictor term
  vector[N] kappa = Intercept_kappa + Xc_kappa * b_kappa;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
    mu[n] = atan2(sin(mu[n]), cos(mu[n])); // Keep mu within (-pi, pi)
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    kappa[n] += r_2_kappa_1[J_2[n]] * Z_2_kappa_1[n] + r_2_kappa_2[J_2[n]] * Z_2_kappa_2[n];
  }
  for (n in 1:N) {
    // apply the inverse link function
    kappa[n] = exp(kappa[n]);
  }
  // priors including all constants
  //prior likelihood is added to total log likelihood
  target += normal_lpdf(Intercept_kappa | 0, 3); //most likely values of log(kappa), 95% mean vectors in (0.001,0.999)
  // A strong prior that the population mean should fall at the
  // circular mean of the dataset aids convergence of SD and Z.
  target += von_mises_unitvector_lpdf(mu_vec| mean_circular(Y,N), 10); //95% prob. at �37�
  target += von_mises_unitvector_lpdf(b_vec| 0, 0.5); //weak bias to small shifts
  //  A prior for indiv. means on a vonmises distribution
  for (n in 1:N_1) {
  target += von_mises_unitvector_lpdf(z_1[n,1] | 0, exp(sd_1[1])); //pop. mean 0, logkappa is sd_1
  target += von_mises_unitvector_lpdf(z_1[n,2] | 0, exp(sd_1[2])); //pop. mean 0, logkappa is sd_1
  }
  // As for the main distribution, expect log kappa in (-6,6)
  target += normal_lpdf(sd_1 | 0, 3);
      // target += uniform_lpdf(z_1[1] | -1, 1); //rng values have 95% probability (-135.6°  135.6°)
      target += normal_lpdf(sd_2 | 0, 3)
        - 1 * normal_lccdf(0 | 0, 3);
      target += normal_lpdf(z_2[1] | 0, 1);
      // target += lkj_corr_cholesky_lpdf(L_1 | 1);
      //for the random effects slope on kappa
      target += normal_lpdf(to_vector(z_2) | 0, 1);
      target += lkj_corr_cholesky_lpdf(L_2 | 1);
  // likelihood including all constants
  if (!prior_only) {
    for (n in 1:N) {
      target += von_mises_real_lpdf(Y[n] | mu[n], kappa[n]);
    }
  }
}
generated quantities {
      // // actual population-level intercept
      // real b_Intercept = Intercept - dot_product(means_X, b);
      // // actual population-level intercept
      // real b_kappa_Intercept = Intercept_kappa - dot_product(means_X_kappa, b_kappa);
      // return population-level unit_vector sin component
      real b_sin1 = mu_vec[2];
      // return population-level unit_vector cos component
      real b_cos2 = mu_vec[1];
      // actual population-level intercept
      real b_Intercept = temp_angle;
      // real b_Intercept = temp_angle - dot_product(means_X, b);//normally this would be temp_Intercept
      // actual population-level intercept
      real b_kappa_Intercept = Intercept_kappa - dot_product(means_X_kappa, b_kappa);
      // group-level correlations
      // corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
      // vector<lower=-1,upper=1>[NC_1] cor_1;
      // group-level correlations
      corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
      vector<lower=-1,upper=1>[NC_2] cor_2;
      // extract upper diagonal of correlation matrix
      // for (k in 1:M_1) {
      //   for (j in 1:(k - 1)) {
      //     cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
      //   }
      // }
      // extract upper diagonal of correlation matrix
      for (k in 1:M_2) {
        for (j in 1:(k - 1)) {
          cor_2[choose(k - 1, 2) + j] = Cor_2[j, k];
        }
      }
}
