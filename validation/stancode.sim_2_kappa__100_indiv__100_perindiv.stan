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
  //brms version of numerically stable log probability density for a von mises distribution
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
   }    // calculate the mean vector length for a circular distribution
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
  int<lower=1> N;// number of observations
  vector[N] Y;// response variable, angles
  int<lower=1> K; 
// X is not assigned
  vector[N] Z_1_1;//offset from population mean in SDs
  int<lower=1> K_kappa; 
// X is not assigned
  vector[N] Z_2_kappa_1;//offset from population mean in SDs
  int<lower=1> J_1[N];//group identity
  int<lower=1> J_2[N];//group identity
  int<lower=1> N_1;//number of groups
  int<lower=1> M_1;//number of random effects SD
  int NC_1;// Mean non-centrality parameter? (not used?)
  int<lower=1> N_2;//number of groups
  int<lower=1> M_2;//number of random effects SD
  int NC_2;// Kappa non-centrality parameter? (not used?)
  int prior_only;// should the likelihood be ignored?

}
transformed data {
}
parameters {
//A unit_vector is used in place of: real temp_Intercept 
// this stops estimates of mu from walking around the circle multiple times 
// which is bad for chain convergence 
unit_vector [2] mu_vec;  //  A unit vector representing the population intercept angle 
// kappa is estimated on a log scale 
// so this intercept is actually log(kappa_intercept) 
real temp_kappa_Intercept;
vector<lower=0>[M_1] sd_1;  // logkappa for group-level means
array[N_1] unit_vector[2] z_1;  // group-level angle to pop mean
vector<lower=0>[M_2] sd_2;  // group-level standard deviations
vector[N_2] z_2[M_2] ;  // unscaled group-level effects
}
transformed parameters {
//The unit_vector is transformed to an angle: real temp_angle 
// this can be used for von mises likelihood calculations 
real temp_angle = atan2(mu_vec[2], mu_vec[1]);  //  An angle representing the population intercept angle 
//Group level effects on mean angle 
vector[N_1] r_1_1; //  The arc between pop. mean and group mean 
for (n in 1:N_1) { 
r_1_1[n] = atan2(z_1[n][2], z_1[n][1]);  //  Convert unit vector to arc 
} 
//Group level effects on log(kappa) 
vector[N_2] r_2_kappa_1 = (sd_2[1] * (z_2[1]));  //  The difference between pop. log kappa and group log kappa 
}
model {
//population level circular mean estimate for mu_vec prior 
real cmean = mean_circular(Y,N);  //  This can be used to place a prior on the intercept mean angle 
vector[N] mu = temp_angle + rep_vector(0, N); // Set up angle estimates for each data point 
vector[N] kappa = temp_kappa_Intercept + rep_vector(0, N); // Set up log kappa estimates for each data point 
//loop through data and transform according to the appropriate coefficients 
for (n in 1:N) { 
mu[n] += r_1_1[J_1[n]] * Z_1_1[n]; // At each data point, mean angle is pop mean + group mean * group index 
mu[n] = atan2(sin(mu[n]), cos(mu[n])); // Keep mu within (-pi, pi) 
kappa[n] += r_2_kappa_1[J_2[n]] * Z_2_kappa_1[n]; //  At each data point, lkappa is pop lkappa + group lkappa * group index 
kappa[n] = exp(kappa[n]);  // kappa is estimated via the log link, avoiding kappa < 0 
} 
//prior likelihood is added to total log likelihood 
target += normal_lpdf(temp_kappa_Intercept | 0, 3); //most likely values of log(kappa), 95% mean vectors in (0.001,0.999) 
// When using the SD*Z formulation for random effects on mean angle 
// A strong prior that the population mean should fall at the 
// circular mean of the dataset aids convergence of SD and Z. 
target += von_mises_unitvector_lpdf(mu_vec| mean_circular(Y,N), 10); //95% prob. at ±37° 
//  A prior for indiv. means on a vonmises distribution 
for (n in 1:N_1) { 
target += von_mises_unitvector_lpdf(z_1[n] | 0, exp(sd_1[1])); //pop. mean 0, logkappa is sd_1 
} 
// As for the main distribution, expect log kappa in (-6,6) 
target += normal_lpdf(sd_1 | 0, 3); 
//  Priors on indiv. differences in log kappa 
target += normal_lpdf(sd_2 | 0, 3) // 95% probability within ±6 log units of population kappa 
- 1 * normal_lccdf(0 | 0, 3); //needed for likelihood scaling (?) 
target += normal_lpdf(z_2[1] | 0, 1); //scaling of individual intercepts predicted on a standard normal 
// likelihood including all constants 
if (!prior_only) { // set prior_only to skip this step and check prior bias 
for (n in 1:N) { // for each data point 
  target += von_mises_real_lpdf(Y[n] | mu[n], kappa[n]); // calculate the von Mises log likelihood 
} 
} 
}
generated quantities {
//Angle values returned to the user 
// return population-level unit_vector sin component 
real b_sin1 = mu_vec[2]; 
// return population-level unit_vector cos component 
real b_cos2 = mu_vec[1]; 
// population-level mean angle 
real b_Intercept = temp_angle; 
//Kappa values returned to the user 
// Population-level log kappa 
real b_kappa_Intercept = temp_kappa_Intercept; 
}

