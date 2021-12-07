#Open file with default program on any OS
# https://stackoverflow.com/a/35044209/3745353
shell.exec.OS  <- function(x){
  # replacement for shell.exec (doesn't exist on MAC)
  if (exists("shell.exec",where = "package:base"))
  {return(base::shell.exec(x))}else
  {comm <- paste0('open "',x,'"')
  return(system(comm))}
}

# Set up circular formats -----------------------------------------------
pipi0 = list(units = 'radians',
             type = 'angles',
             modulo = '2pi',
             zero = 0,
             rotation = 'clock',
             template = 'none')

deg360 = list(units = 'degrees',
              type = 'angles',
              modulo = '2pi',
              zero = 0,
              rotation = 'clock',
              template = 'none')

# Test C++ ----------------------------------------------------------------
CPPtest = function(sys_win = Sys.info()[['sysname']] == 'Windows',
                   hm = ifelse(test = sys_win,
                               yes =  gsub('\\\\', '/', Sys.getenv('USERPROFILE')),#use user profile instead of home directory
                               no =  Sys.getenv('HOME')#Use home directory
                   )
)
{
  hw = file.path(hm,"helloworld.cpp")
  if(file.exists(hw)){file.remove(hw)}else{file.create(hw)}
  cat(paste(
    c('#include <RcppArmadillo.h>',   
      '// [[Rcpp::depends(RcppArmadillo)]]',
      '// [[Rcpp::export]]',
      'void hello_world() {',
      '  Rcpp::Rcout << "C++ is working!" << std::endl;  ',
      '}',
      '// After compile, this function will be immediately called using',
      '// the below snippet and results will be sent to the R console.',
      '/*** R',
      'hello_world() ',
      '*/'),
    "
"),
    file = hw, 
    sep = '\n',
    append = T)
  tryCatch(Rcpp::sourceCpp(hw),
           error = function(e)
           {if(!sys_win)
           {message('R compiler tools need to be reinstalled\n',
                    'see: ', 
                    'thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/',
                    '\nfor more information')}else
                    {message("You need to fix your Windows C++ compiler... ", 'good luck!')}
           }
  )
}



# Generate circular parameters based on input data ------------------------
#Generate parameters for each individual
ParGenerator = function(i,
                        params,
                        nn = 1)
{
  list(mu = switch(EXPR = params$mu,#mu can be "uniform" or "vonmises" (i.e. distribution)
                   uniform = runif(
                     n = nn,
                     min = 0,
                     max = 360-.Machine$double.xmin),
                   vonmises = circular::rvonmises(
                     n = nn, 
                     mu = as.circular(x = 0,
                                      control.circular = deg360),
                     # kappa = (params$sd*pi/180)^-2),#kappa ~= 1/(sd)^2
                     kappa = A1inv(exp(((params$sd*pi/180)^2)/(-2))),
                     control.circular = deg360),#sd = sqrt(-2 * log(r))
                   runif(
                     n = nn,
                     min = 0,
                     max = 360-.Machine$double.xmin)
  ),
  kappa = exp( rep(x = log(params$kappa), times = nn) + 
                 rnorm(n = nn, mean = 0, sd = params$sd_logkappa) )
  )
}

#generate full distribution from those parameters
RVMgenerator = function(param,
                        nn = 20,
                        cc = list(units = 'degrees',
                                  type = 'angles',
                                  modulo = '2pi',
                                  zero = 0,
                                  rotation = 'clock',
                                  template = 'none'))
{
  xx= circular::rvonmises(n = nn,
                          mu = circular::as.circular(x = param$mu,
                                                     control.circular = cc), 
                          kappa = param$kappa)-pi
  return(xx)
}


# Write Stan code ---------------------------------------------------------
#BRMS does not know how to write a von Mises distribution as a unit vector
MakeStan_MEvonMises = function(brms_data,
                               mu_prior = 'vonmises', # 'uniform' # 'uniform' or 'vonmises' 
                               mu_bar_kappa_prior = 10, #a strong prior that the population mean is within 37° of the sample mean
                               mu_raneff = 'uniform' # 'vonmises' # random intercepts in uniform arc or a vonmises distribution
)
{
  #Write the separate components and combine them
  sc <- list()
  # . Functions block -------------------------------------------------------
  
  fun.prior.uv_vonmises <- c(
    "  //prior for a unit_vector estimate of mean angle
   real von_mises_unitvector_lpdf(vector in_vec, real mu, real kappa) {
        real angle = atan2(in_vec[2], in_vec[1]);
        //extreme kappa correction for numeric stability
     if (kappa < 100) {
       return von_mises_lpdf(angle | mu, kappa);
     } else {
       real sd = sqrt(1 / kappa) +1e-16; // avoid sd ~= 0
       return normal_lpdf(angle | mu, sd);
     }
   }\n"
  )
  # Avoids "Error in function 
  # boost::math::cyl_bessel_i<double>(double,double): 
  # numeric overflow"
  
  fun.lpdf.safe_vonmises <- c(
    "  //brms version of numerically stable log probability density for a von mises distribution
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
  "
  )
  
  fun.mang <- c(
    " // calculate the mean angle of a circular distribution (in radians)
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
   }"
  )
  
  fun.mvec <- c(
    "    // calculate the mean vector length for a circular distribution
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
   }"
  )
  
  sc$fun.block <- paste0(
    c("functions {\n",
      fun.prior.uv_vonmises, 
      fun.lpdf.safe_vonmises,
      fun.mang,
      fun.mvec,
      "\n}"),
    collapse = ''
  )
  
  # . Data block ------------------------------------------------------------
  # summary(brms_data)
  Dblock <- function(i,#index of data variable
                     x,#data to search in
                     nms)#names of data variables
  {
    if( is.null(dim(x[[i]])) )#if not a matrix or vector
    {
      if(any(x[[i]]>0))#could be single integer value
      {paste0("  ",
              "int<lower=1> ",
              nms[i],
              ";",
              Cmt(nms[i]),
              "\n")
      }else
      {#could be empty (unused)
        paste0("  ",
               "int ",
               nms[i],
               ";",
               Cmt(nms[i]),
               "\n")
      }
    }else
    {
      if(grepl(pattern = 'J',
               x = nms[i]))#J is a vector of group identities, all of which are integers
      {
        paste0("  ",
               "int<lower=1> ",
               nms[i],
               "[N]",
               ";",
               Cmt(nms[i]),
               "\n")
      }else
      {
        if(grepl(pattern = 'X',
                 x = nms[i]))
        {
          "// X is not assigned\n"
        }else{
          paste0("  ",
                 "vector[N] ",
                 nms[i],
                 ";",
                 Cmt(nms[i]),
                 "\n")
        }#for all variabiles of non-zero length that are not "J" or "X
        
      }
    }
  }
  #Function to add comments for each variable name
  Cmt <- function(nm)
  {switch(EXPR = nm,
          N = "// number of observations",
          Y = "// response variable, angles",
          N_1 = "//number of groups",
          M_1 = "//number of random effects SD",
          J_1 = "//group identity",
          Z_1_1 = "//offset from population mean in SDs",
          N_2 = "//number of groups",
          M_2 = "//number of random effects SD",
          J_2 = "//group identity",
          Z_2_kappa_1 = "//offset from population mean in SDs",
          prior_only = "// should the likelihood be ignored?",
          NC_1 = "// Mean non-centrality parameter? (not used?)",
          NC_2 = "// Kappa non-centrality parameter? (not used?)",
          " "
  )
  }
  #run the data-block construction function with necessary inputs
  dblock <- sapply(X = 1:length(brms_data), #index in data
                   FUN = Dblock, #writing function
                   x = brms_data, #data structure
                   nms = names(brms_data) #names of data variables
  )
  sc$data.block <- paste0(
    c(
      "data {\n",
      dblock,
      "\n}",
      "\ntransformed data {\n}"),#we don't use this block
    collapse = ''
  )
  
  # . Parameters block ------------------------------------------------------
  #Population mean angle
  par.intercept <- paste(
    c("//A unit_vector is used in place of: real temp_Intercept",
      "// this stops estimates of mu from walking around the circle multiple times",
      "// which is bad for chain convergence",
      "unit_vector [2] mu_vec;  //  A unit vector representing the population intercept angle"
    ),
    '\n',
    collapse = ''
  )
  
  #Population mean log(kappa)
  par.kappa <- paste0(
    "// kappa is estimated on a log scale \n",
    "// so this intercept is actually log(kappa_intercept) \n",
    "real temp_kappa_Intercept",
    ";",
    '\n'
  )
  
  #Random effects on mean angle and logkappa
  #Input variable mu_raneff can be: "uniform" or "vonmises"
  #Random mean angles can be drawn from a uniform arc (pop mean +-half arc)
  #Or from a vonmises distribution: vonmises(pop mean, exp(log kappa))
  nm.sd <- names(brms_data)[grep(x = names(brms_data),
                                 pattern = 'M_')
  ]#search for specified numbers of random effects
  par.sd <- paste0(
    sapply(X = 1:length(nm.sd),
           FUN = function(i, mu_raneff)
           {
             switch(EXPR = mu_raneff,
                    uniform = 
                      if(i==1){
                        paste0("vector<lower=0, upper=pi()>",
                               "[", nm.sd[i], "] ",
                               "sd_", i,
                               ";",
                               "  // ±arc containing group-level means", "\n",
                               "vector<lower=-1, upper=1>",
                               "[", "N_",i, "] ",
                               "z_",i,
                               "[", nm.sd[i], "] ",
                               ";",
                               "  // unscaled group-level effects", "\n"
                        )
                      }else
                      {
                        paste0("vector<lower=0>",
                               "[", nm.sd[i], "] ",
                               "sd_", i,
                               ";",
                               "  // group-level standard deviations", "\n",
                               "vector",
                               "[", "N_",i, "] ",
                               "z_",i,
                               "[", nm.sd[i], "] ",
                               ";",
                               "  // unscaled group-level effects", "\n"
                        )
                      },
                    vonmises = 
                      if(i==1){
                        paste0("vector<lower=0>",
                               "[", nm.sd[i], "] ",
                               "sd_", i,
                               ";",
                               "  // logkappa for group-level means", "\n",
                               "array",
                               "[", "N_",i, "] ",
                               "unit_vector[2] ",
                               "z_",i,
                               ";",
                               "  // group-level angle to pop mean", "\n"
                        )
                      }else
                      {
                        paste0("vector<lower=0>",
                               "[", nm.sd[i], "] ",
                               "sd_", i,
                               ";",
                               "  // group-level standard deviations", "\n",
                               "vector",
                               "[", "N_",i, "] ",
                               "z_",i,
                               "[", nm.sd[i], "] ",
                               ";",
                               "  // unscaled group-level effects", "\n"
                        )
                      }
             )
           },
           mu_raneff = mu_raneff
    ),
    collapse = ''
  )
  
  sc$par.block <- paste0(
    c("parameters {\n",
      par.intercept,
      par.kappa,
      par.sd,
      "}"),
    collapse = ''
  )
  
  
  # . Transformed parameters block ------------------------------------------
  
  tpar.intercept <- paste(
    c("//The unit_vector is transformed to an angle: real temp_angle",
      "// this can be used for von mises likelihood calculations",
      "real temp_angle = atan2(mu_vec[2], mu_vec[1]);  //  An angle representing the population intercept angle"
    ),
    '\n',
    collapse = ''
  )
  
  #TODO function for writing indiv. effects from data
  tpar.pm_mu_indiv <- paste(
    switch(EXPR = mu_raneff,
           uniform = 
             c("//Group level effects on mean angle",
               "vector[N_1] r_1_1 = (sd_1[1] * (z_1[1]));  //  The arc between pop. mean and group mean"
             ),
           vonmises = 
             c("//Group level effects on mean angle",
               "vector[N_1] r_1_1; //  The arc between pop. mean and group mean",
               "for (n in 1:N_1) {",
               "r_1_1[n] = atan2(z_1[n][2], z_1[n][1]);  //  Convert unit vector to arc",
               "}"  
             )
    ),
    '\n',
    collapse = ''
  )
  
  tpar.pm_kappa_indiv <- paste(
    c("//Group level effects on log(kappa)",
      "vector[N_2] r_2_kappa_1 = (sd_2[1] * (z_2[1]));  //  The difference between pop. log kappa and group log kappa"
    ),
    '\n',
    collapse = ''
  )
  
  sc$tpar.block <- paste0(
    c("transformed parameters {\n",
      tpar.intercept,
      tpar.pm_mu_indiv,
      tpar.pm_kappa_indiv,
      "}"),
    collapse = ''
  )
  
  # . Model block -----------------------------------------------------------
  
  mod.params <- paste(
    c("//population level circular mean estimate for mu_vec prior",
      "real cmean = mean_circular(Y,N);  //  This can be used to place a prior on the intercept mean angle",
      "vector[N] mu = temp_angle + rep_vector(0, N); // Set up angle estimates for each data point",
      "vector[N] kappa = temp_kappa_Intercept + rep_vector(0, N); // Set up log kappa estimates for each data point"
    ),
    '\n',
    collapse = ''
  )
  
  mod.loop <- paste(
    c("//loop through data and transform according to the appropriate coefficients",
      "for (n in 1:N) {",
      "mu[n] += r_1_1[J_1[n]] * Z_1_1[n]; // At each data point, mean angle is pop mean + group mean * group index",
      "mu[n] = atan2(sin(mu[n]), cos(mu[n])); // Keep mu within (-pi, pi)",
      "kappa[n] += r_2_kappa_1[J_2[n]] * Z_2_kappa_1[n]; //  At each data point, lkappa is pop lkappa + group lkappa * group index",
      "kappa[n] = exp(kappa[n]);  // kappa is estimated via the log link, avoiding kappa < 0",
      "}"
    ),
    '\n',
    collapse = ''
  )
  
  
  mod.priors_intercepts <- paste(
    c("//prior likelihood is added to total log likelihood",
      "target += normal_lpdf(temp_kappa_Intercept | 0, 3); //most likely values of log(kappa), 95% mean vectors in (0.001,0.999)",
      "// When using the SD*Z formulation for random effects on mean angle",
      "// A strong prior that the population mean should fall at the",
      "// circular mean of the dataset aids convergence of SD and Z.",
      switch(EXPR = mu_prior,
             uniform = "target += uniform_lpdf(temp_angle | -pi(), pi()); //any population mean within one circle has equal likelihood",
             # vonmises = "target += von_mises_unitvector_lpdf(mu_vec| mean_circular(Y,N), 10); //95% prob. at ±37°")
             vonmises = "target += von_mises_unitvector_lpdf(mu_vec| mean_circular(Y,N), 5); //95% prob. at ±54°")
    ),
    '\n',
    collapse = ''
  )
  
  mod.priors_indivs_mu <- paste(
    
    switch(EXPR = mu_raneff,
           uniform = c(
             "//  A prior limiting the indiv. means arc to ±half the circle",
             "target += uniform_lpdf(sd_1 | 0, pi())  //any half_arc in one half circle",
             "- 1 * uniform_lccdf(0 | 0, pi()); //needed for likelihood scaling (?)",
             "// An inverted normal. Strong prior that individual means are not at population mean.",
             "target += uniform_lpdf(z_1[1] | -1, 1) - normal_lpdf(z_1[1]  | 0, 0.25);"),
           vonmises = c(
             "//  A prior for indiv. means on a vonmises distribution",
             "for (n in 1:N_1) {",
             "target += von_mises_unitvector_lpdf(z_1[n] | 0, exp(sd_1[1])); //pop. mean 0, logkappa is sd_1",
             "}",
             "// As for the main distribution, expect log kappa in (-6,6)",
             "target += normal_lpdf(sd_1 | 0, 3);")
    ),
    '\n',
    collapse = ''
  )
  
  mod.priors_indivs_kappa <- paste(
    c(
      "//  Priors on indiv. differences in log kappa",
      "target += normal_lpdf(sd_2 | 0, 3) // 95% probability within ±6 log units of population kappa",
      "- 1 * normal_lccdf(0 | 0, 3); //needed for likelihood scaling (?)",
      "target += normal_lpdf(z_2[1] | 0, 1); //scaling of individual intercepts predicted on a standard normal"
    ),
    '\n',
    collapse = ''
  )
  
  mod.constants <- paste(
    c("// likelihood including all constants",
      "if (!prior_only) { // set prior_only to skip this step and check prior bias",
      "for (n in 1:N) { // for each data point",
      "  target += von_mises_real_lpdf(Y[n] | mu[n], kappa[n]); // calculate the von Mises log likelihood",
      "}",
      "}"),
    '\n',
    collapse = ''
  )
  
  sc$model.block <- paste0(
    c("model {\n",
      mod.params,
      mod.loop,
      mod.priors_intercepts,
      mod.priors_indivs_mu,
      mod.priors_indivs_kappa,
      mod.constants,
      "}"),
    collapse = ''
  )
  
  
  # . Generated quantities block --------------------------------------------
  
  gquant.mu <- paste(
    c("//Angle values returned to the user",
      "// return population-level unit_vector sin component",
      "real b_sin1 = mu_vec[2];",
      "// return population-level unit_vector cos component",
      "real b_cos2 = mu_vec[1];",
      "// population-level mean angle",
      "real b_Intercept = temp_angle;"
    ),
    '\n',
    collapse = ''
  )
  
  gquant.kappa <- paste(
    c("//Kappa values returned to the user",
      "// Population-level log kappa",
      "real b_kappa_Intercept = temp_kappa_Intercept;"
    ),
    '\n',
    collapse = ''
  )
  
  #Anything else?
  
  sc$gquant.block <- paste0(
    c("generated quantities {\n",
      gquant.mu,
      gquant.kappa,
      "}"),
    collapse = ''
  )
  
  # . Combine and save ------------------------------------------------------
  sc <- lapply(sc, #add carriage returns after each block
               function(x)
               {paste0(x,'\n')})
  #combine into a single text block
  st.code <- do.call(paste0, sc)
  
  return(st.code)
}
