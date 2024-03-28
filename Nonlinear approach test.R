#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 03 28
#     MODIFIED:	James Foster              DATE: 2024 03 28
#
#  DESCRIPTION: Simulate RE von Mises data and fit using the unit-vector transform
#               through a custom family and non-linear model.
#               
#       INPUTS: 
#               
#      OUTPUTS: A BRMS dataset (.Rdata), custom Stan code (.stan), and a fitted 
#               Stan model (.Rdata).
#
#	   CHANGES: - 
#             - 
#             - 
#
#   REFERENCES: Batschelet E (1981).
#               Graphical presentation, Chap 1.2, p. 4-6
#               Chapter 1: Measures of Location
#               In: Circular Statistics in Biology
#               Academic Press (London)
#               
#               Bürkner, P.-C., 2018. 
#               Advanced Bayesian Multilevel Modeling with the R Package brms.
#               The R Journal 10, 395–411. 
#               https://doi.org/10.32614/RJ-2018-017
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Simulate data   
#- Get modelling consistent
#- Save results  

# Check environment -------------------------------------------------------

#check the operating system and assign a logical flag (TRUE or FALSE)
sys_win <- Sys.info()[['sysname']] == 'Windows'

# Load packages -----------------------------------------------------------
#Package for transferring data from R to Stan
#Relies on System C compiler, may require special installation
library(cmdstanr) # observe startup messages #Use development verion: remotes::install_github("stan-dev/rstan", ref = "develop", subdir = "rstan/rstan")
#Package for managing circular data
library(circular)
#Package for automatic setup of Stan models and data
library(brms)
#For convenience, leave one CPU for the user to continue their work
options(mc.cores = parallel::detectCores()-1)

# Simulate data (not used) ------------------------------------------------
n_indiv = 20
n_trials = 10
mu_offset = circular(x = rad(-30) )
# minimum discriminable angle appears to be approx 35°
kappa_both = A1inv(0.7) #concentration around each trial mean
kappa_indiv = A1inv(0.5) #concentration across individuals (pairs)
#mean angle in trail 1 for each individual (pair)
mu1_sim = rvonmises(n = n_indiv,
                      mu = rcircularuniform(1),#random angle
                      kappa = kappa_indiv#the wider the distribution of individual biases, the greater the influence of pairing
                      )
#simulate the full dataset
sim = data.frame(
                 angle = round(c(suppressWarnings( #rvonmises converts to circular and warns
                   sapply(X = mu1_sim,
                          FUN = rvonmises,
                          n = n_trials,
                           kappa = kappa_both
                   )
                 ))*180/pi),#convert to angles and round to nearest degree
                 angle_2 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                             sapply(X =mu1_sim + mu_offset,# true difference,
                                    FUN = rvonmises,
                                    n = n_trials,
                                    kappa = kappa_both
                             )
                 ))*180/pi), #convert to angles and round to nearest degree
                 ID = sort(rep(sample(LETTERS, 
                                      size = n_indiv, 
                                      replace = FALSE),
                               times = n_trials))
                 )
#save somewhere the user likely keeps data
write.table(x = sim,
            file = file.path(ltp,'Documents', "simulated_REangles.csv"),
            sep = csv_sep,
            row.names = FALSE
            )


# Plot data and ML estimates ----------------------------------------------


par(mfrow = rep(ceiling(sqrt(n_indiv)), times = 2),
    mar = c(0,0,0,0))
for(ii in unique(sim$ID))
{
  cid = circular(subset(sim, ID %in% ii)$angle, template = 'geographics')
  plot.circular(x = cid, 
                stack  = T, 
                sep = 0.1,
                shrink = 1.2,
                bins = 72)
  arrows.circular(x = circular(mean.circular(cid), template = 'geographics'),
                  y = 1-rho.circular(cid))
}

# Define a custom unit vector von Mises family ---------------------------

#no real progress here

      # von_mises_unitvector <- custom_family(name = "von_mises_unitvector", 
      #                                       dpars = c("mu", "kappa"),
      #                                       links = c("identity", "log"),
      #                                       lb = c(-pi, 0), 
      #                                       ub = c(pi, NA),
      #                                       type = "real",
      #                                       vars = "vreal[n]"
      # )
      # # 
      # stan_funs <- stanvar(scode =
      #               "
      #                   //prior for a unit_vector estimate of mean angle
      #                  real von_mises_unitvector_lpdf(vector in_vec, real mu, real kappa) {
      #                       real angle = atan2(in_vec[2], in_vec[1]);
      #                       //extreme kappa correction for numeric stability
      #                    if (kappa < 100) {
      #                      return von_mises_lpdf(angle | mu, kappa);
      #                    } else {
      #                      real sd = sqrt(1 / kappa) +1e-16; // avoid sd ~= 0
      #                      return normal_lpdf(angle | mu, sd);
      #                    }
      #                  }
      #               ", block = "functions")
      # # stan_funs <- stanvar(scode = 
      # #               "
      # #                   //prior for a unit_vector estimate of mean angle
      # #                  real von_mises_unitvector_lpdf(vector in_vec, real mu, real kappa, real indmu, real indkappa) {
      # #                       real angle = atan2(in_vec[2], in_vec[1]);
      # #                       //extreme kappa correction for numeric stability
      # #                    if (kappa < 100) {
      # #                      return von_mises_lpdf(angle | mu + indmu, kappa + indkappa);
      # #                    } else {
      # #                      real sd = sqrt(1 / (kappa + indkappa) +1e-16; // avoid sd ~= 0
      # #                      return normal_lpdf(angle | mu + indmu, sd);
      # #                    }
      # #                  }
      # #               ", block = "functions")
      # 
      # stan_parm <- stanvar(scode = 
      #               "
      #                   //unit_vector estimate of mean angle
      #                  unit_vector[2] muv;
      #                  real<lower=0> kappa;
      #               ", block = "parameters")
      # 
      # stan_tran <- stanvar(scode = 
      #               "
      #                   //mean angle as an angle
      #                  real<lower = -pi(), upper = pi()> mu = atan2(muv[2], muv[1]);
      #               ", block = "tparameters")
      # 
      # stan_vars<- stan_parm + stan_tran
      # # frm_vm = bf(angle ~ mu, 
      # #             kappa ~ 1, 
      # #             family = von_mises(link = 'identity', link_kappa = 'log'))
      # 
      # frm_simple = bf(angle ~1, 
      #                 kappa ~1, 
      #                 family = von_mises(link = 'identity', link_kappa = 'log') )
      # 
      # scode_vm = make_stancode(formula = frm_simple,
      #                          data = sim_rad,
      #                          stanvars = stan_vars,
      #                          family = von_mises(link = 'identity', link_kappa = 'log'))
      # brm_vm = brm(formula = frm_simple,
      #                          data = sim_rad,
      #                          stanvars = stan_vars,
      #                          family = von_mises(link = 'identity', link_kappa = 'log'))
      # 
      # sim_rad = within(sim,
      #                  {
      #                  angle = atan2(cos(rad(angle)), sin(rad(angle))) 
      #                  angle_2 = atan2(cos(rad(angle_2)), sin(rad(angle_2))) 
      #                  }
      #                  )
      # frm_nl = bf(angle ~ mangle, 
      #             kappa ~ 1, 
      #             family = von_mises(link = 'identity', link_kappa = 'log')) + 
      #   nlf(mangle ~ atan2(muv[2], muv[1]),
      #       muv ~ 1)
      # 
      #                 
      # prs = get_prior(frm_nl, data = sim_rad, nl  = TRUE)
      # prs = within(prs, 
      #              {
      #                prior[nlpar == 'kappa'] = "normal(0,3)"
      #                prior[nlpar == 'muv'] = "uniform(-1,1)"
      #                  }
      #              )
      # 
      # scode_nl = make_stancode(formula = frm_nluvvm,
      #                               data = sim_rad,
      #                          family = von_mises(link = 'identity', link_kappa = 'log'))
      # 
      # brm_nl = brm(frm_nl, data = sim, prior = prs)
      # 
      # sc_simple = '
      # data{
      #     int n;
      #     real <lower=-pi(), upper=pi()> angle[n];
      #     }
      # parameters{
      #     unit_vector[2] muangleuv;
      #     real<lower=0> kappa;
      # }
      # transformed parameters{
      #     real<lower = -pi(), upper = pi()> muangle = atan2(muangleuv[2], muangleuv[1]);
      # }
      # model{
      #     kappa~exponential(1);
      #     angle~von_mises(muangle,kappa);
      # }'
      # 
      # fit_simple = brm(formula = bf(angle ~ 1,
      #                               kappa ~ 1),
      #                  data = sim,
      #                  stanvars = stanvar(scode = sc_simple),
      #                  family = von_mises(link = 'identity', link_kappa = 'identity')
      #                  )
      # 
      # 
      # frm_nluvvm = bf(angle ~ 1,
      #                 kappa ~ 1)
      # 
      # fit_nl_uvvm <- brm(formula = frm_nluvvm,
      #                    data = sim,
      #                    family = frm_nluvvm,
      #                    stanvars = stan_vars,
      #                    
      # )
      # 
      # frm_nl = bf(angle ~ 1,
      #             kappa ~ 1,
      #             )
