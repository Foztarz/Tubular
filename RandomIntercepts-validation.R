#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2021 12 02
#     MODIFIED:	James Foster              DATE: 2021 12 08
#
#  DESCRIPTION: Loads functions, generates simulated dance angles and fits a 
#               von Mises model with individual effects. This model is then 
#               extracted and compared to the input parameters.
#               
#       INPUTS: An ".R" file with the required functions.
#               User should specify data parameters (line 100).
#               
#      OUTPUTS: A BRMS dataset (.Rdata), custom Stan code (.stan), and a fitted 
#               Stan model (.Rdata).
#
#	   CHANGES: - waitbar
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
#    EXAMPLES:  Fill out user input (lines 100-105), then press ctrl+shift+s to run
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Simulate data   +
#- Set up input - output comparison +
#- Batch run  +
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
#Package for flexible data reorganisation
library(bayesplot)
#For convenience, leave one CPU for the user to continue their work
options(mc.cores = parallel::detectCores()-1)

# Useful Functions --------------------------------------------------------
fun_path = tryCatch(expr = 
                      {file.path(dirname(sys.frame(1)$ofile), "Validation-functions.R")},
                    error = function(e){'Validation-functions.R'}
)
if(!file.exists(fun_path))
{
  msg <- 'Please select "Validation-functions.R"'
  # ask user to find data
  if(sys_win){#choose.files is only available on Windows
    message('\n\n',msg,'\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    fun_path  <- choose.files(
      default = file.path(gsub('\\\\', '/', Sys.getenv('USERPROFILE')),#user
                          'Documents'),#For some reason this is not possible in the "root" user
      caption = msg
    )
  }else{
    message('\n\n',msg,'\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    fun_path <- file.choose(new=F)
  }
}
#read in relevant functions
source(file = fun_path, 
       encoding = 'UTF-8')

# . Check C++ is working --------------------------------------------------
CPPtest()
cmdstanr::check_cmdstan_toolchain()








# User input --------------------------------------------------------------
no_indiv = c(5, 10, 30)#, 100) #number of individuals
nper_indiv = c(5, 10, 30)#, 100) # observations per individual
kappa = A1inv(c(0.3, 0.5, 0.7))# mean vector lengths
sd_logkappa = log(1+ A1inv(c(0.1, 0.2, 0.3)))# ~sd of mean vector length
sd_mu = c(NA, 15, 30, 90)#~ sd of individual means in ° (NA for uniform)
mu_dist = 'vonmises' # 'uniform' #distribution to use for modelling individual means
n_rep = 5 # number of replicates of each 
# . Batch setup -----------------------------------------------------------
par_in_df = expand.grid(
                        no_indiv = no_indiv,
                        nper_indiv = nper_indiv,
                        kappa = kappa,
                        sd_logkappa = sd_logkappa,
                        sd = sd_mu
                        )
par_in_df = within(par_in_df,
                   {
                     mu = ifelse(test = is.na(sd),
                                 yes = 'uniform',
                                 no = 'vonmises')
                   }
                   )

# Replicate conditions ----------------------------------------------------
par_in_df = do.call(what = rbind, 
                    args = lapply(X = 1:n_rep,
                                  FUN = function(i){return(par_in_df)}
                                  )
)
#set up output data
par_inout_df = par_in_df
progb = txtProgressBar(min = 0,
                    max = dim(par_in_df)[1],
                    initial = 0, 
                    char = "-",
                    width = NA, 
                    style = 3)
# Open the loop -----------------------------------------------------------
for(i in 1:dim(par_in_df)[1])
{
  
  par_in = par_in_df[i,]
  setTxtProgressBar(pb = progb,
                    value = i)
# Derived variables -------------------------------------------------------
#make a list of parameters for each individual
par_lst = lapply(X = 1:par_in$no_indiv,
                 FUN = ParGenerator,
                 params = par_in)
#find a location to save output
save_path = file.path(dirname(fun_path), 'validation')
if(!dir.exists(save_path)){dir.create(save_path)}
#naming based on parameters
cond_nm = paste0(round(par_in$kappa,1),'_kappa__',
                 par_in$no_indiv,'_indiv__', 
                 par_in$nper_indiv    , '_perindiv')

# Simulate data -----------------------------------------------------------
#function for generating data based on input parameters


# . Simulate data ---------------------------------------------------------
#generate bearings for each individual from the parameter list
fake_dta  = data.frame(bearing = c(sapply(X = par_lst,
                                          FUN = RVMgenerator,
                                          nn = par_in$nper_indiv
)
),
#individual IDs appear alongside bearings, in order
indiv = c(sapply(X = 1:par_in$no_indiv,
                    FUN = rep,
                    times = par_in$nper_indiv
)
)
)
#convert to format for Stan model
fake_dta = within(fake_dta,
                  {
                    angle = bearing * pi/180  #convert to radians
                    angle = atan2(sin(angle), cos(angle)) # convert to [-pi,pi]
                  }
)

# Summarise simulated data ------------------------------------------------
MLEparm = function(i, #individual number or name
                   data, #data.frame
                   cc = list(units = 'radians',#circular format
                             type = 'angles',
                             modulo = '2pi',
                             zero = 0,
                             rotation = 'clock',
                             template = 'none')
)
{
  with(subset(data, indiv == i),#for each individual
       {
         prm = circular::mle.vonmises(#maximum likelihood parameters
           x = as.circular(x = angle, control.circular = cc),
           bias = TRUE # correct for bias
         )
       }
  )
}
#calculate maximum likelihood parameters for simulated data
mle_parm = lapply(1:par_in$no_indiv,
                  FUN = MLEparm,
                  data = fake_dta)

# . Plot summary ----------------------------------------------------------
  # PltIndiv = function(i,#individual number or name 
  #                        data)#data.frame
  # {
  #   with(subset(data, indiv == i), #for each individual
  #        {plot.circular(x = as.circular(bearing, control.circular = deg360),
  #                       stack = TRUE,
  #                       sep = -0.04,
  #                       bins = (360-5)/5,
  #                       zero = pi/2)
  #                       })
  # }
  # #set up plot with enough space for all individuals
  # par(mfrow = rep(x = ceiling(sqrt(par_in$no_indiv)), times = 2),
  #     mar = c(0,0,0,0)#narrow margins for circular plots
  # )
  # #plot and hide loop output
  # invisible(
  #   {
  #     lapply(X = 1:par_in$no_indiv,#for each individual
  #            FUN = PltIndiv,#plotting function
  #            data = fake_dta #simulated data
  #     )
  #   }
  # )

# Make Stan data ------------------------------------------------------
# a Stan formula
form1 <- bf( angle ~ (1 |indiv), #angle varies by individual
             kappa ~ (1 |indiv) #concentration varies by individual
)			
#use the formula to organise the data for Stan
br_data <- brms::make_standata(
  formula = form1,
  data = fake_dta,
  family = von_mises
)
#file to save this data in
dt_path <- file.path(save_path,
                     paste0('br.sim_',
                            cond_nm,
                            '.Rdata')
)	
  # #save the data for later
  # save(x = br_data, 
  #      file = dt_path)	
  # #load it if it has been saved
  # load(file = dt_path)

#Make a custom Mixed effects von Mises for Stan
scode = MakeStan_MEvonMises(brms_data = br_data,
                            mu_raneff = mu_dist)
#where to save this code in ('.Stan')
cd_path = file.path(save_path,
                    paste0('stancode.sim_',
                           cond_nm,
                           '.stan')
)	
#save it as an ASCII file
write.table(x = scode,
            file = cd_path,
            col.names = F,
            row.names = F, 
            quote = F)
#I recommend downloading https://atom.io/ and setting it as default for '.stan' files
  # shell.exec.OS(cd_path)


# Run model on simulated data ---------------------------------------------

# . Set up model ----------------------------------------------------------


shortrun = TRUE #FALSE # 

stan_mod = cmdstanr::cmdstan_model(cd_path)



# . Run model -------------------------------------------------------------

message('----  ', 'Running model with parameters: ', cond_nm, '  ----')
#beware, this could take a while!
fit = stan_mod$sample(
  data = br_data,#data object brms made
  chains = 4, 
  parallel_chains = getOption("mc.cores", 1L),
  iter_warmup = ifelse(test = shortrun, yes = 500, no = 2000),
  iter_sampling = ifelse(test = shortrun, yes = 1000, no = 4000),
  refresh = 500
)

  # #file to save model in
  # mod_path = file.path(save_path,
  #                      paste0('model.sim_',
  #                             cond_nm,
  #                             '.Rdata')
  # )
  # #save model
  # save(x = fit,
  #      file = mod_path)
  # #load it
  # load(mod_path)

# . Preliminary inspection of fit -----------------------------------------
  #' parms_of_interest <- c(#'b_vec',
  #'   #'b', 
  #'   'mu_vec', 
  #'   'b_Intercept', 
  #'   'sd_1', 
  #'   'b_kappa_Intercept', 
  #'   #'b_kappa',
  #'   'sd_2')
  #' fit$print(parms_of_interest)
  #' dev.new()
  #' bayesplot::mcmc_trace(fit$draws(parms_of_interest))

  # suppressWarnings({mean.circular(fake_dta$angle)})


# Extract and plot fitted parameters --------------------------------------
# . Stan fit --------------------------------------------------------------

fitted_stan = as.list(fit$draws(format = 'df'))	#extract draws
q_alpha = 0.05
q_probs = c(0,1)+c(1,-1)*q_alpha/2
par_out = with(fitted_stan, 
               {
                 #extract and summarise mu
                 #estimates of mu
                 mu_est = as.circular(x = b_Intercept, 
                                      control.circular = pipi0)
                 #estimates of mu sd
                 mu_sd = as.circular(x = `sd_1[1]`, 
                                     control.circular = pipi0)
                 #estimates of kappa
                 kappa_est = b_kappa_Intercept
                 #estimates of mu sd
                 kappa_sd = `sd_2[1]`
                 #estimates of individual mu
                 r_1_1 = mget(paste0('r_1_1[',1:par_in$no_indiv,']'))
                 mu_ind_est = lapply(X = as.data.frame(r_1_1), 
                                     FUN = as.circular,
                                     control.circular = pipi0)
                 #estimates of individual kappa
                 r_2_kappa_1 = mget(paste0('r_2_kappa_1[',1:par_in$no_indiv,']'))
                 kappa_ind_est = as.data.frame(r_2_kappa_1)
                 mu = list(#list for results
                   #robust estimate of population mean angle
                   med = median.circular(x = mu_est
                   ),
                   #quantiles of population mean angle
                   quant = quantile.circular(x = mu_est,
                                             probs = q_probs),
                   #random effects on mean angle
                   sdev = list(#median
                     med = median.circular(x = mu_sd),
                     #quantiles
                     quant = quantile.circular(x = mu_sd,
                                               probs = q_probs)
                   )
                 )
                 #extract and summarise kappa
                 kappa = list(#list for results
                   #robust estimate of population concentration
                   med = median(x = kappa_est),
                   #quantiles of population concentration
                   quant = quantile(x = kappa_est,
                                    probs = q_probs),
                   #random effects on concentration
                   sdev = list(#median
                     med = median(x = kappa_sd),
                     #quantiles
                     quant = quantile(x = kappa_sd,
                                      probs = q_probs)
                   )
                 )
                 #extract and summarise individual effects
                 mu_ind = list(
                   med = lapply(X = mu_ind_est,
                                FUN = median.circular),
                   quant = lapply(X = mu_ind_est, 
                                  FUN = quantile.circular, 
                                  probs = q_probs)
                 )
                 kappa_ind = list(
                   med = lapply(X = kappa_ind_est,
                                FUN = median),
                   quant = lapply(X = kappa_ind_est, 
                                  FUN = quantile, 
                                  probs = q_probs)
                 )
                 
                 return(
                   list(
                     mu = mu,
                     kappa = kappa,
                     mu_ind = mu_ind,
                     kappa_ind = kappa_ind
                   )
                 )
               }
)

# Extract values and compare ----------------------------------------------
StanVMaComp = function(i, 
                       data,
                       params_out,
                       params_in)
{
  mle = with(subset(data, indiv == i), 
       {
       tmp_ang = as.circular(
         x = angle,
         control.circular = pipi0
       )
       return(
         list(
           mu = deg(mean.circular(tmp_ang)),
           kappa = A1inv(rho.circular(tmp_ang))
         )
       )
       }
  )
  orig = with(params_in[[i]],
       {
         return(
           list(
               mu = mu,
               kappa = kappa
               )
             )
       }
  )
  modstan = with(params_out,
       {
        return(
          list(
            mu = as.circular(deg(mu$med[[1]] + mu_ind$med[i][[1]]),
                                  control.circular = deg360),
            kappa = exp(
                          kappa$med[[1]] + kappa_ind$med[i][[1]]
                        )
              )
            )
       }
       )
  return(
    list(orig = orig,
         mle = mle,
         stan = modstan
          # mle_mu = abs(mle$mle_mu - orig$true_mu),
          # mle_kappa = abs(mle$mle_kappa - orig$true_kappa),
          # stan_mu = abs(modstan$stan_mu - orig$true_mu),
          # stan_kappa = abs(modstan$stan_kappa - orig$true_kappa)
          )
  )
}
preds = lapply(X = 1:par_in$no_indiv,
                 FUN = StanVMaComp,
                 data = fake_dta,
                 params_out = par_out,
                 params_in = par_lst
          )
errors = sapply(X = 1:par_in$no_indiv, #print)
                FUN = function(i, prd){
                  with(prd[[i]],
                       {
                         mm = mle$mu - orig$mu
                         sm = stan$mu - orig$mu
                         mm = deg(atan2(y = sin(rad(mm)), x = cos(rad(mm))))
                         sm = deg(atan2(y = sin(rad(sm)), x = cos(rad(sm))))
                         list(mle_mu = as.numeric(abs(mm)),
                           stan_mu = as.numeric(abs(sm)),
                           mle_kappa = abs(mle$kappa - orig$kappa),
                           stan_kappa = abs(stan$kappa - orig$kappa)
                         )
                       }
                       )
                },
                prd = preds
                )
sm_errors = sapply(X = data.frame(t(errors)),
                  FUN = function(x,...){quantile(unlist(x),...)},
                  probs = c(0.025,0.25,0.5,0.75,0.975)
                )
for(nm in colnames(sm_errors))
{
  for(qn in rownames(sm_errors))
  {
par_inout_df[i,paste0(nm,'_',qn)] = sm_errors[qn, nm]
}
}

# Close loop --------------------------------------------------------------
rm(fit, fake_dta, par_in, par_out, fake_dta)
}


# Save data ---------------------------------------------------------------
#file to save this data in
sv_path <- file.path(save_path,
                     paste0('validation_',
                            Sys.Date(),
                            'UNIF',
                            '.csv')
                            # '.Rdata')
)	

write.csv(x = par_inout_df,
          file = sv_path
          )
# q('no')
# . Plot Mu & Kappa -------------------------------------------------------
  # StanVMaPlot = function(i, 
  #                        data,
  #                        params_out,
  #                        params_in)
  # {
  #   with(subset(data, indiv == i), 
  #        {plot.circular(x = as.circular(bearing, control.circular = deg360),
  #                       stack = TRUE,
  #                       sep = -0.04,
  #                       bins = (360-5)/5,
  #                       zero = pi/2)
  #          tmp_ang = as.circular(
  #            x = angle,
  #            control.circular = pipi0
  #          )
  #          arrows.circular(x = -pi/2+
  #                            mean.circular(tmp_ang),
  #                          shrink = rho.circular(tmp_ang),
  #                          col = 'green3',
  #                          lwd = 2,
  #                          length = 0)
  #        }
  #   )
  #   with(params_in[[i]],
  #        {
  #          arrows.circular(x = 
  #                            as.circular(
  #                              x = mu,
  #                              control.circular = deg360
  #                            ),
  #                          shrink = A1(kappa),
  #                          col = 'blue',
  #                          lwd = 1,
  #                          length = 0,
  #                          zero = pi/2)
  #          lines(x = A1(kappa)*cos(seq(from = -pi, to  = pi, length.out = 1e3)),
  #                y = A1(kappa)*sin(seq(from = -pi, to  = pi, length.out = 1e3)),
  #                lty = 2,
  #                col = 'blue')
  #        }
  #   )
  #   with(params_out,
  #        {
  #          lines.circular(x = -pi/2+
  #                           as.circular(
  #                             x = rep(mu$med[[1]] + mu_ind$med[i][[1]],2),
  #                             control.circular = pipi0
  #                           ),
  #                         y = -(1-A1(exp(
  #                           kappa$med[[1]] + kappa_ind$quant[i][[1]]
  #                         ))),
  #                         col = adjustcolor('orange',alpha.f = 70/255),
  #                         lwd = 10,
  #                         lend = 'butt'
  #          )
  #          arrows.circular(x = -pi/2+
  #                            as.circular(x = mu$med[[1]] + mu_ind$med[i][[1]],
  #                                        control.circular = pipi0
  #                            ),
  #                          shrink = A1(exp(kappa$med + kappa_ind$med[i][[1]])),
  #                          col = adjustcolor('red',alpha.f = 150/255),
  #                          lwd = 5,
  #                          length = 0
  #          )
  #        })
  #   
  # }
  # par(mfrow = rep(x = ceiling(sqrt(par_in$no_indiv)), times = 2),
  #     mar = c(0,0,0,0)
  # )
  # invisible(
  #   {
  #     lapply(X = 1:par_in$no_indiv,
  #            FUN = StanVMaPlot,
  #            data = fake_dta,
  #            params_out = par_out,
  #            params_in = par_lst
  #     )
  #   }
  # )
  # plot.new()
  # legend(x = 'center',#topright',
  #        legend = c('input',
  #                   'Maximum likelihood',
  #                   'Stan fit',
  #                   'Credible interval'),
  #        lty = 1,
  #        col = c('blue',
  #                'green3',
  #                adjustcolor('red',alpha.f = 200/255),
  #                adjustcolor('orange',alpha.f = 100/255)),
  #        lwd = c(1,2,3,5,10),
  #        cex = 0.75,
  #        bg = 'white',
  # )
