#FOR A 'CLEAN' RUN, RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 05 08
#     MODIFIED:	James Foster              DATE: 2024 05 08
#
#  DESCRIPTION: Use an optimiser to search for angles on a toroid constrained
#               close to a unit vector, without forcing a unit vector relationship.
#               
#       INPUTS: 
#               
#      OUTPUTS: 
#	   CHANGES: - 
#
#   REFERENCES: Raoul-Kima
#               Divergence / Treedepth issues with unit vector
#               Stan Forums
#               https://discourse.mc-stan.org/t/divergence-treedepth-issues-with-unit-vector/8059
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- 

# Useful functions --------------------------------------------------------

## Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircStats)#package for circular hypothesis tests
  }
)


## General functions -----------------------------------------------------

Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
}


## Circular statistics ---------------------------------------------------


## Critical value calculation --------------------------------------------

#critical mean vector length for a von Mises distribution
CritVMr = function(n, test = 'Rayleigh')
{
  switch(EXPR = test,
         Rayleigh = sqrt(-log(0.05)/n),
         V = 1.644854/(sqrt(2 * n)),
         sqrt(-log(0.05)/n)
  )
}




# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
# paired_data = TRUE # Are the data in the two columns paired (each from the same animal or group)?
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
angle_rot = "counter" # "clock" or "counter"
paired_data = FALSE
ref_angle = 0

#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#User profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

# # Simulate data ------------------------------------------------
n_angles = 44
# minimum discriminable angle appears to be approx 35°
mu_offset = rad(35)
kappa_both = A1inv(0.7) #concentration around each trial mean
logkappa_var = 1.0 #scale of random variation in concentration (log units)
if(paired_data)
{
kappa_indiv = A1inv(0.98) #concentration across individuals (pairs)
#mean angle in trail 1 for each individual (pair)
mu1_sim = rvonmises(n = n_angles,
                      mu = rcircularuniform(1),#random angle
                      kappa = kappa_indiv#the wider the distribution of individual biases, the greater the influence of pairing
                      )
#simulate the full dataset
sim = data.frame(
                 angle_1 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                   sapply(X = mu1_sim,
                          FUN = rvonmises,
                          n = 1,
                           kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
                   )
                 ))*180/pi),#convert to angles and round to nearest degree
                 angle_2 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                             sapply(X =mu1_sim + mu_offset,# true difference,
                                    FUN = rvonmises,
                                    n = 1,
                                    kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
                             )
                 ))*180/pi) #convert to angles and round to nearest degree
                 )
}else
{
n_angles2 = ceiling(0.75*n_angles)
mu1_sim = rcircularuniform(n = 1,control.circular = list(units = angle_unit))
sim = data.frame(
                 angle_1 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                   rvonmises(n = n_angles,
                             mu = circular(x = mu1_sim, units = angle_unit),
                           kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
                   )
                 ))),#convert to angles and round to nearest degree
                 angle_2 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                   rvonmises(n = n_angles2,
                             mu = circular(x = mu1_sim+deg(mu_offset), units = angle_unit),
                             kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
                   )
                 ),
                 circular(x = rep(x = NA, times = n_angles - n_angles2),
                         units = angle_unit) #convert to angles and round to nearest degree
                 ) )
                )
}
#save somewhere the user likely keeps data
write.table(x = sim,
            file = file.path(ltp,'Documents', "simulated_angles.csv"),
            sep = csv_sep,
            row.names = FALSE
            )

# Organise data -----------------------------------------------------------
adata = sim
dt_dim = dim(adata)
# Plot simulated data -----------------------------------------------------

par(mar =rep(0,4))
plot.circular(x = circular(x = adata$angle_1, 
                           type = 'angles',
                           unit = angle_unit,
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
stack = TRUE,
bins = 360/5,
sep = 0.5/dt_dim[1],
col = 'cyan4'
)
par(new = T)
plot.circular(x = circular(x = adata$angle_2, 
                           type = 'angles',
                           unit = 'degrees',
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
stack = TRUE,
bins = 360/5,
sep = -0.5/dt_dim[1],
col = 'darkblue',
shrink = 1.05,
axes = F
)
arrows.circular(x = circular(x = ref_angle, 
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
),
y = 1,
lwd = 5, 
col = rgb(0,0,0,0.1),
length = 0
)
with(adata,
arrows.circular(x = mean.circular(x = 
                    circular(x = angle_1,
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot),
                    na.rm = TRUE
                              ),
                y = rho.circular(x = 
                    circular(x = angle_1,
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot),
                    na.rm = TRUE),
                  col = 'cyan4',
                lwd = 3)
)

with(adata,
arrows.circular(x = mean.circular(x = 
                    circular(x = angle_2,
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot),
                    na.rm = TRUE),
                y = rho.circular(x = 
                    circular(x = angle_2,
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot),
                    na.rm = TRUE),
                  col = 'darkblue',
                lwd = 3)
)

#Add the critical values
lines.circular(x = circular(x = 
                              seq(from = -180,
                                  to = 180,
                                  length.out = 1e3),
                            units = angle_unit, 
                            rotation = angle_rot),  
               y = rep(x = CritVMr(n = sum(!is.na(adata$angle_1)) ),
                       times = 1e3)-1,
               col = 'cyan4', 
               lty = 2,
               lwd = 1,
               lend = 'butt')
lines.circular(x = circular(x = 
                              seq(from = -180,
                                  to = 180,
                                  length.out = 1e3),
                            units = angle_unit, 
                            rotation = angle_rot),  
               y = rep(x = CritVMr(n = sum(!is.na(adata$angle_2)) ),
                       times = 1e3)-1,
               col = 'darkblue', 
               lty = 2,
               lwd = 1,
               lend = 'butt')


# Set up optimiser --------------------------------------------------------

DoughnutOptim = function(params,
                         angles,
                         prior_uv = 1.0, #gentle encouragement is enough
                         prior_logkappa = c(1.0,3.0), #gentle encouragement is enough
                         au = 'degrees',
                         ar = 'clock')
{
  mu_cos = params[1]
  mu_sin = params[2]
  mu = atan2(y = mu_cos,
             x = mu_sin)
  if(au %in% 'degrees')
  {mu = deg(mu)}
  kappa = exp(params[3])
  vlength = sqrt(mu_cos^2 + mu_sin^2)
  #von Mises likelihood
  lpd = sum( # add together
    dvonmises(x = circular(x = angles[!is.na(angles)],
                           units = au,
                           rotation = ar), # probability density for each observed angle
              mu = circular(x = mu,
                            units = au,
                            rotation = ar),#,rotation = ar),# ML estimated mean
              kappa = kappa, # ML estimated concentration
              log = TRUE) # on a log scale (i.e. add instead of multiplying)
  )
  #unit vector prior
  lpd = lpd + dnorm(x = vlength,
                    mean = 1.0,
                    sd = prior_uv,
                    log = TRUE)  
  #nonzero kappa prior
  lpd = lpd + dnorm(x = params[3],
                    mean = prior_logkappa[1],
                    sd = prior_logkappa[2],
                    log = TRUE)
  return(-lpd)
}

DO_getpar = function(params = c(mu_cos = 0.0,
                                mu_sin = 1.0,
                                log_kappa = 1.0),
                     angles,
                     prior_uv = 1.0, #gentle encouragement is enough
                     au = 'degrees',
                     ar = 'clock',
                     ...)#passed to optim
{
  opar = optim(par = dn_start,
                   fn = DoughnutOptim,
                   angles = angles,
                   prior_uv = 1.0,
                   ...)
  dt_opar = data.frame(t(opar$par))
  dt_opar = within(dt_opar, 
         {
           vlength = sqrt(mu_cos^2 + mu_sin^2)
           mu = atan2(y = mu_cos,
                      x = mu_sin)
           mu_deg = deg(mu)
           kappa = exp(log_kappa)
           rho = A1(kappa)
         }
  )
  return(dt_opar)
}

#compare with mle.vonmises
mlvm = with(adata,
            mle.vonmises(x = circular(x = angle_1,
                                      units = angle_unit,
                                      rotation = angle_rot),
                                      bias = TRUE)
            )
dno = with(adata,
            DO_getpar(angles = angle_1,
                      method = 'BFGS',
                      control = list(maxit = 1e3)
                      )
            )
print(
cbind(MLE = with(mlvm, c(as.numeric(mu), kappa)),
      Optim = with(dno, c(mu_deg, kappa))),
digits = 4)


ldata = with(adata,
             data.frame(
             angle = c(angle_1,
                       angle_2),
             condition = c(rep(1, times = length(angle_1)),
                           rep(2, times = length(angle_2))
                           )
             )
)

two_mod = sapply(X = with(ldata, unique(condition)),
         FUN = function(i)
         {
           with(subset(ldata, condition %in% i) ,
                {
                  DO_getpar(angles = angle)
                }
                )
         }
)

with(data.frame(t(two_mod)),
     print(cbind(mu_deg, rho), digits = 4))
