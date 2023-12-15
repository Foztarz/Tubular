# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2021 12 07
#     MODIFIED:	James Foster              DATE: 2023 04 06
#
#  DESCRIPTION: Loads a data frame containing parameter estimate error 
#               from a batch run of "RandomIntercepts-validation.R" and plots
#                pairs of errors for Maximum Likelihood and Bayesian Estimation.
#               
#       INPUTS: An ".Rdata" file with the required functions.
#               User should specify data parameters (line 100).
#               
#      OUTPUTS: Plots (.pdf; .png).
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
#    EXAMPLES:  Fill out user input (lines 100-105), then press ctrl+shift+s to run
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Load data   
#- General comparison
#- Comparison by parameter
#- Quantile comparison
#- Save results  

# Check environment -------------------------------------------------------
#check the operating system and assign a logical flag (TRUE or FALSE)
sys_win <- Sys.info()[['sysname']] == 'Windows'
# Load packages -----------------------------------------------------------
#Package for handling circular data
library(circular)
#Package for plotting MCMC chains
library(bayesplot)
#Package for handing large datasets
library(data.table)
# Useful Functions --------------------------------------------------------
#Open file with default program on any OS
# https://stackoverflow.com/a/35044209/3745353
shell.exec.OS  <- function(x){
  # replacement for shell.exec (doesn't exist on MAC)
  if (exists("shell.exec",where = "package:base"))
  {return(base::shell.exec(x))}else
  {comm <- paste0('open "',x,'"')
  return(system(comm))}
}
PlotNamer = function(file,
                     namend, 
                     saveas = 'pdf')
{
  fp = file.path(dirname(file), 
                 paste0(basename(file),
                        namend,'.', saveas))
  if(file.exists(fp)) # if this file would overwrite an existing file
  {
    message('A plot called "', basename(fp), '" already exists in this folder.')
    nnm = readline(prompt = 'New plot name: '
                  )#ask the user for a new name
    fp = file.path(dirname(file),
                   ifelse(test = nchar(nnm),#if a new name is supplied
                                 yes = paste0(nnm,'.', save_type),#use it
                                 no = basename(fp))#otherwise overwrite
                   )
  }
  return(fp)
}
















# User input --------------------------------------------------------------
point_col = 'darkblue'#or any of these: https://htmlcolorcodes.com/color-names/
plette_lines = 'Dark3'#'Viridis'#  for options see http://colorspace.r-forge.r-project.org/articles/approximations.html
plette_heat = 'Plasma'#'YlGnBu'#  for options see http://colorspace.r-forge.r-project.org/articles/approximations.html
mgn = c(5,5,3,0)# margin size in lines of text: bottom, left, top, right
save_plt = TRUE # FALSE# should the plot be saved
save_type ="png"# "pdf"# 

# Load the data -----------------------------------------------------------
vdir = dirname(sys.frame(1)$ofile)
dt_path = tryCatch(expr = 
                      {
                        file.path(path = vdir, 
                                 max(list.files(path = vdir, 
                                                pattern = "^validation_.*.csv$")
                                     )#named by date, max should find the most recent
                                 )
                        },
                    error = function(e){paste0('validation_',
                                               Sys.Date(),
                                               '.csv')}
)
if(!file.exists(dt_path))
{
  msg <- 'Please select "validation_YYYY-MM-DD.csv"'
  # ask user to find data
  if(sys_win){#choose.files is only available on Windows
    message('\n\n',msg,'\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    dt_path  <- choose.files(
      default = file.path(gsub('\\\\', '/', Sys.getenv('USERPROFILE')),#user
                          'Documents'),#For some reason this is not possible in the "root" user
      caption = msg
    )
  }else{
    message('\n\n',msg,'\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    dt_path <- file.choose(new=F)
  }
}
#read in relevant functions
plt_dta = fread(file = dt_path)

# Derive variables --------------------------------------------------------
ind_nos = with(plt_dta, unique(no_indiv))
nper_nos = with(plt_dta, unique(nper_indiv))
kappas = with(plt_dta, unique(kappa))
sdmus = with(plt_dta, unique(sd))
sdlkappas = with(plt_dta, unique(sd_logkappa))

# Set up plot to save ---------------------------------------------------
if(save_plt)
{
  plt_path = PlotNamer(file = dt_path,
                       namend = '_overview',
                       saveas = save_type)
switch(save_type,
       pdf = 
         pdf(file = plt_path,
             paper = 'a4',
             height = 10,
             bg = 'white',
             useDingbats = F
         ),
       png = png(file = plt_path,
                 res = 150,
                 width = 210*10,
                 height = 297*10,
                 bg = 'white'
       ),
       jpeg(file = paste0(plt_path,'.jpeg'),
            quality = 100,
            width = 210*10,
            height = 297*10,
            bg = 'white'
       )
)
}

# Set up plot area ------------------------------------------------------
par(mar = mgn,
    mfrow = c(6,2))

# General comparison ------------------------------------------------------

# . Plot mean angle errors ------------------------------------------------
with(plt_dta,
     {
      plot(x = NULL,
          xlim = c(0,180),
          ylim = c(0,180),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )  
       title(main =  'median error in mean angle (°)',
             line = -1)
      axis(side = 1, at = seq(from = 0, to  = 180, by  = 30))
      axis(side = 2, at = seq(from = 0, to  = 180, by  = 30))
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = c(0,180, 180,0), 
              y = c(0,180, 0,0), 
              col = gray(250/255),
              border = NA
              )
      points(x = `mle_mu_50%`,
             y = `stan_mu_50%`,
             pch = 19,
             col=  adjustcolor(col = point_col,
                               alpha.f = 150/255)
             )
     }
     )

# . Plot mean vector length errors ---------------------------------------
with(plt_dta,
     {
      plot(x = NULL,
          xlim = A1inv(c(0,0.7)),
          ylim = A1inv(c(0,0.7)),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )  
       title(main =  'median error in concentration (kappa)',
             line = -1)
      axis(side = 1)
      axis(side = 2)
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = A1inv(c(0,0.7, 0.7,0)), 
              y = A1inv(c(0,0.7, 0,0)), 
              col = gray(250/255),
              border = NA
      )
      points(x = `mle_kappa_50%`,
             y = `stan_kappa_50%`,
             pch = 19,
             col=  adjustcolor(col = point_col,
                               alpha.f = 150/255)
             )
     }
     )

# By number of individuals-------------------------------------------------

# . Plot mean angle errors ------------------------------------------------
cls = hcl.colors(n = length(ind_nos),
                 palette = plette_lines)
with(plt_dta,
     {
      plot(x = NULL,
          xlim = c(0,180),
          ylim = c(0,180),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
      title(main =  'median error in mean angle (°)',
            line = -1)
      axis(side = 1, at = seq(from = 0, to  = 180, by  = 30))
      axis(side = 2, at = seq(from = 0, to  = 180, by  = 30))
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = c(0,180, 180,0), 
              y = c(0,180, 0,0), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(nn in ind_nos)
{
  with(subset(plt_dta, subset = no_indiv %in% nn),
       {
        points(x = `mle_mu_50%`,
               y = `stan_mu_50%`,
               pch = 19,
               col=  adjustcolor(col = cls[which(ind_nos %in% nn)],
                                 alpha.f = 150/255)
        )
       }
       )
}
legend(x = 'bottomright',
       legend = paste(ind_nos, 'indiv.'),
       col = cls,
       pch = 19
       )

# . Plot kappa errors -----------------------------------------------------
with(plt_dta,
     {
      plot(x = NULL,
          xlim = A1inv(c(0,0.7)),
          ylim = A1inv(c(0,0.7)),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
       title(main =  'median error in concentration (kappa)',
             line = -1)
      axis(side = 1)
      axis(side = 2)
      abline(a = 0, b = 1, lty = 2)
     }
     )
for(nn in ind_nos)
{
  with(subset(plt_dta, subset = no_indiv %in% nn),
       {
         points(x = `mle_kappa_50%`,
                y = `stan_kappa_50%`,
                pch = 19,
                col=  adjustcolor(col = cls[which(ind_nos %in% nn)],
                                  alpha.f = 150/255)
         )
       }
  )
}
legend(x = 'bottomright',
       legend = paste(ind_nos, 'indiv.'),
       col = cls,
       pch = 19
)

# By number of observations per individual---------------------------------

# . Plot mean angle errors ------------------------------------------------
cls = hcl.colors(n = length(nper_nos),
                 palette = plette_lines)
with(plt_dta,
     {
      plot(x = NULL,
          xlim = c(0,180),
          ylim = c(0,180),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
      title(main =  'median error in mean angle (°)',
             line = -1) 
      axis(side = 1, at = seq(from = 0, to  = 180, by  = 30))
      axis(side = 2, at = seq(from = 0, to  = 180, by  = 30))
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = c(0,180, 180,0), 
              y = c(0,180, 0,0), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(np in nper_nos)
{
  with(subset(plt_dta, subset = nper_indiv %in% np),
       {
        points(x = `mle_mu_50%`,
               y = `stan_mu_50%`,
               pch = 19,
               col=  adjustcolor(col = cls[which(nper_nos %in% np)],
                                 alpha.f = 150/255)
        )
       }
       )
}
legend(x = 'bottomright',
       legend = paste('n=',nper_nos, 'per indiv.'),
       col = cls,
       pch = 19
       )

# . Plot kappa errors -----------------------------------------------------
with(plt_dta,
     {
      plot(x = NULL,
          xlim = A1inv(c(0,0.7)),
          ylim = A1inv(c(0,0.7)),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
      title(main =  'median error in concentration (kappa)',
             line = -1) 
      axis(side = 1)
      axis(side = 2)
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = A1inv(c(0,0.7, 0.7,0)), 
              y = A1inv(c(0,0.7, 0,0)), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(np in nper_nos)
{
  with(subset(plt_dta, subset = nper_indiv %in% np),
       {
         points(x = `mle_kappa_50%`,
                y = `stan_kappa_50%`,
                pch = 19,
                col=  adjustcolor(col = cls[which(nper_nos %in% np)],
                                  alpha.f = 150/255)
         )
       }
  )
}
legend(x = 'bottomright',
       legend = paste('n=',nper_nos, 'per indiv.'),
       col = cls,
       pch = 19
)

# By kappa ----------------------------------------------------------------

# . Plot mean angle errors ------------------------------------------------
cls = hcl.colors(n = length(kappas),
                 palette = plette_lines)
with(plt_dta,
     {
      plot(x = NULL,
          xlim = c(0,180),
          ylim = c(0,180),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
       title(main =  'median error in mean angle (°)',
             line = -1) 
      axis(side = 1, at = seq(from = 0, to  = 180, by  = 30))
      axis(side = 2, at = seq(from = 0, to  = 180, by  = 30))
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = c(0,180, 180,0), 
              y = c(0,180, 0,0), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(kp in kappas)
{
  with(subset(plt_dta, subset = kappa %in% kp),
       {
        points(x = `mle_mu_50%`,
               y = `stan_mu_50%`,
               pch = 19,
               col=  adjustcolor(col = cls[which(kappas %in% kp)],
                                 alpha.f = 150/255)
        )
       }
       )
}
legend(x = 'bottomright',
       legend = paste('mean vector length=',round(A1(kappas),2)),
       col = cls,
       pch = 19
       )

# . Plot kappa errors -----------------------------------------------------
with(plt_dta,
     {
      plot(x = NULL,
          xlim = A1inv(c(0,0.7)),
          ylim = A1inv(c(0,0.7)),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
       title(main =  'median error in concentration (kappa)',
             line = -1) 
      axis(side = 1)
      axis(side = 2)
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = A1inv(c(0,0.7, 0.7,0)), 
              y = A1inv(c(0,0.7, 0,0)), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(kp in kappas)
{
  with(subset(plt_dta, subset = kappa %in% kp),
       {
         points(x = `mle_kappa_50%`,
                y = `stan_kappa_50%`,
                pch = 19,
                col=  adjustcolor(col = cls[which(kappas %in% kp)],
                                  alpha.f = 150/255)
         )
       }
  )
}
legend(x = 'bottomright',
       legend = paste('mean vector length=',round(A1(kappas),2)),
       col = cls,
       pch = 19
)
# By sd in mean angle -----------------------------------------------------

# . Plot mean angle errors ------------------------------------------------
cls = hcl.colors(n = length(sdmus),
                 palette = plette_lines)
with(plt_dta,
     {
      plot(x = NULL,
          xlim = c(0,180),
          ylim = c(0,180),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
       title(main =  'median error in mean angle (°)',
             line = -1) 
      axis(side = 1, at = seq(from = 0, to  = 180, by  = 30))
      axis(side = 2, at = seq(from = 0, to  = 180, by  = 30))
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = c(0,180, 180,0), 
              y = c(0,180, 0,0), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(sm in sdmus)
{
  with(subset(plt_dta, subset = sd %in% sm),
       {
        points(x = `mle_mu_50%`,
               y = `stan_mu_50%`,
               pch = 19,
               col=  adjustcolor(col = cls[which(sdmus %in% sm)],
                                 alpha.f = 150/255)
        )
       }
       )
}
legend(x = 'bottomright',
       legend = sapply(X = sdmus,
                       FUN = function(x)
                       {
                       if(is.na(x)){'uniform distribution'}else
                       {paste0('sd mean angles = ', x, '°')}
                       }
                     ),
       col = cls,
       pch = 19
       )

# . Plot kappa errors -----------------------------------------------------
with(plt_dta,
     {
      plot(x = NULL,
          xlim = A1inv(c(0,0.7)),
          ylim = A1inv(c(0,0.7)),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
       title(main =  'median error in concentration (kappa)',
             line = -1) 
      axis(side = 1)
      axis(side = 2)
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = A1inv(c(0,0.7, 0.7,0)), 
              y = A1inv(c(0,0.7, 0,0)), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(sm in sdmus)
{
  with(subset(plt_dta, subset = sd %in% sm),
       {
         points(x = `mle_kappa_50%`,
                y = `stan_kappa_50%`,
                pch = 19,
                col=  adjustcolor(col = cls[which(sdmus %in% sm)],
                                  alpha.f = 150/255)
         )
       }
  )
}
legend(x = 'bottomright',
       legend = sapply(X = sdmus,
                       FUN = function(x)
                       {
                         if(is.na(x)){'uniform distribution'}else
                         {paste0('sd mean angles = ', x, '°')}
                       }
       ),
       col = cls,
       pch = 19
)

# By sd in log kappa ------------------------------------------------------

# . Plot mean angle errors ------------------------------------------------
cls = hcl.colors(n = length(sdlkappas),
                 palette = plette_lines)
with(plt_dta,
     {
      plot(x = NULL,
          xlim = c(0,180),
          ylim = c(0,180),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
       title(main =  'median error in mean angle (°)',
             line = -1) 
      axis(side = 1, at = seq(from = 0, to  = 180, by  = 30))
      axis(side = 2, at = seq(from = 0, to  = 180, by  = 30))
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = c(0,180, 180,0), 
              y = c(0,180, 0,0), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(sk in sdlkappas)
{
  with(subset(plt_dta, subset = sd_logkappa %in% sk),
       {
        points(x = `mle_mu_50%`,
               y = `stan_mu_50%`,
               pch = 19,
               col=  adjustcolor(col = cls[which(sdlkappas %in% sk)],
                                 alpha.f = 150/255)
        )
       }
       )
}
legend(x = 'bottomright',
       legend = paste0('sd log(kappa) = ', round(sdlkappas,2)),
       col = cls,
       pch = 19
       )

# . Plot kappa errors -----------------------------------------------------
with(plt_dta,
     {
      plot(x = NULL,
          xlim = A1inv(c(0,0.7)),
          ylim = A1inv(c(0,0.7)),
          xlab = 'Maximum Likelihood',
          ylab = 'Bayesian Estimation',
          axes = F
          )
       title(main =  'median error in concentration (kappa)',
             line = -1) 
      axis(side = 1)
      axis(side = 2)
      abline(a = 0, b = 1, lty = 2)
      abline(h = 0, v = 0, lty = 1)
      polygon(x = A1inv(c(0,0.7, 0.7,0)), 
              y = A1inv(c(0,0.7, 0,0)), 
              col = gray(250/255),
              border = NA
      )
     }
     )
for(sk in sdlkappas)
{
  with(subset(plt_dta, subset = sd_logkappa %in% sk),
       {
         points(x = `mle_kappa_50%`,
                y = `stan_kappa_50%`,
                pch = 19,
                col=  adjustcolor(col = cls[which(sdlkappas %in% sk)],
                                  alpha.f = 150/255)
         )
       }
  )
}
legend(x = 'bottomright',
       legend = paste0('sd log(kappa) = ', round(sdlkappas,2)),
       col = cls,
       pch = 19
)

# Save and open plot ----------------------------------------------------
if(save_plt)
{
  dev.off()
  shell.exec.OS(plt_path)
}


# Set up boxplot to save ------------------------------------------------
if(save_plt)
{
  plt_path = PlotNamer(file = dt_path,
                       namend = '_boxplot',
                       saveas = save_type)
  switch(save_type,
         pdf = 
           pdf(file = plt_path,
               paper = 'a4',
               height = 10,
               bg = 'white',
               useDingbats = F
           ),
         png = png(file = plt_path,
                   res = 150,
                   width = 210*10,
                   height = 297*10,
                   bg = 'white'
         ),
         jpeg(file = paste0(plt_path,'.jpeg'),
              quality = 100,
              width = 210*10,
              height = 297*10,
              bg = 'white'
         )
  )
}

# Set up plot area ------------------------------------------------------
par(mar = mgn,
    mfrow = c(6,2),
    oma = if(save_plt){c(0,5,0,5)}else{c(0,0,0,0)} 
)

# General comparison ------------------------------------------------------

BXplotmeth = function(data,
           palette,
           varnm,
           lbl,
           pltmar = 0.3,
           alpha = 50/255
           )
{
  SdmuRelabel = function(x)
  {
    if(is.na(x)){'uniform'}else
    {paste0(x,'°')}
  }
  unq = sort( unique(with(data, { eval(str2lang(varnm)) })) ) #N.B. boxplot sorts
  lunq = length(unq)
  cll = hcl.colors(n = lunq,
                   palette = palette)
  expmu = str2lang(paste0('(`mle_mu_50%`-`stan_mu_50%`)~',varnm))
  expkappa = str2lang(paste0('(`mle_kappa_50%`-`stan_kappa_50%`)~',varnm))
  
  boxplot(x = NULL,
          data = data,
          ylim = c(-180,180)/1.5,
          xlim = c(1,lunq)+c(-1,1)*pltmar,
          xlab = lbl,
          ylab = 'Bayesian Estimation improvement (°)',
          axes = F
  )  
  polygon(x = c(1,lunq, lunq,1) + c(-1,1,1,-1)*pltmar, 
          y = c(0,0, 180,180), 
          col = gray(250/255),
          border = NA
  )
  title(main =  'median error in mean angle (°)',
        line = -1)
  axis(side = 1, at = 1:lunq, labels = unq, col = NA)
  axis(side = 2, at = seq(from = -180, to  = 180, by  = 30))
  boxplot(formula = eval(expmu),
          data = data,
          col = cll,
          add = T ,
          xlab = '',
          ylab = '',
          axes = F,
          outline = F,
          pars = list(boxwex = diff(range(par('xaxp')))*0.06,
                      staplewex = diff(range(par('xaxp')))*0.10, 
                      outwex = diff(range(par('xaxp')))*0.05)
  )
  stripchart(x = eval(expmu),
             data = data,
             bg = gray(200/255, alpha= alpha),
             col = adjustcolor(col = cll, alpha.f = alpha),
             pch = 21,
             vertical = T,
             jitter = 0.3,
             method = 'jitter',
             add = T ,
             xlab = '',
             ylab = '',
             axes = F
  )
  abline(h = 0, lty = 2)
  
  boxplot(x = NULL,
          data = data,
          ylim = c(-1,1)*A1inv(0.85),
          xlim = c(1,lunq) + c(-1,1)*pltmar,
          xlab = lbl,
          ylab = 'Bayesian Estimation improvement (kappa)',
          axes = F
  )  
  polygon(x = c(1,lunq, lunq,1) + c(-1,1,1,-1)*pltmar, 
          y = A1inv(c(0,0,0.85,0.85)), 
          col = gray(250/255),
          border = NA
  )
  title(main =  'median error in concentration (kappa)',
        line = -1)
  axis(side = 1, at = 1:lunq, labels = unq, col = NA)
  axis(side = 2, at = seq(from = round(-A1inv(0.85)), to  = round(A1inv(0.85)), by  = 0.5))
  boxplot(formula = eval(expkappa),
          data = data,
          col = cll,
          add = T ,
          xlab = '',
          ylab = '',
          axes = F,
          outline = F,
          pars = list(boxwex = diff(range(par('xaxp')))*0.06,
                      staplewex = diff(range(par('xaxp')))*0.10, 
                      outwex = diff(range(par('xaxp')))*0.05)
  )
  stripchart(x = eval(expkappa),
             data = data,
             bg = gray(200/255, alpha= alpha),
             col = adjustcolor(col = cll, alpha.f = alpha),
             pch = 21,
             vertical = T,
             jitter = 0.3,
             method = 'jitter',
             add = T ,
             xlab = '',
             ylab = '',
             axes = F
  )
  abline(h = 0, lty = 2)
}

BXplotmeth(data = plt_dta,
           palette = 'Blues',
           varnm = '(no_indiv>0)',
           lbl = 'All')
BXplotmeth(data = plt_dta,
           palette = plette_lines,
           varnm = 'no_indiv',
           lbl = 'Number of Individuals')
BXplotmeth(data = plt_dta,
           palette = plette_lines,
           varnm = 'nper_indiv',
           lbl = 'Observations per individual')
BXplotmeth(data = plt_dta,
           palette = plette_lines,
           varnm = 'round(A1(kappa),2)',
           lbl = 'Mean vector length')
BXplotmeth(data = plt_dta,
           palette = plette_lines,
           varnm = "sapply(X = sd, FUN = SdmuRelabel)",
           lbl = 'sd in mean angle (°)')
BXplotmeth(data = plt_dta,
           palette = plette_lines,
           varnm = 'round(sd_logkappa,2)',
           lbl = 'sd in log(kappa)')
# Save and open plot ----------------------------------------------------
if(save_plt)
{
  dev.off()
  shell.exec.OS(plt_path)
}
