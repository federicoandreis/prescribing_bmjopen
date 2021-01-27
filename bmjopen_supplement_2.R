#######################################################################
# Code to reproduce the analyses and the plots in the the manuscript  #
# "Prescribing patterns for medical treatment of suspected prostatic  #
#  obstruction: a longitudinal register-based study the Scottish      #
#  Health and Social Care Open Data" (2021)                           #
# Authors: Federico Andreis, Richard Bryant, Emanuele Giorgi,         #
#          Andrea Williamson, Ashleigh Ward                           #
#######################################################################

#############################
# Code by: Federico Andreis #
#############################

#############
# libraries #
#############
library(tidyverse)    # to manipulate data objects

library(rgdal)        # to deal with spatial objects
library(sf)           
library(raster) 

library(caramellar)   # to build the adjacency matrix, not on CRAN
                      # devtools::install_github("barryrowlingson/caramellar")

library(INLA)         # to fit the model
library(brinla)       # to manipulate inla objects

######################################
# load the dataset and the shapefile #
######################################
dataset <- read_csv('dataset.csv')

map_scotland_hb <- readOGR('SG_NHS_HealthBoards_2018.shp')

###########################
# create adjacency matrix #
###########################

# create adjacency matrix and transform it into an inla graph
adjacency_graph <- voronoi_adjacency(dataset %>%
                                       group_by(postcode_n) %>%
                                       slice(1) %>%
                                       dplyr::select(postcode_n,Easting,Northing) %>%
                                       ungroup,
                                     postcode_n~Easting+Northing)$Adjacencies %>% inla.read.graph

#################
# fit the model #
#################

# set inla controls
control <-list(predictor=list(compute=TRUE,link=1), 
               results=list(return.marginals.random=TRUE, 
                            return.marginals.predictor=TRUE), 
               compute=list(hyperpar=TRUE, return.marginals=TRUE, 
                            dic=TRUE, mlik=TRUE, cpo=TRUE, po=TRUE, 
                            waic=TRUE, graph=TRUE, gdensity=TRUE))

# model specification
model <- n_daily_doses~drug_group+gp_run+dispensing+     # fixed effects
                 log(list_size)+prop_m_45p+p_15+ur_code+ # fixed effects
  f(postcode_n,model='besag',graph=adjacency_graph)+     # structured spatial component (iCAR)
  f(time_n,model='ar1',group=hb_drug)+                   # temporal component (grouped AR1)
  f(gp_drug,model='iid')                                 # unstructured spatial component (iid)

# fit the model using a cheap Gaussian approximation to obtain reasonable starting values
# run time ~18 minutes on a AMD Ryzen 7 2700x processor with 32GB DDR4 RAM, under Windows 10
cheap_approximation <- inla(model,
                            family='poisson',
                            data=dataset,
                            control.inla=list(diagonal=100,
                                              strategy="gaussian",
                                              int.strategy="eb"),
                            control.compute=control$compute,
                            control.predictor=control$predictor,
                            control.results=control$results,
                            verbose=TRUE)

# use the command
#
# model_fit <- cheap_approximation
#
# and ignore the following model fit block to be able to use the rest of the 
# code without having to wait for the more accurate approximation to be obtained
# Note: if you use this approach, the resulting estimates and plots will differ
# from those in the paper


# fit the model using the cheap approximation estimates as starting values
# note: this step takes ~7 hours on the machine described earlier
model_fit <- inla(model,
                  family='poisson',
                  data=dataset,
                  control.inla=list(diagonal=10),
                  control.fixed = list(prec.intercept = 0.1),
                  control.compute=control$compute,
                  control.predictor=control$predictor,
                  control.results=control$results,
                  control.mode=list(result=cheap_approximation,
                                    restart=TRUE),
                  verbose=TRUE)

rm(adjacency_graph,control,
   model,cheap_approximation) # clean up

#######################################
# post-processing of the model output #
#######################################

# add the fitted values to the dataset and make health board names into factors
dataset <- dataset %>% 
  mutate(fitted=model_fit$summary.fitted.values$mean,
         hb_name=factor(hb_name))

# obtain predictions by health board and drug group
average_fitted <- dataset %>%
  group_by(time=time,hb_name,drug_group) %>% 
  summarise(avg_fitted=mean(fitted),
            avg_observed=mean(n_daily_doses))

# extract the exponentiated spatial residuals and add them to the dataset
exp_residuals <- numeric(length(model_fit$marginals.random$postcode_n))

for (i in 1:length(exp_residuals)) {
  
  tmp <- model_fit$marginals.random$postcode_n[[i]]
  
  exp_residuals[i] <- inla.emarginal(exp,tmp) # exponentiate the spatial residuals
  
}

rm(tmp,i) # clean up

# make the residuals  into a tibble that also contains the practice postcodes 
# and health boards names
exp_residuals <- data.frame(postcode_n=1:length(exp_residuals),
                            postcode=unique(dataset$postcode),
                            exp_residuals=exp_residuals) %>%
  mutate(postcode=as.character(postcode)) %>%
  left_join(.,dataset %>%
              dplyr::select(postcode,hb_name),
            by=c('postcode','postcode')) %>% 
  group_by(postcode) %>% 
  slice(1)

# add the exponentiated residuals to the dataset
dataset <- left_join(dataset,
                     exp_residuals %>% 
                       dplyr::select(postcode_n,exp_residuals),
                     by='postcode_n')

# set up the shapefile and data needed to to plot the map in Figure 3
geo_df <- dataset
coordinates(geo_df) <- ~Easting+Northing
crop_map_hb <- crop(map_scotland_hb,geo_df)

rm(map_scotland_hb) # clean up

# obtain the quantiles of the exponentiate residuals distribution
# residuals are constant for each GP practice, extract only one line each
dataset_singletons <- dataset %>% 
  group_by(gp_code) %>% 
  slice(1) %>% 
  ungroup

exp_residuals_quantiles <- quantile(dataset_singletons$exp_residuals,
                                    p=seq(0,1,by=.025))

##################################
# produce the plots in the paper #
##################################

# Figure 1
average_fitted %>%
  ggplot(aes(x=time))+
  geom_line(aes(y=avg_fitted,linetype=drug_group))+
  theme_bw()+facet_wrap(~hb_name)+
  xlab('date')+ylab('number of daily doses per month (estimated)')

# Figure 2
exp_residuals %>% 
  group_by(hb_name) %>% 
  mutate(q1=quantile(exp_residuals,.25),
         q2=quantile(exp_residuals,.5),
         q3=quantile(exp_residuals,.75),
         w=1.5*(q3-q1)) %>% 
  ungroup %>% 
  ggplot()+geom_histogram(aes(exp_residuals))+
  geom_vline(xintercept = 1, lty=2)+
  ylim(-3,18)+
  geom_hline(yintercept=0,col='grey')+
  geom_segment(aes(y=-2,yend=-2,x=q1-w,xend=q1))+
  geom_segment(aes(y=-2,yend=-2,x=q1,xend=q3),lwd=2)+
  geom_segment(aes(y=-2,yend=-2,x=q3,xend=q3+w))+
  geom_point(aes(y=-2,x=q2),lwd=2,pch=16,col='white')+
  geom_point(aes(y=-2,x=q2),lwd=2,pch=1)+
  theme_bw()+facet_wrap(~hb_name)

# Figure 3

# change quantiles as needed
lower_quantile <- 2   # corresponding to 2.5% in the quantiles vector
upper_quantile <- 40  # corresponding to 97.5% in the quantiles vector  

# create the proportions of GP practices meeting the quantile requirements
prop_q_by_hb <- dataset_singletons %>% 
  group_by(hb_name) %>% 
  summarise(HBCode=first(hb_code),
            prop_l_q2.5=mean(exp_res<quant_res[2]),
            prop_g_q97.5=mean(exp_res>quant_res[40]))

crop_map_hb@data <- crop_map_hb@data %>% 
  left_join(prop_q_by_hb,by='HBCode')

# define a scaling factor for the grey scale: the max observed proportion
# is slightly less than 0.2, rescaling aids visualisation
scale_factor <- .3

par(mfrow=c(1,2))
plot(crop_map_hb,
     main='high-volumes prescribers',
     #border='black',
     col=grey(1-crop_map_hb$prop_l_q2.5/scale_factor))

for (i in 0:5) {
  
  rect(470000,750000+i*20000,520000,750000+(i+1)*20000,
       col=grey(seq(1,0,length.out=6))[i+1])
  text(x=545000,y=760000+i*20000,
       paste0(round(scale_factor*seq(0,1,length.out=6)[i+1]*100,2),' %'),
       cex=.5)
}

plot(crop_map_hb,
     main='low-volumes prescribers',
     #border='black',
     col=grey(1-crop_map_hb$prop_g_q97.5/scale_factor))

for (i in 0:5) {
  
  rect(470000,750000+i*20000,520000,750000+(i+1)*20000,
       col=grey(seq(1,0,length.out=6))[i+1])
  text(x=545000,y=760000+i*20000,
       paste0(round(scale_factor*seq(0,1,length.out=6)[i+1]*100,2),' %'),
       cex=.5)
}
par(mfrow=c(1,1))

rm(prop_q_by_hb,scale_factor,
   lower_quantile,upper_quantile) # clean up

###################################
# posterior distributions summary #
###################################

# table of fixed effects posterior estimates
model_fit$summary.fixed

# table of random effects posterior estimates
model_fit$summary.random        # expressed in terms of precision
bri.hyperpar.summary(model_fit) # expressed in terms of standard deviation
