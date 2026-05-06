#### --------- description of script --------- ####

# this script calculates FluB COP and creates COP plots
# Author:Lauren Grant
# Editor: Skyler Rogers
# 05/09/2025

# Manuscript colors:
#B/Phuket: #f28e69
#B/Brisbane: "#104E85"

#### --------- load relevant libraries --------- #### 

library(tidyverse)
library(haven)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(knitr)
library(splines)
library(gridExtra)
library(epiDisplay)
library(readxl)##

#clear workspace
rm(list=ls(all=TRUE))

#### --------- read in the data set --------- #### 
secure_data <- Sys.getenv("SECURE_DATA_PATH")

shiri_aim3_coxph_phuk <- read_csv("Datasets/shiri_aim3_coxph_phuk.csv")
shiri_aim3_coxph_bris <- read_csv("Datasets/shiri_aim3_coxph_bris.csv")

#left-truncated survival time
study_start_date <- as.Date("2017-06-01") #first collection (6-26-17)
study_end_date <- as.Date("2018-04-01") #last onset (3-20-2018)

# make sq and cubic follow-up time variable
#shiri_aim3_coxph_phuk = shiri_aim3_coxph_phuk %>%
#  mutate(study_time_n = as.numeric(study_time),
#         study_time_n_sq = study_time*study_time,
#         study_time_n_cub = study_time*study_time*study_time)
#shiri_aim3_coxph_bris = shiri_aim3_coxph_bris%>%
#  mutate(study_time_n = as.numeric(study_time),
#         study_time_n_sq = study_time*study_time,
#         study_time_n_cub = study_time*study_time*study_time)


coxphag_inf_bris_left <- as.formula("Surv(time_start_numeric, time_stop_numeric,infected_long)~ BRIS_BVIC_NAI_LOG2 + cluster(STUDY_ID)")
coxphag_inf_phuk_left <- as.formula("Surv(time_start_numeric, time_stop_numeric,infected_long)~ PHU_BYAM_NAI_LOG2 + cluster(STUDY_ID)")
fit.coxphag_inf_bris_left <- coxph(coxphag_inf_bris_left, data = shiri_aim3_coxph_bris)
fit.coxphag_inf_phuk_left <- coxph(coxphag_inf_phuk_left, data = shiri_aim3_coxph_phuk)

#pull results###
summary(fit.coxphag_inf_bris_left)
summary(fit.coxphag_inf_phuk_left)


####################################################
#Create predicted risk reduction from model results#
####################################################

#If any numeric covariates are used will need their means
#summary(data$numeric_covariate)

### Make a function to summarize model results###

#antigen= name of the variable with the titers/concentration of interest used in the model
#range= values to create the range of titers/concentrations for that antigen. 
    #For 2 fold assay with minimum detection of 10, it will be 0:9, which creates 5,10,40,40,...2560. 
    #Can adjust accordingly but first number should be the assumed 'no titer' level i.e. 5 for HAI 
#modelc= name of the cox regression for that antigen

model_summarize_log2 = function(antigen,range,modelc) {
  
##create new data frame to estimate predicted values
    #Data frame should include all variables used as covariates in the model
    #if any covariates are numeric, set them to their mean
pred<-data.frame(num=log2(5*2^(range))) #Creates titer values and logs them
pred[[antigen]]<-pred$num #Needed to pull in correct variable name 


##use predict function to get predicted values and standard errors
  preds<-predict(modelc, pred,se.fit = TRUE,type = 'lp')   
  beta<-as.vector(summary(modelc)$coef[,1])
  
#create data frame 
  est<-data.frame(titer=log2(5*2^(range)))
  est$titer2<-2^est$titer
  
#get basic estimates and CIs
  est$mu <-beta[1] * est$titer
  est$hr<-exp(preds$fit)
  est$se<-preds$se.fit
  est$ci_lwr <- with(preds, exp(fit + qnorm(0.025)*se.fit))
  est$ci_upr <- with(preds, exp(fit + qnorm(0.975)*se.fit))
  
# Relative risk ratio. 
  est$rrr <-1-(est$hr/est$hr[1])
  est$rrr <-1-(est$hr/est$hr[1]) #this assumes the first number in the range is the reference value
  est$rrr.ci1 <-1-(est$ci_lwr/est$hr[1])
  est$rrr.ci2 <-1-(est$ci_upr/est$hr[1])
  est$rrr.ci2_2 <- ifelse(est$rrr.ci2 <0,0,est$rrr.ci2 ) #setting lower bound to zero if less than zero for geom_ribbon to show up
  
#Calculate 50 and 80% thresholds. Can calculate others by changing the 0.5/0.8
  est$thresh_50_log=log(exp(est$mu[1])*(-0.5+1))/beta[1] 
  est$thresh_50=2^(est$thresh_50_log) 
  est$thresh_80_log=log(exp(est$mu[1])*(-0.8+1))/beta[1]
  est$thresh_80=2^(est$thresh_80_log)
  
  # return the data frame
  return(est)
  
} 

pred_phuk= model_summarize_log2("PHU_BYAM_NAI_LOG2",0:9,fit.coxphag_inf_phuk_left)
pred_bris= model_summarize_log2("BRIS_BVIC_NAI_LOG2",0:9,fit.coxphag_inf_bris_left)

pred_phuk$Virus = rep("PHU_BYAM_NAI_LOG2")
pred_bris$Virus = rep("BRIS_BVIC_NAI_LOG2")

pred_NAI = rbind(pred_phuk,pred_bris)
pred_NAI$Virus = factor(pred_NAI$Virus,levels=c("PHU_BYAM_NAI_LOG2","BRIS_BVIC_NAI_LOG2"),labels=c('B/Phuket','B/Brisbane'))



#print thresholds
thresholds_NAI = pred_NAI %>%
  filter(titer2 < 10) %>%
  dplyr::select(Virus,thresh_50,thresh_80) 

###############################
#Plot predicted risk reduction#
###############################

#NAI plot
correlates_nai_plot = ggplot(data=pred_NAI,aes(x=titer2,y=rrr,color=Virus,fill=Virus)) +
  geom_line(linewidth=1) +
  geom_ribbon(aes(ymin = rrr.ci1, ymax = rrr.ci2_2), alpha = 0.1,linetype = "dotted", show.legend = FALSE) +
  theme_bw() +
  xlab("NAI titer") +
  ylab("Relative reduction in hazard of infection") +
  geom_hline(yintercept = c(0.50,0.80),linetype = "dashed",linewidth=0.3) +
  scale_color_manual(values=c('#f28e69','#1E76C3')) +
  scale_fill_manual(values=c('#f28e69','#1E76C3')) +
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),limits=c(0,1)) +
  scale_x_continuous(breaks=c(5,10,20,40,80,160,320,640,1280,2560),trans='log2',limits=c(5,2560))+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),legend.text = element_text(size=10))
correlates_nai_plot
ggsave(plot = correlates_nai_plot, width = 10, units = "in", dpi = 300, filename = file.path(secure_data, "Tables + Figures/cop_curve_coxph.png"))

