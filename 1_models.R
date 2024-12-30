rm(list=ls())
## open the libraries
library(tidyverse);library(ggplot2);library(cowplot)
library(brms);library(tidybayes)
set.seed(123)
## set up the work directory 
dir.data<-"C:/Users/chqq3/work/NutNet data/"
dir.graphs<-"C:/Users/chqq3/work/homogenization/Nature communications/graphs/"
setwd(dir.graphs)

##########################################################################################
###### write a function to calculate species diversity at alpha and gamma spatial scales 
###########################################################################################
### focusing on sites with 4-year treatments, also focus on control and NPK
load("data for control and NPK.rdata")

source("funtion to calculate diversity metrics.R")

d7.yr4.10.14<-d7%>%filter(year_trt%in% c(4, 10, 14))
alpha.gamma<-calculate.alpha.gamma.diversity(d7.yr4.10.14)

eff.nut.alpha<-alpha.gamma%>%filter(scale %in% c("alpha"))%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.alpha=log(NPK/Control))%>%mutate(scale="dif.alpha", dif.value=dif.alpha)%>%select(site_code, year_trt, functional_group, q, scale, dif.value)
# check.na<-eff.nut.alpha%>% filter(is.na(Control))%>%filter(is.na(NPK))

eff.nut.beta<-alpha.gamma%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.HillDiv=log(NPK/Control))%>%
  group_by(site_code, year_trt, functional_group, q, scale)%>%summarise(dif.HillDiv1=mean(dif.HillDiv))%>%pivot_wider(names_from = "scale", values_from = "dif.HillDiv1")%>%mutate(dif.beta=gamma - alpha)%>%rename(dif.gamma=gamma)
eff.nut.alpha.beta.gamma<-eff.nut.beta%>% filter(!is.na(alpha))%>%filter(!is.na(dif.gamma))%>%select(site_code, year_trt, functional_group, q, dif.gamma, dif.beta)%>%
  pivot_longer(cols = c("dif.gamma", "dif.beta"), names_to = "scale", values_to = "dif.value")%>%select(site_code, year_trt, functional_group, q, scale, dif.value)%>%
  bind_rows(eff.nut.alpha)%>% mutate(variable.id=paste0(functional_group, scale, year_trt, q))
unique(eff.nut.alpha.beta.gamma$variable.id)

eff.nut.alpha.beta.gamma1<-eff.nut.alpha.beta.gamma%>%filter(q==0 | (year_trt==4 &  q%in%c(1, 2) & functional_group=="all"))
unique(eff.nut.alpha.beta.gamma1$variable.id)
variable.id<-eff.nut.alpha.beta.gamma1%>%ungroup()%>%select(variable.id, functional_group, scale, year_trt, q)%>%distinct()
write.csv(variable.id, file="combinations for species groups and diversity metrics and years and hill q.csv")

# show raw data 
i<-0
(p.overall<-eff.nut.alpha.beta.gamma1%>%filter(q==i)%>%
    ggplot(aes(year_trt, dif.value,  shape=scale, linetype=scale))+theme_cowplot(font_size = 30)+panel_border()+
    facet_wrap(~functional_group, scale="free")+
    geom_point(size=2.5, alpha=0.2)+
    # geom_line()+
    geom_smooth(se=F)+
    labs(x="Years after treatments", y=paste0("Hill number (Q = ", i, ")"), color=NULL, shape=NULL, linetype=NULL)+
    scale_y_continuous(name = paste0("Hill number (Q = ", i, ")")))
# ggsave(p.overall, width = 21, height=14.8, file=paste0("overall trends of alpha and gamma diversity based on Q of ", i, ".pdf"))

###########################################################################################
# run linear mixed-effect models 
###########################################################################################
# effects of nutrient addition on diversity 4 years after treatments

for(vi in unique(eff.nut.alpha.beta.gamma1$variable.id)){
 # vi <- "alldif.alpha40"
  diversity.temp<-eff.nut.alpha.beta.gamma1%>%filter(variable.id==vi)
  mod.eff.nut <- brm( dif.value ~ 1  + (1 | site_code), 
                      data = diversity.temp , cores = 6, iter=3000, warmup = 1000, chains = 6, control = list(adapt_delta = 0.99))
  save(mod.eff.nut, file=paste0("model for diversity change for combination of ",vi, ".Rdata"))
  }


##########################################################################################
### robustness test using sites with 4 blocks at year_trt 4
##########################################################################################
data.more.blocks<-read.csv(file = "raw cover data for sites have four blocks.csv")

alpha.gamma<-calculate.alpha.gamma.diversity(data.more.blocks)
eff.nut.alpha<-alpha.gamma%>%filter(scale %in% c("alpha"))%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.alpha=log(NPK/Control))%>%mutate(scale="dif.alpha", dif.value=dif.alpha)%>%select(site_code, year_trt, functional_group, q, scale, dif.value)

eff.nut.beta<-alpha.gamma%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.HillDiv=log(NPK/Control))%>%
  group_by(site_code, year_trt, functional_group, q, scale)%>%summarise(dif.HillDiv1=mean(dif.HillDiv))%>%pivot_wider(names_from = "scale", values_from = "dif.HillDiv1")%>%mutate(dif.beta=gamma - alpha)%>%rename(dif.gamma=gamma)
eff.nut.alpha.beta.gamma<-eff.nut.beta%>% filter(!is.na(alpha))%>%filter(!is.na(dif.gamma))%>%select(site_code, year_trt, functional_group, q, dif.gamma, dif.beta)%>%
  pivot_longer(cols = c("dif.gamma", "dif.beta"), names_to = "scale", values_to = "dif.value")%>%select(site_code, year_trt, functional_group, q, scale, dif.value)%>%
  bind_rows(eff.nut.alpha)%>% mutate(variable.id=paste0(functional_group, scale, year_trt, q))
unique(eff.nut.alpha.beta.gamma$variable.id)

eff.nut.alpha.beta.gamma1<-eff.nut.alpha.beta.gamma%>%filter(q==0 )
unique(eff.nut.alpha.beta.gamma1$variable.id)

# show raw data 
i<-0
(p.overall<-eff.nut.alpha.beta.gamma1%>%filter(q==i)%>%
    ggplot(aes(year_trt, dif.value,  shape=scale, linetype=scale))+theme_cowplot()+panel_border()+
    facet_grid(site_code~functional_group)+
    geom_point(size=2.5, alpha=0.2)+
    # geom_line()+
    geom_smooth(se=F)+
    labs(x="Years after treatments", y=paste0("Hill number (Q = ", i, ")"), color=NULL, shape=NULL, linetype=NULL))

for(vi in unique(eff.nut.alpha.beta.gamma1$variable.id)){
 # vi <- "alldif.alpha40"
  diversity.temp<-eff.nut.alpha.beta.gamma1%>%filter(variable.id==vi)
  mod.eff.nut <- brm( dif.value ~ 1  + (1 | site_code), 
                      data = diversity.temp , cores = 6, iter=3000, warmup = 1000, chains = 6, control = list(adapt_delta = 0.99))
  save(mod.eff.nut, file=paste0("model for diversity change for combination of ",vi, " for sites with four blocks.Rdata"))
}

##########################################################################################
### robustness test using sites with 5 blocks at year_trt 4
##########################################################################################
data.more.blocks<-read.csv(file = "raw cover data for sites have five blocks.csv")

alpha.gamma<-calculate.alpha.gamma.diversity(data.more.blocks)
eff.nut.alpha<-alpha.gamma%>%filter(scale %in% c("alpha"))%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.alpha=log(NPK/Control))%>%mutate(scale="dif.alpha", dif.value=dif.alpha)%>%select(site_code, year_trt, functional_group, q, scale, dif.value)

eff.nut.beta<-alpha.gamma%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.HillDiv=log(NPK/Control))%>%
  group_by(site_code, year_trt, functional_group, q, scale)%>%summarise(dif.HillDiv1=mean(dif.HillDiv))%>%pivot_wider(names_from = "scale", values_from = "dif.HillDiv1")%>%mutate(dif.beta=gamma - alpha)%>%rename(dif.gamma=gamma)
eff.nut.alpha.beta.gamma<-eff.nut.beta%>% filter(!is.na(alpha))%>%filter(!is.na(dif.gamma))%>%select(site_code, year_trt, functional_group, q, dif.gamma, dif.beta)%>%
  pivot_longer(cols = c("dif.gamma", "dif.beta"), names_to = "scale", values_to = "dif.value")%>%select(site_code, year_trt, functional_group, q, scale, dif.value)%>%
  bind_rows(eff.nut.alpha)%>% mutate(variable.id=paste0(functional_group, scale, year_trt, q))
unique(eff.nut.alpha.beta.gamma$variable.id)

eff.nut.alpha.beta.gamma1<-eff.nut.alpha.beta.gamma%>%filter(q==0 )
unique(eff.nut.alpha.beta.gamma1$variable.id)

# show raw data 
i<-0
(p.overall<-eff.nut.alpha.beta.gamma1%>%filter(q==i)%>%
    ggplot(aes(year_trt, dif.value,  shape=scale, linetype=scale))+theme_cowplot()+panel_border()+
    facet_grid(site_code~functional_group)+
    geom_point(size=2.5, alpha=0.2)+
    # geom_line()+
    geom_smooth(se=F)+
    labs(x="Years after treatments", y=paste0("Hill number (Q = ", i, ")"), color=NULL, shape=NULL, linetype=NULL))

for(vi in unique(eff.nut.alpha.beta.gamma1$variable.id)){
  # vi <- "alldif.alpha40"
  diversity.temp<-eff.nut.alpha.beta.gamma1%>%filter(variable.id==vi)
  mod.eff.nut <- brm( dif.value ~ 1  + (1 | site_code), 
                      data = diversity.temp , cores = 6, iter=3000, warmup = 1000, chains = 6, control = list(adapt_delta = 0.99))
  save(mod.eff.nut, file=paste0("model for diversity change for combination of ",vi, " for sites with five blocks.Rdata"))
}

# the end 
