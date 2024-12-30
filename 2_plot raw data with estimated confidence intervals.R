
rm(list=ls())
## open the library
library(tidyverse);library(ggplot2);library(cowplot);library(brms);library(tidybayes)
set.seed(123)
## set up the work directory 
dir.data<-"C:/Users/chqq3/work/NutNet data/"
dir.graphs<-"C:/Users/chqq3/work/homogenization/Nature communications/graphs/"
setwd(dir.graphs)
options(scipen = 999)
###########################################################################################
##predicted intercepts in alpha, gamma, and beta diversity
###########################################################################################
### focusing on sites with 4-year treatments, also focus on control and NPK
load("data for control and NPK.rdata")
d7.yr4.10.14<-d7%>%filter(year_trt%in% c(4, 10, 14))
source("funtion to calculate diversity metrics.R")
alpha.gamma<-calculate.alpha.gamma.diversity(d7.yr4.10.14)

eff.nut<-alpha.gamma%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.HillDiv=log(NPK) - log(Control))%>%
  group_by(site_code, year_trt, functional_group, q, scale)%>%summarise(dif.HillDiv1=mean(dif.HillDiv))%>%
  pivot_wider(names_from = "scale", values_from = "dif.HillDiv1")%>%mutate(beta=gamma - alpha)%>%
  filter(!is.na(alpha))%>%filter(!is.na(gamma))%>%
  mutate(cat.beta=case_when((beta<0 )~"Homogenization", 
                            (beta>0 )~"Differentiation", 
                            TRUE~"Little change in beta diversity"))%>%
  mutate(cat.process=case_when((alpha>gamma & alpha>0 &gamma>0 )~"Gain of widespread species",
                               (alpha>0 &gamma<0 )~"Spatially restricted species replaced by widespread species",
                               (alpha>gamma & alpha<0 &gamma<0)~"Loss of spatially restricted species",
                               (alpha<gamma & alpha<0 &gamma<0)~"Loss of widespread species",
                               (alpha<0 &gamma>0)~"Widespread species replaced by spatially restricted species",
                               (alpha<gamma & alpha>0 &gamma>0 )~"Gain of spatially restricted species" ,
                               (alpha==gamma & alpha>0 & gamma>0 )~"Gain of spatially restricted and widespread species at simimar magnitude" ,
                               (alpha==gamma & alpha<0 &gamma<0 )~"Loss of spatially restricted and widespread species at simimar magnitude",
                               TRUE ~ "Other situations"))%>% mutate(across(c("alpha", "gamma", "beta" ), ~ round(.x, digits = 4)))

eff.nut$cat.beta<-factor(eff.nut$cat.beta, levels = c("Homogenization", "Differentiation", "Little change in beta diversity"))

eff.nut$cat.process<-factor(eff.nut$cat.process, levels = c("Gain of widespread species", "Spatially restricted species replaced by widespread species", "Loss of spatially restricted species",
                                                            "Loss of widespread species",  "Widespread species replaced by spatially restricted species", "Gain of spatially restricted species" ,
                                                            "Other situations"))

unique(eff.nut$functional_group)
check.year.4<- eff.nut%>%filter(year_trt==4 & q==0 & functional_group=="all")
(beta.group<-table(check.year.4$cat.beta)%>%data.frame())
names(beta.group)<-c("Groups by ∆β", "Number of sites")
(process.group<-table(check.year.4$cat.process)%>%data.frame()%>%mutate(`Groups by ∆β` = rep(c("Homogenization", "Differentiation", " "), times=c(3,3,1))))

check.year.4.native.nonnative<- eff.nut%>%filter(year_trt==4 & q==0 )%>%filter(functional_group %in% c("NAT", "INT"))%>%
   arrange(functional_group, cat.beta, cat.process)
table(check.year.4.native.nonnative$functional_group, check.year.4.native.nonnative$cat.beta)
table(check.year.4.native.nonnative$functional_group, check.year.4.native.nonnative$cat.process)

check.year.4.life.forms<- eff.nut%>%filter(year_trt==4 & q==0 )%>%filter(functional_group %in% c("FORB", "GRAMINOID", "WOODY",  "LEGUME"))%>%arrange(functional_group, cat.beta, cat.process)
table(check.year.4.life.forms$functional_group, check.year.4.life.forms$cat.beta)
table(check.year.4.life.forms$functional_group, check.year.4.life.forms$cat.process)

variable.id.combi<- read.csv("combinations for species groups and diversity metrics and years and hill q.csv")%>%mutate(X=NULL)

summary.diversity.across.sites<-c(); summary.diversity.individual.sites<-c(); summary.model.fit<-c()
for(vi in unique(variable.id.combi$variable.id)){
    # vi<-"alldif.beta41"
    load(paste0("model for diversity change for combination of ",vi, ".Rdata") )   
    #plot(mod.eff.nut)
    #pp_check(mod.eff.nut)
 
   # fixed effect coefficients 
   diversity_fixef <- fixef(mod.eff.nut)%>%data.frame()
   # coefficients for site-level (random) effects
   diversity_coef <- coef(mod.eff.nut)
   diversity_coef2 <-  bind_cols(site_code = rownames(diversity_coef$site_code[,,'Intercept']),
                                       diversity_coef$site_code[,,'Intercept'] %>% 
                                          as_tibble() %>% mutate(Intercept = Estimate,
                                                                 Intercept_lower = Q2.5,
                                                                 Intercept_upper = Q97.5) %>% 
                                          select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) 
  
    # save the predicted intercepts and intercepts at global and site level  
    summary.diversity.across.sites<-summary.diversity.across.sites%>%bind_rows(diversity_fixef%>%mutate(variable.id=vi))
    summary.diversity.individual.sites<-summary.diversity.individual.sites%>%bind_rows(diversity_coef2%>%mutate(variable.id=vi))
    
    # check model fit 
    t.model<-summary(mod.eff.nut)
    t.model1<-t.model$fixed%>%mutate(terms=rownames(.), variable.id=vi)
    N.sites<-length(unique(diversity_coef2$site_code))
    summary.model.fit<-summary.model.fit%>%bind_rows(t.model1%>%mutate(N.sites=N.sites))  
  }

colnames(summary.model.fit)
summary.model.fit1<-summary.model.fit%>%merge(variable.id.combi, by=c("variable.id")) %>%
  mutate(across(c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS" ), ~ round(.x, digits = 4)))%>%
  mutate(Diversity.cross.scales=case_when(scale=="dif.alpha"~"∆α", 
                                     scale=="dif.gamma"~"∆γ",
                                     scale=="dif.beta"~"∆β"), 
         Diversity.metrics=case_when(q==0 ~"Species richness", 
                                     q==1 ~"Shannon diversity",
                                     q==2 ~"Simpson diversity"), 
         Species.groups=case_when(functional_group=="all" ~ "All species", 
                                  functional_group=="NAT" ~ "Native", 
                                  functional_group=="INT" ~ "Non-native", TRUE~functional_group))%>%
  select("year_trt", "Diversity.metrics", "Diversity.cross.scales", "Species.groups", "Estimate", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS", "N.sites")%>%arrange(year_trt, Diversity.metrics, Diversity.cross.scales)

write.csv(summary.model.fit1, file="estimated overall effects and model fit.csv")
# check for significant effects 
summary.model.fit1.sig<-summary.model.fit1%>%filter((`l-95% CI` <0 & `u-95% CI` <0)| (`l-95% CI` >0 &`u-95% CI` >0))
# check for sample size 
check.sample.size<-summary.model.fit1%>%select("year_trt", "Diversity.metrics", "Diversity.cross.scales", "Species.groups","N.sites")%>%distinct()%>%arrange(year_trt, Species.groups, Diversity.cross.scales, Diversity.metrics)

summary.diversity.across.sites1<-summary.diversity.across.sites%>%merge(variable.id.combi, by=c("variable.id")) %>%
  mutate(ci_lower=Q2.5 , ci_upper=Q97.5, name=str_sub(scale, 5))%>%select(year_trt, functional_group, name, Estimate,  ci_lower, ci_upper, year_trt, q)

summary.diversity.individual.sites1<-summary.diversity.individual.sites%>%merge(variable.id.combi, by=c("variable.id")) %>%
  mutate(Estimate=Intercept, ci_lower=Intercept_lower, ci_upper=Intercept_upper)%>%
  mutate(Diversity.cross.scales=case_when(scale=="dif.alpha"~"∆α", 
                                          scale=="dif.gamma"~"∆γ",
                                          scale=="dif.beta"~"∆β"), 
         Diversity.metrics=case_when(q==0 ~"Species richness", 
                                     q==1 ~"Shannon diversity",
                                     q==2 ~"Simpson diversity"), 
         Species.groups=case_when(functional_group=="all" ~ "All species", 
                                  functional_group=="NAT" ~ "Native", 
                                  functional_group=="INT" ~ "Non-native", TRUE~functional_group))%>%
  mutate(across(c("Estimate", "ci_lower", "ci_upper" ), ~ round(.x, digits = 4)))%>%
   select(year_trt, Diversity.metrics,  Diversity.cross.scales, Species.groups, site_code, Estimate,  ci_lower, ci_upper)
write.csv(summary.diversity.individual.sites1, file="estimated effects at individual sites.csv")

summary.diversity.individual.sites1.sig<-summary.diversity.individual.sites1%>%filter((ci_lower <0 & ci_upper <0)| (ci_lower >0 & ci_upper>0))%>%filter(year_trt==4)%>%
  filter(Diversity.metrics=="Species richness")%>%filter(Species.groups=="All species")
table(summary.diversity.individual.sites1.sig$Diversity.cross.scales)

eff.nut.l<-eff.nut%>%pivot_longer(cols = c("alpha", "gamma", "beta"), names_to = "scale0", values_to = "raw data")%>%mutate(scale=paste0("dif.", scale0))
  
summary.diversity.individual.sites1.yr4<- summary.diversity.individual.sites%>%merge(variable.id.combi, by=c("variable.id"))  %>%
  filter(year_trt==4)%>%filter(q==0)%>% mutate(across(c("Intercept",  "Intercept_lower", "Intercept_upper"), ~ round(.x, digits = 4)))%>%
  mutate(mean.Intercept= paste0(Intercept, " [", Intercept_lower, ", ", Intercept_upper, "]"))%>%select("year_trt", "scale", "functional_group", "site_code", "Intercept",  "Intercept_lower", "Intercept_upper", "mean.Intercept")%>%
  pivot_wider(names_from = "scale", values_from = c("Intercept",  "Intercept_lower", "Intercept_upper", "mean.Intercept"))%>%
  mutate(cat.beta=case_when((Intercept_dif.beta<0 & Intercept_lower_dif.beta<0 &Intercept_upper_dif.beta<0)~"Homogenization", 
                            (Intercept_dif.beta>0 & Intercept_lower_dif.beta>0 &Intercept_upper_dif.beta>0)~"Differentiation", 
                            TRUE~"Little change in beta diversity"))%>%
  mutate(cat.process.ci=case_when((Intercept_dif.alpha>Intercept_dif.gamma & Intercept_dif.alpha>0 &Intercept_dif.gamma>0 & ((Intercept_lower_dif.alpha > 0 & Intercept_upper_dif.alpha > 0) & (Intercept_lower_dif.gamma > 0 & Intercept_upper_dif.gamma > 0)))~"Gain of widespread species",
                                  (Intercept_dif.alpha>0 &Intercept_dif.gamma<0 & ((Intercept_lower_dif.alpha > 0 & Intercept_upper_dif.alpha > 0) & (Intercept_lower_dif.gamma < 0 & Intercept_upper_dif.gamma < 0)))~"Spatially restricted species replaced by widespread species",
                                  (Intercept_dif.alpha>Intercept_dif.gamma & Intercept_dif.alpha<0 &Intercept_dif.gamma<0 & ((Intercept_lower_dif.alpha < 0 & Intercept_upper_dif.alpha < 0) & (Intercept_lower_dif.gamma < 0 & Intercept_upper_dif.gamma < 0)))~"Loss of spatially restricted species",
                                  (Intercept_dif.alpha<Intercept_dif.gamma & Intercept_dif.alpha<0 &Intercept_dif.gamma<0 & ((Intercept_lower_dif.alpha < 0 & Intercept_upper_dif.alpha < 0) & (Intercept_lower_dif.gamma < 0 & Intercept_upper_dif.gamma < 0)))~"Loss of widespread species",
                                  (Intercept_dif.alpha<0 &Intercept_dif.gamma>0 & ((Intercept_lower_dif.alpha < 0 & Intercept_upper_dif.alpha < 0) & (Intercept_lower_dif.gamma > 0 & Intercept_upper_dif.gamma > 0)))~"Widespread species replaced by spatially restricted species",
                                  (Intercept_dif.alpha<Intercept_dif.gamma & Intercept_dif.alpha>0 &Intercept_dif.gamma>0 & ((Intercept_lower_dif.alpha > 0 & Intercept_upper_dif.alpha > 0) & (Intercept_lower_dif.gamma > 0 & Intercept_upper_dif.gamma > 0)))~"Gain of spatially restricted species" ,
                                  TRUE~"Other situations"))%>%
  mutate(cat.process=case_when((Intercept_dif.alpha>Intercept_dif.gamma & Intercept_dif.alpha>0 &Intercept_dif.gamma>0 )~"Gain of widespread species",
                               (Intercept_dif.alpha>0 &Intercept_dif.gamma<0 )~"Spatially restricted species replaced by widespread species",
                               (Intercept_dif.alpha>Intercept_dif.gamma & Intercept_dif.alpha<0 &Intercept_dif.gamma<0)~"Loss of spatially restricted species",
                               (Intercept_dif.alpha<Intercept_dif.gamma & Intercept_dif.alpha<0 &Intercept_dif.gamma<0)~"Loss of widespread species",
                               (Intercept_dif.alpha<0 &Intercept_dif.gamma>0)~"Widespread species replaced by spatially restricted species",
                               (Intercept_dif.alpha<Intercept_dif.gamma & Intercept_dif.alpha>0 &Intercept_dif.gamma>0 )~"Gain of spatially restricted species" ,
                               TRUE~"Other situations"))%>%
  merge(eff.nut%>%filter(year_trt==4)%>%filter(q==0)%>% mutate(across(c("alpha", "beta",  "gamma"), ~ round(.x, digits = 4)), cat.beta=NULL, cat.process=NULL), by=c("site_code", "year_trt", "functional_group"))%>%
  select("year_trt", "functional_group", "site_code", "alpha", "mean.Intercept_dif.alpha", "beta", "mean.Intercept_dif.beta" , "gamma", "mean.Intercept_dif.gamma" , "cat.beta", "cat.process", "cat.process.ci")%>%
  arrange(year_trt, functional_group)

table(summary.diversity.individual.sites1.yr4$functional_group, summary.diversity.individual.sites1.yr4$cat.beta)
table(summary.diversity.individual.sites1.yr4$functional_group, summary.diversity.individual.sites1.yr4$cat.process)
table(summary.diversity.individual.sites1.yr4$functional_group, summary.diversity.individual.sites1.yr4$cat.process.ci)
write.csv(summary.diversity.individual.sites1.yr4, file="estimated effects at individual sites in year 4 with different diversity change scenarios.csv")

##########################################################################################
###### change in species richness across scales in year 4, 10, and 14  
###########################################################################################
source("function to draw background for main figures.R")
for(yr in c(4, 10, 14)){ 
  # yr<-4
  eff.nut.rich.yr<-eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c("all"))
  range.values<-range(range(eff.nut.rich.yr$alpha), range(eff.nut.rich.yr$gamma))
 
   (beta.group<-table(eff.nut.rich.yr$cat.beta)%>%data.frame())
  names(beta.group)<-c("Groups.beta", "Number.sites")
  (process.group<-table(eff.nut.rich.yr$cat.process)%>%data.frame())
  names(process.group)<-c("Groups.process", "Number.sites")
  
  eff.nut.overall.yr<-summary.diversity.across.sites1%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c("all"))%>%pivot_wider(names_from = "name", values_from = c("Estimate",  "ci_lower", "ci_upper"))%>%
    mutate(dif.positive=ci_upper_beta - Estimate_beta)
  
  ppp1<-plot.six.scenarios(range.values[1]-0.05, range.values[2]+0.05)
  
  angle.positive<- 120; angle.negative<- -60; dif.positive<-eff.nut.overall.yr$dif.positive
  min_value<- range.values[1];  max_value<- ifelse (range.values[2] >0.5, range.values[2], 0.5)
  (ppp2<-ppp1+
      annotate(geom = "text", x = min_value/1.7, y = max_value/1.7, angle = 45, label = "Differentiation", color="white", size =5, fontface = "bold") +  
      annotate(geom = "text", x = min_value/1.85, y = max_value/1.85, angle = 45, label = "Increase in ∆β (LRR)", color="white", size =4) +  
      annotate(geom = "text", x = min_value/2.05, y = max_value/2.05, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Differentiation", ]$Number.sites,")"), color="white", size =4) +  
      annotate(geom = "text", x = max_value/1.7, y = min_value/1.7, angle = 45, label = "Homogenization", color="white", size =5, fontface = "bold") +
      annotate(geom = "text", x = max_value/1.55, y = min_value/1.55, angle = 45, label = "Decrease in ∆β (LRR)", color="white", size =4) +
      annotate(geom = "text", x = max_value/1.45, y = min_value/1.45, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Homogenization", ]$Number.sites,")"),  color="white", size =4) +
 
      annotate(geom = "text", x = 0.9*max_value, y=0.1*max_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Gain of widespread species", ]$Number.sites,")"), size =4) +
      annotate(geom = "text",  x = 0.9*max_value, y =  0.1*min_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Spatially restricted species replaced by widespread species", ]$Number.sites,")"),  size =4) +
      annotate(geom = "text", x = 0.1*min_value, y = 0.9*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of spatially restricted species", ]$Number.sites,")"), size =4) +
     
       annotate(geom = "text", x = 0.1*max_value,  y = 0.9*max_value,angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Gain of spatially restricted species", ]$Number.sites,")"),  size =4) +
      annotate(geom = "text", x =0.9*min_value, y = 0.1*max_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Widespread species replaced by spatially restricted species", ]$Number.sites,")"),  size =4) +
      annotate(geom = "text", x = 0.9*min_value, y = 0.1*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of widespread species", ]$Number.sites,")"),  size =4)+
     
       geom_point(data=eff.nut.rich.yr, aes(alpha, gamma), size=1.5, alpha=0.3)+
      geom_errorbar(data=eff.nut.overall.yr, aes(x=Estimate_alpha, y=Estimate_gamma, xmin=ci_lower_alpha, xmax=ci_upper_alpha), width=0.0001, color="black", linewidth=1, alpha=0.7) +
      geom_errorbar(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, ymin=ci_lower_gamma, ymax=ci_upper_gamma), width=0.0001, color="black", linewidth=1,  alpha=0.7) +
      geom_segment(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.positive * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.positive * pi / 180)),  lineend = "square", color="black", linewidth=1, alpha=0.7) +
      geom_segment(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.negative * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.negative * pi / 180)),   lineend = "square", color="black", linewidth=1, alpha=0.7) +
      geom_point(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma), size=5, stroke = 1, pch=21, alpha=0.7, color="black"))
if(yr==4){
  ppp3<-ppp2
}else{
  (ppp3<-ppp2+labs(title = paste0(nrow(eff.nut.rich.yr), " sites with ", yr, " years after nutrient addition began")))
}
 ggsave(ppp3, width = 15, height = 15, dpi=600, unit="cm",  file=paste0("change in species richness across scales at year ",yr, ".png"))
}

###########################################################################################
####### model fit for change in species richness across scales in year 4
###########################################################################################
variable.id.combi_yr4_all.species<-variable.id.combi%>%filter(year_trt==4 & functional_group=="all" & q==0)

list_plots_check <- vector('list', 3)
for( vi in unique(variable.id.combi_yr4_all.species$variable.id)){
  # vi<-"alldif.beta40"
  diversity.metric.temp<-ifelse(grepl("beta", vi), "beta", ifelse(grepl("alpha", vi), "alpha", "gamma"))
  
  load(paste0("model for diversity change for combination of ",vi, ".Rdata") )   
  # check model fitting 
  check.pred<-pp_check(mod.eff.nut)+ theme_cowplot()+panel_border()+
    labs(x=paste0("Change in ", diversity.metric.temp, " diversity"), y="Density")
  list_plots_check[[vi]]<-check.pred
}
legend<-get_legend(check.pred)

(check.fit<-plot_grid(list_plots_check$alldif.alpha40 +theme(legend.position = "none"), 
                      list_plots_check$alldif.gamma40+theme(legend.position = "none"), 
                      list_plots_check$alldif.beta40+theme(legend.position = "none"),
                      labels = "AUTO", nrow=1))
(check.fit1<-plot_grid(check.fit, 
                       legend, rel_widths = c(9, 1)))
ggsave(check.fit1, width = 21, height=10, unit="cm",  file=paste0("check model fitting for change in species richness across scales in year 4.png"))

##########################################################################################
###### change in species richness for native and non-native across scales in year 4, 10, and 14 
###########################################################################################

list_plots <- vector('list', 20)
for(yr in c(4, 10, 14)){ 
  for(fg in c("FORB", "GRAMINOID",  "LEGUME", "WOODY", "NAT", "INT")){
    # fg <-"WOODY" ; yr<-4
    eff.nut.rich.yr.groups<-eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c(fg)) %>% filter(!is.na(alpha))%>%filter(!is.na(gamma))
    range.values<-range(range(eff.nut.rich.yr.groups$alpha), range(eff.nut.rich.yr.groups$gamma))
    
    (beta.group<-table(eff.nut.rich.yr.groups$cat.beta)%>%data.frame())
    names(beta.group)<-c("Groups.beta", "Number.sites")
    (process.group<-table(eff.nut.rich.yr.groups$cat.process)%>%data.frame())
    names(process.group)<-c("Groups.process", "Number.sites")
    
    eff.nut.overall.yr.groups<-summary.diversity.across.sites1%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c(fg))%>%pivot_wider(names_from = "name", values_from = c("Estimate", "ci_lower", "ci_upper"))%>%
      mutate(dif.positive=ci_upper_beta - Estimate_beta)
    dif.positive<-eff.nut.overall.yr.groups$dif.positive
    
    pp.temp<-plot.six.scenarios(range.values[1]-0.1, range.values[2]+0.1)
    
    min_value<- ifelse (range.values[1] > - 0.5, -0.5, range.values[1]);  max_value<- ifelse (range.values[2] >0.5, range.values[2], 0.5)
    (pp.temp1<-pp.temp+ 
        annotate(geom = "text", x = min_value/1.7, y = max_value/1.7, angle = 45, label = "Differentiation", color="white", size =4, fontface = "bold") +  
        annotate(geom = "text", x = min_value/1.9, y = max_value/1.9, angle = 45, label = "Increase in ∆β (LRR)", color="white", size =3) +  
        annotate(geom = "text", x = min_value/2.15, y = max_value/2.15, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Differentiation", ]$Number.sites,")"), color="white", size =3) +  
        annotate(geom = "text", x = max_value/1.7, y = min_value/1.7, angle = 45, label = "Homogenization", color="white", size =4, fontface = "bold") +
        annotate(geom = "text", x = max_value/1.5, y = min_value/1.5, angle = 45, label = "Decrease in ∆β (LRR)", color="white", size =3) +
        annotate(geom = "text", x = max_value/1.35, y = min_value/1.35, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Homogenization", ]$Number.sites,")"),  color="white", size =3) +
        
        annotate(geom = "text", x = 0.9*max_value, y=0.1*max_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Gain of widespread species", ]$Number.sites,")"), size =3) +
        annotate(geom = "text",  x = 0.9*max_value, y =  0.1*min_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Spatially restricted species replaced by widespread species", ]$Number.sites,")"),  size =3) +
        annotate(geom = "text", x = 0.1*min_value, y = 0.9*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of spatially restricted species", ]$Number.sites,")"), size =3) +
        
        annotate(geom = "text", x = 0.1*max_value,  y = 0.9*max_value,angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Gain of spatially restricted species", ]$Number.sites,")"),  size =3) +
        annotate(geom = "text", x =0.9*min_value, y = 0.1*max_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Widespread species replaced by spatially restricted species", ]$Number.sites,")"),  size =3) +
        annotate(geom = "text", x = 0.9*min_value, y = 0.1*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of widespread species", ]$Number.sites,")"),  size =3)+
        
        geom_point(data=eff.nut.rich.yr.groups, aes(alpha, gamma), size=1.5, alpha=0.3)+
        geom_errorbar(data=eff.nut.overall.yr.groups, aes(x=Estimate_alpha, y=Estimate_gamma, xmin=ci_lower_alpha, xmax=ci_upper_alpha), width=0.0001, color="black", linewidth=1, alpha=0.7) +
        geom_errorbar(data=eff.nut.overall.yr.groups,  aes(x=Estimate_alpha, y=Estimate_gamma, ymin=ci_lower_gamma, ymax=ci_upper_gamma), width=0.0001, color="black", linewidth=1,  alpha=0.7) +
        geom_segment(data=eff.nut.overall.yr.groups,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.positive * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.positive * pi / 180)),  lineend = "square", color="black", linewidth=1, alpha=0.7) +
        geom_segment(data=eff.nut.overall.yr.groups,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.negative * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.negative * pi / 180)),   lineend = "square", color="black", linewidth=1, alpha=0.7) +
        geom_point(data=eff.nut.overall.yr.groups,  aes(x=Estimate_alpha, y=Estimate_gamma), size=5, stroke = 1, pch=21, alpha=0.7, color="black")+
        labs(x=NULL, y=NULL))
    
    list_plots[[fg]]<-pp.temp1
  }
  
  # combine the graphs for life forms
  nsites.FORB<- nrow(eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c("FORB")) %>% filter(!is.na(alpha))%>%filter(!is.na(gamma)))
  nsites.GRAMINOID<- nrow(eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c("GRAMINOID")) %>% filter(!is.na(alpha))%>%filter(!is.na(gamma)))
  nsites.LEGUME<- nrow(eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c("LEGUME")) %>% filter(!is.na(alpha))%>%filter(!is.na(gamma)))
  nsites.WOODY<- nrow(eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c("WOODY")) %>% filter(!is.na(alpha))%>%filter(!is.na(gamma)))
  
  (pp.life.forms1<-plot_grid(list_plots$FORB,
                             list_plots$GRAMINOID,
                             list_plots$LEGUME,
                             list_plots$WOODY,
                             nrow=2, ncol=2, 
                             hjust = -0.5,
                             labels = c(paste0("A (forb; ", nsites.FORB," sites)"), paste0("B (graminoid; ", nsites.GRAMINOID," sites)"), paste0("C (legume; ", nsites.LEGUME," sites)"), paste0("D (woody; ", nsites.WOODY," sites)")), label_fontface = "bold", label_size = 11))
  (pp.life.forms3<-pp.life.forms1+ theme(plot.margin = margin(l = 10, b = 10, t=15))+
      draw_label(expression(bar("∆α") ~ " (LRR)"), x = 0.5, y = 0, size = 11) + draw_label("∆γ (LRR)", x = 0, y = 0.5, angle = 90, size = 11))
  ggsave(pp.life.forms3, width = 20, height = 20, dpi=600, unit="cm",  file=paste0("effects on change in diversity across scales for life forms at year ",yr, ".png"))
  
  # combine the graphs for native and non-native species 
  nsites.NAT<- nrow(eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c("NAT")) %>% filter(!is.na(alpha))%>%filter(!is.na(gamma)))
  nsites.INT<- nrow(eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c("INT")) %>% filter(!is.na(alpha))%>%filter(!is.na(gamma)))
  
  (pp.origins1<-plot_grid(list_plots$NAT,
                          list_plots$INT,
                          nrow=1, ncol=2, 
                          # vjust = -0.01, 
                          hjust=-0.25,
                          labels = c(paste0("A (native species; ", nsites.NAT," sites)"), paste0("B (non-native species; ", nsites.INT," sites)") ), label_fontface = "bold", label_size = 11))
  # adjust the t (top), b (bottom), l (left), and r (right) values for controlling the margins 
  (pp.origins3<-pp.origins1+ theme(plot.margin = margin(l = 10, b = 10, t=15))+
      draw_label(expression(bar("∆α") ~ " (LRR)"), x = 0.5, y = 0, size = 11) + draw_label("∆γ (LRR)", x = 0, y = 0.5, angle = 90, size = 11))
  ggsave(pp.origins3, width = 20, height = 10, dpi=600, unit="cm",  file=paste0("effects on change in diversity across scales for native and non-native species at year ",yr, ".png"))
}

##########################################################################################
######supplementary figures: change in Hill Numbers with q of 1 and 2 in year 4
###########################################################################################
yr<-4
list_plots1 <- vector('list', 4)
for(hill.q in c(1, 2)){ 
  # hill.q<- 1
  eff.nut.rich.yr<-eff.nut%>%filter(year_trt%in% c(4))%>%filter(q%in% c(hill.q))%>%filter(functional_group%in% c("all"))
  range.values<-range(range(eff.nut.rich.yr$alpha), range(eff.nut.rich.yr$gamma))
 
   (beta.group<-table(eff.nut.rich.yr$cat.beta)%>%data.frame())
  names(beta.group)<-c("Groups.beta", "Number.sites")
  (process.group<-table(eff.nut.rich.yr$cat.process)%>%data.frame())
  names(process.group)<-c("Groups.process", "Number.sites")
  
  eff.nut.overall.yr<-summary.diversity.across.sites1%>%filter(year_trt%in% c(4))%>%filter(q%in% c(hill.q))%>%filter(functional_group%in% c("all"))%>%pivot_wider(names_from = "name", values_from = c("Estimate", "ci_lower", "ci_upper"))%>%
    mutate(dif.positive=ci_upper_beta - Estimate_beta)
  
  ppp1<-plot.six.scenarios(range.values[1]-0.05, range.values[2]+0.05)
  
  angle.positive<- 120; angle.negative<- -60; dif.positive<-eff.nut.overall.yr$dif.positive
  min_value<- range.values[1];  max_value<- ifelse (range.values[2] >0.5, range.values[2], 0.5)
  (ppp2<-ppp1+
      annotate(geom = "text", x = min_value/1.7, y = max_value/1.7, angle = 45, label = "Differentiation", color="white", size =4, fontface = "bold") +  
      annotate(geom = "text", x = min_value/1.9, y = max_value/1.9, angle = 45, label = "Increase in ∆β (LRR)", color="white", size =3) +  
      annotate(geom = "text", x = min_value/2.15, y = max_value/2.15, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Differentiation", ]$Number.sites,")"), color="white", size =3) +  
      annotate(geom = "text", x = max_value/1.7, y = min_value/1.7, angle = 45, label = "Homogenization", color="white", size =4, fontface = "bold") +
      annotate(geom = "text", x = max_value/1.5, y = min_value/1.5, angle = 45, label = "Decrease in ∆β (LRR)", color="white", size =3) +
      annotate(geom = "text", x = max_value/1.35, y = min_value/1.35, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Homogenization", ]$Number.sites,")"),  color="white", size =3) +
      
      annotate(geom = "text", x = 0.9*max_value, y=0.1*max_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Gain of widespread species", ]$Number.sites,")"), size =3) +
      annotate(geom = "text",  x = 0.9*max_value, y =  0.1*min_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Spatially restricted species replaced by widespread species", ]$Number.sites,")"),  size =3) +
      annotate(geom = "text", x = 0.1*min_value, y = 0.9*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of spatially restricted species", ]$Number.sites,")"), size =3) +
      
      annotate(geom = "text", x = 0.1*max_value,  y = 0.9*max_value,angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Gain of spatially restricted species", ]$Number.sites,")"),  size =3) +
      annotate(geom = "text", x =0.9*min_value, y = 0.1*max_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Widespread species replaced by spatially restricted species", ]$Number.sites,")"),  size =3) +
      annotate(geom = "text", x = 0.9*min_value, y = 0.1*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of widespread species", ]$Number.sites,")"),  size =3)+
      geom_point(data=eff.nut.rich.yr, aes(alpha, gamma), size=1.5, alpha=0.3)+
      geom_errorbar(data=eff.nut.overall.yr, aes(x=Estimate_alpha, y=Estimate_gamma, xmin=ci_lower_alpha, xmax=ci_upper_alpha), width=0.0001, color="black", linewidth=1, alpha=0.7) +
      geom_errorbar(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, ymin=ci_lower_gamma, ymax=ci_upper_gamma), width=0.0001, color="black", linewidth=1,  alpha=0.7) +
      geom_segment(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.positive * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.positive * pi / 180)),  lineend = "square", color="black", linewidth=1, alpha=0.7) +
      geom_segment(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.negative * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.negative * pi / 180)),   lineend = "square", color="black", linewidth=1, alpha=0.7) +
      geom_point(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma), size=5, stroke = 1, pch=21, alpha=0.7, color="black")+
      labs(x=NULL, y=NULL))
  hill.q1<-paste0("hill", hill.q)
  list_plots1[[hill.q1]]<-ppp2
}

(pp.div1<-plot_grid(list_plots1$hill1,
                    list_plots1$hill2,
                    nrow=1, ncol=2, 
                    # vjust = -0.01, hjust=-0.2,
                    labels = c("A (Hill number, q = 1)", "B (Hill number, q = 2)"), label_fontface = "bold", label_size = 11))
(pp.div3<-pp.div1+ theme(plot.margin = margin(l = 10, b = 10, t=15))+
    draw_label(expression(bar("∆α") ~ " (LRR)"), x = 0.5, y = 0, size = 11) + draw_label("∆γ (LRR)", x = 0, y = 0.5, angle = 90, size = 11))
ggsave(pp.div3, width = 20, height = 10, dpi=600, unit="cm",  file=paste0("effects on change in Shannon diversity and Simpson diversity at year ",yr, ".png"))

# the end 
