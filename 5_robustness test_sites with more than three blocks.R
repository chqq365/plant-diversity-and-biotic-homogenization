
rm(list=ls())
## open the library
library(tidyverse);library(brms);library(tidybayes);library(ggplot2);library(cowplot)
set.seed(123)
## set up the work directory 
dir.data<-"C:/Users/chqq3/work/NutNet data/"
dir.graphs<-"C:/Users/chqq3/work/homogenization/Nature communications/graphs1/"
setwd(dir.graphs)

###########################################################################################
##predicted intercepts in alpha, gamma, and beta diversity
###########################################################################################

num.block<- "five"
  
data.more.blocks<-read.csv(file = paste0("raw cover data for sites have ", num.block ," blocks.csv"))

source("funtion to calculate diversity metrics.R")
alpha.gamma<-calculate.alpha.gamma.diversity(data.more.blocks)

eff.nut<-alpha.gamma%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.HillDiv=log(NPK) - log(Control))%>%
  group_by(site_code, year_trt, functional_group, q, scale)%>%summarise(dif.HillDiv1=mean(dif.HillDiv))%>%pivot_wider(names_from = "scale", values_from = "dif.HillDiv1")%>%mutate(beta=gamma - alpha)%>%
  filter(!is.na(alpha))%>%filter(!is.na(gamma))%>%
  mutate(cat.beta=case_when((beta<0 )~"Homogenization", 
                            (beta>0 )~"Differentiation", 
                            TRUE~"Little change in beta diversity"))%>%
  mutate(cat.process=case_when((alpha>gamma & alpha>0 &gamma>0 )~"Gain of widespread species",
                               (alpha>0 &gamma<0 )~"Spatially restricted replaced by widespread species",
                               (alpha>gamma & alpha<0 &gamma<0)~"Loss of spatially restricted species",
                               (alpha<gamma & alpha<0 &gamma<0)~"Loss of widespread species",
                               (alpha<0 &gamma>0)~"Widespread replaced by spatially restricted species",
                               (alpha<gamma & alpha>0 &gamma>0 )~"Gain of spatially restricted species" ,
                               (alpha==gamma & alpha>0 & gamma>0 )~"Gain of spatially restricted and widespread species at simimar magnitude" ,
                               (alpha==gamma & alpha<0 &gamma<0 )~"Loss of spatially restricted and widespread species at simimar magnitude",
                               TRUE ~ "Other situations"))%>% mutate(across(c("alpha", "gamma", "beta" ), ~ round(.x, digits = 4)))

eff.nut$cat.beta<-factor(eff.nut$cat.beta, levels = c("Homogenization", "Differentiation", "Little change in beta diversity"))

eff.nut$cat.process<-factor(eff.nut$cat.process, levels = c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                                            "Loss of widespread species",  "Widespread replaced by spatially restricted species", "Gain of spatially restricted species" ,
                                                            "Loss of spatially restricted and widespread species at simimar magnitude", "Gain of spatially restricted and widespread species at simimar magnitude",
                                                            "Other situations"))

variable.id.combi<- read.csv("combinations for species groups and diversity metrics and years and hill q.csv")%>%mutate(X=NULL)%>%filter(q==0)%>%filter(year_trt==4)

summary.diversity.across.sites<-c(); summary.diversity.individual.sites<-c(); summary.model.fit<-c()
for(vi in unique(variable.id.combi$variable.id)){
  # vi<-"alldif.gamma40"
  load(paste0("model for diversity change for combination of ",vi, " for sites with ", num.block," blocks.Rdata"))
  #plot(mod.eff.nut)
  #pp_check(mod.eff.nut)
  
  # fixed effect coefficients 
  diversity_fixef <- fixef(mod.eff.nut)%>%data.frame()
  # coefficients for site-level (random) effects
  diversity_coef <- coef(mod.eff.nut)
  extract.site.names<-rownames(diversity_coef$site_code[,,'Intercept'])
  if(is.null(extract.site.names)){ 
    diversity_coef2 <-  t(diversity_coef$site_code[,,'Intercept']) %>%data.frame()%>% 
      mutate(Intercept = Estimate, Intercept_lower = Q2.5,  Intercept_upper = Q97.5) %>% 
      select(-Estimate, -Est.Error, -Q2.5, -Q97.5)%>%mutate(site_code="unknown")
  } else{
    diversity_coef2 <-  diversity_coef$site_code[,,'Intercept'] %>%data.frame()%>% 
      mutate(Intercept = Estimate, Intercept_lower = Q2.5,  Intercept_upper = Q97.5) %>% 
      select(-Estimate, -Est.Error, -Q2.5, -Q97.5)%>%mutate(site_code= rownames(diversity_coef$site_code[,,'Intercept']))
  }
  
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
  mutate(across(c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS" ), ~ round(.x, digits = 2)))%>%
  mutate(Diversity.cross.scales=case_when(scale=="dif.alpha"~"∆α", 
                                          scale=="dif.gamma"~"∆γ",
                                          scale=="dif.beta"~"∆β"), 
         Diversity.metrics=case_when(q==0 ~"Species richness", 
                                     q==1 ~"Shannon diversity",
                                     q==2 ~"Simpson diversity"), 
         Species.groups=case_when(functional_group=="all" ~ "All species", 
                                  functional_group=="NAT" ~ "Native", 
                                  functional_group=="INT" ~ "Non-native", TRUE~functional_group))%>%
  select("year_trt", "Diversity.cross.scales", "Diversity.metrics", "Species.groups", "Estimate", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS", "N.sites")%>%arrange(year_trt, Diversity.metrics, Diversity.cross.scales)

# write.csv(summary.model.fit1, file="estimated overall effects and model fit for sites with ", num.block," blocks.csv")
summary.model.fit1.sig<-summary.model.fit1%>%filter((`l-95% CI` <0 & `u-95% CI` <0)| (`l-95% CI` >0 &`u-95% CI` >0))
# check for sample size 
check.sample.size<-summary.model.fit1%>%select("year_trt", "Diversity.metrics", "Diversity.cross.scales", "Species.groups","N.sites")%>%distinct()%>%arrange(year_trt, Species.groups, Diversity.cross.scales, Diversity.metrics)

summary.diversity.across.sites1<-summary.diversity.across.sites%>%merge(variable.id.combi, by=c("variable.id")) %>%
  mutate(ci_lower=Q2.5 , ci_upper=Q97.5, name=str_sub(scale, 5))%>%select(functional_group, name, Estimate,  ci_lower, ci_upper, year_trt, q)

##########################################################################################
###### change in species richness across scales in year 4
###########################################################################################
source("function to draw background for main figures.R")
yr<-4
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
ppp2<-ppp1+
  annotate(geom = "text", x = min_value/1.7, y = max_value/1.7, angle = 45, label = "Differentiation", color="white", size =6, fontface = "bold") +  
  annotate(geom = "text", x = min_value/1.85, y = max_value/1.85, angle = 45, label = "Increase in ∆β (LRR)", color="white", size =5) +  
  annotate(geom = "text", x = min_value/2.05, y = max_value/2.05, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Differentiation", ]$Number.sites,")"), color="white", size =5) +  
  annotate(geom = "text", x = max_value/1.7, y = min_value/1.7, angle = 45, label = "Homogenization", color="white", size =6, fontface = "bold") +
  annotate(geom = "text", x = max_value/1.55, y = min_value/1.55, angle = 45, label = "Decrease in ∆β (LRR)", color="white", size =5) +
  annotate(geom = "text", x = max_value/1.45, y = min_value/1.45, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Homogenization", ]$Number.sites,")"),  color="white", size =5) +
  
  annotate(geom = "text", x = 0.9*max_value, y=0.1*max_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Gain of widespread species", ]$Number.sites,")"), size =5) +
  annotate(geom = "text",  x = 0.9*max_value, y =  0.1*min_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Spatially restricted replaced by widespread species", ]$Number.sites,")"),  size =5) +
  annotate(geom = "text", x = 0.1*min_value, y = 0.9*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of spatially restricted species", ]$Number.sites,")"), size =5) +
  
  annotate(geom = "text", x = 0.1*max_value,  y = 0.9*max_value,angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Gain of spatially restricted species", ]$Number.sites,")"),  size =5) +
  annotate(geom = "text", x =0.9*min_value, y = 0.1*max_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Widespread replaced by spatially restricted species", ]$Number.sites,")"),  size =5) +
  annotate(geom = "text", x = 0.9*min_value, y = 0.1*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of widespread species", ]$Number.sites,")"),   size =5) +
  
  geom_point(data=eff.nut.rich.yr, aes(alpha, gamma), size=2.5, alpha=0.3)+
  geom_errorbar(data=eff.nut.overall.yr, aes(x=Estimate_alpha, y=Estimate_gamma, xmin=ci_lower_alpha, xmax=ci_upper_alpha), width=0.0001, color="black", linewidth=1.5, alpha=0.7) +
  geom_errorbar(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, ymin=ci_lower_gamma, ymax=ci_upper_gamma), width=0.0001, color="black", linewidth=1,  alpha=0.7) +
  geom_segment(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.positive * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.positive * pi / 180)),  lineend = "square", color="black", linewidth=1.5, alpha=0.7) +
  geom_segment(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.negative * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.negative * pi / 180)),   lineend = "square", color="black", linewidth=1.5, alpha=0.7) +
  geom_point(data=eff.nut.overall.yr,  aes(x=Estimate_alpha, y=Estimate_gamma), size=6, stroke = 1, pch=21, alpha=0.7, color="black")+ labs(title = paste0(nrow(eff.nut.rich.yr), " sites with ", num.block," blocks"))
leg.space_1<-plot_grid(spacer, leg, nrow = 1, rel_widths = c(3, 97)) # from the file function to draw background for main figures.R
(ppp4 <-plot_grid(ppp2 + theme(plot.margin = margin(0, 5, 5, 5)) , 
                  leg.space_1, 
                  nrow=2, align = "v",  rel_heights = c(8, 1))) 
ggsave(ppp4, width = 18, height = 19, dpi=600, unit="cm",  file=paste0("change in species richness across scales at year ",yr, " for sites with ", num.block," blocks.png"))

##########################################################################################
###### change in species richness for native and non-native across scales in year 4
###########################################################################################

list_plots <- vector('list', 20)
for(fg in c("FORB", "GRAMINOID",  "LEGUME", "WOODY", "NAT", "INT")){
  # fg <-"NAT" ; yr<-4
  eff.nut.rich.yr.groups<-eff.nut%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c(fg)) %>% filter(!is.na(alpha))%>%filter(!is.na(gamma))
  range.values<-range(range(eff.nut.rich.yr.groups$alpha), range(eff.nut.rich.yr.groups$gamma))
  
  (beta.group<-table(eff.nut.rich.yr.groups$cat.beta)%>%data.frame())
  names(beta.group)<-c("Groups.beta", "Number.sites")
  (process.group<-table(eff.nut.rich.yr.groups$cat.process)%>%data.frame())
  names(process.group)<-c("Groups.process", "Number.sites")
  
  eff.nut.overall.yr.groups<-summary.diversity.across.sites1%>%filter(year_trt%in% c(yr))%>%filter(q%in% c(0))%>%filter(functional_group%in% c(fg))%>%pivot_wider(names_from = "name", values_from = c("Estimate", "ci_lower", "ci_upper"))%>%
    mutate(dif.positive=ci_upper_beta - Estimate_beta)
  dif.positive<-eff.nut.overall.yr.groups$dif.positive;
  # min_value<- 1.2* range.values[1];  max_value<-  1.2*range.values[2]
  min_value<- ifelse (range.values[1] > - 0.5, -0.5, range.values[1]);  max_value<- ifelse (range.values[2] >0.6, range.values[2], 0.6)
  
  pp.temp<-plot.six.scenarios(min_value, max_value)
  (pp.temp1<-pp.temp+ 
      annotate(geom = "text", x = min_value/1.7, y = max_value/1.7, angle = 45, label = "Differentiation", color="white", size =3.6, fontface = "bold") +  
      annotate(geom = "text", x = min_value/1.9, y = max_value/1.9, angle = 45, label = "Increase in ∆β (LRR)", color="white", size =2.6) +  
      annotate(geom = "text", x = min_value/2.15, y = max_value/2.15, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Differentiation", ]$Number.sites,")"), color="white", size =3) +  
      annotate(geom = "text", x = max_value/1.7, y = min_value/1.7, angle = 45, label = "Homogenization", color="white", size =3.6, fontface = "bold") +
      annotate(geom = "text", x = max_value/1.5, y = min_value/1.5, angle = 45, label = "Decrease in ∆β (LRR)", color="white", size =2.6) +
      annotate(geom = "text", x = max_value/1.35, y = min_value/1.35, angle = 45, label = paste0("(", beta.group[beta.group$Groups.beta=="Homogenization", ]$Number.sites,")"),  color="white", size =3) +
      
      annotate(geom = "text", x = 0.9*max_value, y=0.1*max_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Gain of widespread species", ]$Number.sites,")"), size =3) +
      annotate(geom = "text",  x = 0.9*max_value, y =  0.1*min_value, angle = 0,label = paste0("(", process.group[process.group$Groups.process=="Spatially restricted replaced by widespread species", ]$Number.sites,")"),  size =3) +
      annotate(geom = "text", x = 0.1*min_value, y = 0.9*min_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Loss of spatially restricted species", ]$Number.sites,")"), size =3) +
      
      annotate(geom = "text", x = 0.1*max_value,  y = 0.9*max_value,angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Gain of spatially restricted species", ]$Number.sites,")"),  size =3) +
      annotate(geom = "text", x =0.9*min_value, y = 0.1*max_value, angle = 0, label = paste0("(", process.group[process.group$Groups.process=="Widespread replaced by spatially restricted species", ]$Number.sites,")"),  size =3) +
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
    draw_label( expression(bar("∆α") ~ " (LRR)"), x = 0.5, y = 0, size = 11) + draw_label("∆γ (LRR)", x = 0, y = 0.5, angle = 90, size = 11))
(pp.life.forms4 <-plot_grid(pp.life.forms3 + theme(plot.margin = margin(1, 5, 5, 5)) , 
                            leg.space, 
                            nrow=2, align = "v",  rel_heights = c(8, 1))) 

ggsave(pp.life.forms4, width = 18, height = 20, dpi=600, unit="cm",  file=paste0("effects on change in diversity across scales for life forms at year ",yr, " for sites with ", num.block," blocks.png"))

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
    draw_label( expression(bar("∆α") ~ " (LRR)"), x = 0.5, y = 0, size = 11) + draw_label("∆γ (LRR)", x = 0, y = 0.5, angle = 90, size = 11))
(pp.origins4 <-plot_grid(pp.origins3 + theme(plot.margin = margin(1, 5, 5, 5)) , 
                         leg.space, 
                         nrow=2, align = "v",  rel_heights = c(8, 2))) 
ggsave(pp.origins4, width = 18, height = 10, dpi=600, unit="cm",  file=paste0("effects on change in diversity across scales for native and non-native species at year ",yr, " for sites with ", num.block," blocks.png"))

##########################################################################################
###### source data for figures 
###########################################################################################
colnames(summary.model.fit1)
eff.nut.l<-eff.nut%>%pivot_longer(cols = c("alpha", "gamma", "beta"), names_to = "scale0", values_to = "raw data")%>%mutate(scale=paste0("dif.", scale0))
eff.nut.l1<-eff.nut.l%>% mutate(Diversity.cross.scales=case_when(scale=="dif.alpha"~"∆α", 
                                                                 scale=="dif.gamma"~"∆γ",
                                                                 scale=="dif.beta"~"∆β"), 
                                Diversity.metrics=case_when(q==0 ~"Species richness", 
                                                            q==1 ~"Shannon diversity",
                                                            q==2 ~"Simpson diversity"), 
                                Species.groups=case_when(functional_group=="all" ~ "All species", 
                                                         functional_group=="NAT" ~ "Native", 
                                                         functional_group=="INT" ~ "Non-native", TRUE~functional_group),
                                Values=`raw data`, Data.type="Individual sites", `l-95% CI`=999, `u-95% CI`=999, N.sites=999)%>%
  ungroup()%>%select("Diversity.metrics", "Species.groups", "site_code", "year_trt",  "Diversity.cross.scales", "Data.type", "Values", "l-95% CI", "u-95% CI", "cat.beta", "cat.process" , "N.sites")

summary.model.fit2<-summary.model.fit1%>%mutate(site_code="",  Values=Estimate, Data.type="Estimate across sites", cat.beta="", cat.process="")%>%
  select("Diversity.metrics", "Species.groups", "site_code", "year_trt",  "Diversity.cross.scales", "Data.type", "Values", "l-95% CI", "u-95% CI", "cat.beta", "cat.process" , "N.sites")
# add all data together 
sites.and.estimate<-eff.nut.l1%>%bind_rows(summary.model.fit2) %>% mutate_all(~ replace(., . == 999, ""))%>%arrange(Diversity.metrics, Species.groups, Data.type, Diversity.cross.scales)
colnames(sites.and.estimate)<-c("Diversity metrics", "Species groups", "site_code", "year_trt",   "Diversity across scales", "Data type",  "Values", "l-95% CI" , "u-95% CI",  "Categories in ∆β", "Categories in six scenarios", "Number of sites")
Fig1_data_source<-sites.and.estimate%>%filter(`Diversity metrics` %in% c("Species richness"))%>%filter(`Species groups` %in% c("All species"))%>%filter(year_trt%in% c(4))
Fig2_data_source<-sites.and.estimate%>%filter(`Diversity metrics` %in% c("Species richness"))%>%filter(`Species groups` %in% c("Native",  "Non-native"))%>%filter(year_trt%in% c(4))
Fig3_data_source<-sites.and.estimate%>%filter(`Diversity metrics` %in% c("Species richness"))%>%filter(`Species groups` %in% c("FORB", "GRAMINOID", "LEGUME", "WOODY"))%>%filter(year_trt%in% c(4))

write.csv(Fig1_data_source, file=paste0("SFig1_data_source ", num.block, " blocks.csv"))
write.csv(Fig2_data_source, file=paste0("SFig2_data_source ", num.block, " blocks.csv"))
write.csv(Fig3_data_source, file=paste0("SFig3_data_source ", num.block, " blocks.csv"))

# the end 