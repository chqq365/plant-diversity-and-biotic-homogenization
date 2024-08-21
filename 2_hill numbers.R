
rm(list=ls())
## open the library
library(tidyverse);library(ggplot2);library(cowplot);library(brms);library(tidybayes)
set.seed(123)
## set up the work directory 
dir.data<-"C:/Users/chqq3/work/NutNet data/"
dir.data1<-"C:/Users/chqq3/work/homogenization/raw data/"
dir.graphs<-"C:/Users/chqq3/work/homogenization/graphs2/"
setwd(dir.graphs)

###########################################################################################
##predicted intercepts in alpha, gamma, and beta diversity
###########################################################################################

summary.intercepts.individual.sites<-c(); summary.intercepts.across.sites<-c(); summary.alpha.gamma.individual.sites<-c(); summary.model.fit<-c()

for(i in c(0, 1, 2)){
  for(fg in c("all")){
    # i<-0; fg<-"all"
    load(paste0("model for alpha diversity for ", fg , " species with q of ", i, " for sites with 4-year nutrient addition.Rdata"))    
    #plot(mod.alpha)
    #pp_check(mod.alpha)
    load(paste0("model for gamma diversity for ", fg , " species with q of ", i, " for sites with 4-year nutrient addition.Rdata"))    
    #plot(mod.gamma)
    #pp_check(mod.gamma)
    
    ## alpha scale  
    # fixed effect coefficients 
    alpha.diversity_fixef <- fixef(mod.alpha)
    # coefficients for site-level (random) effects
    alpha.diversity_coef <- coef(mod.alpha)
    alpha.diversity_coef2 <-  bind_cols(site_code = rownames(alpha.diversity_coef$site_code[,,'Intercept']),
                                        alpha.diversity_coef$site_code[,,'Intercept'] %>% 
                                          as_tibble() %>% mutate(Intercept = Estimate,
                                                                 Intercept_lower = Q2.5,
                                                                 Intercept_upper = Q97.5) %>% 
                                          select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                        
                                        alpha.diversity_coef$site_code[,,'trtNPK'] %>% 
                                          as_tibble() %>% 
                                          mutate(TE = Estimate,
                                                 TE_lower = Q2.5,
                                                 TE_upper = Q97.5) %>% 
                                          select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) 
    # gamma scale
    # fixed effect coefficients 
    gamma.diversity_fixef <- fixef(mod.gamma)
    # coefficients for site-level (random) effects
    gamma.diversity_coef <- coef(mod.gamma)
    gamma.diversity_coef2 <-  bind_cols(site_code = rownames(gamma.diversity_coef$site_code[,,'Intercept']),
                                        gamma.diversity_coef$site_code[,,'Intercept'] %>% 
                                          as_tibble() %>% mutate(Intercept = Estimate,
                                                                 Intercept_lower = Q2.5,
                                                                 Intercept_upper = Q97.5) %>% 
                                          select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                        
                                        gamma.diversity_coef$site_code[,,'trtNPK'] %>% 
                                          as_tibble() %>% 
                                          mutate(TE = Estimate,
                                                 TE_lower = Q2.5,
                                                 TE_upper = Q97.5) %>% 
                                          select(-Estimate, -Est.Error, -Q2.5, -Q97.5) )
    
    
    # add intercepts for alpha and gamma diversity for individual sites
    #colnames(alpha.diversity_coef2);colnames(gamma.diversity_coef2)
    alpha.gamma.intercepts.individual.sites0<-alpha.diversity_coef2%>%mutate(scale="alpha")%>%
      bind_rows(gamma.diversity_coef2%>%mutate(scale="gamma"))%>%mutate(q=i, functional_group=fg)
    alpha.gamma.intercepts.individual.sites<-alpha.gamma.intercepts.individual.sites0%>%select(scale, site_code, TE)%>%
      pivot_wider(names_from = "scale", values_from=c("TE"))
    
    # add intercepts for alpha and gamma diversity across all sites 
    #colnames(alpha.diversity_fixef); colnames(gamma.diversity_fixef)
    alpha.gamma.intercepts0<-alpha.diversity_fixef%>%as.data.frame()%>%mutate(terms=rownames(.), scale="alpha")%>%
      bind_rows(gamma.diversity_fixef%>%as.data.frame()%>%mutate(terms=rownames(.), scale="gamma"))%>%   
      select(scale, terms, Estimate, Q2.5, Q97.5)%>%
      mutate(q=i, functional_group=fg)
    
    ############################calculate change in beta diversity#############################
    
    # change in beta diversity per year calculated as the distance from 1:1 line (left = homogenization, right = differentiation) of 
    # alpha and gamma-scale slope posterior distributions at the global (population) level
    beta.across.sites0 <- bind_cols(gather_draws(mod.alpha, `b_trtNPK`) %>%
                                      ungroup() %>%
                                      select(x = .value),
                                    gather_draws(mod.gamma, `b_trtNPK`) %>%
                                      ungroup() %>%
                                      select(y = .value)) %>%
      mutate(.value = y - x)%>%
      summarise(Estimate = mean(.value), Q2.5 = quantile(.value, 0.025), Q97.5 = quantile(.value, 0.975))%>%
      mutate(scale="beta", terms="", q=i, functional_group=fg)%>%select(scale, terms, Estimate, Q2.5,  Q97.5, q, functional_group)
    
    # change in beta at individual sites
    # 1000 draws of 𝛼- and 𝛾-scale slope posterior distributions 
    alpha.gamma.global.1000<-bind_cols(gather_draws(mod.alpha, `b_trtNPK`,  ndraws = 1000, seed = 123) %>%
                                         ungroup() %>%
                                         select(alpha.global = .value),
                                       gather_draws(mod.gamma, `b_trtNPK`, ndraws = 1000, seed = 123) %>%
                                         ungroup() %>%
                                         select(gamma.global = .value))
    
    sites.yr<-mod.alpha$data%>%select(site_code)%>%distinct()
    
    diversity.individual.sites0<-c()
    for(s in unique(sites.yr$site_code)){
      # s<- "smith.us" 
      # check the answer for question "How to dynamically define variable name in gather_draws?"
      # https://github.com/mjskay/tidybayes/issues/166
      
      var<-sym(paste0("r_site_code[", s , ",trtNPK]"))
      diversity.individual.sites.temp<-bind_cols(gather_draws(mod.alpha, !!var, ndraws = 1000, seed = 123) %>% ungroup() %>%select(x = .value),
                                                 gather_draws(mod.gamma, !!var, ndraws = 1000, seed = 123)  %>% ungroup() %>%select(y = .value),
                                                 alpha.gamma.global.1000)%>%
        mutate(true.alpha=x+alpha.global, true.gamma=y+gamma.global)%>%
        mutate(change.beta =true.gamma - true.alpha)%>%
        summarise_at(c("change.beta", "true.alpha", "true.gamma"), list(mean=mean, Q2.5 = ~quantile(., probs=0.025), Q97.5 = ~quantile(., prob=0.975)))%>%mutate(site_code=s)
      
      diversity.individual.sites0<-diversity.individual.sites0%>%rbind(diversity.individual.sites.temp)  
    }
    
    diversity.individual.sites00<-diversity.individual.sites0%>%mutate(q=i, functional_group=fg)  
    
    # save the predicted alpha, gamma, and beta slope at site level  
    summary.intercepts.individual.sites<-summary.intercepts.individual.sites%>%bind_rows(diversity.individual.sites00)
    # save the predicted intercepts and intercepts at global and site level  
    summary.intercepts.across.sites<-summary.intercepts.across.sites%>%bind_rows(alpha.gamma.intercepts0)%>%bind_rows(beta.across.sites0)
    summary.alpha.gamma.individual.sites<-summary.alpha.gamma.individual.sites%>%bind_rows(alpha.gamma.intercepts.individual.sites0)
    
    # check model fit 
    t.alpha<-summary(mod.alpha)
    t.alpha1<-t.alpha$fixed%>%mutate(terms=rownames(.), scale="alpha")
    t.gamma<-summary(mod.gamma)
    t.gamma1<-t.gamma$fixed%>%mutate(terms=rownames(.), scale="gamma")
    t.alpha.gamma<-t.alpha1%>%bind_rows(t.gamma1)%>%mutate(q=i, functional_group=fg)
    N.sites<-length(unique(alpha.gamma.intercepts.individual.sites$site_code))
    summary.model.fit<-summary.model.fit%>%bind_rows(t.alpha.gamma%>%mutate(N.sites=N.sites))  
    
  }
}

summary.intercepts.individual.sites1<-summary.intercepts.individual.sites%>%mutate_if(is.numeric, round, digits=2)
colnames(summary.intercepts.individual.sites1)
length(unique(summary.intercepts.individual.sites1$site_code))
# determine category counts that incorporate uncertainty
categorize.beta<-summary.intercepts.individual.sites1%>%
  mutate(cat.beta=case_when((change.beta_mean<0 & change.beta_Q2.5<0 & change.beta_Q97.5<0)~"Homogenization", 
                            (change.beta_mean>0 & change.beta_Q2.5>0 & change.beta_Q97.5>0)~"Differentiation", 
                            TRUE~"Little change in beta diversity"))%>%
  mutate(cat.process.CI=case_when((true.alpha_mean>true.gamma_mean & true.alpha_mean>0 &true.gamma_mean>0 & ((true.alpha_Q2.5 > 0 & true.alpha_Q97.5 > 0) & (true.gamma_Q2.5 > 0 & true.gamma_Q97.5 > 0)))~"Gain widespread species",
                                  (true.alpha_mean>0 &true.gamma_mean<0 & ((true.alpha_Q2.5 > 0 & true.alpha_Q97.5 > 0) & (true.gamma_Q2.5 < 0 & true.gamma_Q97.5 < 0)))~"Localized species replaced by widespread species",
                                  (true.alpha_mean>true.gamma_mean & true.alpha_mean<0 &true.gamma_mean<0 & ((true.alpha_Q2.5 < 0 & true.alpha_Q97.5 < 0) & (true.gamma_Q2.5 < 0 & true.gamma_Q97.5 < 0)))~"Loss of localized species",
                                  (true.alpha_mean<true.gamma_mean & true.alpha_mean<0 &true.gamma_mean<0 & ((true.alpha_Q2.5 < 0 & true.alpha_Q97.5 < 0) & (true.gamma_Q2.5 < 0 & true.gamma_Q97.5 < 0)))~"Loss of widespread species",
                                  (true.alpha_mean<0 &true.gamma_mean>0 & ((true.alpha_Q2.5 < 0 & true.alpha_Q97.5 < 0) & (true.gamma_Q2.5 > 0 & true.gamma_Q97.5 > 0)))~"Widespread species replaced by localized species",
                                  (true.alpha_mean<true.gamma_mean & true.alpha_mean>0 &true.gamma_mean>0 & ((true.alpha_Q2.5 > 0 & true.alpha_Q97.5 > 0) & (true.gamma_Q2.5 > 0 & true.gamma_Q97.5 > 0)))~"Gain localized species" ,
                                   TRUE~"Other situations"))%>%
  mutate(cat.process=case_when((true.alpha_mean>true.gamma_mean & true.alpha_mean>0 &true.gamma_mean>0 )~"Gain widespread species",
                               (true.alpha_mean>0 &true.gamma_mean<0 )~"Localized species replaced by widespread species",
                               (true.alpha_mean>true.gamma_mean & true.alpha_mean<0 &true.gamma_mean<0)~"Loss of localized species",
                               (true.alpha_mean<true.gamma_mean & true.alpha_mean<0 &true.gamma_mean<0)~"Loss of widespread species",
                               (true.alpha_mean<0 &true.gamma_mean>0)~"Widespread species replaced by localized species",
                               (true.alpha_mean<true.gamma_mean & true.alpha_mean>0 &true.gamma_mean>0 )~"Gain localized species" ,
                               TRUE~"Other situations"))%>%arrange(cat.beta, cat.process.CI)

# save the data 
colnames(categorize.beta)
i<-0; fg<-"all"
categorize.beta1<-categorize.beta%>%select(q, site_code, cat.beta, cat.process, cat.process.CI)%>%distinct()
table(categorize.beta$q, categorize.beta$cat.beta)
che<-as.data.frame(table(categorize.beta$q, categorize.beta$cat.process))%>%arrange(Var1)
table(categorize.beta$q, categorize.beta$cat.process.CI)
table(categorize.beta$q, categorize.beta$cat.process)

colnames(categorize.beta1)<-c("q", "site_code", "Category change in beta diversity", "Category processes_without CI", "Category processes")
# writexl::write_xlsx(categorize.beta1, path="Category change in beta diversity and processes.xlsx")  

# save change in alpha, gamma, and beta diversity and 95% confident intervals 
summary.intercepts.individual.sites2<-summary.intercepts.individual.sites1%>%filter(functional_group==fg)%>%mutate(functional_group=NULL, continent=NULL)%>%
  mutate(Alpha.diversity=paste0(true.alpha_mean, " [", true.alpha_Q2.5, ", ", true.alpha_Q97.5, "]"))%>%
  mutate(Beta.diversity=paste0(change.beta_mean, " [", change.beta_Q2.5, ", ", change.beta_Q97.5, "]"))%>%
  mutate(Gamma.diversity=paste0(true.gamma_mean, " [", true.gamma_Q2.5, ", ", true.gamma_Q97.5, "]"))%>%
  select(q, site_code, Alpha.diversity, Gamma.diversity, Beta.diversity)%>%arrange(q)%>%
  pivot_wider(names_from = "q", values_from = c("Alpha.diversity", "Gamma.diversity", "Beta.diversity"))
# writexl::write_xlsx(summary.intercepts.individual.sites2, path="change in alpha, gamma, and beta diversity for hill number of 0 to 2.xlsx")  
# also save the raw data 
# write.csv(categorize.beta, file="estimated alpha, gamma, and beta diversity for hill number of 0 to 2 at individual sites.csv")

###########################################################################################
################################ make customized figures #############################
###########################################################################################
summary.intercepts.across.sites<-summary.intercepts.across.sites%>%mutate(q1=paste0("hill_", q))%>%mutate_if(is.numeric, round, digits=2)
summary.intercepts.across.sites1<-summary.intercepts.across.sites%>%mutate_if(is.numeric, round, digits=2)%>%mutate(mean.ci= paste0(Estimate, " [", Q2.5, ", ", Q97.5, "]"))%>%
  filter(terms %in% c("trtNPK", ""))%>%arrange(q, functional_group)%>%select(q, scale, mean.ci)
colnames(summary.intercepts.across.sites1)<-c("q", "Spatial scale", "Estimate")
# writexl::write_xlsx(summary.intercepts.across.sites1, path="overall effects of nutrient addition on alpha, gamma, and beta diversity for hill number of 0 to 2.xlsx")  

categorize.beta<-categorize.beta%>%mutate(q1=paste0("hill_", q))
 source("function to draw figures.R")
# for hill number of 0
i<-"hill_0"; fg<-"all"; sc<-"beta"
summary.intercepts.across.sites.w<-summary.intercepts.across.sites%>%filter(terms =="trtNPK")%>%filter(scale!=sc)%>%filter(q1==i & functional_group==fg)%>%
  pivot_wider(names_from = "scale", values_from = c("Estimate", "Q2.5", "Q97.5"))
overall.change.in.beta<- summary.intercepts.across.sites%>%filter(q1==i & functional_group==fg & scale==sc)%>%
  mutate(dif.positive=Q97.5 - Estimate, dif.negative=Estimate - Q2.5)
# site-level data 
summary.intercepts.individual.sites_i_fg<-categorize.beta%>%filter(q1==i & functional_group==fg)
length(unique(summary.intercepts.individual.sites_i_fg$site_code))

# Create a data frame for the diagonal line
df.ticks <- data.frame(x = seq(-3, 3, by = 0.03))%>%mutate(x.position= x + 0.02)
df.ticks$y.position <- -  df.ticks$x + 0.02
df.ticks1<-df.ticks%>%mutate(x1=-round(sqrt(2)*x, 2), x2=ifelse(x1==0, NA, x1))
(pp.hill_0 <- pp.diversity.annotation(data1=summary.intercepts.individual.sites_i_fg, data2=summary.intercepts.across.sites.w, data3=overall.change.in.beta)+
  geom_abline(slope=-1, intercept=0)+
  geom_text(data=df.ticks1, aes(x.position, y.position, label=x2), size = 3, angle=45)) 
ggsave(pp.hill_0, width = 18.4, height = 18.4, dpi=600, unit="cm",  file="effects on change in diversity across scales for hill number 0.png")

# for hill number of 1
i<-"hill_1"; fg<-"all"; sc<-"beta"
summary.intercepts.across.sites.w<-summary.intercepts.across.sites%>%filter(terms =="trtNPK")%>%filter(scale!=sc)%>%filter(q1==i & functional_group==fg)%>%
  pivot_wider(names_from = "scale", values_from = c("Estimate", "Q2.5", "Q97.5"))
overall.change.in.beta<- summary.intercepts.across.sites%>%filter(q1==i & functional_group==fg & scale==sc)%>%
  mutate(dif.positive=Q97.5 - Estimate, dif.negative=Estimate - Q2.5)
# site-level data 
summary.intercepts.individual.sites_i_fg<-categorize.beta%>%filter(q1==i & functional_group==fg)
# position for annotation
label.position<-summary.intercepts.individual.sites_i_fg%>%
  summarise(max.alpha=max(true.alpha_Q97.5), max.gamma=max(true.gamma_Q97.5), min.alpha=min(true.alpha_Q2.5), min.gamma=min(true.gamma_Q2.5))%>%
  rowwise() %>% mutate(min_value = 0.9*min(min.alpha, min.gamma), max_value_0 =  0.8*max(max.alpha, max.gamma), max_value=ifelse(max_value_0 < 0.2, 0.2, max_value_0))

df.ticks <- data.frame(x = seq(-3, 3, by = 0.15))%>%mutate(x.position= x + 0.06)
df.ticks$y.position <- -  df.ticks$x + 0.06
df.ticks1<-df.ticks%>%mutate(x1=-round(sqrt(2)*x, 2), x2=ifelse(x1==0, NA, x1))
(pp.hill_1 <- pp.diversity.simple(data1=summary.intercepts.individual.sites_i_fg, data2=summary.intercepts.across.sites.w, data3=overall.change.in.beta)+
  annotate(geom = "text", x = label.position$min_value/2, y = 0, angle = 45, label = "Differentiation", size = 4, color="grey68") +  
  annotate(geom = "text", x = 0, y = label.position$min_value/2, angle = 45, label = "Homogenization", size = 4, color="grey68") +
  annotate(geom = "text", x = 0.80*label.position$max_value, y=0.25*label.position$max_value, angle = 0, label = "Gain of\n widespread\n species", size = 3, color="grey78") +
  annotate(geom = "text",  x = 0.30*label.position$max_value, y =  0.85*label.position$min_value, angle = 0, label = "Localized\nreplaced by \nwidespread\n species", size = 3, color="grey78") +
  annotate(geom = "text", x = 0.2*label.position$min_value, y = 0.90*label.position$min_value, angle = 0, label = "Loss of\n localized\n species", size = 3, color="grey78") +
  annotate(geom = "text", x = 0.85*label.position$min_value, y = 0.3*label.position$min_value, angle = 0, label = "Loss of\n widespread\n species", size = 3, color="grey78") +
  annotate(geom = "text", x =0.85*label.position$min_value, y = 0.32*label.position$max_value, angle = 0, label = "Widespread\nreplaced by \nlocalized\n species", size = 3, color="grey78") +
  annotate(geom = "text", x = 0.3*label.position$max_value,  y = 0.9*label.position$max_value,  angle = 0, label = "Gain of\n localized\n species", size = 3, color="grey78")+
    geom_abline(slope=-1, intercept=0)+
    geom_text(data=df.ticks1, aes(x.position, y.position, label=x2), size = 2, angle=45)) 

# for simplified figures without 
list_plots <- vector('list', 1)
for(i in c("hill_2")){
  # i<-"hill_0"; fg<-"all";    sc<-"beta"
  summary.intercepts.across.sites.w<-summary.intercepts.across.sites%>%filter(terms =="trtNPK")%>%filter(scale!=sc)%>%filter(q1==i & functional_group==fg)%>%
    pivot_wider(names_from = "scale", values_from = c("Estimate", "Q2.5", "Q97.5"))
  
  overall.change.in.beta<- summary.intercepts.across.sites%>%filter(q1==i & functional_group==fg & scale==sc)%>%
    mutate(dif.positive=Q97.5 - Estimate, dif.negative=Estimate - Q2.5)
  
  # site-level data 
  summary.intercepts.individual.sites_i_fg<-categorize.beta%>%filter(q1==i & functional_group==fg)
  # color for each section 
  # unique(summary.intercepts.individual.sites_i_fg$cat.process)
  # position for labels to show 8 sections and x y limits 
  label.position<-summary.intercepts.individual.sites_i_fg%>%
    summarise(max.alpha=max(true.alpha_Q97.5), max.gamma=max(true.gamma_Q97.5), min.alpha=min(true.alpha_Q2.5), min.gamma=min(true.gamma_Q2.5))%>%
    rowwise() %>% mutate(min_value = min(min.alpha, min.gamma), max_value_0 = max(max.alpha, max.gamma), max_value=ifelse(max_value_0 < 0.2, 0.2, max_value_0))
  
  list_plots[[i]]<- pp.diversity.simple(data1=summary.intercepts.individual.sites_i_fg, data2=summary.intercepts.across.sites.w, data3=overall.change.in.beta)+
    geom_abline(slope=-1, intercept=0)+
    geom_text(data=df.ticks1, aes(x.position, y.position, label=x2), size = 2, angle=45)
}

# combine the graphs 
(pp.div1<-plot_grid(pp.hill_1,
                    list_plots$hill_2,
                    nrow=1, ncol=2, 
                    vjust = -0.01, 
                    labels = c("A (Shannon diversity)", "B (Simpson diversity)"), label_fontface = "bold", label_size = 11))
(pp.div3<-pp.div1+ theme(plot.margin = margin(l = 10, b = 10, t=15))+
    draw_label("Effects of nutrient addition on alpha diversity (LRR)", x = 0.5, y = 0, size = 11) + draw_label("Effects of nutrient addition on gamma diversity (LRR)", x = 0, y = 0.5, angle = 90, size = 11))#ggsave(pp.div3, width = 18.4, height = 12, dpi=600, unit="cm",  file="effects on change in diversity across scales.png")
ggsave(pp.div3, width = 18.4, height = 9.2, dpi=600, unit="cm",  file="effects on change in diversity across scales for hill number 1 to 2.png")

###########################################################################################
####### sort out summary tables and check model fit 
###########################################################################################
i<-0; fg<-"all"
load(paste0("model for alpha diversity for ", fg , " species with q of ", i, " for sites with 4-year nutrient addition.Rdata"))    
load(paste0("model for gamma diversity for ", fg , " species with q of ", i, " for sites with 4-year nutrient addition.Rdata"))    
# check model fitting 
check.pred.alpha<-pp_check(mod.alpha)+ theme_cowplot()+panel_border()+
  labs(x="Change in alpha diversity", y="Density")
check.pred.gamma<-pp_check(mod.gamma)+ theme_cowplot()+panel_border()+
  labs(x="Change in gamma diversity", y="Density")
legend<-get_legend(check.pred.alpha)
(check.fit<-plot_grid(check.pred.alpha +theme(legend.position = "none"), 
                      check.pred.gamma+theme(legend.position = "none"), 
                      labels = "AUTO"))
(check.fit1<-plot_grid(check.fit, 
                       legend, rel_widths = c(9, 1)))
ggsave(check.fit1, width = 21, height=14.8, unit="cm",  file=paste0("check model fitting for alpha and gamma diversity based on Q of ", i, " for ",fg, " species.png"))# summary table 


###########################################################################################
### relationships between change in  diversity and site biotic and abiotic factors 
###########################################################################################
drought<-read.csv(file="site drought index.csv")
site.productivity<-read.csv(file="site productivity.csv")
site.species.pool<-read.csv(file="site species pool.csv")
biomass.removed.by.herbivore<-read.csv(file="site herbivory as biomass removed by herbivores.csv")
site.herbivores<-read.csv(file="site herbivores.csv")

# select years 
d7<-read.csv("cover data for control and NPK at 72 NutNet sites.csv")
year.used.for.climate<-d7%>%filter(year_trt==4)%>%select(site_code, year, year_trt)%>%distinct()
drought1<-drought%>%merge(year.used.for.climate, by=c("site_code"))%>%filter(year.x <= year.y & year.x >= (year.y-4) )%>%
  group_by(site_code)%>%summarise(drought.index1=mean(drought.index))
  
site.block.select<-d7%>%filter(year_trt==4)%>%mutate(site.block=paste(site_code, block, sep="_"))%>%select(site.block)%>%distinct()

site.productivity1<-site.productivity%>%mutate(site.block=paste(site_code, block, sep="_"))%>%
  filter(site.block %in% site.block.select$site.block)%>% filter(trt %in% c("Control"))%>%filter(year_trt<=4)
# check years included for each site 
check.year.biomass<-site.productivity1%>%ungroup()%>%select(site_code, year_trt)%>%distinct()%>%group_by(site_code)%>%summarise(N.years=length(year_trt))
table(check.year.biomass$N.years)# majority of the sites have 5 years' data 
site.productivity2<-site.productivity1%>%group_by(site_code)%>%summarise(site.productivity=mean(live_mass))

biomass.removed.by.herbivore1<-biomass.removed.by.herbivore%>%mutate(site.block=paste(site_code, block, sep="_"))%>%
  filter(site.block %in% site.block.select$site.block)%>% filter(year_trt<=4)
# check years included for each site 
check.year.biomass.removed<-biomass.removed.by.herbivore1%>%ungroup()%>%select(site_code, year_trt)%>%distinct()%>%group_by(site_code)%>%summarise(N.years=length(year_trt))
table(check.year.biomass.removed$N.years)# majority of the sites have 5 years' data 
biomass.removed.by.herbivore2<-biomass.removed.by.herbivore1%>%group_by(site_code)%>%summarise(biomass.removed.by.herbivore=mean(herbivory.biomass.removed))
table(biomass.removed.by.herbivore2$biomass.removed.by.herbivore >0)
# at 45 sites biomass removed is >0, at 17 sites, biomass removed is < 0, not sure this is a good proxy for herbivore intensity 
site.species.pool1<-site.species.pool%>%filter(site_code %in% year.used.for.climate$site_code)%>% filter(year_trt<=4)
# check years included for each site 
check.year.richness<-site.species.pool1%>%ungroup()%>%select(site_code, year_trt)%>%distinct()%>%group_by(site_code)%>%summarise(N.years=length(year_trt))
table(check.year.richness$N.years)# majority of the sites have 5 years' data 
# it is very important that each site should have the same number of years. Sampling more times may find more species 
check.year.richness.4.years<-check.year.richness%>%filter(N.years <5)
# for those sites that only sample for 4 years, add year 5 data if possible 
site.species.pool_sub.sites<-site.species.pool%>%filter(site_code %in% check.year.richness.4.years$site_code)%>% filter(year_trt<=5)
check.year.richness_sub<-site.species.pool_sub.sites%>%ungroup()%>%select(site_code, year_trt)%>%distinct()%>%group_by(site_code)%>%summarise(N.years=length(year_trt))%>% filter(N.years==5)

site.species.pool2<-site.species.pool1%>%filter(!site_code %in%check.year.richness.4.years$site_code)%>%
  bind_rows(site.species.pool1%>%filter(site_code %in%check.year.richness_sub$site_code)%>% filter(year_trt<=5))%>%
  group_by(site_code)%>%summarise(site.species.pool=mean(HillDiv))

# add all data together 
change.in.div.cov<-categorize.beta%>%filter(q==0 )%>%select(site_code, functional_group, true.alpha_mean, true.gamma_mean, change.beta_mean)%>%
  pivot_longer(cols = c("true.alpha_mean",   "true.gamma_mean", "change.beta_mean"), names_to = "diversity.facet", values_to = "diversity.value")%>%
  left_join(drought1, by=c("site_code"))%>%
  left_join(site.species.pool2, by=c("site_code"))%>%
  left_join(site.productivity2, by=c("site_code"))%>%
  left_join(biomass.removed.by.herbivore2, by=c("site_code"))%>%
  left_join(site.herbivores%>%mutate(X=NULL), by=c("site_code"))
# visualize the relationships 
colnames(change.in.div.cov) 
change.in.div.all<-change.in.div.cov%>% mutate( diversity.facet1=case_when(grepl("alpha", diversity.facet) ~ "Change in alpha diversity",
                                                                          grepl("change.beta_mean", diversity.facet) ~ "Change in beta diversity", 
                                                                          TRUE~"Change in gamma diversity"))
# use color for different intercepts from conceptual figure
cat.col<-categorize.beta%>%filter(q==0 & functional_group=="all")%>%select(site_code, functional_group, cat.process)
unique(cat.col$cat.process)
# themes needed 
concept_colour = c(  "Gain widespread species" = '#D9D956',
                     "Localized species replaced by widespread species"  = '#E9AE27',
                     "Loss of widespread species" = '#155F49',
                     "Loss of localized species" = '#D17538',
                     "Widespread species replaced by localized species"  ='#2CABA4', 
                     "Gain localized species" = '#61CDE0',
                     # for points that fall on the boundary
                     "Other situations" = '#f0f0f0')

(pp.relation.raw<-change.in.div.all%>%
    select(site_code, diversity.facet, diversity.facet1, diversity.value, site.species.pool, 
           site.productivity, drought.index1, biomass.removed.by.herbivore, herb.index, herb.rich, grazer.biomass)%>%
    pivot_longer(cols = 5:ncol(.))%>%
    merge(cat.col, by = c("site_code")) %>%mutate(grp=paste0(diversity.facet1, name, set="_"))%>%
    ggplot(aes(value, diversity.value, fill = cat.process, shape = diversity.facet1, group=grp, linetype=diversity.facet1)) + theme_cowplot()+panel_border()+
    geom_point(size=3, alpha=0.3)+
    geom_smooth(color="black")+
    facet_wrap(~name, scales="free", nrow = 3)+
    scale_shape_manual(values = c(21, 24,  23)) +
    scale_fill_manual(values = concept_colour) +
    scale_linetype_manual(values=c("solid", "dashed", "dotted"))+
    guides(fill = "none", alpha = guide_legend(nrow = 2), linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    theme(legend.position = "top", strip.placement = "outside", strip.background = element_blank()) +
    labs(x="Values for covariates", y="Effects of nutrient addition on change in diversity", color=NULL, shape = NULL, fill = NULL, alpha = NULL, linetype = NULL))
# ggsave(pp.relation.raw, unit="cm", width = 21, height=14.9, file="relationships between change in beta diversity and site biotic and abiotic factors.png")

###########################################################################################
# Bivariate analysis for four site covariates that often used in previous literature 
change.in.div.all_l<-change.in.div.all%>%select(site_code, functional_group, diversity.facet1, diversity.value, site.species.pool, site.productivity, biomass.removed.by.herbivore, herb.index, drought.index1)%>%
  pivot_longer(cols = c("site.species.pool", "site.productivity", "herb.index", "biomass.removed.by.herbivore", "drought.index1"))%>% 
  mutate(name1=case_when(name=="site.species.pool" ~"Site species pool", 
                         name=="site.productivity" ~"Site productivity", 
                         name=="herb.index" ~"Grazing intensity", 
                         name=="biomass.removed.by.herbivore" ~"Herbivory", 
                         name=="drought.index1" ~"Drought index"))%>%mutate(combi.id=paste(functional_group, diversity.facet1, name1, sep="_"))              
# 
summary.bivariate.estimate.coef<-c(); summary.residuals<-c()
for(df in unique(change.in.div.all_l$combi.id)){
  #  df<-"all_Change in beta diversity_Site productivity" 
  data.temp<-change.in.div.all_l%>%filter(combi.id==df)%>%na.omit()%>%filter_all(all_vars(!is.infinite(.)))
  #hist(data.temp$diversity.value)
  n.sites<-length(unique(data.temp$site_code))
  # model 
  Bivariate.mod <- brm( diversity.value ~ value, data=data.temp, iter=3000, warmup = 1000, cores = 6)
  # coefficients
  t.temp<-summary(Bivariate.mod)
  t.temp1<-t.temp$fixed
  summary.bivariate.estimate.coef<-summary.bivariate.estimate.coef%>%bind_rows(t.temp1%>%as.data.frame()%>%mutate(terms=row.names(.), combi.id=df, n.sites=n.sites))
  # residuals
  resid<-as.data.frame(residuals(Bivariate.mod))%>%bind_cols(data.temp%>%select(site_code, value)%>%distinct())%>%mutate(combi.id=df)
  summary.residuals<-summary.residuals%>%bind_rows(resid)
  # pp_check(Bivariate.mod)+ theme_cowplot(font_size = 18)+panel_border()+labs(x=paste0(df), y="Density")
  # ggsave(check.pred, file=paste0("check model predict for", combi.id , ".PNG"))
}

colnames(summary.bivariate.estimate.coef)
summary.bivariate.estimate.coef1<-summary.bivariate.estimate.coef%>%mutate_at(vars(c(1:7)), round, digits=4)%>%
  mutate(variable=ifelse(grepl("Intercept", terms), "Intercept", "Slope"))%>%dplyr::rename(Q2.5='l-95% CI', Q97.5='u-95% CI')%>%
  mutate(sig=ifelse(((Q2.5 > 0 & Q97.5>0)|(Q2.5<0 & Q97.5< 0)), "Significant", "Non-significant"))%>%
  merge(summary.bivariate.estimate.coef%>%filter(terms=="Intercept")%>%select(combi.id, Estimate), by=c("combi.id"))%>%
  filter(terms!="Intercept")%>% dplyr::rename(Intercept=Estimate.y, Slope=Estimate.x)%>%mutate(terms=NULL, variable=NULL)%>%
  merge(change.in.div.all_l%>%select(diversity.facet1, name1, combi.id)%>%distinct(), by=c("combi.id"))%>%
  filter(!name1 %in% c("Herbivory"))
# relevel 
table(summary.bivariate.estimate.coef1$sig)
summary.bivariate.estimate.coef1$sig<-factor(summary.bivariate.estimate.coef1$sig, levels = c("Significant", "Non-significant"))
unique(summary.bivariate.estimate.coef1$diversity.facet1 )
summary.bivariate.estimate.coef1$diversity.facet1 <- factor(summary.bivariate.estimate.coef1$diversity.facet1, levels=c("Change in alpha diversity", "Change in gamma diversity", "Change in beta diversity" , "Change in beta_C diversity"))

alpha.crosswalk <- c("Significant" = 0.9, "Non-significant" = 0.1)
change.in.div.all_l$diversity.facet1 <- factor(change.in.div.all_l$diversity.facet1, levels=c("Change in alpha diversity", "Change in gamma diversity", "Change in beta diversity" ))

(pp.relation.pred_1 <- change.in.div.all_l%>% merge(cat.col%>%mutate(functional_group=NULL), by = c("site_code")) %>%filter(!name1 %in% c("Herbivory"))%>%
    ggplot(aes(value, diversity.value, fill = cat.process)) +
    theme_cowplot() +    panel_border() +
    geom_point(size = 3, alpha = 0.5, pch=21) +
    geom_abline(data = summary.bivariate.estimate.coef1, aes(slope = Slope, intercept = Intercept, alpha = sig), linewidth = 1) +
    geom_hline(yintercept = 0, linetype="dotted")+
    facet_grid(diversity.facet1~name1, scales = "free", switch = "both") +
    scale_alpha_manual(values = alpha.crosswalk) +
    # scale_shape_manual(values = c(21, 24,  23)) +
    scale_fill_manual(values = concept_colour) +
    scale_linetype_manual(values=c("solid", "dashed", "dotted"))+
    guides(fill = "none", alpha = "none", shape = "none") +
    theme(legend.position = "top", strip.placement = "outside", strip.background = element_blank()) +
    labs(x = NULL, y = NULL, shape = NULL, fill = NULL, alpha = NULL, linetype = NULL, shape=NULL))
ggsave(pp.relation.pred_1, unit="cm", width = 20, height=20, file="predicted bivariate relationships between change in diversity separately and site biotic and abiotic factors.png")


# the end 