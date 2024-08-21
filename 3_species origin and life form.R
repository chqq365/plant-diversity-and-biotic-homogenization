
rm(list=ls())
## open the library
library(tidyverse);library(ggplot2);library(cowplot);library(brms);library(tidybayes)
set.seed(123)
## set up the work directory 
dir.data<-"C:/Users/chqq3/work/homogenization/raw data/"
dir.graphs<-"C:/Users/chqq3/work/homogenization/graphs2/"
setwd(dir.graphs)

###########################################################################################
##################predicted intercepts in alpha, gamma, and beta diversity#####################
###########################################################################################
summary.intercepts.individual.sites<-c(); summary.intercepts.across.sites<-c(); summary.alpha.gamma.individual.sites<-c(); summary.model.fit<-c()
for(i in c(0)){
  for(fg in c("FORB", "GRAMINOID",  "LEGUME",  "WOODY", "NAT", "INT")){
    # i<-0; fg<-"FORB"
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
                                          select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) 
    
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
    t.alpha.gamma<-t.alpha1%>%bind_rows(t.gamma1)%>%mutate(q=i, functional_group=fg)%>%mutate_at(vars(c(1:7)), round, digits=3)
    N.sites<-length(unique(alpha.gamma.intercepts.individual.sites$site_code))
    summary.model.fit<-summary.model.fit%>%bind_rows(t.alpha.gamma%>%mutate(N.sites=N.sites))  
  }
}
# write.csv(summary.model.fit, file="95% confidence intervals for intercepts and slopes and rhat for funtional groups.csv")  

summary.intercepts.individual.sites1<-summary.intercepts.individual.sites%>%mutate_if(is.numeric, round, digits=3)
colnames(summary.intercepts.individual.sites1)
table(summary.intercepts.individual.sites1$functional_group)

# determine category counts that incorporate uncertainty
categorize.beta<-summary.intercepts.individual.sites1%>%
  mutate(cat.beta=case_when((change.beta_mean<0 & change.beta_Q2.5<0 & change.beta_Q97.5<0)~"homogenization", 
                            (change.beta_mean>0 & change.beta_Q2.5>0 & change.beta_Q97.5>0)~"differentiation", 
                            TRUE~"no change in beta"))%>%
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

# save change in alpha, gamma, and beta diversity and 95% confident intervals 
summary.intercepts.individual.sites2<-summary.intercepts.individual.sites1%>%filter(q %in% c(0))%>%
  mutate(continent=NULL)%>%
  mutate(Alpha.diversity=paste0(true.alpha_mean, " [", true.alpha_Q2.5, ", ", true.alpha_Q97.5, "]"))%>%
  mutate(Beta.diversity=paste0(change.beta_mean, " [", change.beta_Q2.5, ", ", change.beta_Q97.5, "]"))%>%
  mutate(Gamma.diversity=paste0(true.gamma_mean, " [", true.gamma_Q2.5, ", ", true.gamma_Q97.5, "]"))%>%
  mutate(grp=ifelse(functional_group %in% c("INT", "NAT"), "G1", "G2"))%>%arrange(q, site_code, grp, functional_group)%>%
  select(q, site_code, functional_group, Alpha.diversity, Gamma.diversity, Beta.diversity)
# also save the raw data 
# write.csv(categorize.beta, file="estimated alpha, gamma, and beta diversity for functional groups at individual sites.csv")

###########################################################################################
################################ make more customized figures #############################
###########################################################################################
colnames(summary.intercepts.across.sites)
summary.intercepts.across.sites1<-summary.intercepts.across.sites%>%mutate_if(is.numeric, round, digits=2)%>%mutate(mean.ci= paste0(Estimate, " [", Q2.5, ", ", Q97.5, "]"))%>%
  filter(terms %in% c("trtNPK", ""))%>%mutate(grp=ifelse(functional_group %in% c("INT", "NAT"), "G1", "G2"))%>%arrange(q, grp, functional_group)%>%select(q, functional_group, scale, mean.ci)
source("function to draw figures.R")
for(i in c(0)){
  # for native and nonnative 
  fg<-"NAT"; sc<-"beta"
  summary.intercepts.across.sites.w<-summary.intercepts.across.sites%>%filter(terms =="trtNPK")%>%filter(scale!=sc)%>%filter(q==i & functional_group==fg)%>%
    pivot_wider(names_from = "scale", values_from = c("Estimate", "Q2.5", "Q97.5"))
  overall.change.in.beta<- summary.intercepts.across.sites%>%filter(q==i & functional_group==fg & scale==sc)%>%
    mutate(dif.positive=Q97.5 - Estimate, dif.negative=Estimate - Q2.5)
  # site-level data 
  summary.intercepts.individual.sites_i_fg<-categorize.beta%>%filter(q==i & functional_group==fg)
  # position for annotation
  label.position<-summary.intercepts.individual.sites_i_fg%>%
    summarise(max.alpha=max(true.alpha_Q97.5), max.gamma=max(true.gamma_Q97.5), min.alpha=min(true.alpha_Q2.5), min.gamma=min(true.gamma_Q2.5))%>%
    rowwise() %>% mutate(min_value = min(min.alpha, min.gamma), max_value_0 = max(max.alpha, max.gamma), max_value=ifelse(max_value_0 < 0.2, 0.2, max_value_0))
 
  df.ticks <- data.frame(x = seq(-3, 3, by = 0.1))%>%mutate(x.position= x + 0.05)
  df.ticks$y.position <- -  df.ticks$x + 0.05
  df.ticks1<-df.ticks%>%mutate(x1=-round(sqrt(2)*x, 2), x2=ifelse(x1==0, NA, x1))
  (pp.nat <- pp.diversity.simple(data1=summary.intercepts.individual.sites_i_fg, data2=summary.intercepts.across.sites.w, data3=overall.change.in.beta)+
    annotate(geom = "text", x = label.position$min_value/2, y = 0, angle = 45, label = "Differentiation", size = 4, color="grey68") +  
    annotate(geom = "text", x = 0, y = label.position$min_value/2, angle = 45, label = "Homogenization", size = 4, color="grey68") +
      annotate(geom = "text", x = 0.80*label.position$max_value, y=0.3*label.position$max_value, angle = 0, label = "Gain of\n widespread\n species", size = 3, color="grey78") +
      annotate(geom = "text",  x = 0.35*label.position$max_value, y =  0.85*label.position$min_value, angle = 0, label = "Localized\nreplaced by \nwidespread\n species", size = 3, color="grey78") +
      annotate(geom = "text", x = 0.2*label.position$min_value, y = 0.90*label.position$min_value, angle = 0, label = "Loss of\n localized\n species", size = 3, color="grey78") +
      annotate(geom = "text", x = 0.85*label.position$min_value, y = 0.3*label.position$min_value, angle = 0, label = "Loss of\n widespread\n species", size = 3, color="grey78") +
      annotate(geom = "text", x =0.85*label.position$min_value, y = 0.35*label.position$max_value, angle = 0, label = "Widespread\nreplaced by \nlocalized\n species", size = 3, color="grey78") +
      annotate(geom = "text", x = 0.3*label.position$max_value,  y = 0.9*label.position$max_value,  angle = 0, label = "Gain of\n localized\n species", size = 3, color="grey78")+
      geom_abline(slope=-1, intercept=0)+
    geom_text(data=df.ticks1, aes(x.position, y.position, label=x2), size = 2, angle=45)) 
  
  # for native and nonnative 
  fg<-"INT"; sc<-"beta"
  summary.intercepts.across.sites.w<-summary.intercepts.across.sites%>%filter(terms =="trtNPK")%>%filter(scale!=sc)%>%filter(q==i & functional_group==fg)%>%
    pivot_wider(names_from = "scale", values_from = c("Estimate", "Q2.5", "Q97.5"))
  overall.change.in.beta<- summary.intercepts.across.sites%>%filter(q==i & functional_group==fg & scale==sc)%>%
    mutate(dif.positive=Q97.5 - Estimate, dif.negative=Estimate - Q2.5)
  # site-level data 
  summary.intercepts.individual.sites_i_fg<-categorize.beta%>%filter(q==i & functional_group==fg)
  (pp.int <- pp.diversity.simple(data1=summary.intercepts.individual.sites_i_fg, data2=summary.intercepts.across.sites.w, data3=overall.change.in.beta) +
      geom_abline(slope=-1, intercept=0)+
    geom_text(data=df.ticks1, aes(x.position, y.position, label=x2), size = 2, angle=45)) 

  # combine the graphs 
  (pp.div1<-plot_grid(pp.nat,
                      pp.int,
                      nrow=1, ncol=2, 
                      vjust = -0.01, 
                      labels = c("A (native species)", "B (non-native species)"), label_fontface = "bold", label_size = 11))
  # adjust the t (top), b (bottom), l (left), and r (right) values for controlling the margins 
  (pp.div3<-pp.div1+ theme(plot.margin = margin(l = 10, b = 10, t=15))+
      draw_label("Effects of nutrient addition on alpha diversity (LRR)", x = 0.5, y = 0, size = 11) + draw_label("Effects of nutrient addition on gamma diversity (LRR)", x = 0, y = 0.5, angle = 90, size = 11))
  ggsave(pp.div3, width = 20, height = 10, dpi=600, unit="cm",  file=paste0("effects on change in diversity across scales for native and invasive for q of ", i, ".png"))
  
  # for "FORB", "GRAMINOID",  "LEGUME",  "WOODY"
  fg<-"FORB"; sc<-"beta"
  summary.intercepts.across.sites.w<-summary.intercepts.across.sites%>%filter(terms =="trtNPK")%>%filter(scale!=sc)%>%filter(q==i & functional_group==fg)%>%
    pivot_wider(names_from = "scale", values_from = c("Estimate", "Q2.5", "Q97.5"))
  overall.change.in.beta<- summary.intercepts.across.sites%>%filter(q==i & functional_group==fg & scale==sc)%>%
    mutate(dif.positive=Q97.5 - Estimate, dif.negative=Estimate - Q2.5)
  # site-level data 
  summary.intercepts.individual.sites_i_fg<-categorize.beta%>%filter(q==i & functional_group==fg)
  # position for annotation
  label.position<-summary.intercepts.individual.sites_i_fg%>%
    summarise(max.alpha=max(true.alpha_Q97.5), max.gamma=max(true.gamma_Q97.5), min.alpha=min(true.alpha_Q2.5), min.gamma=min(true.gamma_Q2.5))%>%
    rowwise() %>% mutate(min_value = 0.85*min(min.alpha, min.gamma), max_value_0 = 0.85*max(max.alpha, max.gamma), max_value=ifelse(max_value_0 < 0.2, 0.2, max_value_0))
  
  (pp.forb <- pp.diversity.simple(data1=summary.intercepts.individual.sites_i_fg, data2=summary.intercepts.across.sites.w, data3=overall.change.in.beta)+
    annotate(geom = "text", x = label.position$min_value/2, y = 0, angle = 45, label = "Differentiation", size = 4, color="grey68") +  
    annotate(geom = "text", x = 0, y = label.position$min_value/2, angle = 45, label = "Homogenization", size = 4, color="grey68") +
    annotate(geom = "text", x = 0.80*label.position$max_value, y=0.3*label.position$max_value, angle = 0, label = "Gain of\n widespread\n species", size = 3, color="grey78") +
    annotate(geom = "text",  x = 0.35*label.position$max_value, y =  0.85*label.position$min_value, angle = 0, label = "Localized\nreplaced by \nwidespread\n species", size = 3, color="grey78") +
    annotate(geom = "text", x = 0.2*label.position$min_value, y = 0.90*label.position$min_value, angle = 0, label = "Loss of\n localized\n species", size = 3, color="grey78") +
    annotate(geom = "text", x = 0.85*label.position$min_value, y = 0.3*label.position$min_value, angle = 0, label = "Loss of\n widespread\n species", size = 3, color="grey78") +
    annotate(geom = "text", x =0.85*label.position$min_value, y = 0.35*label.position$max_value, angle = 0, label = "Widespread\nreplaced by \nlocalized\n species", size = 3, color="grey78") +
    annotate(geom = "text", x = 0.3*label.position$max_value,  y = 0.9*label.position$max_value,  angle = 0, label = "Gain of\n localized\n species", size = 3, color="grey78")+
    geom_abline(slope=-1, intercept=0)+
    geom_text(data=df.ticks1, aes(x.position, y.position, label=x2), size = 2, angle=45)) 

  # for simplified figures without annotate
  list_plots <- vector('list', 3)
  for(fg in c("GRAMINOID", "LEGUME", "WOODY" )){
    sc<-"beta"; 
    summary.intercepts.across.sites.w<-summary.intercepts.across.sites%>%filter(terms =="trtNPK")%>%filter(scale!=sc)%>%filter(q==i & functional_group==fg)%>%
      pivot_wider(names_from = "scale", values_from = c("Estimate", "Q2.5", "Q97.5"))
    
    overall.change.in.beta<- summary.intercepts.across.sites%>%filter(q==i & functional_group==fg & scale==sc)%>%
      mutate(dif.positive=Q97.5 - Estimate, dif.negative=Estimate - Q2.5)
    
    # site-level data 
    summary.intercepts.individual.sites_i_fg<-categorize.beta%>%filter(q==i & functional_group==fg)
    # color for each section 
    # unique(summary.intercepts.individual.sites_i_fg$cat.process)
    
    list_plots[[fg]]<- pp.diversity.simple(data1=summary.intercepts.individual.sites_i_fg, data2=summary.intercepts.across.sites.w, data3=overall.change.in.beta)+
      geom_abline(slope=-1, intercept=0)+
      geom_text(data=df.ticks1, aes(x.position, y.position, label=x2), size = 2, angle=45)
  }
  # combine the graphs 
  (pp.div1<-plot_grid(pp.forb,
                      list_plots$GRAMINOID,
                      list_plots$LEGUME,
                      list_plots$WOODY,
                      nrow=2, ncol=2, 
                      vjust = -0.01, 
                      labels = c("A (forb)", "B (graminoid)", "C (legume)", "D (woody)"), label_fontface = "bold", label_size = 11))
  (pp.div3<-pp.div1+ theme(plot.margin = margin(l = 10, b = 10, t=15))+
      draw_label("Effects of nutrient addition on alpha diversity (LRR)", x = 0.5, y = 0, size = 11) + draw_label("Effects of nutrient addition on gamma diversity (LRR)", x = 0, y = 0.5, angle = 90, size = 11))
  ggsave(pp.div3, width = 20, height = 20, dpi=600, unit="cm",  file=paste0("effects on change in diversity across scales for life forms for q of ", i, ".png"))
  
}

# the end
