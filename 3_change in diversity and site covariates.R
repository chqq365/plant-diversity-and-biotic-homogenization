rm(list=ls())
## open the library
library(tidyverse);library(ggplot2);library(cowplot);library(brms);library(tidybayes)
set.seed(123)
## set up the work directory 
dir.data<-"C:/Users/chqq3/work/NutNet data/"
dir.graphs<-"C:/Users/chqq3/work/homogenization/Nature communications/graphs/"
setwd(dir.graphs)

###########################################################################################
##effects of nutrient addition on alpha, gamma, and beta diversity
###########################################################################################
### focusing on sites with 4-year treatments, also focus on control and NPK
load("data for control and NPK.rdata")
d7.yr4<-d7%>%filter(year_trt%in% c(4))
source("funtion to calculate diversity metrics.R")
alpha.gamma<-calculate.alpha.gamma.diversity(d7.yr4)

eff.nut<-alpha.gamma%>%group_by(site_code, trt, year_trt, functional_group, q, scale, block)%>%summarise(HillDiv1=mean(HillDiv))%>%
  pivot_wider(names_from = "trt", values_from = "HillDiv1")%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK)) %>% mutate(dif.HillDiv=log(NPK) - log(Control))%>%
  group_by(site_code, year_trt, functional_group, q, scale)%>%summarise(dif.HillDiv1=mean(dif.HillDiv))%>%
  pivot_wider(names_from = "scale", values_from = "dif.HillDiv1")%>%mutate(beta=gamma - alpha) %>%
  filter(!is.na(alpha))%>%filter(!is.na(gamma))

##########################################################################################
###### autocorrelation for change in species richness across scales in year 4
###########################################################################################
# add geolocation data and continents 
sites<-read.csv(paste0(dir.data, "comb-by-plot-clim-soil-diversity_2023-11-07.csv"), header=T)
alpha.gamma.beta1<-eff.nut%>%filter(year_trt%in% c(4))%>%filter(q==0 & functional_group=="all")%>%pivot_longer(cols = c("alpha", "gamma", "beta"))
alpha.gamma.beta2<-sites%>%select(site_code, continent, latitude, longitude, country, habitat)%>%distinct()%>%merge(alpha.gamma.beta1, by=c("site_code"))%>%
  mutate(variable.id=paste0(name, functional_group, year_trt))
alpha.gamma.beta2_summary<-alpha.gamma.beta2 %>%group_by(name, continent)%>%summarise(avg=mean(value), sd=sd(value))

library(sf);library(spdep);
spatial.autocorrelation<-data.frame()
for(vi in unique(alpha.gamma.beta2$variable.id)){
  # vi<-"alphaall4"
  temp.data<-alpha.gamma.beta2%>%filter(variable.id== vi)%>%filter(!is.na(value))%>%arrange(continent)
  n.sites<-length(unique(temp.data$site_code))
  
  df_sf <- st_as_sf(temp.data, coords = c("latitude", "longitude"), crs = 4326)
  # Convert the sf object to a SpatialPointsDataFrame for spdep
  df_spatial <- as(df_sf, "Spatial")
  # Calculate the spatial neighbors (3 nearest neighbors)
  k_neigh <- knn2nb(knearneigh(df_spatial, k = 3))
  # Create a spatial weights object
  w <- nb2listw(k_neigh, style = "W")
  # Calculate Moran's I statistic 
  moran_result <- moran.test(temp.data$value, w)
  auto.cor<-data.frame(t(moran_result$estimate))%>%mutate(variable.id=vi, standard.deviate= moran_result$statistic, p=moran_result$p.value )
  spatial.autocorrelation<-spatial.autocorrelation%>%bind_rows(auto.cor)
}
spatial.autocorrelation.sig<-spatial.autocorrelation%>%filter(p <= 0.05)

spatial.autocorrelation1<- spatial.autocorrelation %>% mutate(across(c("Moran.I.statistic", "Expectation", "standard.deviate", "p" ), ~ round(.x, digits = 4)))%>%
  mutate(Diversity.facets=case_when(grepl("alpha", variable.id)~ "∆α", 
                                    grepl("beta", variable.id)~ "∆β", 
                                    grepl("gamma", variable.id)~ "∆γ"))%>%
  select("Diversity.facets", "Moran.I.statistic", "standard.deviate",  "p" )
write.csv(spatial.autocorrelation1, file="check for spatial autocorrelation for change in diversity across scales.csv")

###########################################################################################
#######  autocorrelation for change in residual of species richness across scales in year 4
###########################################################################################
variable.id.combi<- read.csv("combinations for species groups and diversity metrics and years and hill q.csv")%>%mutate(X=NULL)
variable.id.combi_yr4_all.species<-variable.id.combi%>%filter(year_trt==4 & functional_group=="all" & q==0)

spatial.autocorrelation.residual<-c()
for( vi in unique(variable.id.combi_yr4_all.species$variable.id)){
  # vi<- "alldif.alpha40"
 
  load(paste0("model for diversity change for combination of ",vi, ".Rdata") ) 
  diversity_coef <- coef(mod.eff.nut)
  site_code <- rownames(diversity_coef$site_code[,,'Intercept'])
  if(grepl("alpha", vi)){
    site_code.list <- rep(site_code, each=3)
  }else{
    site_code.list <- site_code}
  
  temp.data<-residuals(mod.eff.nut)%>%data.frame()%>%mutate(site_code=site_code.list)%>%group_by(site_code)%>%summarise(Estimate=mean(Estimate))%>%
    merge(sites%>%select(site_code, continent, latitude, longitude, country, habitat)%>%distinct(), by=c("site_code"))
  n.sites<-length(unique(temp.data$site_code))
  
  df_sf <- st_as_sf(temp.data, coords = c("latitude", "longitude"), crs = 4326)
  df_spatial <- as(df_sf, "Spatial")
  # Calculate the spatial neighbors 
  k_neigh <- knn2nb(knearneigh(df_spatial, k = 3))
  w <- nb2listw(k_neigh, style = "W")
  # Calculate Moran's I statistic for residuals
  moran_result <- moran.test(temp.data$Estimate, w)
  auto.cor<-data.frame(t(moran_result$estimate))%>%mutate(variable.id=vi, standard.deviate= moran_result$statistic, p=moran_result$p.value )
  spatial.autocorrelation.residual<-spatial.autocorrelation.residual%>%bind_rows(auto.cor)
}
spatial.autocorrelation.residual.sig<-spatial.autocorrelation.residual%>%filter(p <= 0.05)

##########################################################################################
######supplementary figures: change in species richness across spatial scales and site covariates at year 4
###########################################################################################
'
To be consistent for all sites,
We quantified the drought index as the sum of annual evapotranspiration/precipitation, and averaged it across from year 0 (the start of the experiment) to year 4 at each site.
We quantified site species pool as the total number of species occurring from year 0 (the start of the experiment) to year 4 from at each site and
site productivity as the average biomass from year 0 to 4 from in the ambient conditions.
'
drought<-read.csv(file="site drought index.csv")
site.productivity<-read.csv(file="site productivity.csv")
site.species.pool<-read.csv(file="site species pool.csv")
site.herbivores<-read.csv(file="site herbivores.csv")
block.distance<-read.csv(file="minimum area and distance among blocks.csv")%>%filter(avg.min.distance1< 2000)

# years should be used for climate 
year.used.for.climate<-d7%>%filter(year_trt== 4)%>%select(site_code, year, year_trt)%>%distinct()
# sites and blocks used in richness data 
site.block.select<-d7%>%filter(year_trt== 4)%>%mutate(site.block=paste(site_code, block, sep="_"))%>%select(site.block)%>%distinct()

drought1<-drought%>%merge(year.used.for.climate, by=c("site_code"))%>%filter(year.x <= year.y & year.x >= (year.y- 4) )%>%
  group_by(site_code)%>%summarise(drought.index1=mean(drought.index))

site.productivity1<-site.productivity%>%mutate(site.block=paste(site_code, block, sep="_"))%>%
  filter(site.block %in% site.block.select$site.block)%>% filter(trt %in% c("Control"))%>%filter(year_trt<= 4)
# check years included for each site 
check.year.biomass<-site.productivity1%>%ungroup()%>%select(site_code, year_trt)%>%distinct()%>%group_by(site_code)%>%summarise(N.years=length(year_trt))
table(check.year.biomass$N.years)#
site.productivity2<-site.productivity1%>%group_by(site_code)%>%summarise(site.productivity=mean(live_mass))

# it is very important that each site should have the same number of years. Sampling more times may find more species.
site.species.pool1<-site.species.pool%>%filter(q==0 )%>%filter(site_code %in% year.used.for.climate$site_code)%>% filter(year_trt<= 4)
# check years included for each site 
check.year.richness<-site.species.pool1%>%ungroup()%>%select(site_code, year_trt)%>%distinct()%>%group_by(site_code)%>%summarise(N.years=length(year_trt))
table(check.year.richness$N.years)# 
check.year.richness.yr.years<-check.year.richness%>%filter(N.years == 5)
site.species.pool2<-site.species.pool1%>%filter(site_code %in%check.year.richness.yr.years$site_code)%>%
  group_by(site_code)%>%summarise(site.species.pool=mean(HillDiv))

# add all site covariates 
all.site.covariates<-drought1%>%
  left_join(block.distance, by=c("site_code"))%>%
  left_join(site.species.pool2, by=c("site_code"))%>%
  left_join(site.productivity2, by=c("site_code"))%>%
  left_join(site.herbivores%>%mutate(X=NULL), by=c("site_code"))%>%mutate(X=NULL)%>%
  select(site_code, site.species.pool, site.productivity, herb.index, drought.index1, avg.min.distance1)%>%
  pivot_longer(cols = c("site.species.pool", "site.productivity", "herb.index",  "drought.index1", "avg.min.distance1"))%>% 
  mutate(name1=case_when(name=="site.species.pool" ~"Site species pool", 
                         name=="site.productivity" ~"Site productivity", 
                         name=="herb.index" ~"Grazing intensity", 
                         name=="drought.index1" ~"Drought index",
                         name=="drought.index1" ~"Drought index",
                         name=="avg.min.distance1" ~"Block distance"))

eff.nut1<- eff.nut%>%filter(q==0 & functional_group=="all")%>%ungroup()%>%
  mutate(cat.beta=case_when((beta<0 )~"homogenization", 
                            (beta>0 )~"differentiation", 
                            TRUE~"no change in beta"))%>%
  mutate(cat.process=case_when((alpha>gamma & alpha>0 &gamma>0 )~"Gain of widespread species",
                               (alpha>0 &gamma<0 )~"Spatially restricted replaced by widespread species",
                               (alpha>gamma & alpha<0 &gamma<0)~"Loss of spatially restricted species",
                               (alpha<gamma & alpha<0 &gamma<0)~"Loss of widespread species",
                               (alpha<0 &gamma>0)~"Widespread replaced by spatially restricted species",
                               (alpha<gamma & alpha>0 &gamma>0 )~"Gain of spatially restricted species" ,
                               TRUE~"Other situations"))%>%  select(year_trt, site_code, alpha, gamma, beta, cat.process)%>%
  pivot_longer(cols = c("beta",  "alpha",   "gamma"), names_to = "diversity.facet", values_to = "diversity.value")%>% 
  mutate( diversity.facet1=case_when(grepl("alpha", diversity.facet) ~ "Change in alpha diversity",
                                     grepl("beta", diversity.facet) ~ "Change in beta diversity", 
                                     TRUE~"Change in gamma diversity"))
# check site sier.us
sier.us<-eff.nut1%>%filter(site_code=="sier.us")

change.in.div.cov<-eff.nut1%>% left_join(all.site.covariates, by=c("site_code"))%>%mutate(combi.id=paste(diversity.facet1, name1, sep="_"))     

# use color for different intercepts from conceptual figure
concept_colour = c("Gain of widespread species" =  "#F0E442",
                   "Spatially restricted replaced by widespread species" = "#E69F00",
                   "Loss of spatially restricted species" = "#D55E00",
                   "Loss of widespread species" ="#009E73",
                   "Widespread replaced by spatially restricted species"  =  "#0072B2", 
                   "Gain of spatially restricted species"  = "#56B4E9",
                   "Other situations" = '#f0f0f0')

change.in.div.cov$cat.process<-factor(change.in.div.cov$cat.process, levels=c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                                                              "Loss of widespread species",   "Widespread replaced by spatially restricted species", "Gain of spatially restricted species",  "Other situations"))
# check raw data and relationships
 yr<-4
  # Bivariate relationships for site covariates that often used in previous literature 
  change.in.div.cov_temp<-change.in.div.cov%>%ungroup()%>%filter(year_trt==yr)%>% mutate(grp=paste0(diversity.facet1, name1, set="_"))
  (pp.relation.raw<-change.in.div.cov_temp%>%
      ggplot(aes(value, diversity.value, fill = cat.process, group=grp)) + theme_cowplot()+panel_border()+
      geom_hline(yintercept = 0, linetype="dotted")+
      geom_point(size=3, pch=21, alpha=0.5)+
      geom_smooth(color="black")+
      facet_grid(diversity.facet1~name1, scales = "free", switch = "both") +
      scale_fill_manual(values = concept_colour) +
      guides(fill = "none") +
      theme(legend.position = "top", strip.placement = "outside", strip.background = element_blank()) +
      labs(x=NULL, y=NULL, fill = NULL))
  ggsave(pp.relation.raw, width = 20, height=20, unit="cm", file=paste0("relationships between change in diversity and environmental factors for all species in year ", yr, ".png"))


summary.bivariate.estimate.coef<-c(); pred.slopes<-c()
# Bivariate relationships for site covariates that often used in previous literature 
  change.in.div.cov_temp<-change.in.div.cov%>%ungroup()%>%filter(year_trt==yr)

    # statistic with simple linear regression 
  for(df in unique(change.in.div.cov_temp$combi.id)){
    #  df<-"Change in beta diversity_Site species pool" 
    data.temp<-change.in.div.cov_temp%>%filter(combi.id==df)%>%na.omit()%>%filter_all(all_vars(!is.infinite(.)))
    #hist(data.temp$diversity.value)
    n.sites<-length(unique(data.temp$site_code))
    # model 
    Bivariate.mod <- brm( diversity.value ~ value, data=data.temp, iter=3000, warmup = 1000, cores = 6)
    #pp_check(Bivariate.mod)
    # coefficients
    t.temp<-summary(Bivariate.mod)
    t.temp1<-t.temp$fixed
    summary.bivariate.estimate.coef<-summary.bivariate.estimate.coef%>%bind_rows(t.temp1%>%as.data.frame()%>%mutate(terms=row.names(.), combi.id=df, n.sites=n.sites, year_trt=yr))
    
    # Generate conditional effects
    effects <- conditional_effects(Bivariate.mod, effects = "value", prob = 0.95)
    # Extract conditional effects data
    effects_data <- effects$value
    pred.slopes<-pred.slopes%>%bind_rows(effects_data%>%data.frame()%>%mutate(combi.id=df, year_trt=yr))
   }

summary.bivariate.estimate.coef1<-summary.bivariate.estimate.coef%>%mutate_at(vars(c(1:7)), round, digits=4)%>%
  mutate(variable=ifelse(grepl("Intercept", terms), "Intercept", "Slope"))%>%dplyr::rename(Q2.5='l-95% CI', Q97.5='u-95% CI')%>%
  mutate(sig=ifelse(((Q2.5 > 0 & Q97.5>0)|(Q2.5<0 & Q97.5< 0)), "Significant", "Non-significant"))%>%
  merge(summary.bivariate.estimate.coef%>%filter(terms=="Intercept")%>%select(combi.id, Estimate), by=c("combi.id"))%>%
  filter(terms!="Intercept")%>% dplyr::rename(Intercept=Estimate.y, Slope=Estimate.x)%>%mutate(terms=NULL, variable=NULL)%>%
  merge(change.in.div.cov%>%select(diversity.facet1, name1, combi.id)%>%distinct(), by=c("combi.id"))%>%
  filter(!name1 %in% c("Herbivory"))
# relevel 
table(summary.bivariate.estimate.coef1$sig)
summary.bivariate.estimate.coef1$sig<-factor(summary.bivariate.estimate.coef1$sig, levels = c("Significant", "Non-significant"))
unique(summary.bivariate.estimate.coef1$diversity.facet1 )
summary.bivariate.estimate.coef1$diversity.facet1 <- factor(summary.bivariate.estimate.coef1$diversity.facet1, levels=c("Change in alpha diversity", "Change in gamma diversity", "Change in beta diversity" , "Change in beta_C diversity"))
# Save the table 
Table.S<-summary.bivariate.estimate.coef%>%mutate_at(vars(c(1:7)), round, digits=4)%>% 
  merge(change.in.div.cov%>%select(diversity.facet1, name1, combi.id)%>%distinct(), by=c("combi.id"))%>%
  select("year_trt", "diversity.facet1", "name1", "terms", "Estimate", "l-95% CI", "u-95% CI",   "Rhat" ,  "Bulk_ESS", "Tail_ESS",  "n.sites" )%>%arrange(year_trt, diversity.facet1, name1, terms)
colnames(Table.S)<-c("year_trt", "Diversity facets", "Site covarites", "terms", "Estimate", "l-95% CI", "u-95% CI",   "Rhat" ,  "Bulk_ESS", "Tail_ESS",  "Number of sites")
write.csv(Table.S, file="model output for change in diversity and site covariates.csv")

alpha.crosswalk <- c("Significant" = 0.9, "Non-significant" = 0.1)
change.in.div.cov$diversity.facet1 <- factor(change.in.div.cov$diversity.facet1, levels=c("Change in alpha diversity", "Change in gamma diversity", "Change in beta diversity" ))

pred.slopes1<-pred.slopes%>%  merge(change.in.div.cov%>%select(diversity.facet1, name1, combi.id)%>%distinct(), by=c("combi.id"))%>% mutate(grp=paste0(diversity.facet1, name1, set="_"))%>%
  mutate(diversity.facet2=gsub("alpha", "average alpha", diversity.facet1))
pred.slopes1$diversity.facet1 <- factor(pred.slopes1$diversity.facet1, levels=c("Change in alpha diversity", "Change in gamma diversity", "Change in beta diversity" ))

# plot predicted slopes 
  # predicted Bivariate relationships for site covariates that often used in previous literature 
  pred.slopes1_temp<-pred.slopes1%>%filter(year_trt==yr)
  change.in.div.cov_temp<-change.in.div.cov%>%ungroup()%>%filter(year_trt==yr)%>% mutate(grp=paste0(diversity.facet1, name1, set="_"))%>%
    mutate(diversity.facet2=gsub("alpha", "average alpha", diversity.facet1))
  n.sites<-change.in.div.cov_temp%>%filter(!is.na(value))%>%group_by(grp, name1, diversity.facet2)%>%summarise(N.sites=n(), max.x=0.9*max(value), max.y=max(diversity.value))
  
  (pp.relation.raw<-change.in.div.cov_temp%>%
      ggplot(aes(value, diversity.value)) + theme_cowplot()+panel_border()+
      geom_hline(yintercept = 0, linetype="dotted")+
      geom_point(aes(fill = cat.process, group=grp), size=3, pch=21, alpha=0.5)+
      geom_line(data=pred.slopes1_temp, aes(value, estimate__), color = "black") +  # Regression line
      geom_ribbon(data=pred.slopes1_temp, aes(ymin = lower__, ymax = upper__), alpha = 0.2) +  # CI band
      facet_grid(diversity.facet2~name1, scales = "free", switch = "both") +
      geom_text(data=n.sites, aes(max.x, max.y, label =paste0("(", N.sites, ")")),  size =4) +  
      scale_fill_manual(values = concept_colour) +  guides(fill = guide_legend(nrow = 3)) + 
      theme(legend.position = "top",  legend.text = element_text(size =8), strip.placement = "outside", strip.background = element_blank(), axis.text.x = element_text(angle = 30)) +
      labs(x=NULL, y=NULL, fill = NULL))
  ggsave(pp.relation.raw, width = 20, height=25, unit="cm", file=paste0("predicted relationships between change in diversity and environmental factors for all species in year ", yr, ".png"))
 
# the end 