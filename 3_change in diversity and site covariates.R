
rm(list=ls())
## open the library
library(tidyverse);library(brms);library(tidybayes);library(ggplot2);library(cowplot)
library(sf);library(spdep); # spatial analyses 
set.seed(123)

# use color for different intercepts from conceptual figure
concept_colour = c("Gain of widespread species" =  "#F0E442",
                   "Spatially restricted replaced by widespread species" = "#E69F00",
                   "Loss of spatially restricted species" = "#D55E00",
                   "Loss of widespread species" ="#009E73",
                   "Widespread replaced by spatially restricted species"  = "#56B4E9", 
                   "Gain of spatially restricted species" = "#0072B2", 
                   "Other situations" = '#f0f0f0')

## set up the work directory 
dir.data<-"C:/Users/chqq3/work/NutNet data/"
dir.graphs<-"C:/Users/chqq3/work/homogenization/Nature communications/graphs1/"
setwd(dir.graphs)

###########################################################################################
##effects of nutrient addition on alpha, gamma, and beta diversity
###########################################################################################
### focusing on sites with 4-year treatments
yr<- 4
eff.nut<- read.csv("effects of nutrient addition on diversity across scales.csv")
eff.nut1<- eff.nut%>%filter(q==0 & functional_group=="all")%>%filter(year_trt%in% c(4))

##########################################################################################
###### autocorrelation for change in species richness across scales in year 4
###########################################################################################
# add geolocation data and Continents 
sites<-read.csv("sites with geolocation and experimental years used.csv")
alpha.gamma.beta1<-eff.nut1%>%filter(q==0 & functional_group=="all")%>%pivot_longer(cols = c("alpha", "gamma", "beta"))
alpha.gamma.beta2<-sites%>%select(site_code,Latitude, Longitude, Continent, Habitat)%>%distinct()%>%merge(alpha.gamma.beta1, by=c("site_code"))%>%
  mutate(variable.id=paste0(name, functional_group, year_trt))
alpha.gamma.beta2_summary<-alpha.gamma.beta2 %>%group_by(name, Continent)%>%summarise(avg=mean(value), sd=sd(value))

spatial.autocorrelation<-c()
for(vi in unique(alpha.gamma.beta2$variable.id)){
  # vi<-"alphaall4"
  temp.data<-alpha.gamma.beta2%>%filter(variable.id== vi)%>%filter(!is.na(value))%>%arrange(Continent)
  n.sites<-length(unique(temp.data$site_code))
  
  df_sf <- st_as_sf(temp.data, coords = c("Latitude", "Longitude"), crs = 4326)
  # Convert the sf object to a SpatialPointsDataFrame for spdep
  df_spatial <- as(df_sf, "Spatial")
  # Calculate the spatial neighbors (2 nearest neighbors)
  k_neigh <- knn2nb(knearneigh(df_spatial, k = 2))
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
site.species.pool<-read.csv(file="site species pool.csv")%>%mutate(site.species.pool=species.pool)
site.herbivores<-read.csv(file="site herbivores.csv")
block.distance<-read.csv(file="minimum area and distance among blocks.csv")%>%filter(avg.min.distance1< 2000)

# add all site covariates 
all.site.covariates<-drought%>%
  left_join(block.distance, by=c("site_code"))%>%
  left_join(site.species.pool, by=c("site_code"))%>%
  left_join(site.productivity, by=c("site_code"))%>%
  left_join(site.herbivores%>%mutate(X=NULL), by=c("site_code"))%>%mutate(X=NULL)%>%
  select(site_code, site.species.pool, site.productivity, herb.index, drought.index1, avg.min.distance1)%>%
  pivot_longer(cols = c("site.species.pool", "site.productivity", "herb.index",  "drought.index1", "avg.min.distance1"))%>% 
  mutate(name1=case_when(name=="site.species.pool" ~"Site species pool", 
                         name=="site.productivity" ~"Site productivity", 
                         name=="herb.index" ~"Grazing intensity", 
                         name=="drought.index1" ~"Drought index",
                         name=="drought.index1" ~"Drought index",
                         name=="avg.min.distance1" ~"Block distance"))

unique(eff.nut1$cat.process)
table(eff.nut1$cat.process)
eff.nut1$cat.process<-factor(eff.nut1$cat.process, levels = c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                                            "Loss of widespread species",  "Widespread replaced by spatially restricted species", "Gain of spatially restricted species" ,
                                                         #   "Loss of spatially restricted and widespread species at simimar magnitude", "Gain of spatially restricted and widespread species at simimar magnitude",
                                                            "Other situations"))

change.in.div.cov<-eff.nut1%>% select(year_trt, site_code, alpha, gamma, beta, cat.process)%>%
  pivot_longer(cols = c("beta",  "alpha",   "gamma"), names_to = "diversity.facet", values_to = "diversity.value")%>% 
  mutate( diversity.facet1=case_when(grepl("alpha", diversity.facet) ~ "Change in alpha diversity",
                                     grepl("beta", diversity.facet) ~ "Change in beta diversity", 
                                     TRUE~"Change in gamma diversity"))

change.in.div.cov_1<-change.in.div.cov%>%ungroup()%>% left_join(all.site.covariates, by=c("site_code"))%>%mutate(combi.id=paste(diversity.facet1, name1, sep="_"))%>%
  mutate(grp=paste0(diversity.facet1, name1, set="_"))

# check raw data and relationships
  (pp.relation.raw<-change.in.div.cov_1%>%
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
   # statistic with simple linear regression 
  for(df in unique(change.in.div.cov_1$combi.id)){
    #  df<- "Change in beta diversity_Site species pool"
    data.temp<-change.in.div.cov_1%>%filter(combi.id==df)%>%na.omit()%>%filter_all(all_vars(!is.infinite(.)))
    #hist(data.temp$diversity.value)
    n.sites<-length(unique(data.temp$site_code))
    # model 
    Bivariate.mod <- brm( diversity.value ~ value, data=data.temp, iter=3000, warmup = 1000, cores = 6)
    # pp_check(Bivariate.mod)
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
  merge(change.in.div.cov_1%>%select(diversity.facet1, name1, combi.id)%>%distinct(), by=c("combi.id"))
# relevel 
table(summary.bivariate.estimate.coef1$sig)
summary.bivariate.estimate.coef1.sig <- summary.bivariate.estimate.coef1%>%filter(sig %in% c("Significant"))
unique(summary.bivariate.estimate.coef1$diversity.facet1 )
summary.bivariate.estimate.coef1$diversity.facet1 <- factor(summary.bivariate.estimate.coef1$diversity.facet1, levels=c("Change in alpha diversity", "Change in gamma diversity", "Change in beta diversity" , "Change in beta_C diversity"))
# Save the table 
Table.S<-summary.bivariate.estimate.coef%>%mutate_at(vars(c(1:7)), round, digits=4)%>% 
  merge(change.in.div.cov_1%>%select(diversity.facet1, name1, combi.id)%>%distinct(), by=c("combi.id"))%>%
  select("year_trt", "diversity.facet1", "name1", "terms", "Estimate", "l-95% CI", "u-95% CI",   "Rhat" ,  "Bulk_ESS", "Tail_ESS",  "n.sites" )%>%arrange(year_trt, diversity.facet1, name1, terms)
colnames(Table.S)<-c("year_trt", "Diversity facets", "Site covarites", "terms", "Estimate", "l-95% CI", "u-95% CI",   "Rhat" ,  "Bulk_ESS", "Tail_ESS",  "Number of sites")
write.csv(Table.S, file="model output for change in diversity and site covariates.csv")

# plot predicted slopes 
# predicted Bivariate relationships for site covariates that often used in previous literature 
pred.slopes1<-pred.slopes%>%  merge(change.in.div.cov_1%>%select(diversity.facet1, name1, combi.id)%>%distinct(), by=c("combi.id"))%>% 
  mutate(diversity.facet2=gsub("alpha", "average alpha", diversity.facet1))%>%
  mutate(sig=ifelse(combi.id %in%summary.bivariate.estimate.coef1.sig$combi.id, "Significant", "Non-significant" ))
pred.slopes1$sig<-factor(pred.slopes1$sig, levels = c("Significant", "Non-significant"))

pred.slopes1$diversity.facet1 <- factor(pred.slopes1$diversity.facet1, levels=c("Change in alpha diversity", "Change in gamma diversity", "Change in beta diversity" ))

change.in.div.cov_1$diversity.facet1 <- factor(change.in.div.cov_1$diversity.facet1, levels=c("Change in alpha diversity", "Change in gamma diversity", "Change in beta diversity" ))
change.in.div.cov_2<-change.in.div.cov_1%>%ungroup()%>%  mutate(diversity.facet2=gsub("alpha", "average alpha", diversity.facet1))

n.sites<-change.in.div.cov_2%>%filter(!is.na(value))%>%group_by(combi.id, name1, diversity.facet2)%>%mutate(N.sites=n(), max.x=0.9*max(value))%>%
    group_by(diversity.facet2)%>%mutate(max.y=max(diversity.value))%>%arrange(combi.id)
  
# check.data<-change.in.div.cov_2%>%filter(diversity.facet2=="Change in gamma diversity" & name1=="Site species pool")%>%filter(!is.na(value))
  
(pp.relation.raw<-change.in.div.cov_2%>%
      ggplot(aes(value, diversity.value)) + theme_cowplot(font_size=11)+panel_border()+
      geom_hline(yintercept = 0, linetype="dotted")+
      geom_point(aes(fill = cat.process, group=combi.id), size=3, pch=21, alpha=0.4)+
      geom_line(data=pred.slopes1, aes(value, estimate__, linetype = sig), color = "black") +  # Regression line
      geom_ribbon(data=pred.slopes1, aes(ymin = lower__, ymax = upper__), alpha = 0.2) +  # CI band
      facet_grid(diversity.facet2~name1, scales = "free", switch = "both") +
      geom_text(data=n.sites, aes(max.x, max.y, label =paste0("(", N.sites, ")")),  size =4) +  
      scale_linetype_manual(values = c("solid", "dashed"))+
      scale_fill_manual(values = concept_colour) +  guides(fill = guide_legend(nrow = 3), linetype="none") + 
      theme(legend.position = "top",  legend.text = element_text(size =8), strip.placement = "outside", strip.background = element_blank(), axis.text.x = element_text(angle = 30)) +
    theme(legend.text = element_text(margin = margin(l = -0.2, unit = "pt")), legend.spacing = unit(5, "mm"))+
    labs(x=NULL, y=NULL, fill = NULL, linetype=NULL))
ggsave(pp.relation.raw, width = 18, height=20, unit="cm", file=paste0("predicted relationships between change in diversity and environmental factors for all species in year ", yr, ".png"))
 
##########################################################################################
######source data file 
###########################################################################################
colnames(pred.slopes1)

source.data<-change.in.div.cov_2%>%mutate(data.type="Observed", lower__ = 999, upper__=999)%>%select("data.type", "combi.id" ,  "site_code", "diversity.facet2", "lower__", "upper__",  "cat.process" , "diversity.value", "name1",  "value")%>%
  bind_rows(pred.slopes1%>%mutate(data.type="Predicted", diversity.value=estimate__, site_code="", cat.process="")%>%select("data.type", "combi.id" ,  "site_code", "diversity.facet2", "lower__", "upper__",  "cat.process" , "diversity.value", "name1",  "value"))%>%
  arrange(combi.id)%>%mutate(across(c("lower__", "upper__",  "diversity.value",  "value"), ~ round(.x, digits = 4)))%>% mutate_all(~ replace(., . == 999.0000, ""))
colnames(source.data)<-c("Data type", "Combi.id" ,  "site_code", "Diversity across scales", "l-95% CI" , "u-95% CI", "Categories in six scenarios", "Diversity value", "Site covariates",  "Covariates value")
write.csv(source.data, file="SFig3_data_source.csv")
# the end 