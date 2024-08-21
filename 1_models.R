rm(list=ls())
## open the libraries
library(tidyverse);library(ggplot2);library(cowplot)
library(chemodiv)
library(brms);library(tidybayes)
set.seed(123)

## set up the work directory 
dir.data<-"C:/Users/chqq3/work/homogenization/raw data/"
dir.graphs<-"C:/Users/chqq3/work/homogenization/graphs2/"
setwd(dir.graphs)

##########################################################################################
###### write a function to calculate species diversity at alpha and gamma spatial scales 
###########################################################################################
### focusing on sites with 4-year treatments, also focus on control and NPK
d7<-read.csv("cover data for control and NPK at 72 NutNet sites.csv")

calculate.alpha.gamma.diversity<- function(cover.data){
  # cover.data<-d7
  d8_1<- cover.data%>%mutate(functional_group_1="all")
  table(d7$local_provenance)
  d8_2<- cover.data%>%mutate(functional_group_1=local_provenance)%>%filter(functional_group_1 %in% c("INT", "NAT"))
  
  table(cover.data$functional_group)
  d8_3<- cover.data%>%mutate(functional_group_1=ifelse(functional_group %in% c("GRAMINOID", "GRASS"), "GRAMINOID", functional_group))%>%
    filter(functional_group_1 %in% c("GRAMINOID", "FORB", "LEGUME", "WOODY"))
  # add all functional groups together 
  colnames(d8_1)
  d8_all_groups<-d8_1%>%ungroup()%>%select(site_code, block, plot, trt, year_trt, functional_group_1, standard_taxon, max_cover)%>%
    bind_rows(d8_2%>%ungroup()%>%select(site_code, block, plot, trt, year_trt, functional_group_1, standard_taxon, max_cover))%>%
    bind_rows(d8_3%>%ungroup()%>%select(site_code, block, plot, trt, year_trt, functional_group_1, standard_taxon, max_cover))
  alpha<-d8_all_groups %>%
    pivot_wider(names_from = "standard_taxon", values_from = "max_cover")%>%distinct()
  table(alpha$functional_group_1)
  
  cover<-alpha[,7:ncol(alpha)]
  cover[is.na(cover)]<-0
  
  diversity.q.alpha<-c()
  for (i in 0:2){
    # i<-0
    div_temp <-calcDiv(cover, type="HillDiv", q=i)%>%bind_cols(alpha%>%select(1:6)%>%distinct())%>%mutate(q=i) 
    diversity.q.alpha<-rbind(diversity.q.alpha, div_temp)
  }
  
  # gamma diversity 
  gamma<-d8_all_groups%>%
    group_by(site_code, trt, year_trt, functional_group_1, standard_taxon)%>%
    summarise(sum_max_cover1=sum(max_cover))%>%
    pivot_wider(names_from = "standard_taxon", values_from = "sum_max_cover1")%>%distinct()
  table(gamma$functional_group_1)
  
  cover1<-gamma[,5:ncol(gamma)]
  cover1[is.na(cover1)]<-0
  
  diversity.q.gamma<-c()
  for (i in 0:2){
    # i<-0
    div_temp <-calcDiv(cover1, type="HillDiv", q=i)%>%bind_cols(gamma%>%select(1:4))%>%
      mutate(q=i)
    diversity.q.gamma<-rbind(diversity.q.gamma, div_temp)
  }
  
  # add them together 
  #colnames(diversity.q.alpha)
  #colnames(diversity.q.gamma)
  alpha.gamma<-diversity.q.alpha%>%mutate(scale="alpha")%>%select(HillDiv, site_code, trt, year_trt, functional_group_1, q, scale, block, plot)%>%
    bind_rows(diversity.q.gamma%>%mutate(scale="gamma", block=0, plot=0))%>%dplyr::rename(functional_group=functional_group_1)
  # return(alpha.gamma)
}

alpha.gamma<-calculate.alpha.gamma.diversity(d7)
check.sites<-alpha.gamma%>%filter(site_code %in% c("gilb.za", "pinj.au"))%>%select(site_code, year_trt)%>%distinct()
# these two sites do not have cover data in year 4. 

# check range for different functional groups 
check.range<-alpha.gamma%>%filter(scale=="gamma" & q==0)%>%filter(functional_group%in% c("FORB", "GRAMINOID",  "LEGUME",  "WOODY"))%>%
  group_by(functional_group)%>%summarise(ran=range(HillDiv))

# show raw data 
i<-0
(p.overall<-alpha.gamma%>%filter(q==i)%>%
    ggplot(aes(year_trt, HillDiv,  color=trt, shape=scale, linetype=scale))+theme_cowplot(font_size = 30)+panel_border()+
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
alpha.gamma1<-alpha.gamma%>%mutate(HillDiv.log=log(HillDiv))%>%filter(year_trt==4)
length(unique(alpha.gamma1$site_code))

for(i in c(0, 1, 2)){
  for(fg in unique(alpha.gamma1$functional_group)){
    diversity.alpha.temp<-alpha.gamma1%>%filter(q==i & scale=="alpha" & functional_group==fg)
    mod.alpha <- brm( HillDiv.log ~ trt  + (trt  | site_code/block), 
                      data = diversity.alpha.temp , cores = 6, iter=3000, warmup = 1000, chains = 6, control = list(adapt_delta = 0.99))
    
    diversity.gamma.temp<-alpha.gamma1%>%filter(q==i & scale=="gamma" & functional_group==fg)
    mod.gamma <- brm( HillDiv.log ~ trt  + (trt | site_code), 
                      data = diversity.gamma.temp , cores = 6, iter=3000, warmup = 1000, chains = 6, control = list(adapt_delta = 0.99))
    
    save(mod.alpha, file=paste0("model for alpha diversity for ", fg , " species with q of ", i, " for sites with 4-year nutrient addition.Rdata"))
    save(mod.gamma, file=paste0("model for gamma diversity for ", fg , " species with q of ", i, " for sites with 4-year nutrient addition.Rdata"))
  }
}

##########################################################################################
### robustness test using sites with long-term treatments 
##########################################################################################

for (yr in c(10, 14)){ 
  # yr<-10
  alpha.gamma1<-alpha.gamma%>% mutate(HillDiv.log=log(HillDiv))%>%filter(year_trt== yr)
  #length(unique(alpha.gamma1$site_code))

  # alpha scale 
  for(i in c(0, 2)){
    for(fg in unique(alpha.gamma1$functional_group)){
      # i<-0; fg<-"all"
      diversity.alpha.temp<-alpha.gamma1%>%filter(q==i & scale=="alpha" & functional_group==fg)
      
      mod.alpha <- brm( HillDiv.log ~ trt + (trt | site_code/block), 
                        data = diversity.alpha.temp , cores = 6, iter=3000, warmup = 1000, chains = 6)
      #
      diversity.gamma.temp<-alpha.gamma1%>%filter(q==i & scale=="gamma" & functional_group==fg)
      
      mod.gamma <- brm( HillDiv.log ~ trt + (trt | site_code), 
                        data = diversity.gamma.temp , cores = 6, iter=3000, warmup = 1000, chains = 6)
      
      save(mod.alpha, file=paste0("model for alpha diversity for ", fg , " species with q of ", i, " for sites with ", yr, "-year nutrient addition.Rdata"))
      save(mod.gamma, file=paste0("model for gamma diversity for ", fg , " species with q of ", i, " for sites with ", yr, "-year nutrient addition.Rdata"))
    }
  }
}

##########################################################################################
### robustness test using sites with 4 blocks at year_trt 4
##########################################################################################
coverDF3<-read.csv("cover data for control and NPK at NutNet sites with more than 3 blocks.csv")
#Jon Bakker: sevi.us: plots completely randomly assigned to treatments. Verified with S. Collins on 181024.
#Assigning plots to 5 'pseudo-blocks' (as contiguous as feasible) for analysis.
coverDF3$block[coverDF3$site_code == "sevi.us" & coverDF3$plot %in% c(2, 3, 4, 12, 13, 17, 21, 26)] <- 2
coverDF3$block[coverDF3$site_code == "sevi.us" & coverDF3$plot %in% c(5, 8, 9, 10, 15, 19, 20, 30)] <- 3
coverDF3$block[coverDF3$site_code == "sevi.us" & coverDF3$plot %in% c(14, 18, 23, 24, 25, 31, 34, 35)] <- 4
coverDF3$block[coverDF3$site_code == "sevi.us" & coverDF3$plot %in% c(28, 29, 33, 36, 37, 38, 39, 40)] <- 5

# check number of blocks per year for each site 
site.year.select<-d7%>%mutate(site.year=paste(site_code, year_trt, sep="_"))%>%select(site.year)%>%distinct()
ddd<-coverDF3%>%mutate(site.year=paste(site_code, year_trt, sep="_"))%>%filter(trt %in% c("Control", "NPK"))%>%
  filter(site.year %in% site.year.select$site.year)
# check number of blocks per year for each site 
che.ddd.block.year<-ddd%>%ungroup()%>%select(site_code, year_trt, block)%>%distinct()%>%
  group_by(site_code, year_trt)%>%summarise(N=length(block))%>%filter(N>3)
unique(che.ddd.block.year$site_code)
table(che.ddd.block.year$N)

# check whether the block number is the same over years at each site
che.ddd.block.year2<-che.ddd.block.year%>%select(site_code, N)%>%distinct()%>%group_by(site_code)%>%summarise(N.unique=length(site_code))%>%filter(N.unique>1)
che.ddd.block.year2$site_code
# select those sites with 4 blocks
ddd.more.than.3.blocks<-ddd%>%filter(site_code %in% che.ddd.block.year$site_code)
# check number of total plots per year for each site 
che.ddd.plot.year<-ddd.more.than.3.blocks%>%ungroup()%>%select(site_code, year_trt,block,plot, trt)%>%distinct()%>%group_by(site_code, year_trt)%>%summarise(N=length(trt))
# check sites where plots are different in different years 
che.ddd.plot.year1<-che.ddd.plot.year%>%select(site_code, N)%>%distinct()%>%group_by(site_code)%>%summarise(N.unique=length(site_code))%>%filter(N.unique>1)

# delete some blocks to keep 4 blocks per site according to Jon Bakker 
che.ddd.block.year%>%select(site_code, N)%>%distinct()
#azitwo.cn     5 (delete 5)
#cbgb.us: six blocks (delete 5, 6)
#cdcr.us: five blocks (delete 5)
#cdpt.us: six blocks. On 170802, J. Knops suggested grouping blocks 1, 5, 6 (different) and 2, 3, 4 (similar). Delete blocks 5, 6
#kbs.us: five blocks (delete  5)
#msla.us: five blocks (delete  5)
#msla_2.us: five blocks (delete  5)
#msla_3.us: five blocks (delete  5)
#sevi.us: five blocks (delete  2) 
#sier.us: five blocks (delete 5); no Yr0 coverDF3 for blocks 4-5
#veluwe.nl: five blocks (delete  5)  
ddd.more.than.3.blocks_1<-ddd.more.than.3.blocks%>%filter(!(site_code == "cbgb.us" & block %in% c(5, 6)))%>%
  filter(!(site_code == "azitwo.cn" & block %in% c(5)))%>%
  filter(!(site_code == "cdcr.us" & block %in% c(5)))%>%
  filter(!(site_code == "cdpt.us" & block %in% c(5, 6)))%>%
  filter(!(site_code == "kbs.us" & block %in% c(5)))%>%
  filter(!(site_code == "msla.us" & block %in% c(5)))%>%
  filter(!(site_code == "msla_2.us" & block %in% c(5)))%>%
  filter(!(site_code == "msla_3.us" & block %in% c(5)))%>%
  filter(!(site_code == "sevi.us" & block %in% c(2)))%>%
  filter(!(site_code == "sier.us" & block %in% c(5)))%>%
  filter(!(site_code == "veluwe.nl" & block %in% c(5)))%>%filter(year_trt==4)
# check number of replicates for treatments for each site at year_trt 4
che.ddd.more.than.3.blocks.year<-ddd.more.than.3.blocks_1%>%ungroup()%>%select(site_code, year_trt,block,plot, trt)%>%distinct()%>%group_by(site_code, year_trt, trt)%>%
  summarise(N=length(trt))
range(che.ddd.more.than.3.blocks.year$N)
length(unique(che.ddd.more.than.3.blocks.year$site_code))

alpha.gamma<-calculate.alpha.gamma.diversity(ddd.more.than.3.blocks_1)
# show raw data 
i<-0
(p.overall<-alpha.gamma%>%filter(q==i)%>%
    ggplot(aes(year_trt, HillDiv,  color=trt, shape=scale, linetype=scale))+theme_cowplot()+panel_border()+
    facet_grid(site_code~functional_group)+
    geom_point(size=2.5, alpha=0.2)+
    # geom_line()+
    geom_smooth(se=F)+
    labs(x="Years after treatments", y=paste0("Hill number (Q = ", i, ")"), color=NULL, shape=NULL, linetype=NULL))

alpha.gamma1<-alpha.gamma%>%mutate(HillDiv.log=log(HillDiv))%>%filter(year_trt== 4)

for(i in c(0, 2)){
  for(fg in unique(alpha.gamma1$functional_group)){
    # i<-0; fg<-"all"
    diversity.alpha.temp<-alpha.gamma1%>%filter(q==i & scale=="alpha" & functional_group==fg)
    mod.alpha <- brm( HillDiv.log ~ trt + (trt | site_code/block), 
                      data = diversity.alpha.temp , cores = 6, iter=3000, warmup = 1000, chains = 6)
    
    diversity.gamma.temp<-alpha.gamma1%>%filter(q==i & scale=="gamma" & functional_group==fg)
    mod.gamma <- brm( HillDiv.log ~ trt + (trt | site_code), 
                      data = diversity.gamma.temp , cores = 6, iter=3000, warmup = 1000, chains = 6)
    
    save(mod.alpha, file=paste0("model for alpha diversity for ", fg , " species with q of ", i, " for sites with 4 blocks.Rdata"))
    save(mod.gamma, file=paste0("model for gamma diversity for ", fg , " species with q of ", i, " for sites with 4 blocks.Rdata"))
  }
}

##########################################################################################
### robustness test using sites with 5 blocks at year_trt 4
##########################################################################################

# check number of blocks per year for each site 
che.ddd.block.year<-ddd%>%ungroup()%>%select(site_code, year_trt, block)%>%distinct()%>%
  group_by(site_code, year_trt)%>%summarise(N=length(block))%>%filter(N>4)
unique(che.ddd.block.year$site_code)
table(che.ddd.block.year$N)
# select those sites with 5 blocks
ddd.more.than.3.blocks<-ddd%>%filter(site_code %in% che.ddd.block.year$site_code)
# check number of total plots per year for each site 
che.ddd.plot.year<-ddd.more.than.3.blocks%>%ungroup()%>%select(site_code, year_trt,block,plot, trt)%>%distinct()%>%group_by(site_code, year_trt)%>%summarise(N=length(trt))
# check sites where plots are different in different years 
che.ddd.plot.year1<-che.ddd.plot.year%>%select(site_code, N)%>%distinct()%>%group_by(site_code)%>%summarise(N.unique=length(site_code))%>%filter(N.unique>1)

# delete some blocks to keep 5 blocks per site according to Jon Bakker 
che.ddd.block.year%>%select(site_code, N)%>%distinct()
#cbgb.us: six blocks (delete 6)
#cdpt.us: six blocks. On 170802, J. Knops suggested grouping blocks 1, 5, 6 (different) and 2, 3, 4 (similar). Delete blocks 6

ddd.more.than.3.blocks_1<-ddd.more.than.3.blocks%>%filter(!(site_code == "cbgb.us" & block %in% c(6)))%>%
  filter(!(site_code == "cdpt.us" & block %in% c(6)))%>%filter(year_trt==4)

# check number of replicates for treatments per year for each site 
che.ddd.more.than.3.blocks.year<-ddd.more.than.3.blocks_1%>%ungroup()%>%select(site_code, year_trt,block,plot, trt)%>%distinct()%>%
  group_by(site_code, year_trt, trt)%>%summarise(N=length(trt))
range(che.ddd.more.than.3.blocks.year$N)
length(unique(che.ddd.more.than.3.blocks.year$site_code))

alpha.gamma<-calculate.alpha.gamma.diversity(ddd.more.than.3.blocks_1)
# show raw data 
i<-0
(p.overall<-alpha.gamma%>%filter(q==i & functional_group=="all")%>%
    ggplot(aes(year_trt, HillDiv,  color=trt, shape=scale, linetype=scale))+theme_cowplot()+panel_border()+
    facet_grid(site_code~functional_group)+
    geom_point(size=2.5, alpha=0.2)+
    # geom_line()+
    geom_smooth(se=F)+
    labs(x="Years after treatments", y=paste0("Hill number (Q = ", i, ")"), color=NULL, shape=NULL, linetype=NULL))

alpha.gamma1<-alpha.gamma%>%mutate(HillDiv.log=log(HillDiv))%>%filter(year_trt==4)

for(i in c(0, 2)){
  for(fg in unique(alpha.gamma1$functional_group)){
    # i<-0; fg<-"LEGUME"
    diversity.alpha.temp<-alpha.gamma1%>%filter(q==i & scale=="alpha" & functional_group==fg)
    mod.alpha <- brm( HillDiv.log ~ trt + (trt | site_code/block), 
                      data = diversity.alpha.temp , cores = 6, iter=3000, warmup = 1000, chains = 6)
    
    diversity.gamma.temp<-alpha.gamma1%>%filter(q==i & scale=="gamma" & functional_group==fg)
    mod.gamma <- brm( HillDiv.log ~ trt + (trt | site_code), 
                      data = diversity.gamma.temp , cores = 6, iter=3000, warmup = 1000, chains = 6)
    
    save(mod.alpha, file=paste0("model for alpha diversity for ", fg , " species with q of ", i, " for sites with 5 blocks.Rdata"))
    save(mod.gamma, file=paste0("model for gamma diversity for ", fg , " species with q of ", i, " for sites with 5 blocks.Rdata"))
  }
}

# the end 
