
library(vegan)
'
We calculated α and γ diversity using Hill numbers with q ranging from 0 to 2 (an increase in q indicating greater weights of abundant species). 
Here, q is a continuous variable, we use q of 0, 1, 2, corresponding to species richness, Shannon diversity, and Simpson diversity. 
'
hill_numbers <- function(abundance, q) {
  if (q == 0) {
    # Species richness
    return(sum(abundance > 0))
  } else if (q == 1) {
    # Exponential of Shannon entropy
    shannon_entropy <- -sum((abundance / sum(abundance)) * log(abundance / sum(abundance)), na.rm = TRUE)
    return(exp(shannon_entropy))
  } else if (q == 2) {
    # Inverse Simpson index
    simpson_index <- sum((abundance / sum(abundance))^2)
    return(1 / simpson_index)
  } else {
    # General case for q
    return((sum((abundance / sum(abundance))^q))^(1 / (1 - q)))
  }
}

calculate.alpha.gamma.diversity<- function(cover.data){
  # cover.data<-d7
  d8_1<- cover.data%>%mutate(functional_group_1="all")
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
  
  # Apply the function to each site for different q values
  results.alpha <- sapply(0:2, function(q) apply(cover, 1, hill_numbers, q = q))
  diversity.q.alpha<-results.alpha%>%data.frame()%>%bind_cols(alpha%>%select(1:6)%>%distinct())%>%pivot_longer(cols = c("X1", "X2", "X3"), names_to = "q0", values_to = "HillDiv")%>%
    mutate(q=case_when(q0=="X1"~ 0, q0=="X2"~ 1, q0=="X3"~2), q0=NULL)

  # gamma diversity 
  gamma<-d8_all_groups%>%
    group_by(site_code, trt, year_trt, functional_group_1, standard_taxon)%>%
    summarise(sum_max_cover1=sum(max_cover))%>%
    pivot_wider(names_from = "standard_taxon", values_from = "sum_max_cover1")%>%distinct()
  table(gamma$functional_group_1)
  
  cover1<-gamma[,5:ncol(gamma)]
  cover1[is.na(cover1)]<-0
  
  # Apply the function to each site for different q values
  results.gamma <- sapply(0:2, function(q) apply(cover1, 1, hill_numbers, q = q))
  diversity.q.gamma<-results.gamma%>%data.frame()%>%bind_cols(gamma%>%select(1:4)%>%distinct())%>%pivot_longer(cols = c("X1", "X2", "X3"), names_to = "q0", values_to = "HillDiv")%>%
    mutate(q=case_when(q0=="X1"~ 0, q0=="X2"~ 1, q0=="X3"~2), q0=NULL)

  # add them together 
  alpha.gamma<-diversity.q.alpha%>%mutate(scale="alpha")%>%select(HillDiv, site_code, trt, year_trt, functional_group_1, q, scale, block, plot)%>%
    bind_rows(diversity.q.gamma%>%mutate(scale="gamma", block=0, plot=0))%>%dplyr::rename(functional_group=functional_group_1)
  # return(alpha.gamma)
}