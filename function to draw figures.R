
# 3 datasets are needed for this function 
# data1: site-level mean and sd in change of alpha and gamma diversity
# data2: global mean and sd in change of alpha, gamma, and beta diversity
# data3: position and angles for global sd in change of beta diversity

pp.diversity.annotation<-function(data1, data2, data3){
  concept_colour = c(  "Gain widespread species" = '#D9D956',
                       "Localized species replaced by widespread species"  = '#E9AE27',
                       "Loss of widespread species" = '#155F49',
                       "Loss of localized species" = '#D17538',
                       "Widespread species replaced by localized species"  ='#2CABA4', 
                       "Gain localized species" = '#61CDE0',
                        "Other situations" = '#f0f0f0')
  #Our transformation function
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  # position for CI for change in beta diversity 
  angle.positive<- 135; angle.negative<- -45; dif.positive<-data3$dif.positive; dif.negative<-data3$dif.negative
  # position for labels to show 8 sections and x y limits 
  label.position<-data1%>%
    summarise(max.alpha=max(true.alpha_mean), max.gamma=max(true.gamma_mean), min.alpha=min(true.alpha_mean), min.gamma=min(true.gamma_mean))%>%
    rowwise() %>% mutate(min_value = 1.5*min(min.alpha, min.gamma), max_value_0 = 1.5*max(max.alpha, max.gamma), max_value=ifelse(max_value_0 < 0.2, 0.2, max_value_0))
  
  ggplot(data=data1)+ theme_cowplot()+panel_border()+
    geom_point(aes(x=true.alpha_mean, y=true.gamma_mean, color=cat.process), size=2, alpha=0.5)+
    geom_errorbar(aes(x=true.alpha_mean, y=true.gamma_mean, xmin=true.alpha_Q2.5, xmax=true.alpha_Q97.5, color=cat.process), width=0.0001, linewidth=0.2, alpha=0.15) +
    geom_errorbar(aes(x=true.alpha_mean, y=true.gamma_mean, ymin=true.gamma_Q2.5, ymax=true.gamma_Q97.5, color=cat.process), width=0.0001, linewidth=0.2, alpha=0.15) +
    geom_segment(aes(x=true.alpha_mean, y=true.gamma_mean, xend = true.alpha_mean + (change.beta_Q97.5 - change.beta_mean) * cos(angle.positive * pi / 180), yend = true.gamma_mean + (change.beta_Q97.5 - change.beta_mean) * sin(angle.positive * pi / 180), color=cat.process), lineend = "square",  linewidth=0.2, alpha=0.15) +
    geom_segment(aes(x=true.alpha_mean, y=true.gamma_mean, xend = true.alpha_mean + abs(change.beta_Q2.5 - change.beta_mean) * cos(angle.negative * pi / 180), yend = true.gamma_mean + abs(change.beta_Q2.5 - change.beta_mean) * sin(angle.negative * pi / 180), color=cat.process),  lineend = "square",  linewidth=0.2, alpha=0.15) +
    
    geom_errorbar(data=data2, aes(x=Estimate_alpha, y=Estimate_gamma, xmin=Q2.5_alpha, xmax=Q97.5_alpha), width=0.0001, color="black", linewidth=1, alpha=0.7) +
    geom_errorbar(data=data2,  aes(x=Estimate_alpha, y=Estimate_gamma, ymin=Q2.5_gamma, ymax=Q97.5_gamma), width=0.0001, color="black", linewidth=1,  alpha=0.7) +
    geom_segment(data=data2,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.positive * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.positive * pi / 180)),  lineend = "square", color="black", linewidth=1, alpha=0.7) +
    geom_segment(data=data2,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.negative * cos(angle.negative * pi / 180), yend = Estimate_gamma + dif.negative * sin(angle.negative * pi / 180)),   lineend = "square", color="black", linewidth=1, alpha=0.7) +
    geom_point(data=data2,  aes(x=Estimate_alpha, y=Estimate_gamma), size=5, stroke = 1, pch=21, alpha=0.7, color="black")+
    
    geom_abline(slope=1, intercept=0, linetype="dotted")+
     geom_hline(yintercept = 0, linetype="dotted")+
    geom_vline(xintercept = 0, linetype="dotted")+
    scale_x_continuous(limits = c(label.position$min_value, label.position$max_value), labels=scaleFUN)+
    scale_y_continuous(limits = c(label.position$min_value, label.position$max_value), labels=scaleFUN)+
    # show line y=-x and add scales for guiding readers to know the magnitude of beta diversity 
    #geom_point(data = data.frame(x = seq(label.position$min_value, label.position$max_value, by = 0.01), y = seq(- (label.position$min_value), - (label.position$max_value), by = -0.01)),
    #           aes(x = x, y = y), alpha=0.5,  size = 1) +  # Ticks on the diagonal line
      scale_color_manual(values = concept_colour)+
    annotate(geom = "text", x = label.position$min_value/2, y = 0, angle = 45, label = "Differentiation", size = 5, color="grey68") +  
    annotate(geom = "text", x = 0, y = label.position$min_value/2, angle = 45, label = "Homogenization", size = 5, color="grey68") +
    annotate(geom = "text", x = 0.83*label.position$max_value, y=0.25*label.position$max_value, angle = 0, label = "Gain of\n widespread\n species", size = 4, color="grey78") +
    annotate(geom = "text",  x = 0.85*label.position$max_value, y =  0.85*label.position$min_value, angle = 0, label = "Localized\nreplaced by \nwidespread\n species", size = 4, color="grey78") +
    annotate(geom = "text", x = 0.15*label.position$min_value, y = 0.90*label.position$min_value, angle = 0, label = "Loss of\n localized\n species", size = 4, color="grey78") +
    annotate(geom = "text", x = 0.92*label.position$min_value, y = 0.15*label.position$min_value, angle = 0, label = "Loss of\n widespread\n species", size = 4, color="grey78") +
    annotate(geom = "text", x =0.9*label.position$min_value, y = 0.82*label.position$max_value, angle = 0, label = "Widespread\nreplaced by \nlocalized\n species", size = 4, color="grey78") +
    annotate(geom = "text", x = 0.3*label.position$max_value,  y = 0.9*label.position$max_value,  angle = 0, label = "Gain of\n localized\n species", size = 4, color="grey78") +
    theme(legend.position = "none")+
    labs(x="Effects of nutrient addition on alpha diversity (LRR)", y="Effects of nutrient addition on gamma diversity (LRR)")
}

# without annotation
pp.diversity.simple<-function(data1, data2, data3){
  concept_colour = c(  "Gain widespread species" = '#D9D956',
                       "Localized species replaced by widespread species"  = '#E9AE27',
                       "Loss of widespread species" = '#155F49',
                       "Loss of localized species" = '#D17538',
                       "Widespread species replaced by localized species"  ='#2CABA4', 
                       "Gain localized species" = '#61CDE0',
                        # for points that fall on the boundary
                       "Other situations" = '#f0f0f0')
  #Our transformation function
  scaleFUN <- function(x) sprintf("%.1f", x)
 
  # position for CI for change in beta diversity 
  angle.positive<- 135; angle.negative<- -45; dif.positive<-data3$dif.positive; dif.negative<-data3$dif.negative
  
  ggplot(data=data1)+ theme_cowplot()+panel_border()+
    geom_point(aes(x=true.alpha_mean, y=true.gamma_mean, color=cat.process), size=2, alpha=0.5)+
    geom_errorbar(aes(x=true.alpha_mean, y=true.gamma_mean, xmin=true.alpha_Q2.5, xmax=true.alpha_Q97.5, color=cat.process), width=0.0001, linewidth=0.2, alpha=0.15) +
    geom_errorbar(aes(x=true.alpha_mean, y=true.gamma_mean, ymin=true.gamma_Q2.5, ymax=true.gamma_Q97.5, color=cat.process), width=0.0001, linewidth=0.2, alpha=0.15) +
    geom_segment(aes(x=true.alpha_mean, y=true.gamma_mean, xend = true.alpha_mean + (change.beta_Q97.5 - change.beta_mean) * cos(angle.positive * pi / 180), yend = true.gamma_mean + (change.beta_Q97.5 - change.beta_mean) * sin(angle.positive * pi / 180), color=cat.process), lineend = "square",  linewidth=0.2, alpha=0.15) +
    geom_segment(aes(x=true.alpha_mean, y=true.gamma_mean, xend = true.alpha_mean + abs(change.beta_Q2.5 - change.beta_mean) * cos(angle.negative * pi / 180), yend = true.gamma_mean + abs(change.beta_Q2.5 - change.beta_mean) * sin(angle.negative * pi / 180), color=cat.process),  lineend = "square",  linewidth=0.2, alpha=0.15) +
    
    geom_errorbar(data=data2, aes(x=Estimate_alpha, y=Estimate_gamma, xmin=Q2.5_alpha, xmax=Q97.5_alpha), width=0.0001, color="black", linewidth=1, alpha=0.7) +
    geom_errorbar(data=data2,  aes(x=Estimate_alpha, y=Estimate_gamma, ymin=Q2.5_gamma, ymax=Q97.5_gamma), width=0.0001, color="black", linewidth=1, alpha=0.7) +
    geom_segment(data=data2,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.positive * cos(angle.positive * pi / 180), yend = Estimate_gamma + dif.positive * sin(angle.positive * pi / 180)), lineend = "square", color="black", linewidth=1, alpha=0.7) +
    geom_segment(data=data2,  aes(x=Estimate_alpha, y=Estimate_gamma, xend = Estimate_alpha + dif.negative * cos(angle.negative * pi / 180), yend = Estimate_gamma + dif.negative * sin(angle.negative * pi / 180)), lineend = "square", color="black", linewidth=1, alpha=0.7) +
    geom_point(data=data2,  aes(x=Estimate_alpha, y=Estimate_gamma), pch=21, size=5, stroke = 1,  alpha=0.7, color="black")+
    
    geom_abline(slope=1, intercept=0, linetype="dotted")+
    geom_hline(yintercept = 0, linetype="dotted")+
    geom_vline(xintercept = 0, linetype="dotted")+
    scale_x_continuous(limits = c(label.position$min_value, label.position$max_value), labels=scaleFUN)+
    scale_y_continuous(limits = c(label.position$min_value, label.position$max_value), labels=scaleFUN)+
    scale_color_manual(values = concept_colour)+
    theme(legend.position = "none")+labs(x=NULL, y=NULL)
}

# the end 





