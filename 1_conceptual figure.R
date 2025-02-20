rm(list=ls())
## set up the work directory 
dir.graphs<-"C:/Users/chqq3/work/homogenization/Nature communications/graphs1/"
setwd(dir.graphs)

library(cowplot);library(scales);library(ggthemes)
## color needed 
colorblind_pal()(8)
## "#000000" "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2" "#D55E00" "#CC79A7"
show_col(colorblind_pal()(8))
concept_colour = c("Gain of widespread species" =  "#F0E442",
                   "Spatially restricted replaced by widespread species" = "#E69F00",
                   "Loss of spatially restricted species" = "#D55E00",
                   "Loss of widespread species" ="#009E73",
                   "Widespread replaced by spatially restricted species"  = "#56B4E9", 
                   "Gain of spatially restricted species" = "#0072B2", 
                   "Other situations" = '#f0f0f0')

min_value<- -10
max_value<- 10

# Create a data frame with polygons for each region
regions1<-data.frame(
  x=c(0,10,10,    # I
      0,10,10,0,
      0,-10,0,    # III
      0,-10,-10,
      0,-10,-10, 0,
      0,0,10),     # VI
  
  y=c(0,10,0,        # I
      0,0,-10,-10,
      0,-10,-10,    # III
      0,0,-10,
      0,0,10,10,
      0,10,10)      # VI
)
# Define the regions and colors
regions1$region <- factor(rep(1:6, c(3,4,3,3,4,3)))  # Six regions, each with four coordinates
regions1$Scenarios <- factor(rep(c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                   "Loss of widespread species", "Widespread replaced by spatially restricted species", "Gain of spatially restricted species"), c(3,4,3,3,4,3)))  #
regions1$Scenarios<-factor(regions1$Scenarios,levels=c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                                       "Loss of widespread species" ,  "Widespread replaced by spatially restricted species", "Gain of spatially restricted species"))

pp<-ggplot() +theme_cowplot(font_size = 11)+  # Apply cowplot theme
  theme(panel.border = element_blank(),   # Remove panel border
        axis.line = element_blank() ,
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        legend.text = element_text(size = 9) )+
  geom_polygon(data = regions1, aes(x = x, y = y, fill = Scenarios), color = "black", alpha = 0.5) +
  guides(fill = guide_legend(nrow = 3)) +  # Specify two rows for the legend
  geom_segment(aes(x = -10, y = -10, xend = 10, yend = 10), linewidth=1, color = "white") + 
  scale_fill_manual(values = concept_colour) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) + 
  scale_x_continuous(limits = c(min_value, max_value))+
  scale_y_continuous(limits = c(min_value, max_value))+
  annotate(geom = "text", x = min_value/1.5, y = max_value/1.5, angle = 45, label = "Differentiation", color="white", size =6, fontface = "bold") +  
  annotate(geom = "text", x = min_value/1.65, y = max_value/1.65, angle = 45, label = "Increase in ∆β (LRR)", color="white", size =5) +  
  annotate(geom = "text", x = max_value/1.5, y = min_value/1.5, angle = 45, label = "Homogenization", color="white", size =6, fontface = "bold") +  
  annotate(geom = "text", x = max_value/1.35, y = min_value/1.35, angle = 45, label = "Decrease in ∆β (LRR)", color="white", size =5) + 
  annotate("label", x = 0, y = 0, angle = 45, label = "No change in beta diversity; ∆β (LRR)", fill="white", size =6, fontface = "bold") +
  
  annotate(geom = "text", x = 0.40*max_value, y=0.2*max_value, angle = 0, label = "I", size =5) +
  annotate(geom = "text",  x = 0.10*max_value, y =  0.2*min_value, angle = 0, label = "II", size =5) +
  annotate(geom = "text", x = 0.6*min_value, y = 0.80*min_value, angle = 0, label = "III", size =5) +
  annotate(geom = "text", x = 0.9*min_value, y = 0.2*min_value, angle = 0, label = "IV", size =5) +
  annotate(geom = "text", x =0.9*min_value, y = 0.2*max_value, angle = 0, label = "V", size =5) +
  annotate(geom = "text", x = 0.1*max_value,  y = 0.8*max_value,  angle = 0, label = "VI", size =5) +
  annotate(geom = "text", x = 0.60*max_value, y=0.2*max_value, angle = 0, label = "Gain of\n widespread\n species", size =5) +
  annotate(geom = "text",  x = 0.55*max_value, y =  0.2*min_value, angle = 0, label = "Spatially restricted replaced \n by widespread species", size =5) +
  annotate(geom = "text", x = 0.3*min_value, y = 0.80*min_value, angle = 0, label = "Loss of spatially \nrestricted species", size =5) +
  annotate(geom = "text", x = 0.7*min_value, y = 0.2*min_value, angle = 0, label = "Loss of\n widespread\n species", size =5) +
  annotate(geom = "text", x =0.5*min_value, y = 0.2*max_value, angle = 0, label = "Widespread replaced by\n spatially restricted species", size =5) +
  annotate(geom = "text", x = 0.4*max_value,  y = 0.8*max_value,  angle = 0, label = "Gain of spatially\n restricted species", size =5) +
  labs( x = expression("-          Change in alpha diversity; " ~ bar("∆α") ~ " (LRR)         +"),  # Combine overbar and other text
       y="-          Change in gamma diversity; ∆γ (LRR)          +",  fill=NULL) 
(pp1<-pp+theme(legend.position = "none"))

ggsave(pp1, width = 18, height = 18, dpi=600, unit="cm",  file="conceptual figure.png")
ggsave(pp1, width = 18, height = 18, unit="cm",  file="conceptual figure.pdf", device = cairo_pdf)


