
options(scipen = 999)
library(scales)  # For formatting functions

concept_colour = c("Gain of widespread species" =  "#F0E442",
                   "Spatially restricted replaced by widespread species" = "#E69F00",
                   "Loss of spatially restricted species" = "#D55E00",
                   "Loss of widespread species" ="#009E73",
                   "Widespread replaced by spatially restricted species"  = "#56B4E9", 
                   "Gain of spatially restricted species" = "#0072B2", 
                   "Other situations" = '#f0f0f0')

# function to draw background for figures 
plot.six.scenarios<-function(xx, yy){
  min_value<- ifelse (xx > - 0.5, -0.5, xx)
  max_value<- ifelse(yy>0.5, yy, 0.5)
  
  
  regions1<-data.frame(
    x=c(0,max_value,max_value,    # I
        0,max_value,max_value,0,
        0,min_value,0,    # III
        0,min_value,min_value,
        0,min_value,min_value, 0,
        0,0,max_value),     # IV
    
    y=c(0,max_value,0,        # I
        0,0,min_value,min_value,
        0,min_value,min_value,    # III
        0,0,min_value,
        0,0,max_value,max_value,
        0,max_value,max_value)      # IV
  )
  # Define the regions and colors
  regions1$Scenarios <- factor(rep(c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                     "Loss of widespread species", "Widespread replaced by spatially restricted species", "Gain of spatially restricted species"), c(3,4,3,3,4,3)))  #
  regions1$Scenarios<-factor(regions1$Scenarios, levels=c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                                          "Loss of widespread species" ,  "Widespread replaced by spatially restricted species", "Gain of spatially restricted species"))
  
  pp<-ggplot() +theme_cowplot(font_size = 11)+ theme(panel.border = element_blank(), legend.position = "none")+
    geom_polygon(data = regions1, aes(x = x, y = y, fill = Scenarios), color = "black", alpha = 0.5) +
    scale_fill_manual(values =concept_colour ) +
    scale_x_continuous(limits = c(min_value, max_value), labels = number_format(accuracy = 0.1))+
    scale_y_continuous(limits = c(min_value, max_value), labels = number_format(accuracy = 0.1))+
    #annotate(geom = "text", x = min_value/1.5, y = max_value/1.5, angle = 45, label = "Differentiation\n (increase in ∆β)", color="white", size =5, fontface = "bold") +  
    #annotate(geom = "text", x = max_value/1.5, y = min_value/1.5, angle = 45, label = "Homogenization\n (decrease in ∆β)", color="white", size =5, fontface = "bold") +
    geom_segment(aes(x = min_value, y = min_value, xend = max_value, yend = max_value), linewidth= 1 , color="white")+
    # annotate("label", x = 0, y = 0, angle = 45, label = "No change in beta diversity (∆β)", fill="white", size =5, fontface = "bold") +
    labs( x = expression(bar("∆α") ~ " (LRR)"),  # Combine overbar and other text
         y="∆γ (LRR)")
  pp
}
##########################################################################################
###### produce a legend 
###########################################################################################

# Create a simple dataframe with 6 groups
set.seed(123)  # For reproducibility

df <- data.frame(
  Group = rep(LETTERS[1:6], each = 5),  # 6 groups (A to F), each with 5 observations
  Value = rnorm(30, mean = 50, sd = 10)  # Random values
)
# Define the regions and colors
df$Scenarios <- factor(rep(c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                   "Loss of widespread species", "Widespread replaced by spatially restricted species", "Gain of spatially restricted species"), each=5))  #
df$Scenarios<-factor(df$Scenarios,levels=c("Gain of widespread species", "Spatially restricted replaced by widespread species", "Loss of spatially restricted species",
                                                       "Loss of widespread species" ,  "Widespread replaced by spatially restricted species", "Gain of spatially restricted species"))

pp.for.legend<- ggplot(data=df, aes(Group, Value, fill=Scenarios)) +theme_cowplot(font_size = 10)+
  geom_point(pch=22, alpha = 0.5)+
  scale_fill_manual(values =concept_colour )+
  guides(fill = guide_legend(nrow = 3, override.aes = list(size = 5))) +
  theme(legend.text = element_text(margin = margin(l = -0.2, unit = "pt")), legend.spacing = unit(5, "mm"))+
  labs(fill=NULL)
leg<-get_legend(pp.for.legend)

# Create a blank spacer plot
spacer <- ggplot() + 
  theme_void() + 
  theme(plot.margin = margin(0, 0, 0, 0))  # Adjust left margin for alignment

leg.space<-plot_grid(spacer, leg, nrow = 1, rel_widths = c(0.08, 0.92))

