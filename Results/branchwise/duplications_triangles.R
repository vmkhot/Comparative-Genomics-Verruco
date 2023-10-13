library(ggplot2)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(colorRamps)

setwd("~/GitHub/Comparative Genomics/Results/target_nodes/")
df <- read.csv("consolidated.csv",sep=",",header=TRUE)

df_m <- melt(df,id=c("Node","Duplications","Loss"))


brks <- c(1,10,100,400)

p1 <-ggplot(df_m,aes(y=Node,x=variable)) +
  geom_point(aes(size=value),shape=17,alpha=0.7)+
  scale_size_continuous(breaks=brks,limits=c(1,400),range = c(0, 20))+
  #scale_x_discrete(position="top")+
  #scale_fill_brewer(palette = "Set3")+
  guides(fill = guide_legend(override.aes = list(size = 5)))+
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(colour="black", vjust=0.5,size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.title=element_text(colour="black",face="bold",size=14),
        strip.text = element_text(colour="black",size=14),
        legend.text = element_text(size = 14,colour ="black"), 
        legend.title = element_text(size = 16, face = "bold"),
        legend.position = "right",
        panel.background = element_blank()) 
        #panel.border = element_rect(colour = "grey95", fill = NA, size = 0.1),
        #panel.grid.major.y = element_line(colour = "grey95"),
        #panel.grid.major.x = element_line(colour="#fafafaff"))
p1
