rm(list = ls())

library(ggplot2)
library(RColorBrewer)

###################################################
#load file

micropeptide_relation$cancer_type = factor(micropeptide_relation$cancer_type,levels = c("BRCA","ccRCC","COAD","LUAD","OV","PDAC","UCEC"))

P1 <- ggplot(micropeptide_relation, aes(x = cancer_type, y = rvalue, fill = cancer_type)) +
  geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width", width = 0.5) +
  #geom_boxplot(width = 0.1, fill = "white", outlier.fill = "white", outlier.color = "white", color = "black") +
  #geom_jitter(aes(fill = cancer_type), width = 0.2, shape = 21, size = 2, color = "gray") +
  scale_fill_manual(values = c("#ADD26A", "#F4AF63", "#7FABCB", "#EE7C6F", "#BCB7D6", "#F8F5B3", "#89CCC0")) +
  theme_bw() +
  coord_flip()+
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black", family = "Arial", size = 12),
        axis.text.y = element_text(family = "Arial", size = 12, face = "plain"),
        axis.title.y = element_text(family = "Arial", size = 12, face = "plain"),
        axis.title.x = element_text(family = "Arial", size = 12, face = "plain"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #geom_hline(yintercept = 0.3, linetype = "dashed", color = "red")

print(P1)

###################################################
library(dplyr)
library(tidyr)
library(ggdist)
library(ggplot2)
library(patchwork)
p1 = micropeptide_relation %>%
  ggplot(aes(y = cancer_type, x = rvalue, fill = cancer_type)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), 
            scale = 0.5) +
  stat_dotsinterval(side = "bottom", 
                    scale = 0.7, 
                    slab_linewidth = NA, 
                    slab_color = NA) +
  scale_fill_brewer(palette = "Set2") +
  labs(title="slab+dotsinterval", 
       fill = "cancer_type", 
       x = "rvalue", 
       y = "cancer_type") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12),
    axis.text = element_text(colour='black', size=12),
    legend.text = element_text(colour='black', size=12)
  )

p2 = micropeptide_relation %>%
  ggplot(aes(y = rvalue, x = cancer_type, fill = cancer_type)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), 
            scale = 0.6) +
  stat_dotsinterval(side = "left", 
                    scale = 0.6, 
                    slab_linewidth = NA, 
                    slab_color = NA) +
  #scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#ADD26A", "#F4AF63", "#7FABCB", "#EE7C6F", "#BCB7D6", "#F8F5B3", "#89CCC0"))+
  labs(x = "Disease", 
       y = "Pearson correlation r value") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12,family = "Arial"),
    axis.text = element_text(colour='black', size=12,family = "Arial"),
    legend.text = element_text(colour='black', size=12,family = "Arial")
  )

p2

