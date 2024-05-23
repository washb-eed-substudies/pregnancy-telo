rm(list=ls())
source(here::here("0-config.R"))
library(cowplot)
library(patchwork)

#load spline data
H1_spline <- readRDS(here("figure-data/H1_adj_spline_data.RDS"))

#load results for quartiles
H1_quartiles <- readRDS(here("results/adjusted/H1_adj_res.RDS"))

d_for_plot <- function(x_name, y_name, x_var, y_var, spline, quart){
  d <- NULL
  for (i in 1: length(x_var)) {
    for (j in 1:length(y_var)){
      exists <- (quart%>%filter(X==x_var[i], Y==y_var[j]) %>% nrow()) != 0
      if (exists){
        new <- data.frame(x=x_name[i], y=y_name[j], quart%>%filter(X==x_var[i], Y==y_var[j]))
        d <- rbind(d, new)
      }
    }
  }
  d
}


d1 <- d_for_plot(c("Vitamin A deficiency"),
                 c("Telomere length Year 2"), 
                 c("vit_A_def"),   
                 c("TS_t3_Z"),
                 H1_spline, H1_quartiles) #%>% left_join(colors,"x")

d1$x <- factor(d1$x,levels=c("Vitamin A deficiency"))
#d1$x <- factor(d1$x, levels=rev(levels(d1$x)))
d1$y <- factor(d1$y)
d1$`Time of outcome measurement` <- factor(ifelse(grepl("t3", d1$Y), "Age 28 months", "Age 14 months"))



ggplot(d1) + geom_col(aes(x=0.5, y=pred.q1)) + 
  geom_col(aes(x=1.5, y=pred.q3)) + 
  scale_y_continuous("Predicted child telomere length\nat Year 2 (Z-score)") + 
  scale_x_continuous("Maternal vitamin A", breaks=c(0.5,1.5), label=c("Normal", "Deficient")) + 
  scale_color_manual(values = tableau10[c(7:1)]) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        axis.title.x = element_text(size=11),
        strip.text = element_text(hjust=0.5, size=10),
        strip.placement = "outside",
        axis.text.y = element_text(hjust = 1, size=9),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "bottom")      


p %>% ggsave(filename="figures/pregnancy-telo_point_diff_micronutrients.jpg", width=10, height=7)

d_for_plot <- function(x_name, y_name, x_var, y_var, spline, quart){
  d <- NULL
  for (i in 1: length(x_var)) {
    for (j in 1:length(y_var)){
      exists <- (quart%>%filter(X==x_var[i], Y==y_var[j]) %>% nrow()) != 0
      if (exists){
        new <- data.frame(x=x_name[i], y=y_name[j], spline%>%filter(Xvar==x_var[i], Yvar==y_var[j]), quart%>%filter(X==x_var[i], Y==y_var[j])%>%select(q1, q3))
        d <- rbind(d, new)
      }
    }
  }
  d
}

d1 <- d_for_plot(c("Vitamin A deficiency"),
                 c("Telomere length Year 2"), 
                 c("vit_A_def"),   
                 c("TS_t3_Z"),
                 H1_spline, H1_quartiles) #%>% left_join(colors,"x")

d1$x <- factor(d1$x,levels=c("Vitamin D", "Vitamin D deficiency", "Ln RBP", "Vitamin A deficiency","Ln ferritin", "Ln sTfR", "Iron deficiency"))
#d1$x <- factor(d1$x, levels=rev(levels(d1$x)))
d1$y <- factor(d1$y)

# d1 <- d1 %>% filter(!grepl("def", Xvar))

d1 %>% filter(grepl("t3", Yvar)) %>% ggplot(aes(x=X))+
  geom_smooth(aes(y = fit, color=x), se = F) +
  geom_vline(aes(xintercept=q1), size=.5, color="grey30", linetype="dashed") +
  geom_vline(aes(xintercept=q3), size=.5, color="grey30", linetype="dashed") +
  #geom_point(aes(y=Y), alpha=0.5) +
  geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=x, color=x), alpha=0.5) + 
  scale_x_continuous("Vitamin A deficiency", breaks=c(0,1)) + 
  scale_y_continuous("Telomere length year 2") + 
  scale_colour_manual(values=tableau10[c(1,5,6,4)]) + 
  scale_fill_manual(values=tableau10[c(1,5,6,4)]) + 
  theme_ki() 

d1 <- d_for_plot(c("Vitamin D"),
                 c("Change in telomere length"), 
                 c("vitD_nmol_per_L"),   
                 c("delta_TS_Z"),
                 H1_spline, H1_quartiles) #%>% left_join(colors,"x")

d1$x <- factor(d1$x,levels=c("Vitamin D", "Vitamin D deficiency", "Ln RBP", "Vitamin A deficiency","Ln ferritin", "Ln sTfR", "Iron deficiency"))
#d1$x <- factor(d1$x, levels=rev(levels(d1$x)))
d1$y <- factor(d1$y)

# d1 <- d1 %>% filter(!grepl("def", Xvar))

d1 %>% ggplot(aes(x=X))+
  geom_smooth(aes(y = fit, color=x), se = F) +
  geom_vline(aes(xintercept=q1), size=.5, color="grey30", linetype="dashed") +
  geom_vline(aes(xintercept=q3), size=.5, color="grey30", linetype="dashed") +
  #geom_point(aes(y=Y), alpha=0.5) +
  geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=x, color=x), alpha=0.5) + 
  scale_x_continuous("Vitamin D (nmol/L)") + 
  scale_y_continuous("Change in telomere length\nbetweens years 1 and 2") + 
  scale_colour_manual(values=tableau10[c(1,5,6,4)]) + 
  scale_fill_manual(values=tableau10[c(1,5,6,4)]) + 
  theme_ki() 
