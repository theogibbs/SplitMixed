library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(reshape2)
library(ggpubr)


bh_loo <- read.csv("BevertonHoltLOO.csv")
lw_loo <- read.csv("LawWatkinsonLOO.csv")
ri_loo <- read.csv("RickerLOO.csv")

loo_comp <- rbind(bh_loo, lw_loo, ri_loo)

loo_comp$Focal <- paste("Focal:", loo_comp$Focal)

plModelComp <- ggplot(loo_comp, aes(x=Model, y=Estimate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0.2) +
  theme_classic() +
  facet_wrap(~ Focal, scales = "free", nrow = 1) +
  labs(x = "Model", y = "") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank())
plModelComp

jpeg("./figs/SIFigModelComp.jpeg",
     width = 3000, height = 1200, res = 300)
plModelComp
dev.off()