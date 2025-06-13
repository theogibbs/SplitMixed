library(dplyr)
library(ggplot2)
library(brms)
library(tidyverse)
library(tidybayes)

chem_data <- read.csv("NC_data.csv")
plots <- read.csv("plots.csv")

chem_data$species <- substr(chem_data$ID, 1, 2)
chem_data$plot <- substr(chem_data$ID, 3, nchar(chem_data$ID))


chem_data$sp1 <- plots$sp1[match(chem_data$plot, plots$plot)]
chem_data$sp2 <- plots$sp2[match(chem_data$plot, plots$plot)]
chem_data$otherspecies <- ifelse(chem_data$species == chem_data$sp1, chem_data$sp2, chem_data$sp1) # other species is the background competitor that wasn't measured in that datapoint
chem_data$otherspecies <- ifelse(chem_data$species == "FE", substr(chem_data$plot, 2, 3), chem_data$otherspecies)

chem_data$treatment <- plots$treatment[match(chem_data$plot, plots$plot)]
chem_data$treatment <- ifelse(chem_data$species == "FE", substr(chem_data$plot, 4, nchar(chem_data$plot) - 1), chem_data$treatment)

chem_data$grp <- ifelse(chem_data$treatment == "split",  "clustered", paste("mixed", chem_data$otherspecies))
chem_data$grp <- as.factor(chem_data$grp)

chem_data$treatment <- as.factor(chem_data$treatment)
chem_data <- chem_data %>% mutate(treatment = forcats::fct_recode(treatment,
                                                            "clustered" = "split",
                                                            "mixed" = "mix"))
chem_data$treatment <- factor(chem_data$treatment, levels = c("clustered", "mixed"))

draws_treatment <- data.frame()
draws_grp <- data.frame()
for (sp in c("AC","FE", "PL", "SA", "UR")) {
  tmp <- chem_data %>% filter(species == sp)
  
  for (trait in c("N_perc", "C_perc", "C.N", "d15N", "d13C")) {  # add phenology?
    
    formula <- as.formula(paste(trait, "~ treatment"))
    mod_treatment <- brm(formula,
                         data = tmp, backend = "cmdstanr",
                         family = gaussian, iter = 100000,
                         control = list(adapt_delta = 0.99))
    
    cur_draws <- mod_treatment %>%
      spread_draws(b_Intercept,
                   b_treatmentmixed)
    cur_draws$species <- sp
    cur_draws$trait <- trait
    
    draws_treatment <- rbind(draws_treatment, cur_draws)
    
    formula <- as.formula(paste(trait, "~ grp"))
    mod <- brm(formula,
               data = tmp, backend = "cmdstanr",
               family = gaussian, iter = 100000,
               control = list(adapt_delta = 0.99))
    
    if (sp == "AC") {
      cur_draws <- mod %>% spread_draws(b_Intercept,
                                        b_grpmixedFE, b_grpmixedPL, b_grpmixedSA, b_grpmixedUR)
      cur_draws$b_grpmixedAC <- NA
    } else if (sp == "FE") {
      cur_draws <- mod %>% spread_draws(b_Intercept,
                                        b_grpmixedAC, b_grpmixedPL, b_grpmixedSA, b_grpmixedUR)
      cur_draws$b_grpmixedFE <- NA
    } else if (sp == "PL") {
      cur_draws <- mod %>% spread_draws(b_Intercept,
                                        b_grpmixedAC, b_grpmixedFE, b_grpmixedSA, b_grpmixedUR)
      cur_draws$b_grpmixedPL <- NA
    } else if (sp == "SA") {
      cur_draws <- mod %>% spread_draws(b_Intercept,
                                        b_grpmixedAC, b_grpmixedFE, b_grpmixedPL, b_grpmixedUR)
      cur_draws$b_grpmixedSA <- NA
    } else if (sp == "UR") {
      cur_draws <- mod %>% spread_draws(b_Intercept,
                                        b_grpmixedAC, b_grpmixedFE, b_grpmixedPL, b_grpmixedSA)
      cur_draws$b_grpmixedUR <- NA
    }
    
    cur_draws$species <- sp
    cur_draws$trait <- trait
    
    draws_grp <- rbind(draws_grp, cur_draws)
  }
}

clustered_summ <- draws_treatment %>%
  group_by(species, trait) %>%
  point_interval(b_Intercept, .width = 0.89) %>%
  rename(value = b_Intercept,
         lower = .lower,
         upper= .upper) %>%
  mutate(treatment = "clustered")

mixed_summ <- draws_treatment %>%
  group_by(species, trait) %>%
  point_interval(b_treatmentmixed, .width = 0.89) %>%
  rename(value = b_treatmentmixed,
         lower = .lower,
         upper = .upper) 

mixed_summ$lowerdiff <- mixed_summ$value - mixed_summ$lower
mixed_summ$upperdiff <- mixed_summ$upper - mixed_summ$value

mixed_summ2 <- mixed_summ

mixed_summ2$value <- mixed_summ$value + clustered_summ$value
mixed_summ2$lower <- mixed_summ2$value - mixed_summ2$lowerdiff
mixed_summ2$upper <- mixed_summ2$value + mixed_summ2$upperdiff
mixed_summ2$treatment <- "mixed"

trait_draws_summ <- bind_rows(clustered_summ, mixed_summ2)


sp_labels = c("AC" = "<i>A. wrangelianus</i>", "FE" = "<i>F. microstachys</i>", 
              "PL" = "<i>P. erecta</i>", "SA" = "<i>S. columbariae</i>", "UR" = "<i>U. lindleyi</i>")

y_labels <- c(N_perc = "foliar % N", C_perc = "foliar % C", 
              C.N = "C:N", d15N = "δ15N", d13C = "δ13C")

chem_data_long <- chem_data %>% select(!c(ID, WT, N_mass, C_mass, sp1, sp2)) %>% pivot_longer( !c(plot, species, otherspecies, treatment, grp), names_to = "trait", values_to = "value")
chem_data_long$trait <- factor(chem_data_long$trait, levels = c("C_perc", "N_perc", "C.N", "d13C", "d15N"))
trait_draws_summ$trait <- factor(trait_draws_summ$trait, levels = c("C_perc", "N_perc", "C.N", "d13C", "d15N"))

ggplot() +  
  geom_point(data = chem_data_long, aes(y = value, x = treatment, color = treatment), 
             position = "jitter", alpha = 0.2) +
  geom_point(data = trait_draws_summ, aes(y = value, x = treatment)) +
  geom_errorbar(data = trait_draws_summ, aes(ymin = lower, ymax = upper, x = treatment), width = 0.2) +
  theme_classic() +
  ggh4x::facet_grid2(trait ~ species, scales = "free", independent = "y", labeller = labeller(trait = y_labels, species = sp_labels),) + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = ggtext::element_markdown(size = 15, face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values=c("#332288", "#117733"))

ggsave("Figures/Supplementary/SuppChemTraits_splitmixed_rawdata.jpeg", width = 15, height = 10, units = "in")

ggplot() +  
  geom_point(data = trait_draws_summ, aes(y = value, x = treatment, color = treatment)) +
  geom_errorbar(data = trait_draws_summ, aes(ymin = lower, ymax = upper, x = treatment, color = treatment), width = 0.2) +
  theme_classic() +
  ggh4x::facet_grid2(trait ~ species, scales = "free", independent = "y", labeller = labeller(trait = y_labels, species = sp_labels),) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = ggtext::element_markdown(size = 15, face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values=c("#332288", "#117733"))

ggsave("Figures/Supplementary/SuppChemTraits_splitmixed.jpeg", width = 15, height = 10, units = "in")

clustered_summ_grp <- draws_grp %>%
  group_by(species, trait) %>%
  point_interval(b_Intercept, .width = 0.89) %>%
  rename(value = b_Intercept,
         lower = .lower,
         upper= .upper) %>%
  mutate(grp = "clustered")

mixed_summ_grp <- draws_grp %>%
  pivot_longer(cols = c(b_grpmixedAC, b_grpmixedFE, b_grpmixedPL, 
                        b_grpmixedSA, b_grpmixedUR),
               names_to = "grp",
               values_to = "estimate") %>%
  group_by(species, trait, grp) %>%
  point_interval(estimate, .width = 0.89) %>%
  rename(lower = .lower,
         upper= .upper)

mixed_summ_grp$lowerdiff <- mixed_summ_grp$estimate - mixed_summ_grp$lower
mixed_summ_grp$upperdiff <- mixed_summ_grp$upper - mixed_summ_grp$estimate

mixed_summ_grp2 <- mixed_summ_grp
mixed_summ_grp2 <- mixed_summ_grp %>%
  left_join(clustered_summ_grp %>% select(species, trait, value), by = c("species", "trait")) %>%
  mutate(value = estimate + value) %>% select(-estimate)

mixed_summ_grp2$lower <- mixed_summ_grp2$value - mixed_summ_grp2$lowerdiff
mixed_summ_grp2$upper <- mixed_summ_grp2$value + mixed_summ_grp2$upperdiff

trait_draws_grp_summ <- bind_rows(clustered_summ_grp, mixed_summ_grp2)

trait_draws_grp_summ$grp <- as.factor(trait_draws_grp_summ$grp)

levels(trait_draws_grp_summ$grp) <- c("mixed AC", "mixed FE",
                                      "mixed PL", "mixed SA",
                                      "mixed UR", "clustered")

levels(chem_data$grp) <- c("clustered", "mixed AC", "mixed FE",
                        "mixed PL", "mixed SA", "mixed UR")

trait_draws_grp_summ <- trait_draws_grp_summ %>% select(!c(lowerdiff, upperdiff))
trait_draws_grp_summ <- na.omit(trait_draws_grp_summ)

trait_draws_grp_summ$trait <- factor(trait_draws_grp_summ$trait, levels = c("C_perc", "N_perc", "C.N", "d13C", "d15N"))

ggplot() +  
  geom_point(data = chem_data_long, aes(y = value, x = grp, color = treatment), 
                  position = "jitter", alpha = 0.2) +
  geom_point(data = trait_draws_grp_summ, aes(y = value, x = grp)) +
  geom_errorbar(data = trait_draws_grp_summ, aes(ymin = lower, ymax = upper, x = grp), width = 0.2) +
  theme_classic() +
  ggh4x::facet_grid2(trait ~ species, scales = "free", independent = "y", labeller = labeller(trait = y_labels, species = sp_labels),) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = ggtext::element_markdown(size = 15, face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = ggtext::element_markdown(angle = 45, hjust = 0.95, vjust = 1)) +
  scale_color_manual(values=c("#332288", "#117733")) +
  scale_x_discrete(labels = c("mixed AC" = 'mixed with<br><i>A. wrangelianus</i>',
                              "mixed FE" = 'mixed with<br><i>F. michrostachys</i>',
                              "mixed PL" = 'mixed with<br><i>P. erecta</i>',
                              "mixed SA" = 'mixed with<br><i>S. columbariae</i>',
                              "mixed UR" = 'mixed with<br><i>U. lindleyi</i>'
  ))

ggsave("Figures/Supplementary/SuppChemTraits_spcombo_rawdata.jpeg", width = 15, height = 10, units = "in")


trait_draws_grp_summ$treatment <- ifelse(trait_draws_grp_summ$grp == "clustered", "clustered", "mixed")
trait_draws_grp_summ$grp <- factor(trait_draws_grp_summ$grp, levels = c("clustered", "mixed AC", "mixed FE",
                                                                        "mixed PL", "mixed SA",
                                                                        "mixed UR"))

ggplot() +  
  geom_point(data = trait_draws_grp_summ, aes(y = value, x = grp, color = treatment)) +
  geom_errorbar(data = trait_draws_grp_summ, aes(ymin = lower, ymax = upper, x = grp, color = treatment), width = 0.2) +
  theme_classic() +
  ggh4x::facet_grid2(trait ~ species, scales = "free", independent = "y", labeller = labeller(trait = y_labels, species = sp_labels),) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = ggtext::element_markdown(size = 15, face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = ggtext::element_markdown(angle = 45, hjust = 0.95, vjust = 1)) +
  scale_color_manual(values=c("#332288", "#117733")) +
  scale_x_discrete(labels = c("mixed AC" = 'mixed with<br><i>A. wrangelianus</i>',
                              "mixed FE" = 'mixed with<br><i>F. michrostachys</i>',
                              "mixed PL" = 'mixed with<br><i>P. erecta</i>',
                              "mixed SA" = 'mixed with<br><i>S. columbariae</i>',
                              "mixed UR" = 'mixed with<br><i>U. lindleyi</i>'
  ))

ggsave("Figures/Supplementary/SuppChemTraits_splcombo.jpeg", width = 15, height = 10, units = "in")




## Canopy index ## 

