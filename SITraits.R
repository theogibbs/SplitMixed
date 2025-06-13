library(dplyr)
library(ggplot2)
library(brms)
library(tidyverse)
library(tidybayes)

traits <- read.csv("traits.csv")
plots <- read.csv("plots.csv")

traits$sp1 <- plots$sp1[match(traits$plot, plots$plot)]
traits$sp2 <- plots$sp2[match(traits$plot, plots$plot)]
traits$otherspecies <- ifelse(traits$species == traits$sp1, traits$sp2, traits$sp1) # other species is the background competitor that wasn't measured in that datapoint

traits$treatment <- plots$treatment[match(traits$plot, plots$plot)]
traits$grp <- ifelse(traits$treatment == "split",  "clustered", paste("mixed", traits$otherspecies))
traits$grp <- as.factor(traits$grp)

traits$treatment <- as.factor(traits$treatment)
traits <- traits %>% mutate(treatment = forcats::fct_recode(treatment,
                                "clustered" = "split",
                                "mixed" = "mix"))
traits$treatment <- factor(traits$treatment, levels = c("clustered", "mixed"))

traits$canopy <- ((traits$axis1 + traits$axis2)/2)/(traits$height)
traits$LDMC <- traits$LDMC * 100 # to convert from g/g to mg/g


# UR heights were measured end of season, URLI heights were measured before maximum
traits[traits$species == "URLI",]$height <- NA
traits[traits$species == "URLI",]$species <- "UR"

draws_treatment <- data.frame()
draws_grp <- data.frame()
for (sp in c("AC", "FE", "PL", "SA", "UR")) {
  tmp <- traits %>% filter(species == sp)
  
  for (trait in c("height", "biomass", "SLA", "LDMC")) {  # add phenology?
    
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

trait <- "canopy"
for (sp in c("AC", "PL")) {
  tmp <- traits %>% filter(species == sp)
    
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
    } else if (sp == "PL") {
      cur_draws <- mod %>% spread_draws(b_Intercept,
                                        b_grpmixedAC, b_grpmixedFE, b_grpmixedSA, b_grpmixedUR)
      cur_draws$b_grpmixedPL <- NA
    }
    
    cur_draws$species <- sp
    cur_draws$trait <- trait
    
    draws_grp <- rbind(draws_grp, cur_draws)
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

y_labels <- c(height = "height (cm)", biomass = "biomass (g)", 
              SLA = "SLA (cm²/g)", LDMC = "LDMC (mg/g)")

traits_long <- traits %>% select(c(species, otherspecies, plot, treatment, grp, height, biomass, SLA, LDMC)) %>% pivot_longer( !c(plot, species, otherspecies, treatment, grp), names_to = "trait", values_to = "value")


ggplot() +  
  geom_point(data = traits_long, aes(y = value, x = treatment, color = treatment), 
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

# ggsave("Figures/Supplementary/SuppTraits_splitmixed_rawdata.jpeg", width = 15, height = 10, units = "in")

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

# ggsave("Figures/Supplementary/SuppTraits_splitmixed.jpeg", width = 15, height = 10, units = "in")



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

levels(traits$grp) <- c("clustered", "mixed AC", "mixed FE",
                        "mixed PL", "mixed SA", "mixed UR")

trait_draws_grp_summ <- trait_draws_grp_summ %>% select(!c(lowerdiff, upperdiff))
trait_draws_grp_summ <- na.omit(trait_draws_grp_summ)

ggplot() +  
  geom_point(data = traits_long, aes(y = value, x = grp, color = treatment), 
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

# ggsave("Figures/Supplementary/SuppTraits_spcombo_rawdata.jpeg", width = 15, height = 10, units = "in")


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

# ggsave("Figures/SuppTraits.jpeg", width = 15, height = 10, units = "in")


## Canopy index ##

ggplot() +  
  geom_point(data = trait_draws_summ %>% filter(trait == "canopy"), aes(y = value, x = treatment, color = treatment)) +
  geom_errorbar(data = trait_draws_summ %>% filter(trait == "canopy"), aes(ymin = lower, ymax = upper, x = treatment, color = treatment), width = 0.2) +
  theme_classic() +
  facet_wrap(species ~ ., scales = "free") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = ggtext::element_markdown(size = 15, face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values=c("#332288", "#117733"))

