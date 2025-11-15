library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(reshape2)
library(ggpubr)
library(ggfortify)
library(gridExtra)
library(grid)

# reading in seed calibration data
seed_head_cal <- read.csv("seed_head_calibration.csv")

# extracting AC seed calibration data
AC_seed_head <- seed_head_cal[,6]
AC_seed_head <- AC_seed_head[2:33]
AC_seed_head <- as.numeric(AC_seed_head)

# computing AC seed calibration statistics
AC_seed_per_head_mean <- mean(AC_seed_head)
AC_seed_per_head_sd   <- sd(AC_seed_head)

# histogram of AC seed calibration data
ggplot(data.frame(SeedsPerHead = AC_seed_head), aes(x = SeedsPerHead)) +
  geom_histogram(bins = 4, size = 2, fill = "white", color = "darkblue") +
  geom_vline(linetype = "dashed", linewidth = 2, color = "gray",
             xintercept = c(AC_seed_per_head_mean - AC_seed_per_head_sd,
                            AC_seed_per_head_mean,
                            AC_seed_per_head_mean + AC_seed_per_head_sd)) +
  ggtitle("AC seeds per head")

# extracting UR seed calibration data
UR_seed_head <- seed_head_cal[,1]
UR_seed_head <- UR_seed_head[2:31]
UR_seed_head <- as.numeric(UR_seed_head)

# computing UR seed calibration statistics
UR_seed_per_head_mean <- mean(UR_seed_head)
UR_seed_per_head_sd   <- sd(UR_seed_head)

# histogram of UR seed calibration data
ggplot(data.frame(SeedsPerHead = UR_seed_head), aes(x = SeedsPerHead)) +
  geom_histogram(size = 2, fill = "white", color = "darkblue") +
  geom_vline(linetype = "dashed", size = 2,, color = "gray",
             xintercept = c(UR_seed_per_head_mean - UR_seed_per_head_sd,
                            UR_seed_per_head_mean,
                            UR_seed_per_head_mean + UR_seed_per_head_sd)) +
  ggtitle("UR seeds per head")

# extracting SA seed calibration data, including measured diameters
SA_seed_head <- seed_head_cal[,c(3, 4)]
SA_seed_head <- SA_seed_head[-c(1, nrow(SA_seed_head)),]
colnames(SA_seed_head) <- c("Diameter_CM", "Seeds")
SA_seed_head$Diameter_CM <- as.numeric(SA_seed_head$Diameter_CM)
SA_seed_head$Seeds <- as.numeric(SA_seed_head$Seeds)

# plotting the relationship between seeds and diameter or diameter^2
ggplot(SA_seed_head, aes(x = Diameter_CM, y = Seeds)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm")

ggplot(SA_seed_head, aes(x = Diameter_CM^2, y = Seeds)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm")

# creating a new variable for diameter^2
SA_seed_head$Diameter_CM_Squared <- SA_seed_head$Diameter_CM^2

# linear regression between seeds and diameter^2
lm_obj <- summary(lm(Seeds ~ Diameter_CM_Squared, data = SA_seed_head))
lm_int <- lm_obj$coefficients[1,1]
lm_slope <- lm_obj$coefficients[2,1]

# function that takes in a diameter and returns a number of seeds
# based on an input regression
DiameterToSeeds_SA <- function(diam, lm_int, lm_slope) {
  ret <- diam^2 * lm_slope + lm_int
  return(ret)
}

# reading in fecundity data
fec_data <- read.csv("fecundity_data.csv")

# converting AC and UR heads to seeds
ur_ac_fec_data <- fec_data[,1:3]
ur_ac_fec_data$Uropappus <- UR_seed_per_head_mean * ur_ac_fec_data$Uropappus
ur_ac_fec_data$Acmispon <- AC_seed_per_head_mean * ur_ac_fec_data$Acmispon

# extracting SA data
sa_fec_data <- fec_data[,c(1, 4, 6:ncol(fec_data))]

melt_sa_data <- sa_fec_data %>%
  melt(id.vars = c("Plot")) %>% # reshaping on the basis of plot
  na.omit() %>% # removing NAs to only keep heads that are actually there
  mutate(Seeds = DiameterToSeeds_SA(value, lm_int = lm_int, lm_slope = lm_slope)) %>% # replacing heads with seeds
  group_by(Plot) %>% # grouping by plot
  dplyr::summarise(Salvia = sum(Seeds)) # adding up the number of seeds from heads on the same plant

# plots that are missing from the above data
# either because they are bare
# or because they are missing a focal
plots_without_salvia_data <- setdiff(500:814, melt_sa_data$Plot)
print(plots_without_salvia_data)

# reading in Salvia bare plot fecundity data
sa_bare_fec <- read.csv("salvia_bare_fecundity.csv")

# renaming columns to be the plot numbers
colnames(sa_bare_fec) <- sa_bare_fec[1,]

sa_bare_fec <- sa_bare_fec[-1,] %>% # removing plot numbers from data
  melt() %>% # reshaping data
  na.omit() %>% # removing missing data
  mutate(Seeds = DiameterToSeeds_SA(value, lm_int = lm_int, lm_slope = lm_slope)) %>% # converting to seeds
  group_by(variable) %>% # grouping by plot
  dplyr::summarise(Salvia = sum(Seeds)) # adding up the number of seeds from heads on the same plant

# renaming columns
colnames(sa_bare_fec) <- c("Plot", "Salvia")

# merging with the previous data from non bare plots
sa_fec_data <- rbind(melt_sa_data, sa_bare_fec) %>%
  arrange(Plot) # reordering by plot

# plots that are missing from the above data
# because they are missing a focal
plots_without_salvia_data <- setdiff(500:814, sa_fec_data$Plot)
print(plots_without_salvia_data)

missing_sa_data <- data.frame(Plot = plots_without_salvia_data, Salvia = NA)

# merging with the previous data from non bare plots
sa_fec_data <- rbind(sa_fec_data, missing_sa_data) %>%
  arrange(Plot) # reordering by plot

# checking to see if plot numbers are the same
print(prod(ur_ac_fec_data$Plot == sa_fec_data$Plot))

# putting the fecundity together for the difference focals together
fec_data <- ur_ac_fec_data
fec_data$Salvia <- sa_fec_data$Salvia

# reshaping data and renaming variables to species name abbreviations
fec_data <- fec_data %>%
  melt(id.vars = c("Plot")) %>%
  mutate(focal = case_when(variable == "Acmispon" ~ "AC",
                           variable == "Salvia" ~ "SA",
                           variable == "Uropappus" ~ "UR")) %>%
  select(-c("variable")) %>% # dropping old variables column
  arrange(Plot, focal) # sorting

# renaming columns
colnames(fec_data) <- c("plot", "Seeds", "focal")

# reading in census data for competitor densities
cen_data <- read.csv("official_census.csv")

# merging census and fecundity data
fec_data <- merge(cen_data, fec_data, by = c("plot", "focal"))
fec_data$Seeds <- as.numeric(fec_data$Seeds)
# making some t-test figures

# creating factors
fecundity <- fec_data
fecundity$treatment <- factor(fecundity$treatment, levels = c("split", "bare", "mix"))
fecundity$combo <- paste(fecundity$sp1, fecundity$sp2)

# for(cur_foc in c("AC", "SA", "UR")) {
#   cur_fec <- fecundity %>%
#     filter(focal == cur_foc) %>%
#     filter(treatment != "bare")
#   
#   show(ggboxplot(cur_fec, x = "combo", y = "Seeds", add = "jitter",
#                  color = "treatment", group="treatment") + 
#          stat_compare_means(aes(group=treatment), method="t.test", vjust = 1)  + 
#          xlab("competitor combination") +
#          ylab("Seeds") + ggtitle(paste(cur_foc, "focals")) + 
#          facet_wrap(~combo, scales = "free") + 
#          scale_color_hue(direction = -1))
#   
#   ## pooled across background competitor combinations
#   show(ggboxplot(cur_fec %>% filter(treatment != "bare"), x = "treatment", y = "Seeds", add = "jitter",
#                 color = "treatment") +
#         stat_compare_means(aes(group=treatment), method="t.test") +
#         ylab("Seeds") + ggtitle(cur_foc) +
#         scale_color_hue(direction = -1))
#   
# }


# removing data for missing focals
rem_bare_na_fec_data <- fecundity %>%
  filter(treatment != "bare") %>%
  mutate(sp1_count = ifelse(treatment == "bare", 0, sp1_count)) %>%
  mutate(sp2_count = ifelse(treatment == "bare", 0, sp2_count)) %>%
  filter(!(is.na(sp1_count & is.na(sp2_count)))) %>%
  filter(!is.na(Seeds))

sp_combos <- unique(rem_bare_na_fec_data$combo)
rem_bare_na_fec_data$treatment <- factor(rem_bare_na_fec_data$treatment, 
                                         levels = c("mix", "split"))

loo_data <- data.frame()
fec_draws <- data.frame()
for(cur_combo in sp_combos) {
  
  cur_fec <- rem_bare_na_fec_data %>%
    filter(combo == cur_combo) %>%
    arrange(treatment)
  
  cur_fit <- brm(formula = Seeds ~ 1 + treatment + focal 
                 + sp1_count + sp2_count
                 ,
                 data = cur_fec, backend = "cmdstanr", family = lognormal)
  
  cur_loo <- loo(cur_fit)$estimates[3, 1:2]
  cur_loo <- as.data.frame(t(cur_loo))
  cur_loo$SpCombo <- cur_combo
  loo_data <- rbind(loo_data, cur_loo)
  
  cur_draws <- cur_fit %>%
    spread_draws(b_Intercept,
                 b_treatmentsplit,
                 b_focalSA,
                 b_focalUR
                 #b_sp1_count,
                 #b_sp2_count
    )
  
  cur_draws$sp1 <- unique(cur_fec$sp1)
  cur_draws$sp2 <- unique(cur_fec$sp2)
  fec_draws <- rbind(fec_draws, cur_draws)
  
}

fec_draws$SpCombo <- paste(fec_draws$sp1, fec_draws$sp2)

# ggplot(fec_draws, aes(y = b_treatmentsplit, x = SpCombo, ,
#                         fill = after_stat(y < 0))) +
#   stat_slab(aes(alpha = (after_stat(level))), .width = c(.95, 1)) +
#   scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
#   geom_hline(yintercept = 0, alpha = 0.75, linetype = "dashed") +
#   theme_classic() +
#   theme(text = element_text(size=10),
#         legend.position = "none",
#         axis.text.x = element_text(size = 10),
#         legend.text=element_text(size = 15),
#         strip.background = element_blank(),
#         plot.title.position = "plot",
#         plot.caption.position =  "plot") +
#   labs(x = "Species Combination",
#        y = "Difference in Focal Fecundity\n in Clustered Relative to Mixed Competitors")


summ_stats <- fec_draws %>% group_by(SpCombo) %>% point_interval(b_treatmentsplit, .width = 0.89)

# reverse species combos
summ_stats <- summ_stats %>% separate(SpCombo, c("sp1", "sp2"), remove = F)
summ_stats2 <- summ_stats %>% mutate(SpCombo = paste(sp2, sp1))
summ_stats <- rbind(summ_stats, summ_stats2) %>% arrange(SpCombo)

traits <- read.csv("traits.csv")
plots <- read.csv("plots.csv")

NC_data <- read.csv("NC_data.csv")
NC_data$species <- substr(NC_data$ID, 1, 2)
NC_data$plot <- substr(NC_data$ID, 3, nchar(NC_data$ID))
NC_data <- NC_data %>% filter(species != "FE")


traits$treatment <- plots$treatment[match(traits$plot, plots$plot)]
traits$treatment <- factor(traits$treatment, levels = c("mix", "split"))
traits$sp1 <- plots$sp1[match(traits$plot, plots$plot)]
traits$sp2 <- plots$sp2[match(traits$plot, plots$plot)]
traits$otherspecies <- ifelse(traits$species == traits$sp1, traits$sp2, traits$sp1) # other species is the background competitor that wasn't measured in that datapoint

# traits$C_perc <- NC_data$C_perc[match(traits$plot, NC_data$plot)]
# traits$N_perc <- NC_data$N_perc[match(traits$plot, NC_data$plot)]
# traits$d15N <- NC_data$d15N[match(traits$plot, NC_data$plot)]

traits$CNratio <- NC_data$C.N[match(traits$plot, NC_data$plot)]
traits$d13C <- NC_data$d13C[match(traits$plot, NC_data$plot)]

traits2 <- traits %>% mutate(canopy = ((axis1 + axis2)/2)/(height))

full_pca <- prcomp(traits2[c(3,9, 10,12)] %>% na.omit(), scale = T)

full_pca_plot <- autoplot(full_pca, traits2[c(1,3,9, 10,12, 14)] %>% na.omit() %>% 
                            mutate(species = case_when(species == "AC" ~ "A. wrangelianus",
                                                       species == "FE" ~ "F. microstachys",
                                                       species == "PL" ~ "P. erecta",
                                                       species == "SA" ~ "S. columbariae",
                                                       species == "URLI" ~ "U. lindleyi")) %>% 
                            mutate(treatment = case_when(treatment == "mix" ~ "mixed",
                                                       treatment == "split" ~ "clustered")), 
                          scale = T,
                          colour = 'species', shape = "treatment", alpha = 0.7,
                          frame = T, frame.type = 'norm',
                          loadings = TRUE, loadings.color = "black", 
                          loadings.label = TRUE, loadings.label.color = "black", loadings.label.vjust = -.25, loadings.label.hjust = -.25) + 
  theme_classic()

full_pca_plot

jpeg("../SplitMixed_Figures/SIFigfullPCA.jpeg",
     width = 2000, height = 1500, res = 300)
print(full_pca_plot)
dev.off()


traits_pca <- data.frame()
plot_list <- list()
for (sp in c("SA", "URLI", "FE", "PL")) {
  
  sp_traits <- traits2 %>% filter(species == sp)
  sp_traits <- sp_traits[,c(2,3,9, 10,12, 14, 17:19)] %>% na.omit()
  sp_traits <- rename(sp_traits, "C:N ratio" = CNratio)
  
  pca <- prcomp(sp_traits[c(2:5, 8:9)], scale = T)
  
  loadings_sp <- pca$rotation
  axes_sp <- predict(pca, newdata = sp_traits)
  
  sp_frame <- cbind(sp_traits, axes_sp)
  sp_frame$sp <- sp
  traits_pca <- rbind(traits_pca, sp_frame)
  
  p <- autoplot(pca, sp_frame[c(2:6, 8:9)], scale = T,
                colour = 'treatment', alpha = 0.7,
                loadings = TRUE, loadings.color = "black", 
                loadings.label = TRUE, loadings.label.color = "black", loadings.label.vjust = -.25, loadings.label.hjust = -.25) + 
         theme_classic() +
         scale_color_manual(values=c("#332288", "#117733")) +
    ggtitle(case_when(sp == "FE" ~ "F. microstachys",
                      sp == "PL" ~ "P. erecta",
                      sp == "SA" ~ "S. columbariae",
                      sp == "URLI" ~ "U. lindleyi")) + 
    theme(plot.title = element_text(face = "italic"))
  
  plot_list[[sp]] <- p
}

## AC
acmispon_traits <- traits2 %>% filter(species == "AC")
acmispon_traits <- acmispon_traits[,c(2, 3,9, 10,12, 14, 17:20)] %>% na.omit()
acmispon_traits <- rename(acmispon_traits, "canopy index" = canopy, "C:N ratio" = CNratio)

pca_ac <- prcomp(acmispon_traits[c(2:5, 8:10)], scale = T)

AC_plot <- autoplot(pca_ac, acmispon_traits[c(2:6, 8:10)], scale = T,
         colour = 'treatment', alpha = 0.7,
         loadings = TRUE, loadings.color = "black", 
         loadings.label = TRUE, loadings.label.color = "black", loadings.label.vjust = -.25, loadings.label.hjust = -.25) + 
  theme_classic() +
  scale_color_manual(values=c("#332288", "#117733")) +
  ggtitle("A. wrangelianus") + theme(plot.title = element_text(face = "italic"))

plot_list[["AC"]] <- AC_plot

legend <- gtable::gtable_filter(ggplotGrob(plot_list[[1]]), "guide-box")

# Remove the legend from each plot (to avoid duplication)
for (i in 1:length(plot_list)) {
  plot_list[[i]] <- plot_list[[i]] + theme(legend.position="none")
}

jpeg("SIFigPCA.jpeg",
     width = 2500, height = 2500, res = 300)
PCA_all <- grid.arrange(grobs = c(plot_list, list(legend)),
                        ncol = 2)
dev.off()

loadings_ac <- pca_ac$rotation
axes_ac <- predict(pca_ac, newdata = acmispon_traits)

ac_frame <- cbind(acmispon_traits, axes_ac)
ac_frame$sp <- "AC"

traits_pca <- rbind(traits_pca, ac_frame[c(1:9, 11:16, 18)])


traits_pca <- traits_pca %>% mutate(sp = recode(sp, URLI = 'UR'), SpCombo = paste(sp, otherspecies))

sp_combos <- unique(traits_pca$SpCombo)

loo_data <- data.frame()
trait_draws <- data.frame()
for(cur_combo in sp_combos) {
  
  cur_trait <- traits_pca %>%
    filter(SpCombo == cur_combo) %>%
    arrange(treatment)
  
  cur_fit <- brm(formula = PC1 ~ 1 + treatment,
                 data = cur_trait, backend = "cmdstanr")
  
  cur_draws <- cur_fit %>%
    spread_draws(b_Intercept,
                 b_treatmentsplit)
  
  cur_draws$SpCombo <- unique(cur_trait$SpCombo)
  trait_draws <- rbind(trait_draws, cur_draws)
  
}

PC1_fig <- ggplot(trait_draws, aes(y = b_treatmentsplit, x = SpCombo, ,
                                   fill = after_stat(y < 0))) +
  stat_slab(aes(alpha = (after_stat(level))), .width = c(.95, 1)) +
  scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
  geom_hline(yintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size=12),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Species Combination",
       y = "Difference in PC1\n Clustered Relative to Mixed")
PC1_fig


summ_traits <- trait_draws %>% group_by(SpCombo) %>% point_interval(b_treatmentsplit, .width = 0.89) %>%
  rename(b_treatmentsplit_PC1 = b_treatmentsplit,
         lower_PC1 = .lower,
         upper_PC1 = .upper,
         width = .width,
         point = .point,
         interval = .interval)

summ_stats2 <- cbind(summ_stats, summ_traits)
summ_stats2$SpCombo2 <- paste(summ_stats2$sp1, summ_stats2$sp2)

summ_stats2 <- summ_stats2[-c(10)]

summ_stats2$abs_b_treatmentsplit_PC1 <- abs(summ_stats2$b_treatmentsplit_PC1)
summ_stats2$abs_lower_PC1 <- summ_stats2$lower_PC1 + (summ_stats2$abs_b_treatmentsplit_PC1 - summ_stats2$b_treatmentsplit_PC1)
summ_stats2$abs_upper_PC1 <- summ_stats2$upper_PC1 + (summ_stats2$abs_b_treatmentsplit_PC1 - summ_stats2$b_treatmentsplit_PC1)

summ_stats_max <- summ_stats2 %>% group_by(SpCombo2) %>% slice_max(abs_b_treatmentsplit_PC1)


max_combos <- unique(summ_stats_max$SpCombo)
trait_draws2 <- trait_draws %>% filter(SpCombo %in% max_combos) %>% rename(traits_b_treatmentsplit = b_treatmentsplit)

fec_draws$SpCombo2 <- paste(fec_draws$sp2, fec_draws$sp1)
fec_draws$SpCombo <- ifelse(fec_draws$SpCombo %in% max_combos, fec_draws$SpCombo, fec_draws$SpCombo2)
combined_draws <- merge(trait_draws2 %>% select(-c(b_Intercept)), fec_draws %>% select(c(b_treatmentsplit, SpCombo, .chain, .iteration, .draw)),
                        by = c("SpCombo", ".chain", ".iteration", ".draw"))

# summ_combined_draws <- combined_draws %>% group_by(.chain, .iteration, .draw) %>%
#   summarize(slope = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[2,1],
#             intercept = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[1,1],
#             p = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[2,4])



summ_fec_vs_traits <- combined_draws %>%
  group_by(SpCombo) %>%
  dplyr::summarise(MeanTraits = abs(mean(traits_b_treatmentsplit)),
                   SdTraits = abs(sd(traits_b_treatmentsplit)),
                   SeTraits = SdTraits/(sqrt(15)),
                   MeanFec = mean(b_treatmentsplit),
                   SdFec = sd(b_treatmentsplit),
                   SeFec = SdFec/(sqrt(15))
  )

meta_fit_traits <- brm(formula = MeanFec | se(SdFec) ~ me(MeanTraits, SdTraits),
                       data = summ_fec_vs_traits, backend = "cmdstanr",
                       family = gaussian, iter = 100000,
                       control = list(adapt_delta = 0.99),
                       prior = set_prior("normal(0,0.5)", class = "Intercept")
)
fix_ef_reg <- fixef(meta_fit_traits)

summ_fec_vs_traits <- summ_fec_vs_traits %>% separate(SpCombo, c("sp1", "sp2")) %>%
  mutate(sp1 = case_when(sp1 == "AC" ~ "A. wrangelianus",
                         sp1 == "FE" ~ "F. microstachys",
                         sp1 == "PL" ~ "P. erecta",
                         sp1 == "SA" ~ "S. columbariae",
                         sp1 == "UR" ~ "U. lindleyi")) %>%
  mutate(sp2 = case_when(sp2 == "AC" ~ "A. wrangelianus",
                         sp2 == "FE" ~ "F. microstachys",
                         sp2 == "PL" ~ "P. erecta",
                         sp2 == "SA" ~ "S. columbariae",
                         sp2 == "UR" ~ "U. lindleyi")) %>%
  mutate(SpCombo = ifelse(sp1 < sp2, paste0(sp1, "\n", sp2), paste0(sp2, "\n", sp1)))

plTraitsFecCorr <- ggplot() +
  geom_abline(intercept = fix_ef_reg[1,1],
              slope = fix_ef_reg[2,1],
              size = 1) +
  geom_point(data = summ_fec_vs_traits,
             aes(x = MeanTraits, y = MeanFec, color = SpCombo),
             size = 3) +
  geom_errorbarh(data = summ_fec_vs_traits,
                 height = 0.01,
                 aes(y = MeanFec,
                     xmin = MeanTraits - SdTraits,
                     xmax = MeanTraits + SdTraits,
                     color = SpCombo)) +
  geom_errorbar(data = summ_fec_vs_traits,
                aes(x = MeanTraits,
                    y = MeanFec,
                    ymin = MeanFec - SdFec,
                    ymax = MeanFec + SdFec,
                    color = SpCombo)) +
  theme_classic() +
  theme(text = element_text(size=15),
        axis.text.x = element_text(size = 15),
        legend.text=element_text(size = 10, face = "italic"),
        strip.background = element_blank(),
        legend.key.spacing.y = unit(0.25, "cm"),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Difference in Competitor Size and Traits (PC1)\n Clustered Relative to Mixed",
       y = "Difference in Focal Seed Production\n in Clustered Relative to Mixed Competitors",
       color = "Species\nCombination")
plot(plTraitsFecCorr)

slope_data <- meta_fit_traits %>%
  spread_draws(b_Intercept, bsp_meMeanTraitsSdTraits)

plSlopeHist <- ggplot(slope_data, aes(x = bsp_meMeanTraitsSdTraits,
                                      fill = after_stat(x < 0))) +
  stat_slab(aes(alpha = (after_stat(level))), .width = c(.95, 1)) +
  scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
  scale_alpha_manual(values = c(0.25, 1)) +
  geom_vline(xintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 10),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot",
        plot.background = element_rect(colour = "black", fill=NA, linewidth=1)) +
  xlim(-2.5, 4) +
  labs(x = "Slope",
       y = "")
plSlopeHist


vp <- viewport(width = 0.25, height = 0.25, x = 0.66, y = 0.3)

print(plTraitsFecCorr)
print(plSlopeHist, vp = vp)

jpeg("../SplitMixed_Figures/FigTraits.jpeg",
     width = 2500, height = 1500, res = 300)
print(plTraitsFecCorr)
print(plSlopeHist, vp = vp)
dev.off()

summ_slopes <- sum(slope_data$bsp_meMeanTraitsSdTraits > 0) / nrow(slope_data)
summ_slopes


## AVG INSTEAD OF MAX ##

trait_draws_avg <- trait_draws %>% separate_wider_delim(SpCombo, delim = " ", names = c("sp1", "sp2"), cols_remove = F) %>% 
  mutate(SpCombo3 = ifelse(sp1 < sp2, paste(sp1, sp2), paste(sp2, sp1)))

trait_draws_avg <- trait_draws_avg %>% group_by(.chain, .iteration, .draw, SpCombo3) %>% 
  summarize(traits_b_treatmentsplit = mean(abs(b_treatmentsplit)), traits_b_Intercept = mean(b_Intercept))

fec_draws$SpCombo3 <- ifelse(fec_draws$sp1 < fec_draws$sp2, paste(fec_draws$sp1, fec_draws$sp2), paste(fec_draws$sp2, fec_draws$sp1))

combined_draws_avg <- merge(trait_draws_avg %>% select(-c(traits_b_Intercept)), fec_draws %>% select(c(b_treatmentsplit, SpCombo3, .chain, .iteration, .draw)), 
                            by = c("SpCombo3", ".chain", ".iteration", ".draw"))

summ_combined_draws_avg <- combined_draws_avg %>% group_by(.chain, .iteration, .draw) %>% 
  summarize(slope = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[2,1], 
            intercept = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[1,1],
            p = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[2,4])

summ_fec_vs_traits_avg <- combined_draws_avg %>%
  group_by(SpCombo3) %>%
  dplyr::summarise(MeanTraits = abs(mean(traits_b_treatmentsplit)),
                   SdTraits = abs(sd(traits_b_treatmentsplit)),
                   MeanFec = mean(b_treatmentsplit),
                   SdFec = sd(b_treatmentsplit)
  )

meta_fit_traits <- brm(formula = MeanFec | se(SdFec) ~ me(MeanTraits, SdTraits),
                       data = summ_fec_vs_traits_avg, backend = "cmdstanr",
                       family = gaussian, iter = 100000,
                       control = list(adapt_delta = 0.99),
                       prior = set_prior("normal(0,0.5)", class = "Intercept")
                       )
fix_ef_reg <- fixef(meta_fit_traits)

summ_fec_vs_traits_avg <- summ_fec_vs_traits_avg %>% separate(SpCombo3, c("sp1", "sp2")) %>% 
  mutate(sp1 = case_when(sp1 == "AC" ~ "A. wrangelianus",
                         sp1 == "FE" ~ "F. microstachys",
                         sp1 == "PL" ~ "P. erecta",
                         sp1 == "SA" ~ "S. columbariae",
                         sp1 == "UR" ~ "U. lindleyi")) %>%
  mutate(sp2 = case_when(sp2 == "AC" ~ "A. wrangelianus",
                         sp2 == "FE" ~ "F. microstachys",
                         sp2 == "PL" ~ "P. erecta",
                         sp2 == "SA" ~ "S. columbariae",
                         sp2 == "UR" ~ "U. lindleyi")) %>%
  mutate(SpCombo = ifelse(sp1 < sp2, paste0(sp1, "\n", sp2), paste0(sp2, "\n", sp1)))

plTraitsFecCorr <- ggplot() +
  geom_abline(intercept = fix_ef_reg[1,1],
              slope = fix_ef_reg[2,1],
              size = 1) +
  geom_point(data = summ_fec_vs_traits_avg,
             aes(x = MeanTraits, y = MeanFec, color = SpCombo),
             size = 3) +
  geom_errorbarh(data = summ_fec_vs_traits_avg,
                 height = 0.01,
                 aes(y = MeanFec,
                     xmin = MeanTraits - SdTraits,
                     xmax = MeanTraits + SdTraits,
                     color = SpCombo)) +
  geom_errorbar(data = summ_fec_vs_traits_avg,
                aes(x = MeanTraits,
                    y = MeanFec,
                    ymin = MeanFec - SdFec,
                    ymax = MeanFec + SdFec,
                    color = SpCombo)) +
  theme_classic() +
  theme(text = element_text(size=15),
        axis.text.x = element_text(size = 15),
        legend.text=element_text(size = 10, face = "italic"),
        strip.background = element_blank(),
        legend.key.spacing.y = unit(0.25, "cm"),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Average Difference in Competitor Performance and Traits\n Clustered Relative to Mixed",
       y = "Difference in Focal Fecundity\n in Clustered Relative to Mixed Competitors",
       color = "Species\nCombination")
plTraitsFecCorr

slope_data <- meta_fit_traits %>%
  spread_draws(b_Intercept, bsp_meMeanTraitsSdTraits)

plSlopeHist <- ggplot(slope_data, aes(x = bsp_meMeanTraitsSdTraits,
                                      fill = after_stat(x < 0))) +
  stat_slab(aes(alpha = (after_stat(level))), .width = c(.95, 1)) +
  scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
  scale_alpha_manual(values = c(0.25, 1)) +
  geom_vline(xintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 10),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot",
        plot.background = element_rect(colour = "black", fill=NA, linewidth=1)) +
  xlim(-2.5, 4) +
  labs(x = "Slope",
       y = "")
plSlopeHist


vp <- viewport(width = 0.25, height = 0.25, x = 0.65, y = 0.3)

print(plTraitsFecCorr)
print(plSlopeHist, vp = vp)

summ_slopes <- sum(slope_data$bsp_meMeanTraitsSdTraits > 0) / nrow(slope_data)
summ_slopes

jpeg("SIFigTraits_avg.jpeg",
     width = 2500, height = 1500, res = 300)
print(plTraitsFecCorr)
print(plSlopeHist, vp = vp)
dev.off()


## Check using PC2 ##

loo_data <- data.frame()
trait_drawsPC2 <- data.frame()
for(cur_combo in sp_combos) {
  
  cur_trait <- traits_pca %>%
    filter(SpCombo == cur_combo) %>%
    arrange(treatment)
  
  cur_fit <- brm(formula = PC2 ~ 1 + treatment, # this is the only difference
                 data = cur_trait, backend = "cmdstanr")
  
  cur_draws <- cur_fit %>%
    spread_draws(b_Intercept,
                 b_treatmentsplit)
  
  cur_draws$SpCombo <- unique(cur_trait$SpCombo)
  trait_drawsPC2 <- rbind(trait_drawsPC2, cur_draws)
  
}

PC2_fig <- ggplot(trait_drawsPC2, aes(y = b_treatmentsplit, x = SpCombo, ,
                        fill = after_stat(y < 0))) +
  stat_slab(aes(alpha = (after_stat(level))), .width = c(.95, 1)) +
  scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
  geom_hline(yintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size=12),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Species Combination",
       y = "Difference in PC2\n Clustered Relative to Mixed")
PC2_fig

jpeg("../SplitMixed_Figures/PC_fig.jpeg",
     width = 3500, height = 2000, res = 300)
grid.arrange(PC1_fig, PC2_fig, ncol = 1)
dev.off()

summ_traits <- trait_drawsPC2 %>% group_by(SpCombo) %>% point_interval(b_treatmentsplit, .width = 0.89) %>%
  rename(b_treatmentsplit_PC2 = b_treatmentsplit,
         lower_PC2 = .lower,
         upper_PC2 = .upper,
         width = .width,
         point = .point,
         interval = .interval)

summ_stats2 <- cbind(summ_stats, summ_traits)
summ_stats2$SpCombo2 <- paste(summ_stats2$sp1, summ_stats2$sp2)

summ_stats2 <- summ_stats2[-c(10)]

summ_stats2$abs_b_treatmentsplit_PC2 <- abs(summ_stats2$b_treatmentsplit_PC2)
summ_stats2$abs_lower_PC2 <- summ_stats2$lower_PC2 + (summ_stats2$abs_b_treatmentsplit_PC2 - summ_stats2$b_treatmentsplit_PC2)
summ_stats2$abs_upper_PC2 <- summ_stats2$upper_PC2 + (summ_stats2$abs_b_treatmentsplit_PC2 - summ_stats2$b_treatmentsplit_PC2)

summ_stats_max <- summ_stats2 %>% group_by(SpCombo2) %>% slice_max(abs_b_treatmentsplit_PC2)


max_combos <- unique(summ_stats_max$SpCombo)
trait_draws2 <- trait_draws %>% filter(SpCombo %in% max_combos) %>% rename(traits_b_treatmentsplit = b_treatmentsplit)

fec_draws$SpCombo2 <- paste(fec_draws$sp2, fec_draws$sp1)
fec_draws$SpCombo <- ifelse(fec_draws$SpCombo %in% max_combos, fec_draws$SpCombo, fec_draws$SpCombo2)
combined_draws <- merge(trait_draws2 %>% select(-c(b_Intercept)), fec_draws %>% select(c(b_treatmentsplit, SpCombo, .chain, .iteration, .draw)),
                        by = c("SpCombo", ".chain", ".iteration", ".draw"))

# summ_combined_draws <- combined_draws %>% group_by(.chain, .iteration, .draw) %>%
#   summarize(slope = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[2,1],
#             intercept = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[1,1],
#             p = coef(summary(lm(b_treatmentsplit ~ abs(traits_b_treatmentsplit))))[2,4])



summ_fec_vs_traits <- combined_draws %>%
  group_by(SpCombo) %>%
  dplyr::summarise(MeanTraits = abs(mean(traits_b_treatmentsplit)),
                   SdTraits = abs(sd(traits_b_treatmentsplit)),
                   SeTraits = SdTraits/(sqrt(15)),
                   MeanFec = mean(b_treatmentsplit),
                   SdFec = sd(b_treatmentsplit),
                   SeFec = SdFec/(sqrt(15))
  )

meta_fit_traits <- brm(formula = MeanFec | se(SdFec) ~ me(MeanTraits, SdTraits),
                       data = summ_fec_vs_traits, backend = "cmdstanr",
                       family = gaussian, iter = 100000,
                       control = list(adapt_delta = 0.99),
                       prior = set_prior("normal(0,0.5)", class = "Intercept")
)
fix_ef_reg <- fixef(meta_fit_traits)

summ_fec_vs_traits <- summ_fec_vs_traits %>% separate(SpCombo, c("sp1", "sp2")) %>%
  mutate(sp1 = case_when(sp1 == "AC" ~ "A. wrangelianus",
                         sp1 == "FE" ~ "F. microstachys",
                         sp1 == "PL" ~ "P. erecta",
                         sp1 == "SA" ~ "S. columbariae",
                         sp1 == "UR" ~ "U. lindleyi")) %>%
  mutate(sp2 = case_when(sp2 == "AC" ~ "A. wrangelianus",
                         sp2 == "FE" ~ "F. microstachys",
                         sp2 == "PL" ~ "P. erecta",
                         sp2 == "SA" ~ "S. columbariae",
                         sp2 == "UR" ~ "U. lindleyi")) %>%
  mutate(SpCombo = ifelse(sp1 < sp2, paste0(sp1, "\n", sp2), paste0(sp2, "\n", sp1)))

plTraitsFecCorr <- ggplot() +
  # geom_abline(intercept = fix_ef_reg[1,1],
  #             slope = fix_ef_reg[2,1],
  #             size = 1) +
  geom_point(data = summ_fec_vs_traits,
             aes(x = MeanTraits, y = MeanFec, color = SpCombo),
             size = 3) +
  geom_errorbarh(data = summ_fec_vs_traits,
                 height = 0.01,
                 aes(y = MeanFec,
                     xmin = MeanTraits - SdTraits,
                     xmax = MeanTraits + SdTraits,
                     color = SpCombo)) +
  geom_errorbar(data = summ_fec_vs_traits,
                aes(x = MeanTraits,
                    y = MeanFec,
                    ymin = MeanFec - SdFec,
                    ymax = MeanFec + SdFec,
                    color = SpCombo)) +
  theme_classic() +
  theme(text = element_text(size=15),
        axis.text.x = element_text(size = 15),
        legend.text=element_text(size = 10, face = "italic"),
        strip.background = element_blank(),
        legend.key.spacing.y = unit(0.25, "cm"),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Difference in PC2\n Clustered Relative to Mixed",
       y = "Difference in Focal Fecundity\n in Clustered Relative to Mixed Competitors",
       color = "Species\nCombination")
plot(plTraitsFecCorr)

slope_data <- meta_fit_traits %>%
  spread_draws(b_Intercept, bsp_meMeanTraitsSdTraits)

plSlopeHist <- ggplot(slope_data, aes(x = bsp_meMeanTraitsSdTraits,
                                      fill = after_stat(x < 0))) +
  stat_slab(aes(alpha = (after_stat(level))), .width = c(.95, 1)) +
  scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
  scale_alpha_manual(values = c(0.25, 1)) +
  geom_vline(xintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 10),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot",
        plot.background = element_rect(colour = "black", fill=NA, linewidth=1)) +
  xlim(-2.5, 4) +
  labs(x = "Slope",
       y = "")
plSlopeHist


vp <- viewport(width = 0.25, height = 0.25, x = 0.66, y = 0.3)

print(plTraitsFecCorr)
print(plSlopeHist, vp = vp)

jpeg("../SplitMixed_Figures/FigTraits_PC2.jpeg",
     width = 2500, height = 1500, res = 300)
print(plTraitsFecCorr)
print(plSlopeHist, vp = vp)
dev.off()

summ_slopes <- sum(slope_data$bsp_meMeanTraitsSdTraits > 0) / nrow(slope_data)
summ_slopes



