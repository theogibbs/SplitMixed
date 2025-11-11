library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(reshape2)
library(ggpubr)
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
  geom_vline(linetype = "dashed", size = 2,, color = "gray",
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

for(cur_foc in c("AC", "SA", "UR")) {
  cur_fec <- fecundity %>%
    filter(focal == cur_foc) %>%
    filter(treatment != "bare")
  
  show(ggboxplot(cur_fec, x = "combo", y = "Seeds", add = "jitter",
                 color = "treatment", group="treatment") + 
         stat_compare_means(aes(group=treatment), method="t.test", vjust = 1)  + 
         xlab("competitor combination") +
         ylab("Seeds") + ggtitle(paste(cur_foc, "focals")) + 
         facet_wrap(~combo, scales = "free") + 
         scale_color_hue(direction = -1))
  
  ## pooled across background competitor combinations
  show(ggboxplot(cur_fec %>% filter(treatment != "bare"), x = "treatment", y = "Seeds", add = "jitter",
                 color = "treatment") + 
         stat_compare_means(aes(group=treatment), method="t.test") + 
         ylab("Seeds") + ggtitle(cur_foc) +
         scale_color_hue(direction = -1))
  
}


# removing data for missing focals
rem_na_fec_data <- fec_data %>%
  mutate(sp1_count = ifelse(treatment == "bare", 0, sp1_count)) %>%
  mutate(sp2_count = ifelse(treatment == "bare", 0, sp2_count)) %>%
  filter(!(is.na(sp1_count & is.na(sp2_count)))) %>%
  filter(!is.na(Seeds))


# relabeling data to record which species are present and
# at what densities in five new columns so
# for example "DensFE" is zero
# when neither sp1 or sp2 is FE but is the value of sp1_count
# when FE is sp1 (and similarly when it is sp2)
dens_data <- rem_na_fec_data %>%
  mutate(DensAC = case_when(sp1 == "AC" ~ sp1_count,
                            sp2 == "AC" ~ sp2_count,
                            .default = 0)) %>%
  mutate(DensFE = case_when(sp1 == "FE" ~ sp1_count,
                            sp2 == "FE" ~ sp2_count,
                            .default = 0)) %>%
  mutate(DensPL = case_when(sp1 == "PL" ~ sp1_count,
                            sp2 == "PL" ~ sp2_count,
                            .default = 0)) %>%
  mutate(DensSA = case_when(sp1 == "SA" ~ sp1_count,
                            sp2 == "SA" ~ sp2_count,
                            .default = 0)) %>%
  mutate(DensUR = case_when(sp1 == "UR" ~ sp1_count,
                            sp2 == "UR" ~ sp2_count,
                            .default = 0)) %>%
  mutate(MixLog = treatment == "mix") %>%
  mutate(SplLog = treatment == "split") %>%
  mutate(SplDensAC = as.numeric(DensAC > 0)) %>% # only taking the presence / absence of species
  mutate(SplDensFE = as.numeric(DensFE > 0)) %>% # rather than their actual density / abundance
  mutate(SplDensPL = as.numeric(DensPL > 0)) %>%
  mutate(SplDensSA = as.numeric(DensSA > 0)) %>%
  mutate(SplDensUR = as.numeric(DensUR > 0)) %>%
  mutate(MixDensAC = SplDensAC * MixLog) %>% 
  mutate(MixDensFE = SplDensFE * MixLog) %>% 
  mutate(MixDensPL = SplDensPL * MixLog) %>%
  mutate(MixDensSA = SplDensSA * MixLog) %>%
  mutate(MixDensUR = SplDensUR * MixLog) %>%
  mutate(SplDensAC = SplDensAC * SplLog) %>% 
  mutate(SplDensFE = SplDensFE * SplLog) %>% 
  mutate(SplDensPL = SplDensPL * SplLog) %>%
  mutate(SplDensSA = SplDensSA * SplLog) %>%
  mutate(SplDensUR = SplDensUR * SplLog)

# summing competitor densities for plotting
plot_data <- dens_data %>%
  mutate(SumDens = DensAC + DensFE + DensPL + DensSA + DensUR)

# seeds versus density plot
ggplot(plot_data, aes(x = SumDens, y = Seeds)) +
  geom_point() + scale_y_log10() +
  facet_wrap(~focal, scales = "free") + theme_classic()


# the following all pwi + hois model fits poorly because of the degenaracy of the terms
#all_pwi_hois_fit <- brm(formula = Seeds ~ 1 + focal + SplDensAC + SplDensFE + SplDensPL + SplDensSA + SplDensUR + MixDensAC + MixDensFE + MixDensPL + MixDensSA + MixDensUR + MixDensACFE + MixDensACPL + MixDensACSA + MixDensACUR + MixDensFEPL + MixDensFESA + MixDensFEUR + MixDensPLSA + MixDensPLUR + MixDensSAUR,
#                    data = dens_data, backend = "cmdstanr", family = lognormal)

loo_data <- data.frame()

all_pwis_fit <- brm(formula = Seeds ~ 1 + focal + SplDensAC + SplDensFE + SplDensPL + SplDensSA + SplDensUR + MixDensAC + MixDensFE + MixDensPL + MixDensSA + MixDensUR,
                    data = dens_data, backend = "cmdstanr", family = lognormal, iter = 10000)

mod_draws <- all_pwis_fit %>%
  spread_draws(b_Intercept,
               b_focalSA,
               b_focalUR,
               b_SplDensAC,
               b_SplDensFE,
               b_SplDensPL,
               b_SplDensSA,
               b_SplDensUR,
               b_MixDensAC,
               b_MixDensFE,
               b_MixDensPL,
               b_MixDensSA,
               b_MixDensUR)

melt_draws <- mod_draws %>%
  melt(id.vars = c(".chain", ".iteration", ".draw"))

ggplot(melt_draws %>%
         filter(variable == "b_Intercept"), aes(x = value)) +
  geom_density() +
  facet_grid( ~ variable) +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank())


ggplot(melt_draws %>%
         filter(substr(variable, 3, 5) == "foc"), aes(x = value)) +
  geom_density() +
  facet_wrap( ~ variable, nrow = 5) +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank())


alpha_draws <- melt_draws %>%
  filter(variable != "b_Intercept") %>%
  filter(substr(variable, 3, 5) != "foc")


ggplot(alpha_draws %>%
         filter(substr(variable, 3, 5) == "Spl"), aes(x = value)) +
  geom_density() +
  facet_grid( ~ variable) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank())

ggplot(alpha_draws %>%
         filter(substr(variable, 3, 5) == "Mix"), aes(x = value)) +
  geom_density() +
  facet_grid( ~ variable) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank())

diff_draws <- mod_draws %>%
  mutate(DiffAC = b_SplDensAC - b_MixDensAC) %>%
  mutate(DiffFE = b_SplDensFE - b_MixDensFE) %>%
  mutate(DiffPL = b_SplDensPL - b_MixDensPL) %>%
  mutate(DiffSA = b_SplDensSA - b_MixDensSA) %>%
  mutate(DiffUR = b_SplDensUR - b_MixDensUR) %>%
  select(c("DiffAC", "DiffFE", "DiffPL", "DiffSA", "DiffUR"))

ggplot(melt(diff_draws), aes(x = value)) +
  geom_density() +
  facet_wrap(~variable) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank())

comp_ab_draws <- mod_draws %>%
  mutate(AC = b_SplDensAC + b_MixDensAC) %>%
  mutate(FE = b_SplDensFE + b_MixDensFE) %>%
  mutate(PL = b_SplDensPL + b_MixDensPL) %>%
  mutate(SA = b_SplDensSA + b_MixDensSA) %>%
  mutate(UR = b_SplDensUR + b_MixDensUR) %>%
  select(c("AC", "FE", "PL", "SA", "UR"))

plCompAb <- ggplot(melt(comp_ab_draws), aes(x = value, color = variable)) +
  geom_density() +
  labs(x = "Competitive Effect",
       y = "Density", color = "Species") +
  theme_classic() +
  theme(text = element_text(size=15),
        axis.text.x = element_text(size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.caption.position =  "plot")
plCompAb

jpeg("./figs/SIFigCompAb.jpeg",
     width = 3000, height = 1000, res = 300)
plCompAb
dev.off()

sp_combos <- paste(dens_data$sp1, dens_data$sp2)
sp_combos <- unique(sp_combos)
sp_combos <- sort(sp_combos)
sp_combos <- sp_combos[sp_combos != "NA NA"]
sp_combos <- data.frame(sp1 = substr(sp_combos, 1, 2), sp2 = substr(sp_combos, 4, 5))

comp_ab <- data.frame()
for(i in 1:nrow(sp_combos)) {
  cur_sp1 <- sp_combos$sp1[i]
  cur_sp2 <- sp_combos$sp2[i]
  
  cur_sp1_draws <- comp_ab_draws[,colnames(comp_ab_draws) == cur_sp1]
  cur_sp2_draws <- comp_ab_draws[,colnames(comp_ab_draws) == cur_sp2]
  
  cur_comp_ab <- data.frame(sp1 = cur_sp1, sp2 = cur_sp2,
                            sp1_comp_ab = cur_sp1_draws,
                            sp2_comp_ab = cur_sp2_draws,
                            Diff = cur_sp1_draws - cur_sp2_draws)
  colnames(cur_comp_ab) <- c("sp1", "sp2", "sp1_comp_ab", "sp2_comp_ab", "Diff")
  
  comp_ab <- rbind(comp_ab, cur_comp_ab)
}

comp_ab_diffs <- comp_ab %>%
  mutate(SpCombo = paste(sp1, sp2))

ggplot(comp_ab_diffs, aes(x = Diff, color = SpCombo)) +
  geom_density() +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank())


# removing data for missing focals
rem_bare_na_fec_data <- fecundity %>%
  filter(treatment != "bare") %>%
  mutate(treatment = factor(treatment, levels = c("mix", "split"))) %>%
  mutate(sp1_count = ifelse(treatment == "bare", 0, sp1_count)) %>%
  mutate(sp2_count = ifelse(treatment == "bare", 0, sp2_count)) %>%
  filter(!(is.na(sp1_count & is.na(sp2_count)))) %>%
  filter(!is.na(Seeds))

sp_combos <- unique(rem_bare_na_fec_data$combo)

loo_data <- data.frame()
fec_draws <- data.frame()
for(cur_combo in sp_combos) {
  
  cur_fec <- rem_bare_na_fec_data %>%
    filter(combo == cur_combo) %>%
    arrange(treatment)
  
  cur_fit <- brm(formula = Seeds ~ 1 + treatment + focal + sp1_count + sp2_count,
                 data = cur_fec, backend = "cmdstanr", family = lognormal, iter = 10000)
  
  cur_loo <- loo(cur_fit)$estimates[3, 1:2]
  cur_loo <- as.data.frame(t(cur_loo))
  cur_loo$SpCombo <- cur_combo
  loo_data <- rbind(loo_data, cur_loo)
  
  cur_draws <- cur_fit %>%
    spread_draws(b_Intercept,
                 b_treatmentsplit,
                 b_focalSA,
                 b_focalUR,
                 b_sp1_count,
                 b_sp2_count)
  
  cur_draws$sp1 <- unique(cur_fec$sp1)
  cur_draws$sp2 <- unique(cur_fec$sp2)
  fec_draws <- rbind(fec_draws, cur_draws)
  
}

fec_draws$SpCombo <- paste(fec_draws$sp1, fec_draws$sp2)


ggplot(fec_draws, aes(y = b_Intercept, x = SpCombo, color = rev(SpCombo))) +
  stat_pointinterval(.width = c(0.89, 0.95)) +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Species Combination", y = "Mean Number of Seeds in Split Treatment")

ggplot(fec_draws, aes(y = b_treatmentsplit, x = SpCombo, color = rev(SpCombo))) +
  stat_pointinterval(.width = c(0.89, 0.95)) +
  #facet_grid( ~ Focal, scales = "free_x") +
  geom_hline(yintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Species Combination", y = "Change in Mean Number of Seeds from Mixed Treatment")


plFecDiff <- ggplot(fec_draws, aes(y = b_treatmentsplit, x = SpCombo, ,
                                   fill = after_stat(y < 0))) +
  stat_slab(aes(alpha = (after_stat(level))), .width = c(.95, 1)) +
  scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
  geom_hline(yintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Species Combination",
       y = "Estimated Change in Mean\nNumber of Seeds from Mixed Treatment")
plFecDiff


merge_fec_draws <- fec_draws %>%
  select("SpCombo", "b_treatmentsplit", ".chain", ".iteration", ".draw") %>%
  arrange(SpCombo)

fec_vs_comp_ab <- comp_ab_diffs %>%
  select(c("Diff", "SpCombo", "sp1", "sp2")) %>%
  mutate(AbsDiff = abs(Diff))

print(prod(merge_fec_draws$SpCombo == fec_vs_comp_ab$SpCombo))

fec_vs_comp_ab$FecDiff <- merge_fec_draws$b_treatmentsplit
fec_vs_comp_ab$.chain <- merge_fec_draws$.chain
fec_vs_comp_ab$.iteration <- merge_fec_draws$.iteration
fec_vs_comp_ab$.draw <- merge_fec_draws$.draw


fec_vs_comp_ab <- fec_vs_comp_ab %>%
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
  mutate(SpCombo = paste0(sp1, "\n", sp2))

summ_fec_vs_comp_ab <- fec_vs_comp_ab %>%
  group_by(SpCombo) %>%
  dplyr::summarise(MeanCompAb = mean(AbsDiff),
                   SdCompAb = sd(AbsDiff),
                   MeanFec = mean(FecDiff),
                   SdFec = sd(FecDiff))


#meta_fit <- brm(formula = MeanFec | mi(SdFec) ~ me(MeanCompAb, SdCompAb),
#                data = summ_fec_vs_comp_ab, backend = "cmdstanr",
#                family = gaussian, iter = 20000, control = list(adapt_delta = 0.9))


meta_fit <- brm(formula = MeanFec | se(SdFec) ~ me(MeanCompAb, SdCompAb),
                data = summ_fec_vs_comp_ab, backend = "cmdstanr",
                family = gaussian, iter = 100000,
                control = list(adapt_delta = 0.99))


#meta_fit <- brm(formula = MeanFec | mi(SdFec) ~ MeanCompAb,
#                data = summ_fec_vs_comp_ab, backend = "cmdstanr",
#                family = gaussian, iter = 10000)


#meta_fit <- brm(formula = MeanFec ~ me(MeanCompAb, SdCompAb),
#                data = summ_fec_vs_comp_ab, backend = "cmdstanr",
#                family = gaussian, iter = 10000)


#meta_fit <- brm(formula = MeanFec ~ MeanCompAb,
#                data = summ_fec_vs_comp_ab, backend = "cmdstanr",
#                family = gaussian, iter = 10000)

fix_ef_reg <- fixef(meta_fit)

plCompFecCorr <- ggplot() +
  geom_abline(intercept = fix_ef_reg[1,1],
              slope = fix_ef_reg[2,1],
              size = 1) +
  geom_point(data = summ_fec_vs_comp_ab,
             aes(x = MeanCompAb, y = MeanFec, color = SpCombo),
             size = 3) +
  geom_errorbarh(data = summ_fec_vs_comp_ab,
                 height = 0.01,
                 aes(y = MeanFec,
                     xmin = MeanCompAb - SdCompAb,
                     xmax = MeanCompAb + SdCompAb,
                     color = SpCombo)) +
  geom_errorbar(data = summ_fec_vs_comp_ab,
                aes(x = MeanCompAb,
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
  labs(x = "Absolute Value of Difference in Competitive Effect",
       y = "Difference in Focal Seed Production\n in Clustered Relative to Mixed Competitors",
       color = "Species\nCombination")
plCompFecCorr

slope_data <- meta_fit %>%
  spread_draws(b_Intercept, bsp_meMeanCompAbSdCompAb)

plSlopeHist <- ggplot(slope_data, aes(x = bsp_meMeanCompAbSdCompAb,
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
  labs(x = "Slope",
       y = "")
plSlopeHist

vp <- viewport(width = 0.25, height = 0.25, x = 0.65, y = 0.3)

print(plCompFecCorr)
print(plSlopeHist, vp = vp)

jpeg("./figs/FigCompAb.jpeg",
     width = 2500, height = 1500, res = 300)
print(plCompFecCorr)
print(plSlopeHist, vp = vp)
dev.off()

summ_slopes <- sum(slope_data$bsp_meMeanCompAbSdCompAb > 0) / nrow(slope_data)
summ_slopes


