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

mod_draws <- data.frame()
for(cur_foc in focal_IDs) {
  cur_data <- dens_data %>%
    filter(focal == cur_foc)
  
  cur_pwis_fit <- brm(formula = Seeds ~ 1 + SplDensAC + SplDensFE + SplDensPL + SplDensSA + SplDensUR + MixDensAC + MixDensFE + MixDensPL + MixDensSA + MixDensUR,
                      data = cur_data, backend = "cmdstanr", family = lognormal, iter = 10000)
  
  cur_draws <- cur_pwis_fit %>%
    spread_draws(b_Intercept,
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
  
  cur_draws$focal <- cur_foc
  mod_draws <- rbind(mod_draws, cur_draws)
}

melt_draws <- mod_draws %>%
  melt(id.vars = c("focal", ".chain", ".iteration", ".draw"))

ggplot(melt_draws %>%
         filter(variable == "b_Intercept"), aes(x = value)) +
  geom_density() +
  facet_wrap( ~ focal, scales = "free") +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank())

alpha_draws <- melt_draws %>%
  filter(variable != "b_Intercept")

alpha_draws <- alpha_draws %>%
  mutate(resident = substr(variable, 10, 12)) %>%
  mutate(treatment = substr(variable, 3, 5))

ggplot(alpha_draws, aes(x = value, color = treatment)) +
  geom_density() +
  facet_grid(focal ~ resident) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank())


intra_alphas <- alpha_draws %>%
  filter(resident == focal)

res_IDs <- unique(alpha_draws$resident)

inter_vs_intra <- data.frame()

for(cur_foc in focal_IDs) {
  cur_intra <- intra_alphas %>%
    filter(focal == cur_foc)
  
  for(cur_res in res_IDs) {
    if(cur_foc == cur_res) {
      next
    }
    cur_inter <- alpha_draws %>%
      filter(focal == cur_foc & resident == cur_res) %>%
      mutate(InterIntra = (value - cur_intra$value))
    
    inter_vs_intra <- rbind(inter_vs_intra, cur_inter)
    
  }
  
}

inter_vs_intra <- inter_vs_intra %>%
  mutate(resident = paste("Resident:", resident)) %>%
  mutate(focal = paste("Focal:", focal)) %>%
  mutate(treatment = case_when(treatment == "Spl" ~ "Clustered",
                               treatment == "Mix" ~ "Mixed"))

# plotting posterior densities of the alphas
plInterIntra <- ggplot(inter_vs_intra, aes(x = InterIntra, color = treatment)) +
  geom_density(size = 1) + theme_classic() +
  facet_wrap(focal ~ resident, scales = "free", nrow = 3) +
  labs(x = "Difference of Interspecific and Intraspecific Competition Effect", y = "") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank()) +
  scale_color_manual(values=c("#332288", "#117733"))
plInterIntra

jpeg("./figs/SIFigFocalCompEff.jpeg",
     width = 3250, height = 1600, res = 300)
plInterIntra
dev.off()
