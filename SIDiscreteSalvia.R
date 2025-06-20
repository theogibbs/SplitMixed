library(tidyverse)
library(reshape2)
library(brms)
library(cmdstanr)
library(tidybayes)
library(reshape2)
library(ggpubr)
library(gridExtra)

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
  mutate(Seeds = 1) %>% # replacing heads with seeds
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
  mutate(Seeds = 1) %>% # converting to seeds
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

# removing data for missing focals
rem_bare_na_fec_data <- fecundity %>%
  filter(treatment != "bare") %>%
  mutate(sp1_count = ifelse(treatment == "bare", 0, sp1_count)) %>%
  mutate(sp2_count = ifelse(treatment == "bare", 0, sp2_count)) %>%
  filter(!(is.na(sp1_count & is.na(sp2_count)))) %>%
  filter(!is.na(Seeds))

sp_combos <- unique(rem_bare_na_fec_data$combo)

rem_bare_na_fec_data$treatment <- fct_rev(rem_bare_na_fec_data$treatment)

loo_data <- data.frame()
fec_draws <- data.frame()
for(cur_foc in c("SA")) {
  for(cur_combo in c("FE SA", "PL SA", "PL UR", "SA UR")) {
    
    cur_fec <- rem_bare_na_fec_data %>%
      filter(focal == cur_foc) %>%
      filter(combo == cur_combo) %>%
      arrange(treatment)
    
    
    cur_fit <- brm(formula = Seeds ~ 1 + treatment + sp1_count + sp2_count,
                   data = cur_fec, backend = "cmdstanr", family = poisson)
    
    cur_loo <- loo(cur_fit)$estimates[3, 1:2]
    cur_loo <- as.data.frame(t(cur_loo))
    cur_loo$Focal <- cur_foc
    cur_loo$SpCombo <- cur_combo
    loo_data <- rbind(loo_data, cur_loo)
    
    cur_draws <- cur_fit %>%
      spread_draws(b_Intercept,
                   b_treatmentsplit,
                   b_sp1_count,
                   b_sp2_count)
    
    colnames(cur_draws) <- c(".chain",
                             ".iteration",
                             ".draw",
                             "Intercept",
                             "treatmentsplit",
                             "sp1_slope",
                             "sp2_slope")
    
    cur_draws$Focal <- cur_foc
    cur_draws$sp1 <- unique(cur_fec$sp1)
    cur_draws$sp2 <- unique(cur_fec$sp2)
    fec_draws <- rbind(fec_draws, cur_draws)
    
  }
}

fec_draws <- fec_draws %>%
  mutate(sp1 = case_when(sp1 == "AC" ~ "A. wrangelianus",
                         sp1 == "FE" ~ "F. microstachys",
                         sp1 == "PL" ~ "P. erecta",
                         sp1 == "SA" ~ "S. columbariae",
                         sp1 == "UR" ~ "U. lindleyi")) %>%
  mutate(sp2 = case_when(sp2 == "AC" ~ "A. wrangelianus",
                         sp2 == "FE" ~ "F. microstachys",
                         sp2 == "PL" ~ "P. erecta",
                         sp2 == "SA" ~ "S. columbariae",
                         sp2 == "UR" ~ "U. lindleyi"))


fec_draws$SpCombo <- paste0(fec_draws$sp1, "\n", fec_draws$sp2)

ggplot(fec_draws, aes(y = Intercept, x = SpCombo, color = rev(SpCombo))) +
  stat_pointinterval(.width = c(0.89, 0.95)) +
  facet_wrap(~Focal, scales = "free_y", nrow = 3) +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = "Species Combination", y = "Mean Number of Seeds in Mixed Treatment")

plFecDiff <- ggplot(fec_draws, aes(y = treatmentsplit, x = SpCombo, ,
                                   fill = after_stat(y < 0))) +
  stat_slab(aes(alpha = (after_stat(level))), .width = c(.95, 1)) +
  scale_alpha_manual(values = c(0.25, 1)) +
  scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
  geom_hline(yintercept = 0, alpha = 0.75, linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position = "none",
        axis.text.x = element_text(size = 10, face = "italic"),,
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot",
        ) +
  labs(x = "Species Combination",
       y = "Estimated Change in Mean\nNumber of Seeds from Clustered Competitors")
plFecDiff

jpeg("SIFigDiscreteSalvia.jpeg",
     width = 3100, height = 1600, res = 300)
plFecDiff
dev.off()
