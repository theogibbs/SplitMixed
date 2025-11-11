library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(reshape2)
library(ggpubr)

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
  mutate(MixLog = treatment %in% c("bare", "mix")) %>%
  mutate(SplitLog = treatment %in% c("bare", "split"))

# summing competitor densities for plotting
plot_data <- dens_data %>%
  mutate(SumDens = DensAC + DensFE + DensPL + DensSA + DensUR)

# seeds versus density plot
ggplot(plot_data, aes(x = SumDens, y = Seeds)) +
  geom_point() + scale_y_log10() +
  facet_wrap(~focal, scales = "free") + theme_classic()

# Law-Watkinson model formula for fitting
lw_formula <- bf(log(Seeds + 1) ~ log(lambda / (1 + MixLog * (DensAC^mixalphaAC + DensFE^mixalphaFE + DensPL^mixalphaPL + DensSA^mixalphaSA + DensUR^mixalphaUR) + SplitLog * (DensAC^splalphaAC + DensFE^splalphaFE + DensPL^splalphaPL + DensSA^splalphaSA + DensUR^splalphaUR))),
                 lambda ~ 1,
                 mixalphaAC ~ 1, mixalphaFE ~ 1, mixalphaPL ~ 1, mixalphaSA ~ 1, mixalphaUR ~ 1,
                 splalphaAC ~ 1, splalphaFE ~ 1, splalphaPL ~ 1, splalphaSA ~ 1, splalphaUR ~ 1,
                 nl = TRUE)

# looping over focal species to fit parameters
loo_data <- data.frame()
focal_IDs <- c("AC", "SA", "UR")
draw_data <- data.frame()
for(cur_foc in focal_IDs) {
  
  # printing current progress
  print(paste0("FOCAL SPECIES:  ",
               cur_foc,
               "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  
  # extracting data for current focal species
  cur_data <- dens_data %>%
    filter(focal == cur_foc)
  
  lambda_data <- cur_data %>%
    filter(treatment == "bare")
  lambda_data <- lambda_data$Seeds
  
  mean_lambda <- mean(lambda_data)
  sd_lambda <- sd(lambda_data) / sqrt(length(lambda_data))
  min_lambda <- min(lambda_data)
  max_lambda <- max(lambda_data)
  
  # setting priors for all the parameters
  cur_priors <- c(set_prior(paste0("normal(", mean_lambda, ", ", sd_lambda, ")"), nlpar = "lambda", lb = min_lambda, ub = max_lambda),
                  set_prior("normal(0.05, 5)", nlpar = "mixalphaAC", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "mixalphaFE", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "mixalphaPL", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "mixalphaSA", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "mixalphaUR", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "splalphaAC", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "splalphaFE", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "splalphaPL", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "splalphaSA", lb = .001, ub = 20),
                  set_prior("normal(0.05, 5)", nlpar = "splalphaUR", lb = .001, ub = 20))
  
  # fitting the model
  cur_fit <- brm(lw_formula,
                 data = cur_data,
                 prior = cur_priors,
                 iter = 10000,
                 backend = "cmdstanr")
  
  # computing cross validation
  cur_loo <- loo(cur_fit)$estimates[3, 1:2]
  cur_loo <- as.data.frame(t(cur_loo))
  cur_loo$Focal <- cur_foc
  loo_data <- rbind(loo_data, cur_loo)
  
  # extracting draws from the posterior distribution
  cur_draw_data <- cur_fit %>%
    spread_draws(b_lambda_Intercept,
                 b_mixalphaAC_Intercept,
                 b_mixalphaFE_Intercept,
                 b_mixalphaPL_Intercept,
                 b_mixalphaSA_Intercept,
                 b_mixalphaUR_Intercept,
                 b_splalphaAC_Intercept,
                 b_splalphaFE_Intercept,
                 b_splalphaPL_Intercept,
                 b_splalphaSA_Intercept,
                 b_splalphaUR_Intercept)
  
  # naming the draw data
  colnames(cur_draw_data) <- c(".chain",
                               ".iteration",
                               ".draw",
                               "lambda",
                               "mixalphaAC",
                               "mixalphaFE",
                               "mixalphaPL",
                               "mixalphaSA",
                               "mixalphaUR",
                               "splalphaAC",
                               "splalphaFE",
                               "splalphaPL",
                               "splalphaSA",
                               "splalphaUR")
  
  # recoding the focal species and binding data together
  cur_draw_data$Focal <- cur_foc
  draw_data <- rbind(draw_data, cur_draw_data)
  
}

# reshaping data for plotting
plot_draws <- draw_data %>%
  melt(id.vars = c(".chain", ".iteration", ".draw", "Focal"))

# extracting growth rate draws and labeling focals
lambda_draws <- plot_draws %>%
  filter(variable == "lambda") %>%
  mutate(Focal = paste("Focal:", Focal))

inferred_lambdas <- plot_draws %>%
  filter(variable == "lambda") %>%
  group_by(Focal) %>%
  dplyr::summarise(value = mean(value))

# plotting posterior densities
plLambdas <- ggplot(lambda_draws, aes(x = value)) +
  geom_density(size = 1, color = "darkblue") + theme_classic() +
  facet_wrap(~Focal, scales = "free") +
  labs(x = "Growth Rate", y = "") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 10),
        legend.text=element_text(size = 15),
        strip.background = element_blank())
plLambdas

# extracting the interaction parameter draws and labeling them
alpha_draws <- plot_draws %>%
  filter(variable != "lambda") %>%
  mutate(Resident = substr(variable, 9, 10)) %>%
  mutate(Treatment = substr(variable, 1, 3)) %>%
  mutate(Resident = paste("Resident:", Resident))

inferred_alphas <- alpha_draws %>%
  filter(variable != "lambda") %>%
  group_by(Focal, Resident, Treatment) %>%
  dplyr::summarise(value = mean(value))

alpha_draws <- alpha_draws %>%
  mutate(Focal = paste("Focal:", Focal))

# plotting posterior densities of the alphas
plAlphas <- ggplot(alpha_draws, aes(x = value, color = Treatment)) +
  geom_density(size = 1) + theme_classic() +
  facet_wrap(Focal ~ Resident, scales = "free", nrow = 3) +
  labs(x = "Competition Coefficient", y = "") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank())
plAlphas

mix_alpha_draws <- alpha_draws %>%
  filter(Treatment == "mix")

split_alpha_draws <- alpha_draws %>%
  filter(Treatment == "spl")

prod(substr(mix_alpha_draws$variable, 4, 10) == substr(split_alpha_draws$variable, 4, 10))

diff_alpha_draws <- merge(mix_alpha_draws[,-c(5, 8)], split_alpha_draws[,-c(5, 8)], by = c(".chain", ".iteration", ".draw", "Focal", "Resident"))
diff_alpha_draws$LogRatio <- log(diff_alpha_draws$value.x / diff_alpha_draws$value.y)

plDiffs <- ggplot(diff_alpha_draws, aes(x = LogRatio)) +
  geom_density(size = 1) + theme_classic() +
  facet_wrap(Focal ~ Resident, scales = "free", nrow = 3) +
  labs(x = expression("Log Ratio of Mixed to Clustered Competition Coefficient:" ~ log(frac(alpha[ij]^mixed, alpha[ij]^clustered))), y = "Posterior Density") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank())
plDiffs

jpeg("./figs/SIFigLWDiffs.jpeg",
     width = 3200, height = 1600, res = 300)
plDiffs
dev.off()

agg_diff_draws <- diff_alpha_draws %>%
  group_by(.chain, .iteration, .draw, Focal) %>%
  summarise(CumLogRatio = sum(LogRatio))

plAggDiffs <- ggplot(agg_diff_draws, aes(x = CumLogRatio)) +
  geom_density(size = 1) + theme_classic() +
  facet_wrap(~ Focal, scales = "free", nrow = 1) +
  labs(x = expression("Cumulative Log Ratio of Mixed to\nClustered Competition Coefficient:" ~ Sigma[j] ~ log(frac(alpha[ij]^mixed, alpha[ij]^clustered))),
       y = "Posterior Density") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        strip.background = element_blank())
plAggDiffs

jpeg("./figs/SIFigLWAggDiffs.jpeg",
     width = 2200, height = 900, res = 300)
plAggDiffs
dev.off()

agg_diff_draws %>%
  group_by(Focal) %>%
  dplyr::summarise(positive = sum(CumLogRatio > 0) / n())

loo_data$Model <- "Law-Watkinson"
write.csv(loo_data, "LawWatkinsonLOO.csv", row.names = FALSE)

