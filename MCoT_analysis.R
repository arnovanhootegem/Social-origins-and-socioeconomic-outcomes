library(OpenMx)
library(here)
library(tidyverse)
library(data.table)
library(haven)
library(patchwork)
mxOption(NULL, "Default optimizer", "SLSQP")

### Read in prepared pedigree data
pedigrees_unstd <- data.table::fread(here("pedigrees-unstd.csv"))
pedigrees <- data.table::fread(here("pedigrees.csv"))


# Descriptive statistics:
#
# Twin pair cohort
# Twin pair zygosity
# Twin SES
# Offspring birth order
# Offspring cohort
# Offspring SES

mz_dz_split <- pedigrees_unstd %>% mutate(Zygosity = if_else(Zygo==1,"Monozygotic","Dizygotic")) %>% split(.,.$Zygosity)

descs <- function(df) {
    twinpair_stats <- df %>% 
      select(Zygo, cohort) %>% 
      summarize(N_twinpairs = n(),
              twin_cohort_mean = mean(cohort),
              twin_cohort_sd = sd(cohort))
    twin_stats <- df %>% 
      select(contains("_father_"),contains("_mother_")) %>% 
      select(-ends_with("sex"),-ends_with("_cohort")) %>% 
      mutate(rn = row_number()) %>% 
      gather(var,value, -rn) %>% 
      mutate(outcome = str_sub(var,9), ped = str_sub(var,1,7)) %>% 
      select(-var) %>% 
      group_by(rn, ped) %>% 
      spread(outcome,value) %>% 
      ungroup %>% 
      select(-rn, -ped) %>% 
      summarize(father_education_mean = mean(father_education,na.rm=TRUE),
                  mother_education_mean = mean(mother_education,na.rm=TRUE),
                  father_education_sd = sd(father_education,na.rm=TRUE),
                  mother_education_sd = sd(mother_education,na.rm=TRUE),
                  father_trei_mean = mean(father_trei,na.rm=TRUE),
                  mother_trei_mean = mean(mother_trei,na.rm=TRUE),
                  father_trei_sd = sd(father_trei,na.rm=TRUE),
                  mother_trei_sd = sd(mother_trei,na.rm=TRUE),
                  father_income_rank_mean = mean(father_income_rank,na.rm=TRUE),
                  mother_income_rank_mean = mean(mother_income_rank,na.rm=TRUE),
                  father_income_rank_sd = sd(father_income_rank,na.rm=TRUE),
                  mother_income_rank_sd = sd(mother_income_rank,na.rm=TRUE),
                  father_wealth_rank_mean = mean(father_wealth_rank,na.rm=TRUE),
                  mother_wealth_rank_mean = mean(mother_wealth_rank,na.rm=TRUE),
                  father_wealth_rank_sd = sd(father_wealth_rank,na.rm=TRUE),
                  mother_wealth_rank_sd = sd(mother_wealth_rank,na.rm=TRUE),
                  father_gross_wealth_rank_mean = mean(father_gross_wealth_rank,na.rm=TRUE),
                  mother_gross_wealth_rank_mean = mean(mother_gross_wealth_rank,na.rm=TRUE),
                  father_gross_wealth_rank_sd = sd(father_gross_wealth_rank,na.rm=TRUE),
                  mother_gross_wealth_rank_sd = sd(mother_gross_wealth_rank,na.rm=TRUE))

    offspring_stats <- df %>% 
      select(contains("_child1_"),contains("_child2_")) %>% 
      mutate(rn = row_number()) %>% 
      gather(var,value, -rn) %>% 
      mutate(outcome = str_sub(var,16), family = str_sub(var,7,7), birth_order = str_sub(var,14,14)) %>% 
      select(-var) %>% 
      group_by(rn, family, birth_order) %>% 
      spread(outcome,value) %>% 
      ungroup %>% 
      mutate(female = sex-1) %>% 
      select(-rn, -family) %>% 
      summarize(female_mean = mean(female,na.rm=TRUE),
                female_sd = sd(female,na.rm=TRUE),
                cohort_mean = mean(cohort,na.rm=TRUE),
                cohort_sd = sd(cohort,na.rm=TRUE),
                education_mean = mean(education,na.rm=TRUE),
                education_sd = sd(education,na.rm=TRUE),
                trei_mean = mean(trei,na.rm=TRUE),
                trei_sd = sd(trei,na.rm=TRUE),
                income_rank_mean = mean(income_rank,na.rm=TRUE),
                income_rank_sd = sd(income_rank,na.rm=TRUE),
                wealth_rank_mean = mean(wealth_rank,na.rm=TRUE),
                wealth_rank_sd = sd(wealth_rank,na.rm=TRUE),
                gross_wealth_rank_mean = mean(gross_wealth_rank,na.rm=TRUE),
                gross_wealth_rank_sd = sd(gross_wealth_rank,na.rm=TRUE),
                N_offspring = n())

    
    
    
    stats <- as_tibble(c(twinpair_stats,twin_stats,offspring_stats)) %>% gather(variable, statistic)
    return(stats)
}

descriptives <- mz_dz_split %>% map(descs) %>% bind_rows(.id="Zygosity") %>% 
  mutate(measure = if_else(str_detect(string = variable,pattern = "mean"),"Mean","SD")) %>% 
  mutate(variable = str_replace(string = variable, pattern = "_mean", replacement = "")) %>% 
  mutate(variable = str_replace(string = variable, pattern = "_sd", replacement = "")) %>% 
  group_by(Zygosity,variable) %>% 
  spread(measure,statistic)

knitr::kable(descriptives)


###            
### CoT analysis
###
makeCoTdf <- function(outcome, predictor, pedigrees) {
  df <- select(pedigrees, 
               m1 = str_c("family1_mother_",predictor),
               m2 = str_c("family2_mother_",predictor),
               f1 = str_c("family1_father_",predictor),
               f2 = str_c("family2_father_",predictor),
               o11 = str_c("family1_child1_",outcome),
               o12 = str_c("family2_child1_",outcome),
               o21 = str_c("family1_child2_",outcome),
               o22 = str_c("family2_child2_",outcome),
               cohort_m1 = family1_mother_cohort,
               cohort_m2 = family2_mother_cohort,
               cohort_f1 = family1_father_cohort,
               cohort_f2 = family2_father_cohort,
               cohort_o11 = family1_child1_cohort,
               cohort_o12 = family2_child1_cohort,
               cohort_o21 = family1_child2_cohort,
               cohort_o22 = family2_child2_cohort,
               sex_o11 = family1_child1_sex,
               sex_o12 = family2_child1_sex,
               sex_o21 = family1_child2_sex,
               sex_o22 = family2_child2_sex,
               relation, Zygo, cohort)
    df$mm = ifelse(df$relation == "m1m2", 1, 0)
    df$ff = ifelse(df$relation == "f1f2", 1, 0)
    df$mf = ifelse(df$relation == "m1f2", 1, 0)
    df$fm = ifelse(df$relation == "m2f1", 1, 0)
    df$mz = ifelse(df$Zygo == "1", 1, 0)
    df$dz = 1 - df$mz
    
  return(df)
}

source("adapted_mod_ACEpACEo_rAC.R")

outcomes <- c("education","trei","income_rank","wealth_rank","gross_wealth_rank")
predictors <- c("education","trei","income_rank","wealth_rank","gross_wealth_rank")

runMCoT <- function(amodel) {
  parts <- str_split(amodel, pattern = ":",n=2)[[1]]
  outcome <- parts[1]
  predictor <- parts[2]
  df <- makeCoTdf(outcome = outcome, predictor = predictor, pedigrees=pedigrees)
  message("Running MCoT for y=", outcome," x=", predictor)
  x <- estimateMCoT(df) %>% omxRunCI
}

mcot_models <- expand.grid(outcomes, predictors) %>% unite(col="model",sep=":", remove=TRUE) %>% 
  pull %>% setNames(nm=.)

mcot_results <- map(mcot_models,runMCoT)
saveRDS(mcot_results, here("mcot_results.RDS"))

# Results reporting
mcot_results <- readRDS(here("mcot_results.RDS"))

reportMCoT <- function (fit) {
  est = coef(fit)
  thecoefs <- est %>% as.data.frame %>% t
  #  se_est = sqrt(diag(vcov(fit)))
  # Parents
  #  par_p = c("a1p", "cp", "ep", "d")
  #  estp = data.frame(est = est[par_p],
  #                    se = se_est[par_p],
  #                    zest = mxEval(rbind(za1p, zcp, ep, u), fit))
  
  # Parents
  VAp = mxEval(a1p^2, fit)
  VCp = mxEval(cp^2, fit)
  VEp = mxEval(ep^2, fit)
  VACp = mxEval(2 * a1p * cp * w, fit)
  VAp + VCp + VEp + VACp
  v2p = mxEval(v2p, fit)
  # z
  zVAp = mxEval(za1p^2, fit)
  zVCp = mxEval(zcp^2, fit)
  zVEp = mxEval(zep^2, fit)
  zVACp = mxEval(2 * za1p * zcp * w, fit)
  zVAp + zVCp + zVEp + zVACp
  
  # Children
  VAo = mxEval(a1o^2, fit)
  VApo = mxEval(a2o^2, fit)
  VFo = mxEval(m^2 *v2p + p^2 * v2p + 2 * m * p * d, fit)
  VAFo = mxEval(2 * 1 * a1o * CovACo, fit)
  VCpo = mxEval(co^2, fit)
  VEo = mxEval(eo^2, fit)
  v2o = VAo + VApo + VFo + VAFo + VCpo + VEo
  mxEval(v2o, fit)
  
  zVAo = VAo / v2o
  zVApo = VApo / v2o
  zVFo = VFo / v2o
  zVAFo = VAFo / v2o
  zVCpo = VCpo / v2o
  zVEo = VEo / v2o
  zVAo + zVApo + zVFo + zVAFo + zVCpo + zVEo
  
  items <- c("VAp","VCp","VEp","VACp","v2p","zVAp","zVCp","zVEp","zVACp",
             "VAo","VApo","VFo","VAFo","VCpo","VEo","v2o", "zVAo","zVApo","zVFo","zVCpo", "zVEo")
  names(items) <- items
  stats <-map(items, function(x) get(x)[1,1]) %>% as.data.frame

  thesum <- summary(fit)
  other <- data.frame("M2LL" = thesum$Minus2LogLikelihood, 
             "AIC" = thesum$AIC.Mx,
             "N_pedigrees" = thesum$numObs,
             "N_observations" = thesum$observedStatistics)

  # intervals for path coefs
  intervals <- summary(fit) %>% pluck("CIdetail") %>% select(parameter, side, value) %>% 
    unite(col=key, parameter, side, sep="_", remove=T) %>% spread(key,value)
  
  all_results <- bind_cols(thecoefs, intervals, stats, other)
  return(all_results)
}

mcot_estimates <- map(mcot_results, reportMCoT) %>% 
  bind_rows(.id="model") %>% 
  separate(col=model, into=c("outcome","predictor"), sep=":")
fwrite(mcot_estimates, here("mcot_estimates.csv"))


# Example figures

# Figure of path coefficients
figurea <- mcot_estimates %>% select(m_upper, m_lower, p_upper,p_lower,  m,p, model) %>% 
  gather(key,value,-model) %>% mutate(parent=str_sub(key,1,1)) %>% 
  mutate(parm = case_when(str_detect(key,"upper") ~ "UpperCI", str_detect(key,"lower")~"LowerCI", TRUE ~"Point")) %>% 
  select(-key) %>% spread(parm,value) %>% 
  ggplot(aes(x=model, color=parent,ymin=LowerCI, ymax=UpperCI, y=Point)) + 
  geom_pointrange(position=position_dodge(width = 0.3)) +
  xlab("Model") +
  ylab("Std. regression coefficients") + 
  coord_flip() +
  theme_minimal() + theme(legend.position="bottom")

# Figure on explained family variance
figureb <- mcot_estimates %>% select(model, zVFo) %>% 
  ggplot(aes(x=model, fill=zVFo, y=zVFo)) +
  geom_bar(stat="identity")+
  xlab("Model") +
  ylab("Variance component F") + 
  coord_flip()+
  theme_minimal() + theme(legend.position="none")

mcotplot <- figurea + figureb + plot_layout(axes = "collect_x", axis_titles="collect_x")
ggsave(filename=here("mcotplot.png"), plot = mcotplot, device="png")
