library(tidyverse)
library(data.table)
library(here)
library(broom)
library(lm.beta)
library(ggpubr)
library(here)
library(haven)

# Load in educ, occ, inc, wealth -----------------------------------

# Faste

faste <- fread("N:/durable/data/registers/original/csv/w19_0634_faste_oppl_ut.csv")
faste$foedselsaar <- as.numeric(substr(faste$foedsels_aar_mnd, 1, 4))


# Education per age

educ <- fread("N:/durable/data/registers/original/csv/w19_0634_utd_1970_2018_ut.csv")
educ <- select(educ, w19_0634_lnr, contains("BU"))
educ <- distinct(educ)
educ <- gather(educ, year, BU, BU_1970:BU_2018, factor_key=TRUE)
educ$year <- substr(educ$year, 4, 7)
educ <- left_join(educ, select(faste, w19_0634_lnr, foedselsaar))
educ$year <- as.numeric(educ$year)
educ$age <- educ$year - educ$foedselsaar
educ <- filter(educ, age < 66 & age > 20)
educ <- filter(educ, !is.na(age))
educ <- select(educ, w19_0634_lnr, BU, age)
educ$BU <- as.numeric(substr(as.character(educ$BU), 1,1))
educ <- mutate(educ, BU = case_when(BU==2 ~ 9,
                                    BU==3 ~ 11,
                                    BU==4 ~ 13,
                                    BU==5 ~ 14,
                                    BU==6 ~ 16,
                                    BU==7 ~ 18,
                                    BU==8 ~ 21))

library(dtplyr)

educ <- lazy_dt(educ)

educ <- educ %>% 
  group_by(w19_0634_lnr, age) %>% 
  arrange(w19_0634_lnr, age, desc(BU)) %>%
  slice(1) %>% 
  ungroup()

educ <- as_tibble(educ)

unloadNamespace("dtplyr")

# Income per age

adult_income_ranks <- fread("N:/durable/projects/openflux/data/income-within-age-cohort-sex-ranks.csv")
adult_income_ranks <- filter(adult_income_ranks, age < 66 & age > 20)
adult_income_ranks <- filter(adult_income_ranks, !is.na(age))
adult_income_ranks <- select(adult_income_ranks, -income_age_cohort_rank)

library(dtplyr)

adult_income_ranks <- lazy_dt(adult_income_ranks)

adult_income_ranks <- adult_income_ranks %>% 
  group_by(w19_0634_lnr, age) %>% 
  arrange(w19_0634_lnr, age, desc(income_age_cohort_sex_rank)) %>%
  slice(1) %>% 
  ungroup()

adult_income_ranks <- as_tibble(adult_income_ranks)

unloadNamespace("dtplyr")

# Occupation

occupation <- fread("N:/durable/projects/joakim/inheritance-of-inequality/data/occ_trei.csv")
occupation$w19_0634_lnr <- occupation$lnr
occupation <- left_join(occupation, faste, by = "w19_0634_lnr")
occupation$age <- occupation$year - occupation$foedselsaar
occupation <- filter(occupation, age < 66 & age > 20)
occupation <- filter(occupation, year > 2002)
occupation <- filter(occupation, !is.na(age))
occupation <- select(occupation, w19_0634_lnr, trei, age)
occupation <- distinct(occupation)

library(dtplyr)

occupation <- lazy_dt(occupation)

occupation <- occupation %>% 
  group_by(w19_0634_lnr, age) %>% 
  arrange(w19_0634_lnr, age, desc(trei)) %>%
  slice(1) %>% 
  ungroup()

occupation <- as_tibble(occupation)

unloadNamespace("dtplyr")

# Wealth

wealth <- fread("N:/durable/projects/openflux/data/wealth-within-age-cohort-sex-ranks.csv")
wealth <- filter(wealth, age < 66 & age > 20)
wealth <- filter(wealth, !is.na(age))
wealth <- select(wealth, -wealth_age_cohort_rank, -gross_wealth_age_cohort_rank)

library(dtplyr)

wealth <- lazy_dt(wealth)

wealth <- wealth %>% 
  group_by(w19_0634_lnr, age) %>% 
  arrange(w19_0634_lnr, age, desc(wealth_age_cohort_sex_rank)) %>%
  slice(1) %>% 
  ungroup()

wealth <- as_tibble(wealth)

unloadNamespace("dtplyr")


# Load in twins data ------------------------------------------------------

twins <- fread("N:/durable/users/torkildl/mcot/persons-involved.csv")
twins <- distinct(twins)
#twins <- twins %>% group_by(w19_0634_lnr) %>% slice(1)

incomerankdata <- twins %>% 
  left_join(adult_income_ranks) %>%
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("income_rank",age)) %>% 
  spread(age,income_age_cohort_sex_rank) %>% 
  ungroup()

educdata <- twins %>% 
  left_join(educ) %>%
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("education",age)) %>% 
  spread(age,BU) %>% 
  ungroup()

occdata <- twins %>% 
  left_join(occupation) %>% 
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("trei",age)) %>% 
  relocate(trei, .after = age) %>% 
  spread(age,trei) %>% 
  ungroup()

wealthdata <- twins %>% 
  left_join(wealth) %>% 
  select(-gross_wealth_age_cohort_sex_rank) %>% 
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("wealth_rank",age)) %>% 
  spread(age,wealth_age_cohort_sex_rank) %>% 
  ungroup()

gwealthdata <- twins %>% 
  left_join(wealth) %>% 
  select(-wealth_age_cohort_sex_rank) %>% 
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("gross_wealth_rank",age)) %>% 
  spread(age,gross_wealth_age_cohort_sex_rank) %>% 
  ungroup()

variables <- incomerankdata %>% 
  ungroup %>% 
  left_join(educdata, by = "w19_0634_lnr") %>% 
  left_join(occdata, by = "w19_0634_lnr") %>% 
  left_join(wealthdata, by = "w19_0634_lnr") %>% 
  left_join(gwealthdata, by = "w19_0634_lnr") %>% 
  left_join(select(faste, w19_0634_lnr, foedselsaar))

variables <- select(variables, -"<NA>", - "<NA>.x.x", - "<NA>.y.y", - "<NA>.y", - "<NA>.x")

write.csv(variables, "N:/durable/users/arnovh/Adoptees/Data/MCOT_socioeco.csv")
saveRDS(variables, "N:/durable/users/arnovh/Adoptees/Data/MCOT_socioeco.rds")
