library(tidyverse)
library(data.table)
library(here)
library(broom)
library(lm.beta)
library(ggpubr)
library(here)
library(haven)
library(MatchIt)
library(hrbrthemes)
library(estimatr)

# Load in educ, occ, inc, wealth -----------------------------------

# Faste - Fixed information

faste <- fread("N:/durable/data/registers/original/csv/w19_0634_faste_oppl_ut.csv")
faste$foedselsaar <- as.numeric(substr(faste$foedsels_aar_mnd, 1, 4))

faste <- faste %>% 
  group_by(mor_lnr, far_lnr) %>% 
  arrange(foedselsaar) %>% 
  mutate(familysize = n()) %>% 
  mutate(birthorder = row_number()) %>% 
  mutate(familysize = ifelse(familysize > 15, NA, familysize)) %>% 
  mutate(birthorder = ifelse(birthorder > 15, NA, birthorder)) %>% 
  ungroup()

# Regstat - Finding first date people appear in registers

regstat <- fread("N:/durable/data/registers/original/csv/w19_0634_regstat_1986_2018_ut.csv")
regstat <- left_join(regstat, select(faste, foedselsaar, w19_0634_lnr))
regstat <- select(regstat, w19_0634_lnr, foedselsaar, contains("dat"))
regstat <- distinct(regstat)
regstat <- gather(regstat, year, regstat, regstat_86_dat:regstat_18_dat, factor_key=TRUE)
regstat$year <- substr(regstat$year, 9, 10)
regstat$year <- ifelse(nchar(regstat$year) == 2, ifelse(as.numeric(regstat$year) > 50, paste0("19", regstat$year), paste0("20", regstat$year)), regstat$year)
regstat$year <- as.numeric(regstat$year)
regstat <- mutate(regstat, regstat = ifelse(regstat == 0, NA, regstat))
regstat$regstat <- as.numeric(substr(as.character(regstat$regstat), 1,4))
regstat <- filter(regstat, foedselsaar > 1899)
regstat <- filter(regstat, regstat > (foedselsaar - 2))

library(dtplyr)

regstat <- lazy_dt(regstat)

regstat <- regstat %>% 
  group_by(w19_0634_lnr) %>% 
  arrange(w19_0634_lnr, year) %>%
  mutate(regdato = first(na.omit(regstat))) %>%
  ungroup()

regstat <- as_tibble(regstat)

unloadNamespace("dtplyr")

regstat$regage <- regstat$regdato - regstat$foedselsaar

regstat <- select(regstat, w19_0634_lnr, regdato, regage)
regstat <- distinct(regstat)
faste <- left_join(faste, regstat)
rm(regstat)

# Education per age

educ <- fread("N:/durable/data/registers/original/csv/w19_0634_utd_1970_2018_ut.csv")
educ <- select(educ, w19_0634_lnr, contains("BU"))
educ <- distinct(educ)
educ <- gather(educ, year, BU, BU_1970:BU_2018, factor_key=TRUE)
educ$year <- substr(educ$year, 4, 7)
educ <- left_join(educ, select(faste, w19_0634_lnr, foedselsaar, kjoenn))
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


# Selecting adoptees, siblings and other Norwegians -----------------------

names(faste)
morfarland <- select(faste, w19_0634_lnr,fodeland)

korea <- 492
kohorter <- 1965:1985

alle <- faste %>% 
  select(contains("lnr"), fodeland, foedselsaar, foedsels_aar_mnd, regage, regdato, familysize, birthorder, fodeland, kjoenn) %>% 
  filter(foedselsaar %in% kohorter) %>% 
  left_join(rename(morfarland,mors_fodeland = fodeland), by=c(mor_lnr = "w19_0634_lnr")) %>% 
  left_join(rename(morfarland,fars_fodeland = fodeland), by=c(far_lnr = "w19_0634_lnr"))

norskbakgrunn <- alle %>% 
  filter(mors_fodeland==0, fars_fodeland==0) 

basic_samples <- norskbakgrunn %>% 
  mutate(koreaner = if_else(fodeland==korea, 1, 0)) %>% 
  mutate(nordmann = if_else(fodeland==0, 1, 0)) %>% 
  mutate(pct = runif(n = nrow(.),min = 0,max = 1)) %>% 
  filter(koreaner==1 | (nordmann==1 & pct<0.05))

norskesosken <- basic_samples %>% 
  filter(koreaner==1) %>% 
  select(mor_lnr, far_lnr) %>% 
  distinct %>% 
  left_join(faste, by=c("mor_lnr","far_lnr")) %>% 
  filter(fodeland==0) %>% 
  mutate(norsksosken = 1) %>% 
  mutate(nordmann = 1)

samples <- bind_rows(basic_samples, norskesosken) %>% 
  replace_na(replace = list(koreaner = 0, norsksosken = 0, nordmann = 0)) %>% 
  group_by(mor_lnr,far_lnr) %>% 
  mutate(sibfamily = max(norsksosken)) %>% 
  ungroup()

persons <- select(samples, w19_0634_lnr, mor_lnr, far_lnr) %>% 
  gather(role, person) %>% 
  select(-role) %>% 
  arrange(person) %>% 
  distinct %>% 
  rename(w19_0634_lnr = person)


# Merging adoptees with socio-economic data -------------------------------

incomerankdata <- persons %>% 
  left_join(adult_income_ranks) %>% 
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("income_rank",age)) %>% 
  spread(age,income_age_cohort_sex_rank) %>% 
  ungroup()

educdata <- persons %>% 
  left_join(educ) %>% 
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("education",age)) %>% 
  spread(age,BU) %>% 
  ungroup()

occdata <- persons %>% 
  left_join(occupation) %>% 
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("trei",age)) %>% 
  relocate(trei, .after = age) %>% 
  spread(age,trei) %>% 
  ungroup()

wealthdata <- persons %>% 
  left_join(wealth) %>% 
  select(-gross_wealth_age_cohort_sex_rank) %>% 
  group_by(w19_0634_lnr) %>% 
  mutate(age = str_c("wealth_rank",age)) %>% 
  spread(age,wealth_age_cohort_sex_rank) %>% 
  ungroup()
  
gwealthdata <- persons %>% 
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
  left_join(select(faste, w19_0634_lnr, foedselsaar, fodeland, kjoenn, regage, regdato, familysize, birthorder))

prefixit <- function(df, prefix) {
  nm <- names(df)
  new_nm <- str_c(prefix, nm)
  setNames(df, new_nm)
}

linked <-  samples %>% 
  left_join(variables) %>% 
  left_join(prefixit(variables,"father_"), by = c(far_lnr = "father_w19_0634_lnr")) %>% 
  left_join(prefixit(variables,"mother_"), by = c(mor_lnr = "mother_w19_0634_lnr"))
nrow(linked)
names(linked)

saveRDS(linked, "N:/durable/users/arnovh/Adoptees/Data/fulldata_adoptees.rds")

analytical <- linked %>%  
  select(w19_0634_lnr, mor_lnr, far_lnr, foedselsaar, regage, regdato, familysize, birthorder, koreaner, norsksosken, sibfamily, nordmann, kjoenn, 
         education30:education55, income_rank30:income_rank45,
         trei30:trei45, wealth_rank30:wealth_rank45, gross_wealth_rank30:gross_wealth_rank45,
         father_foedselsaar, 
         mother_foedselsaar, father_education30:father_education65, 
         contains("father_income_rank"),
         mother_education30:mother_education65, 
         contains("mother_income_rank"), 
         contains("mother_trei"),
         contains("father_trei"),
         contains("mother_wealth_rank"),
         contains("father_wealth_rank"), 
         contains("mother_gross_wealth_rank"),
         contains("father_gross_wealth_rank")) %>% 
  mutate(child_ed = rowMeans(across(education30:education45), na.rm = TRUE)) %>% 
  mutate(father_ed = rowMeans(across(father_education45:father_education60), na.rm = TRUE))%>% 
  mutate(mother_ed = rowMeans(across(mother_education45:mother_education60), na.rm = TRUE)) %>% 
  mutate(parent_ed = case_when(is.na(father_ed) ~ mother_ed, is.na(mother_ed) ~ father_ed, 
                               father_ed>mother_ed ~ father_ed, TRUE ~ mother_ed)) %>% 
  mutate(child_irank = 100*(rowMeans(across(income_rank30:income_rank45), na.rm = TRUE))) %>%
  mutate(father_irank = 100*(rowMeans(across(father_income_rank45:father_income_rank60), na.rm = TRUE))) %>%
  mutate(mother_irank = 100*(rowMeans(across(mother_income_rank45:mother_income_rank60), na.rm = TRUE))) %>%
  mutate(parent_irank = case_when(is.na(father_irank) ~ mother_irank, is.na(mother_irank) ~ father_irank, 
                                  father_irank>mother_irank ~ father_irank, TRUE ~ mother_irank)) %>% 
  mutate(child_occ = rowMeans(across(trei30:trei45), na.rm = TRUE)) %>%
  mutate(father_occ = rowMeans(across(father_trei50:father_trei65), na.rm = TRUE)) %>%
  mutate(mother_occ = rowMeans(across(mother_trei50:mother_trei65), na.rm = TRUE)) %>%
  mutate(parent_occ = case_when(is.na(father_occ) ~ mother_occ, is.na(mother_occ) ~ father_occ, 
                                  father_occ>mother_occ ~ father_occ, TRUE ~ mother_occ)) %>% 
  mutate(child_wealth = 100*(rowMeans(across(wealth_rank30:wealth_rank45), na.rm = TRUE))) %>%
  mutate(father_wealth = 100*(rowMeans(across(father_wealth_rank45:father_wealth_rank60), na.rm = TRUE))) %>%
  mutate(mother_wealth = 100*(rowMeans(across(mother_wealth_rank45:mother_wealth_rank60), na.rm = TRUE))) %>%
  mutate(parent_wealth = case_when(is.na(father_wealth) ~ mother_wealth, is.na(mother_wealth) ~ father_wealth, 
                                   father_wealth>mother_wealth ~ father_wealth, TRUE ~ mother_wealth)) %>% 
  mutate(child_gwealth = 100*(rowMeans(across(gross_wealth_rank30:gross_wealth_rank45), na.rm = TRUE))) %>%
  mutate(father_gwealth = 100*(rowMeans(across(father_gross_wealth_rank45:father_gross_wealth_rank60), na.rm = TRUE))) %>%
  mutate(mother_gwealth = 100*(rowMeans(across(mother_gross_wealth_rank45:mother_gross_wealth_rank60), na.rm = TRUE))) %>%
  mutate(parent_gwealth = case_when(is.na(father_gwealth) ~ mother_gwealth, is.na(mother_gwealth) ~ father_gwealth, 
                                   father_gwealth>mother_gwealth ~ father_gwealth, TRUE ~ mother_gwealth)) 


saveRDS(analytical, "N:/durable/users/arnovh/Adoptees/Data/socioeco_adoptees.rds")

# Add marital status parents at age 16

status <- fread("N:/durable/data/registers/original/csv/w19_0634_sivilstand_1975_2018_ut.csv")

status <- gather(status, year, sivilstand, sivilstand_1975:sivilstand_2018, factor_key=TRUE)
status$year <- substr(status$year, 12, 15)
status$year <- as.numeric(as.character(status$year))

analytical$year16 <- analytical$foedselsaar + 16

analytical <- left_join(analytical, status, by = c("mor_lnr" = "w19_0634_lnr", "year16" = "year"))

analytical <- mutate(analytical, married = ifelse(sivilstand == 2, 1, 0))
analytical <- mutate(analytical, married = ifelse(is.na(married), 0, married))

# Recode birth order extra variable

analytical <- analytical %>% 
  mutate(birth_ordered = case_when(
    birthorder == 1 & familysize == 1 ~ "One of one", 
    birthorder == 1 & familysize == 2 ~ "One of two", 
    birthorder == 1 & familysize == 3 ~ "One of three", 
    birthorder == 1 & familysize == 4 ~ "One of four", 
    birthorder == 1 & familysize > 4 ~ "One of many", 
    birthorder == 2 & familysize == 2 ~ "Two of two", 
    birthorder == 2 & familysize == 3 ~ "Two of three", 
    birthorder == 2 & familysize == 4 ~ "Two of four", 
    birthorder == 2 & familysize > 4 ~ "Two of many", 
    birthorder == 3 & familysize == 3 ~ "Three of three", 
    birthorder == 3 & familysize == 4 ~ "Three of four",
    birthorder == 3 & familysize > 4 ~ "Three of many",
    birthorder == 4 & familysize == 4 ~ "Four of four", 
    birthorder == 4 & familysize > 4 ~ "Four of many", 
    birthorder > 4 ~ "Later born"))

analytical$birth_ordered <- as.factor(analytical$birth_ordered)

# Define models -----------------------------------------------------------

analytical <- analytical %>% 
  mutate(child_eds = scale.default(child_ed)) %>% 
  mutate(mother_eds = scale.default(mother_ed)) %>% 
  mutate(father_eds = scale.default(father_ed)) %>% 
  mutate(child_occs = scale.default(child_occ)) %>% 
  mutate(mother_occs = scale.default(mother_occ)) %>% 
  mutate(father_occs = scale.default(father_occ)) %>% 
  mutate(child_iranks = scale.default(child_irank)) %>% 
  mutate(mother_iranks = scale.default(mother_irank)) %>% 
  mutate(father_iranks = scale.default(father_irank)) %>% 
  mutate(child_gwealths = scale.default(child_gwealth)) %>% 
  mutate(mother_gwealths = scale.default(mother_gwealth)) %>% 
  mutate(father_gwealths = scale.default(father_gwealth)) %>% 
  mutate(child_wealths = scale.default(child_wealth)) %>% 
  mutate(mother_wealths = scale.default(mother_wealth)) %>% 
  mutate(father_wealths = scale.default(father_wealth))

analytical <- analytical %>% group_by(mor_lnr, far_lnr) %>% mutate(family_id = cur_group_id()) %>% ungroup()

outcomes <- c("Education" = "child_ed", "Income rank" = "child_irank", "Occupation" = "child_occ", "Wealth" = "child_wealth", "Elite" = "child_elite")
regressors <- c("Education" = "parent_ed", "Income rank" = "parent_irank", "Occupation" = "parent_occ",  "Wealth" = "parent_wealth", "Elite" = "parent_elite")

samples = c("Biological" = 1,
            "Adoptees" = 2,
            "Adoptees, <2yrs" = 3,
            "Adoptee's siblings" = 4)

model_defs <- list(child_eds ~ mother_eds + father_eds,
                   child_eds ~ mother_occs + father_occs,
                   child_eds ~ mother_iranks + father_iranks,
                   child_eds ~ mother_wealths + father_wealths,
                   child_eds ~ mother_gwealths + father_gwealths,
                   child_occs ~ mother_eds + father_eds,
                   child_occs ~ mother_occs + father_occs,
                   child_occs ~ mother_iranks + father_iranks,
                   child_occs ~ mother_wealths + father_wealths,
                   child_occs ~ mother_gwealths + father_gwealths,
                   child_iranks ~ mother_eds + father_eds,
                   child_iranks ~ mother_occs + father_occs,
                   child_iranks ~ mother_iranks + father_iranks,
                   child_iranks ~ mother_wealths + father_wealths,
                   child_iranks ~ mother_gwealths + father_gwealths,
                   child_wealths ~ mother_eds + father_eds,
                   child_wealths ~ mother_occs + father_occs,
                   child_wealths ~ mother_iranks + father_iranks,
                   child_wealths ~ mother_wealths + father_wealths,
                   child_wealths ~ mother_gwealths + father_gwealths,
                   child_gwealths ~ mother_eds + father_eds,
                   child_gwealths ~ mother_occs + father_occs,
                   child_gwealths ~ mother_iranks + father_iranks,
                   child_gwealths ~ mother_wealths + father_wealths,
                   child_gwealths ~ mother_gwealths + father_gwealths,
                   child_eds ~ birth_ordered,
                   child_occs ~ birth_ordered,
                   child_iranks ~ birth_ordered,
                   child_wealths ~ birth_ordered,
                   child_gwealths ~ birth_ordered,
                   child_eds ~ kjoenn,
                   child_occs ~ kjoenn,
                   child_iranks ~ kjoenn,
                   child_wealths ~ kjoenn,
                   child_gwealths ~ kjoenn,
                   child_eds ~ mother_eds,
                   child_eds ~ mother_occs,
                   child_eds ~ mother_iranks,
                   child_eds ~ mother_wealths,
                   child_eds ~ mother_gwealths,
                   child_occs ~ mother_eds,
                   child_occs ~ mother_occs,
                   child_occs ~ mother_iranks,
                   child_occs ~ mother_wealths,
                   child_occs ~ mother_gwealths,
                   child_iranks ~ mother_eds,
                   child_iranks ~ mother_occs,
                   child_iranks ~ mother_iranks,
                   child_iranks ~ mother_wealths,
                   child_iranks ~ mother_gwealths,
                   child_wealths ~ mother_eds,
                   child_wealths ~ mother_occs,
                   child_wealths ~ mother_iranks,
                   child_wealths ~ mother_wealths,
                   child_wealths ~ mother_gwealths,
                   child_gwealths ~ mother_eds,
                   child_gwealths ~ mother_occs,
                   child_gwealths ~ mother_iranks,
                   child_gwealths ~ mother_wealths,
                   child_gwealths ~ mother_gwealths,
                   child_eds ~ father_eds,
                   child_eds ~ father_occs,
                   child_eds ~ father_iranks,
                   child_eds ~ father_wealths,
                   child_eds ~ father_gwealths,
                   child_occs ~ father_eds,
                   child_occs ~ father_occs,
                   child_occs ~ father_iranks,
                   child_occs ~ father_wealths,
                   child_occs ~ father_gwealths,
                   child_iranks ~ father_eds,
                   child_iranks ~ father_occs,
                   child_iranks ~ father_iranks,
                   child_iranks ~ father_wealths,
                   child_iranks ~ father_gwealths,
                   child_wealths ~ father_eds,
                   child_wealths ~ father_occs,
                   child_wealths ~ father_iranks,
                   child_wealths ~ father_wealths,
                   child_wealths ~ father_gwealths,
                   child_gwealths ~ father_eds,
                   child_gwealths ~ father_occs,
                   child_gwealths ~ father_iranks,
                   child_gwealths ~ father_wealths,
                   child_gwealths ~ father_gwealths)


biol_models <- map(model_defs, ~lm_robust(data=analytical, formula=.x, clusters = family_id, se_type = "stata", subset=(nordmann==1))) 
adopt_models <- map(model_defs, ~lm_robust(data=analytical, formula=.x, clusters = family_id, se_type = "stata", subset=(koreaner==1 & regage<3))) 
adoptsib_models <- map(model_defs, ~lm_robust(data=analytical, formula=.x, clusters = family_id, se_type = "stata", subset=(koreaner==1 & sibfamily == 1 & regage<3))) 
sib_models <- map(model_defs, ~lm_robust(data=analytical, formula=.x, clusters = family_id, se_type = "stata", subset=(norsksosken == 1))) 


all_models <- list("1. Birth" = biol_models,"2. Adoptees" =  adopt_models,
                   "3. Adoptees with siblings" = adoptsib_models, "4. Siblings" = sib_models)
all_results <- map(all_models, ~bind_rows(map(.x, ~tidy((.x), conf.int = TRUE)), .id="Model")) %>% bind_rows(.id="Sample")

fit_statistics <- map(all_models, ~bind_rows(map(.x, ~ cbind(glance(.x))), .id="Model")) %>% bind_rows(.id="Sample")

all_results <- all_results %>% group_by(Sample, Model) %>% mutate(nocoef = n()) %>% ungroup()

all_results <- all_results %>% 
  mutate(independent = case_when(
    str_detect(term, "ed") ~ "Education",
    str_detect(term, "gwealth") ~ "Gross wealth",
    str_detect(term, "occ") ~ "Occupation",
    str_detect(term, "_wealths") ~ "Wealth",
    str_detect(term, "irank") ~ "Income",
    str_detect(term, "birth") ~ "Birthorder",
    str_detect(term, "kjoenn") ~ "Gender"))

all_results <- all_results %>% 
  mutate(dependent = case_when(
    str_detect(outcome, "ed") ~ "Education",
    str_detect(outcome, "gwealth") ~ "Gross Wealth",
    str_detect(outcome, "occ") ~ "Occupation",
    str_detect(outcome, "_wealths") ~ "Wealth",
    str_detect(outcome, "irank") ~ "Income"))

all_results <- all_results %>% 
  mutate(Model = paste(dependent, independent, sep = ":"))

all_results <- mutate(all_results, dependent = sub(":.*", "", all_results$Model)) 
all_results <- mutate(all_results, independent = sub(".*:", "", all_results$Model)) 

saveRDS(all_results, "N:/durable/users/arnovh/Adoptees/Data/all_results.rds")
saveRDS(fit_statistics, "N:/durable/users/arnovh/Adoptees/Data/fit_statistics.rds")

# Matching ----------------------------------------------------------------

df <- filter(analytical, norsksosken==0)
df <- mutate(df, excl = ifelse((koreaner == 1 & regage > 2), 1, 0))
df <- filter(df, excl == 0)

df <- df %>% filter(!is.na(kjoenn), !is.na(parent_ed), !is.na(parent_irank), !is.na(parent_occ), !is.na(parent_gwealth), !is.na(parent_wealth), !is.na(mother_foedselsaar), !is.na(father_foedselsaar), !is.na(birthorder), !is.na(familysize))
df <- matchit(koreaner ~ kjoenn + foedselsaar + familysize + birthorder + parent_ed + parent_irank + parent_occ + parent_gwealth + parent_wealth + mother_foedselsaar + father_foedselsaar + married, data = df, method = "nearest")
summary(df)

plot(summary(df))
plot(df, type = "jitter")
plot(df, type = "density", which.xs = ~foedselsaar)

df <- match.data(df)

saveRDS(df, "N:/durable/users/arnovh/Adoptees/Data/socioeco_adoptees_matchedsamples.rds")

biol_models2 <- map(model_defs, ~lm_robust(data=df, formula=.x,clusters = family_id, se_type = "stata", subset=(nordmann==1))) 
adopt_models2 <- map(model_defs, ~lm_robust(data=df, formula=.x, clusters = family_id, se_type = "stata", subset=(koreaner==1))) 

all_models2 <- list("1. Birth matched" = biol_models2,"2. Adoptees matched" =  adopt_models2)
all_results2 <- map(all_models2, ~bind_rows(map(.x, ~tidy((.x),conf.int=TRUE)), .id="Model")) %>% bind_rows(.id="Sample")

fit_statistics2 <- map(all_models2, ~bind_rows(map(.x, ~ cbind(glance(.x))), .id="Model")) %>% bind_rows(.id="Sample")

all_results2 <- all_results2 %>% group_by(Sample, Model) %>% mutate(nocoef = n()) %>% ungroup()

all_results2 <- all_results2 %>% 
  mutate(independent = case_when(
    str_detect(term, "ed") ~ "Education",
    str_detect(term, "gwealth") ~ "Gross Wealth",
    str_detect(term, "occ") ~ "Occupation",
    str_detect(term, "_wealths") ~ "Wealth",
    str_detect(term, "irank") ~ "Income",
    str_detect(term, "birth") ~ "Birthorder",
    str_detect(term, "kjoenn") ~ "Gender"))

all_results2 <- all_results2 %>% 
  mutate(dependent = case_when(
    str_detect(outcome, "ed") ~ "Education",
    str_detect(outcome, "gwealth") ~ "Gross Wealth",
    str_detect(outcome, "occ") ~ "Occupation",
    str_detect(outcome, "_wealths") ~ "Wealth",
    str_detect(outcome, "irank") ~ "Income"))

all_results2 <- all_results2 %>% 
  mutate(Model = paste(dependent, independent, sep = ":"))

all_results2 <- mutate(all_results2, dependent = sub(":.*", "", all_results2$Model)) 
all_results2 <- mutate(all_results2, independent = sub(".*:", "", all_results2$Model)) 

saveRDS(all_results2, "N:/durable/users/arnovh/Adoptees/Data/all_results2.rds")
saveRDS(fit_statistics2, "N:/durable/users/arnovh/Adoptees/Data/fit_statistics2.rds")

# Plotting distributions --------------------------------------------------

# Distributions unmatched samples

educ <- analytical %>% 
  filter(norsksosken == 0) %>% 
  filter(!(koreaner == 1 & regage > 2)) %>% 
  ggplot(aes(x=child_ed, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(name = "Group", values = c("#0072B2", "#F0E442"), labels = c("Norwegian", "Adopted")) +
  labs(x = "Child education",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

inc <- analytical %>% 
  filter(norsksosken == 0) %>% 
  filter(!(koreaner == 1 & regage > 2)) %>% 
  ggplot(aes(x=child_irank, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(name = "Group", values = c("#0072B2", "#F0E442"), labels = c("Norwegian", "Adopted")) +
  labs(x = "Child income rank",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

occ <- analytical %>% 
  filter(norsksosken == 0) %>% 
  filter(!(koreaner == 1 & regage > 2)) %>% 
  ggplot(aes(x=child_occ, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(name = "Group", values = c("#0072B2", "#F0E442"), labels = c("Norwegian", "Adopted")) +
  labs(x = "Child occupation",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

wealth <- analytical %>% 
  filter(norsksosken == 0) %>% 
  filter(!(koreaner == 1 & regage > 2)) %>% 
  ggplot(aes(x=child_wealth, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(name = "Group", values = c("#0072B2", "#F0E442"), labels = c("Norwegian", "Adopted")) +
  labs(x = "Child wealth",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

ggarrange(educ, inc, occ, wealth, nrow = 2, ncol = 2, legend = "bottom", common.legend = TRUE)
ggsave("N:/durable/users/arnovh/Adoptees/Results/Distribution unmatched.jpeg", width = 30, height = 15, units = "cm", dpi = 300)


# Distribution matched samples

educ <- df %>% 
  ggplot(aes(x=child_ed, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(name = "Group", values = c("#0072B2", "#F0E442"), labels = c("Norwegian", "Adopted")) +
  labs(x = "Child education",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

inc <- df %>% 
  ggplot(aes(x=child_irank, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(name = "Group", values = c("#0072B2", "#F0E442"), labels = c("Norwegian", "Adopted")) +
  labs(x = "Child income rank",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

occ <- df %>% 
  ggplot(aes(x=child_occ, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(name = "Group", values = c("#0072B2", "#F0E442"), labels = c("Norwegian", "Adopted")) +
  labs(x = "Child occupation",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

wealth <- df %>% 
  ggplot(aes(x=child_wealth, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(name = "Group", values = c("#0072B2", "#F0E442"), labels = c("Norwegian", "Adopted")) +
  labs(x = "Child wealth",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

ggarrange(educ, inc, occ, wealth, nrow = 2, ncol = 2, legend = "bottom", common.legend = TRUE)
ggsave("N:/durable/users/arnovh/Adoptees/Results/Distribution matched.jpeg", width = 30, height = 15, units = "cm", dpi = 300)

# Distribution siblings vs adoptees

educ <- analytical %>% 
  filter(norsksosken == 1 | koreaner == 1) %>% 
  filter(!(koreaner == 1 & sibfamily == 0)) %>% 
  filter(!(koreaner == 1 & regage > 2)) %>% 
  ggplot(aes(x=child_ed, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_discrete(name = "Group", labels = c("Sibling", "Adopted")) +
  labs(x = "Child education",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

inc <- analytical %>% 
  filter(norsksosken == 1 | koreaner == 1) %>% 
  filter(!(koreaner == 1 & sibfamily == 0)) %>% 
  filter(!(koreaner == 1 & regage > 2)) %>% 
  ggplot(aes(x=child_irank, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_discrete(name = "Group", labels = c("Norwegian", "Adopted")) +
  labs(x = "Child income rank",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

occ <- analytical %>% 
  filter(norsksosken == 1 | koreaner == 1) %>% 
  filter(!(koreaner == 1 & sibfamily == 0)) %>% 
  filter(!(koreaner == 1 & regage > 2)) %>% 
  ggplot(aes(x=child_occ, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_discrete(name = "Group", labels = c("Norwegian", "Adopted")) +
  labs(x = "Child occupation",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

wealth <- analytical %>% 
  filter(norsksosken == 1 | koreaner == 1) %>% 
  filter(!(koreaner == 1 & sibfamily == 0)) %>% 
  filter(!(koreaner == 1 & regage > 2)) %>% 
  ggplot(aes(x=child_wealth, group=factor(koreaner), fill=factor(koreaner))) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_discrete(name = "Group", labels = c("Norwegian", "Adopted")) +
  labs(x = "Child wealth",
       y = "Density") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme_classic()

ggarrange(educ, inc, occ, wealth, nrow = 2, ncol = 2, legend = "bottom", common.legend = TRUE)
ggsave("N:/durable/users/arnovh/Adoptees/Results/Distribution siblings.jpeg", width = 30, height = 15, units = "cm", dpi = 300)


# Making one dataset ------------------------------------------------------

all_results3 <- rbind(all_results, all_results2)
fit_statistics3 <- rbind(fit_statistics, fit_statistics2)

write.csv(all_results3, "N:/durable/users/arnovh/Adoptees/Data/all_model_results.csv")
write.csv(fit_statistics3, "N:/durable/users/arnovh/Adoptees/Data/all_model_fit.csv")

# Checking differences between age of adoption  -----------------------------------------------------------

adoptold_models <- map(model_defs, ~lm_robust(data=df, formula=.x,clusters = family_id, se_type = "stata", subset=(koreaner==1 & regage > 0))) 
adoptyoung_models <- map(model_defs, ~lm_robust(data=df, formula=.x, clusters = family_id, se_type = "stata", subset=(koreaner==1 & regage == 0))) 

all_models4 <- list("1. After 1" = adoptold_models,"2. Before 1" =  adoptyoung_models)
all_results4 <- map(all_models4, ~bind_rows(map(.x, ~tidy((.x),conf.int=TRUE)), .id="Model")) %>% bind_rows(.id="Sample")

fit_statistics4 <- map(all_models4, ~bind_rows(map(.x, ~ cbind(glance(.x))), .id="Model")) %>% bind_rows(.id="Sample")

all_results4 <- all_results4 %>% group_by(Sample, Model) %>% mutate(nocoef = n()) %>% ungroup()

all_results4 <- all_results4 %>% 
  mutate(independent = case_when(
    str_detect(term, "ed") ~ "Education",
    str_detect(term, "gwealth") ~ "Gross Wealth",
    str_detect(term, "occ") ~ "Occupation",
    str_detect(term, "_wealths") ~ "Wealth",
    str_detect(term, "irank") ~ "Income",
    str_detect(term, "birth") ~ "Birthorder",
    str_detect(term, "kjoenn") ~ "Gender"))

all_results4 <- all_results4 %>% 
  mutate(dependent = case_when(
    str_detect(outcome, "ed") ~ "Education",
    str_detect(outcome, "gwealth") ~ "Gross Wealth",
    str_detect(outcome, "occ") ~ "Occupation",
    str_detect(outcome, "_wealths") ~ "Wealth",
    str_detect(outcome, "irank") ~ "Income"))

all_results4 <- all_results4 %>% 
  mutate(Model = paste(dependent, independent, sep = ":"))

all_results4 <- mutate(all_results4, dependent = sub(":.*", "", all_results4$Model)) 
all_results4 <- mutate(all_results4, independent = sub(".*:", "", all_results4$Model)) 

write.csv(all_results3, "N:/durable/users/arnovh/Adoptees/Data/all_model_results_age.csv")
write.csv(fit_statistics3, "N:/durable/users/arnovh/Adoptees/Data/all_model_fit_age.csv")


# Checking differences between adoptees of different genders -------------------------------------------

bioboy_models <- map(model_defs, ~lm_robust(data=df, formula=.x,clusters = family_id, se_type = "stata", subset=(koreaner==0 & kjoenn == 1))) 
biogirl_models <- map(model_defs, ~lm_robust(data=df, formula=.x, clusters = family_id, se_type = "stata", subset=(koreaner==0 & kjoenn == 2))) 
adoptboy_models <- map(model_defs, ~lm_robust(data=df, formula=.x,clusters = family_id, se_type = "stata", subset=(koreaner==1 & kjoenn == 1))) 
adoptgirl_models <- map(model_defs, ~lm_robust(data=df, formula=.x, clusters = family_id, se_type = "stata", subset=(koreaner==1 & kjoenn == 2))) 


all_models5 <- list("1. Birth boys" = bioboy_models,"2. Birth girls" =  biogirl_models, 
                    "3. Adopted boys" = adoptboy_models,"4. Adopted girls" =  adoptgirl_models)
all_results5 <- map(all_models5, ~bind_rows(map(.x, ~tidy((.x),conf.int=TRUE)), .id="Model")) %>% bind_rows(.id="Sample")

fit_statistics5 <- map(all_models5, ~bind_rows(map(.x, ~ cbind(glance(.x))), .id="Model")) %>% bind_rows(.id="Sample")

all_results5 <- all_results5 %>% group_by(Sample, Model) %>% mutate(nocoef = n()) %>% ungroup()

all_results5 <- all_results5 %>% 
  mutate(independent = case_when(
    str_detect(term, "ed") ~ "Education",
    str_detect(term, "gwealth") ~ "Gross Wealth",
    str_detect(term, "occ") ~ "Occupation",
    str_detect(term, "_wealths") ~ "Wealth",
    str_detect(term, "irank") ~ "Income",
    str_detect(term, "birth") ~ "Birthorder",
    str_detect(term, "kjoenn") ~ "Gender"))

all_results5 <- all_results5 %>% 
  mutate(dependent = case_when(
    str_detect(outcome, "ed") ~ "Education",
    str_detect(outcome, "gwealth") ~ "Gross Wealth",
    str_detect(outcome, "occ") ~ "Occupation",
    str_detect(outcome, "_wealths") ~ "Wealth",
    str_detect(outcome, "irank") ~ "Income"))

all_results5 <- all_results5 %>% 
  mutate(Model = paste(dependent, independent, sep = ":"))

all_results5 <- mutate(all_results5, dependent = sub(":.*", "", all_results5$Model)) 
all_results5 <- mutate(all_results5, independent = sub(".*:", "", all_results5$Model)) 

write.csv(all_results5, "N:/durable/users/arnovh/Adoptees/Data/all_model_results_gender.csv")
write.csv(fit_statistics5, "N:/durable/users/arnovh/Adoptees/Data/all_model_fit_gender.csv")

