library(tidyverse)
library(data.table)
library(here)
library(broom)
library(performance)

# regdata = "/tsd/p805/data/durable/data/registers/original/csv/"
# mypath <- function(x) {str_c(regdata,x)}

### Standard demographic information incl mother/father links
faste <- fread("/tsd/p805/data/durable/data/registers/SSB/01_data/data_v1.0/csv/POPULATION_FASTE_OPPL.csv")

### Twin registry data
twins <- fread("/tsd/p805/data/durable/data/registers/twin registry/PDB2601_NTR_2023.csv") %>% 
  setNames(nm=c("twinpair","twinid","twinno","w19_0634_lnr","Zygo","Group","DNAZygo"))

### Fertility histories
fertility <- fread("/tsd/p805/data/durable/projects/openflux/data/fertility-histories.csv")

### All individuals
persons <- twins %>% 
  left_join(faste) %>% 
  left_join(fertility) %>% 
  select(contains("_lnr"), contains("_otherparent")) %>%
  rename(twin = w19_0634_lnr) %>% 
  gather(key=role, value=w19_0634_lnr) %>% 
  filter(w19_0634_lnr!="") %>% 
  select(-role)

# Write out file with ids for persons involved
fwrite(persons, here("persons-involved.csv"))

# Get outcomes from file
outcomes <-fread("/tsd/p805/data/durable/users/arnovh/Adoptees/Data/MCOT_socioeco.csv") %>% select(-V1, -foedselsaar)
names(outcomes)

individual_outcomes <- outcomes %>% 
  gather(key,value, -w19_0634_lnr) %>% 
  filter(!is.na(value)) %>% 
  mutate(age = str_sub(key,start=-2)) %>% 
  mutate(outcome = str_sub(key,end=-3)) %>% 
  select(-key) %>% 
  group_by(w19_0634_lnr, outcome) %>% 
  summarize(value = mean(value)) %>% 
  spread(outcome, value)

### Which cohorts?
twin_cohort_range <- 1940:1960

### Build basic MCoT dataset, with w19lnrs and indicators of relationships
our_twins <- twins %>% left_join(faste) %>% 
  mutate(cohort = as.integer(foedsels_aar_mnd/100)) %>% 
  filter(cohort %in% twin_cohort_range)

all_children <- our_twins %>% select(w19_0634_lnr) %>% 
  left_join(fertility) %>% gather(key,value,-w19_0634_lnr, -numkids) %>% 
  mutate(childno = str_sub(key,start=-2)) %>% 
  mutate(info = str_sub(key,end=-3)) %>% 
  select(-key) %>% 
  group_by(w19_0634_lnr, childno) %>% 
  spread(info, value) %>% 
  filter(!is.na(child_bdate))

# Keep only maximum 2 children per twin, with the same otherparent, lowest birth orders
# where children have known ID and known parents
our_children <- all_children %>% 
  filter(child_otherparent!="") %>% 
  filter(child_lnr!="") %>% 
  group_by(w19_0634_lnr, child_otherparent) %>% 
  mutate(child_byear = as.numeric(str_sub(child_bdate,1,4))) %>% 
  mutate(target_diff = abs(1975-mean(child_byear))) %>% 
  arrange(w19_0634_lnr, child_otherparent, target_diff) %>% 
  slice(1:2) %>% 
  select(-child_bdate, -childno) %>% 
  arrange(child_byear) %>% 
  mutate(total_kids = n()) %>% 
  mutate(childno = row_number()) %>%
  gather(info, value, child_lnr, child_byear, child_sex) %>% 
  mutate(info = str_c(info, childno)) %>% 
  select(-childno) %>% 
  spread(info, value) %>% 
  ungroup %>% 
  group_by(w19_0634_lnr) %>% 
  arrange(w19_0634_lnr, target_diff)  %>% 
  slice(1) %>% 
  rename(partner_lnr = child_otherparent)

twinpair_vars <- select(our_twins, twinpair, Zygo, Group, cob = fodeland, cohort) %>% distinct %>% 
  drop_na

twins_with_families <- our_twins %>%
  left_join(our_children) %>% 
  select(twinpair, twinid, twinno, twin = w19_0634_lnr, partner = partner_lnr, child1 = child_lnr1, child2 = child_lnr2)

common_vars <- select(twins_with_families, twin, partner, child1, child2) %>% 
  gather(role, w19_0634_lnr) %>% 
  select(-role) %>% 
  distinct %>% 
  left_join(select(faste, w19_0634_lnr, sex = kjoenn, cohort = foedsels_aar_mnd)) %>% 
  mutate(cohort = str_sub(cohort,1,4)) %>% 
  left_join(individual_outcomes)

complete_families <- twins_with_families %>% 
  left_join(prefixit(common_vars,"twin_"), by=c(twin = "twin_w19_0634_lnr")) %>% 
  left_join(prefixit(common_vars,"partner_"), by=c(partner = "partner_w19_0634_lnr")) %>% 
  left_join(prefixit(common_vars,"child1_"), by=c(child1 = "child1_w19_0634_lnr")) %>% 
  left_join(prefixit(common_vars,"child2_"), by=c(child2 = "child2_w19_0634_lnr"))

twinpairs_with_families <- complete_families %>% gather(key,value, -twinpair, -twinid, -twinno) %>% 
  mutate(family = str_c("family", twinno,"_",key)) %>% 
  select(-twinno, -key, -twinid) %>% 
  spread(family,value) %>% 
  left_join(twinpair_vars) %>% 
  filter(!is.na(Zygo)) %>% 
  filter(Zygo %in% c(1,2)) %>% 
  mutate(relation = case_when(Group==1 ~ "m1m2",
                              Group==3 ~ "f1f2",
                              Group==2 ~ "m1m2",
                              Group==4 ~ "f1f2"))

theNames <- names(twinpairs_with_families)
ffNames <- theNames %>% 
  gsub(pattern = "_twin",replacement = "_mother") %>% 
  gsub(pattern = "_partner",replacement = "_father")

mmNames <- theNames %>% 
  gsub(pattern = "_twin",replacement = "_father") %>% 
  gsub(pattern = "_partner",replacement = "_mother")

dfs <- arrange(twinpairs_with_families, relation) %>% 
  group_split(relation) %>% 
  setNames(c("ff","mm"))

pedigrees_names <- bind_rows(setNames(dfs$ff, ffNames),
  setNames(dfs$mm, mmNames))

pedigrees <- pedigrees_names %>% 
  select(Zygo, relation, cohort,
         ends_with("_sex"),
         ends_with("_cohort"),
         ends_with("_education"),
         ends_with("_income_rank"),
         ends_with("_wealth_rank"),
         ends_with("_trei"),
  ) %>% 
  mutate(across(-relation, as.numeric)) %>% 
  mutate(across(c(family1_child1_education,family1_child2_education,
                  family2_child1_education,family2_child2_education), ~as.vector(scale.default(.x)))) %>% 
  mutate(across(c(family1_child1_income_rank,family1_child2_income_rank,
                  family2_child1_income_rank,family2_child2_income_rank), ~as.vector(scale.default(.x)))) %>% 
  mutate(across(c(family1_child1_trei,family1_child2_trei,
                  family2_child1_trei,family2_child2_trei), ~as.vector(scale.default(.x)))) %>% 
  mutate(across(c(family1_child1_wealth_rank,family1_child2_wealth_rank,
                  family2_child1_wealth_rank,family2_child2_wealth_rank), ~as.vector(scale.default(.x)))) %>% 
  mutate(across(c(family1_father_education,family1_mother_education,
                  family2_father_education,family2_mother_education), ~as.vector(scale.default(.x)))) %>% 
  mutate(across(c(family1_father_income_rank,family1_mother_income_rank,
                  family2_father_income_rank,family2_mother_income_rank), ~as.vector(scale.default(.x)))) %>% 
  mutate(across(c(family1_father_trei,family1_mother_trei,
                  family2_father_trei,family2_mother_trei), ~as.vector(scale.default(.x)))) %>% 
  mutate(across(c(family1_father_wealth_rank,family1_mother_wealth_rank,
                  family2_father_wealth_rank,family2_mother_wealth_rank), ~as.vector(scale.default(.x))))
  

fwrite(pedigrees, here("./pedigrees.csv"))

### Unstandardized data for descriptive tables
pedigrees_unstd <- pedigrees_names %>% 
  select(Zygo, relation, cohort,
         ends_with("_sex"),
         ends_with("_cohort"),
         ends_with("_education"),
         ends_with("_income_rank"),
         ends_with("_wealth_rank"),
         ends_with("_trei"),
  ) %>% 
  mutate(across(-relation, as.numeric))

fwrite(pedigrees_unstd, here("./pedigrees-unstd.csv"))

