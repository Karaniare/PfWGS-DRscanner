#######################################################################################
################ COHEN's EFFECT SIZE CALCULATION ######################################
#######################################################################################

library(lme4)        # for mixed-effects models
library(lmerTest)    # adds p-values to lmer()
library(emmeans)     # for estimated marginal means
library(performance) # for R², ICC, diagnostics
library(effectsize)  # for effect sizes (Cohen's d, etc.)
library(tidyverse)   # for data wrangling/plotting

## Running Mixed-effects regression model #

exvivo<-read.table("ex_vivo_phenotype_data_long_table_current.csv", header = T, sep = "\t")
ps<-read.table("parasitemia.csv", header = T, sep = "\t")
coi<-read.table("COI_calls_all_uganda_samples.csv", sep = "\t", header = T)
coi<-coi %>% select(sample_name, coi_median)

exvivo1<-exvivo %>%
  filter(wsaf>0.9,assay%in%c("Lumefantrine","Mefloquine","DHA","RSA"), px1_haplo%in%c("PIN","LMD"),
         k13_genotype%in%c("WT","A675V","C469Y","C469F","R561H","T348I","A675V_T348I"))#,
exvivo2<-merge(x=exvivo1, y=ps, by="sample_name")
exvivo2<-merge(x=exvivo2, y=coi, by="sample_name", all = T)
exvivo2<- exvivo2 %>% filter(assay!="NA")


lum<- exvivo2 %>% 
  # filter(assay=="Lumefantrine") %>%
  mutate(subject=px1_haplo)

lum$subject[lum$subject=="PIN"]<-1
lum$subject[lum$subject=="LMD"]<-0
lum$subject<-as.numeric(lum$subject)

mem_total<-data.frame()

for (i in unique(lum$assay)){
  dat1<- lum %>%
    filter(assay==i)
  ## fitting mixed-effect model
  fit<- lm(phenotype_score ~ px1_haplo + site + parasitemia +coi_median+ year+k13_genotype,
           data = dat1)
  
  ### Extracting model parameters
  emm <- emmeans(fit, ~ px1_haplo)
  dat_emm <- as.data.frame(emm)
  dat_pval<-as.data.frame(pairs(emm)) 
  
  ## calculating Cohen's effect size
  r2<-summary(fit)$r.squared
  r2_adj<-summary(fit)$adj.r.squared
  f2<-r2 / (1 - r2)
  
  
  ### combining results
  mem<-data.frame(assay=i,dat_emm,pvalue=dat_pval$p.value,effect_size=f2)
  mem_total<-rbind(mem,mem_total)
  
}


write_tsv(x=mem_total,file="results_mixed_effect_model_cohen.csv")

##############################################################################################
######################## MARGINAL AND CONDITIONAL R2 CALCULATIONS ###########################
#############################################################################################

library(lme4)
library(performance)

library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(effectsize)
library(tidyverse)

exvivo <- read.table("ex_vivo_phenotype_data_long_table_current.csv", header = TRUE, sep = "\t")
ps <- read.table("parasitemia.csv", header = TRUE, sep = "\t")
exvivo<-exvivo %>% filter(px1_haplo%in%c("PIN","LMD"))
exvivo$px1_haplo[exvivo$px1_haplo!="PIN"]<-"Others"

exvivo1 <- exvivo %>%
  filter(
    wsaf > 0.9,
    assay %in% c("Lumefantrine", "Mefloquine", "DHA", "RSA"),
    px1_haplo %in% c("PIN", "Others"),
    ! k13_genotype %in% c( "A675V", "C469Y", "C469F", "R561H", "A675V_T348I")
  )

exvivo2 <- merge(x = exvivo1, y = ps, by = "sample_name")
exvivo2 <- merge(x = exvivo2, y = coi, by = "sample_name", all = TRUE)
exvivo2 <- exvivo2 %>% 
  filter(assay != "NA")

lum <- exvivo2 %>%
  mutate(
    px1_haplo = factor(px1_haplo, levels = c("Others", "PIN")),
    site = factor(site),
    year = factor(year),
    k13_genotype = factor(k13_genotype)
  )


mem_total <- data.frame()

for (i in unique(lum$assay)) {
  dat1 <- lum %>% filter(assay == i)
  
  ################ FULL MODEL #############################
  
  # Fixed: px1_haplo, parasitemia, coi_median
  # Random: site, year, k13_genotype
  fit_full <- lmer(phenotype_score ~ px1_haplo + parasitemia+ coi_median +#k13_genotype+
                   (1 | site) + (1 | year), #+ (1 | k13_genotype), 
                   data = dat1)
  
  ################ REDUCED MODEL (Remove only PIN) #######################
  # Keeps other fixed and all random effects
  fit_red <- lmer(phenotype_score ~ parasitemia + coi_median +# k13_genotype+
                  (1 | site) + (1 | year),# + (1 | k13_genotype), 
                  data = dat1)
  
  ################# R2 CALCULATION #########################################
  r2_full <- r2_nakagawa(fit_full)
  r2_red  <- r2_nakagawa(fit_red)
  
  # Marginal R2: Variance from PIN + Parasitemia + COI
  # Conditional R2: Variance from PIN + Parasitemia + COI + Random Effects
  # Incremental R2: The unique contribution of PIN over the other fixed/random effects
  m_r2 <- r2_full$R2_marginal
  c_r2 <- r2_full$R2_conditional
  inc_r2 <- m_r2 - r2_red$R2_marginal
  
  #####################  EXTRACT RESULTS ####################
  emm_obj <- emmeans(fit_full, ~ px1_haplo)
  emm_df <- as.data.frame(emm_obj)
  pval <- as.data.frame(pairs(emm_obj))$p.value
  
  mem <- data.frame(assay = i, 
                    emm_df, 
                    pvalue = pval,
                    Marginal_R2 = m_r2, 
                    Conditional_R2 = c_r2,
                    PIN_Incremental_R2 = inc_r2)
  
  mem_total <- rbind(mem_total, mem)
}

write_tsv(mem_total, "results_PIN_and_Continuous_Fixed.csv")


##################################################################################
############## INTERACTION TERM BETWE PX1 and K13 GENOTYPES ######################
##################################################################################


################
library(lme4)
library(lmerTest) # Ensures p-values are calculated for lmer objects
library(emmeans)
library(performance)
library(effectsize)
library(tidyverse)

############ Data Loading and Cleaning ###################

exvivo <- read.table("ex_vivo_phenotype_data_long_table_current.csv", header = TRUE, sep = "\t")
ps <- read.table("parasitemia.csv", header = TRUE, sep = "\t")
coi <- read.table("COI_calls_all_uganda_samples.csv", sep = "\t", header = TRUE)

exvivo <- exvivo %>% filter(px1_haplo %in% c("PIN", "LMD"))
exvivo$px1_haplo[exvivo$px1_haplo != "PIN"] <- "LMD"
coi <- coi %>% select(sample_name, coi_median)

exvivo1 <- exvivo %>%
  filter(
    wsaf > 0.9,
    assay %in% c("Lumefantrine", "Mefloquine", "DHA", "RSA"),
    px1_haplo %in% c("PIN", "LMD"),
    k13_genotype %in% c("WT","A675V", "C469Y", "C469F", "R561H", "A675V_T348I")
  )

exvivo2 <- merge(x = exvivo1, y = ps, by = "sample_name")
exvivo2 <- merge(x = exvivo2, y = coi, by = "sample_name", all = TRUE)
exvivo2 <- exvivo2 %>% filter(assay != "NA")

lum <- exvivo2 %>%
  mutate(
    px1_haplo = factor(px1_haplo, levels = c("LMD", "PIN")),
    site = factor(site),
    year = factor(year),
    k13_genotype = factor(k13_genotype)
  )

lum$k13_genotype<-as.character(lum$k13_genotype)
lum$k13_genotype[lum$k13_genotype=="A675V"]<-"Mutant"
lum$k13_genotype[lum$k13_genotype=="C469Y"]<-"Mutant"


############### Loop for Interaction and Sample Sizes #################################

mem_total <- data.frame()

for (i in unique(lum$assay)) {
  dat1 <- lum %>% filter(assay == i)
  
  # Calculate sample sizes for this specific assay
  sample_counts <- dat1 %>%
    group_by(px1_haplo, k13_genotype) %>%
    summarise(n_group = n(), .groups = 'drop')
  
  n_total_assay <- nrow(dat1)
  
  #####  FULL MODEL WITH INTERACTION ############
  
  fit_full <- lmer(phenotype_score ~ px1_haplo * k13_genotype + parasitemia + coi_median + 
                     (1 | site) + (1 | year), 
                   data = dat1)
  
  ######## REDUCED MODEL (No Interaction) ###############
  
  fit_no_inter <- lmer(phenotype_score ~ px1_haplo + k13_genotype + parasitemia + coi_median + 
                         (1 | site) + (1 | year), 
                       data = dat1)
  
  ############ EXTRACT INTERACTION SIGNIFICANCE #############
  
  anova_table <- as.data.frame(anova(fit_full, type = "3"))
  p_interaction <- anova_table["px1_haplo:k13_genotype", "Pr(>F)"]
  
  #### R2 CALCULATION ##############
  
  r2_full <- r2_nakagawa(fit_full)
  r2_red  <- r2_nakagawa(fit_no_inter)
  
  m_r2 <- r2_full$R2_marginal
  c_r2 <- r2_full$R2_conditional
  inter_inc_r2 <- m_r2 - r2_red$R2_marginal 
  
  ########## POST-HOC TESTING ###################
  
  emm_obj <- emmeans(fit_full, ~ px1_haplo | k13_genotype)
  emm_df <- as.data.frame(emm_obj)
  
  ########### MERGE SAMPLE SIZES INTO RESULTS ###############
  
  # Joins the n_group count to the emmeans dataframe
  emm_df <- left_join(emm_df, sample_counts, by = c("px1_haplo", "k13_genotype"))
  
  ############ BUILD FINAL RESULTS ROWS ################
  mem <- emm_df %>%
    mutate(assay = i, 
           n_total_assay = n_total_assay,
           p_interaction = p_interaction,
           Marginal_R2 = m_r2, 
           Conditional_R2 = c_r2,
           Interaction_Incremental_R2 = inter_inc_r2)
  
  mem_total <- rbind(mem_total, mem)
}

write_tsv(mem_total, "results_Interaction_with_SampleSizes.csv")


#######################################

mem_total <- data.frame()

for (i in unique(lum$assay)) {
  dat1 <- lum %>% filter(assay == i)
  
  # 1. CALCULATE SAMPLE SIZE PER GROUP
  # This creates a small table with px1_haplo and the count (n)
  sample_sizes <- dat1 %>%
    group_by(px1_haplo) %>%
    summarise(n = n(), .groups = "drop")
  
  ## fitting model
  fit <- lm(phenotype_score ~ px1_haplo + site + parasitemia + coi_median + year + k13_genotype,
            data = dat1)
  
  ### Extracting model parameters
  emm <- emmeans(fit, ~ px1_haplo)
  dat_emm <- as.data.frame(emm)
  dat_pval <- as.data.frame(pairs(emm)) 
  
  # 2. MERGE SAMPLE SIZE INTO EMM RESULTS
  # Matches the 'n' to the correct 'px1_haplo' row
  dat_emm <- left_join(dat_emm, sample_sizes, by = "px1_haplo")
  
  ## calculating Cohen's effect size
  r2 <- summary(fit)$r.squared
  f2 <- r2 / (1 - r2)
  
  ### combining results
  # Added 'n' to the final data frame
  mem <- data.frame(assay = i, 
                    dat_emm, 
                    pvalue = dat_pval$p.value, 
                    effect_size = f2)
  
  mem_total <- rbind(mem, mem_total)
}
