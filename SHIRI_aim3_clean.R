#_____________________________________________________
#SHIRI: Influenza B Neuraminidase Correlates of Protection Analysis
#R file name: SHIRI_aim3_clean.R
#Author: Skyler Rogers
#Email: skyroger@umich.edu

##Last updated: 05/21/2026
#_____________________________________________________

#Notes
#For manuscript: 'Protective effects of influenza B neuraminidase antibodies against symptomatic influenza virus infection'
#Code was run using RStudio version 4.5.2

#manuscript colors:
#B/Phuket: #f28e69
#B/Brisbane: "#1E76C3"

# Packages ----------------------------------------------------------------
pacman::p_load(
  readr,
  haven,
  dplyr,
  tidyr,
  writexl,
  gtsummary,
  gt,
  ggplot2,
  patchwork, 
  tidyverse,
  survival,
  lubridate,
  Hmisc) 

# Functions ---------------------------------------------------------------
#calculate geometric mean titer (GMT) values
geometric_summary_log2 <- function(titer_val) {
  data_log2 <- titer_val[!is.na(titer_val)]
  gm_mean <- mean(data_log2)
  se_log2 <- sd(data_log2) / sqrt(length(data_log2))
  ci_lower <- gm_mean - 1.96 * se_log2
  ci_upper <- gm_mean + 1.96 * se_log2
  sample_size <- length(data_log2)
  return(c(geom_mean =gm_mean, geom_mean_log2 = 2^(gm_mean), se_log2 = se_log2, ci_lower_log2 = 2^(ci_lower), ci_upper_log2 = 2^(ci_upper), N = sample_size))
} 

#calculate geometric mean fold ratio of a titer value between two groups
geometric_mean_fold_ratio <- function(vax_gmt, unvax_gmt) {
  gm_fold_rise <- vax_gmt[["geom_mean_log2"]]/unvax_gmt[["geom_mean_log2"]]
  CI_lower <- 2^((vax_gmt[["geom_mean"]]-unvax_gmt[["geom_mean"]]) - 1.96*sqrt((vax_gmt[["se_log2"]])^2+(unvax_gmt[["se_log2"]])^2))
  CI_upper <- 2^((vax_gmt[["geom_mean"]]-unvax_gmt[["geom_mean"]]) + 1.96*sqrt((vax_gmt[["se_log2"]])^2+(unvax_gmt[["se_log2"]])^2))
  return(c(gm_fold_rise =gm_fold_rise, CI_lower = CI_lower, CI_upper = CI_upper))
} 

#create indicator variables to identify 2- and 4-fold titer increases and decreases from vaccination
fold_change_ind <- function(df, s1_val, s2_val, dif_val, twof_name, fourf_name, twof_name_d, fourf_name_d) {
  df %>% mutate({{twof_name}} := case_when(
    {{dif_val}}>=1 & !is.na({{s1_val}}) ~ "1", 
    is.na({{s1_val}}) | is.na({{s2_val}}) ~ "NA",
    TRUE ~ '0'),
    {{fourf_name}} := case_when(
      {{dif_val}}>=2 & !is.na({{s1_val}}) ~ "1", 
      is.na({{s1_val}}) | is.na({{s2_val}}) ~ "NA",
      TRUE ~ '0'),
    {{twof_name_d}} := case_when(
      {{dif_val}}<=-1 & !is.na({{s1_val}}) ~ "1", 
      is.na({{s1_val}}) | is.na({{s2_val}}) ~ "NA",
      TRUE ~ '0'),
    {{fourf_name_d}} := case_when(
      {{dif_val}}<=-2 & !is.na({{s1_val}}) ~ "1", 
      is.na({{s1_val}}) | is.na({{s2_val}}) ~ "NA",
      TRUE ~ '0'))
} 

#calculate geometric mean fold ratio (GMFR) values from pre- and post-vaccination titers
gmfr <- function(s1_val, s2_val) {
  d <- s2_val[!is.na(s2_val)] - s1_val[!is.na(s1_val)]
  mean_d <- mean(d)
  n <- length(d)
  se_d  <- sd(d) / sqrt(n)
  ci_lower <- mean_d - 1.96 * se_d
  ci_upper <- mean_d + 1.96 * se_d
  GMFR     <- 2^(mean_d) 
  CI_lower <- 2^(ci_lower) 
  CI_upper <- 2^(ci_upper) 
  return(c(GMFR= GMFR, CI_lower = CI_lower, CI_upper= CI_upper))
} 

# Data cleaning / setup -----------------------------------------------------------

#set working dictionary
secure_data <- Sys.getenv("SECURE_DATA_PATH")
setwd(secure_data)

#Import data
shiri_aim3_cleaned <- read_sas("Datasets/shiri_aim3_final.sas7bdat")

#check variable names
ls(shiri_aim3_cleaned) 

#select relevant variables
shiri_aim3_cleaned <- dplyr::select(shiri_aim3_cleaned, STUDY_ID, COLLECTION_DATE, 
                                    COLLECTION_POINT, ONSET_DATE, TYPE_FINAL, 
                                    SUBTYPE_FINAL, CASE_CONTROL, SPECIMEN_ID, 
                                    ALIQUOT_ID, AGE_ENROLL, SEX, STUDY_YEAR, 
                                    VACC_DATE_2017, VACC_STAT_2013, VACC_STAT_2014, 
                                    VACC_STAT_2015, VACC_STAT_2016, VACC_STAT_2017, 
                                    BRIS_BVIC_NAI, PHU_BYAM_NAI, VICTORIA_NAI,
                                    YAMAGATA_NAI, HAI_AVAILABLE, BRIS_BVIC_HAI, 
                                    PHU_BYAM_HAI, BRIS_BVIC_LOG_HAI, 
                                    PHU_BYAM_LOG_HAI)


#check for duplicate observations
duplicated(shiri_aim3_cleaned)

#titer variable data cleaning
#NAI
#create consistent naming structure
shiri_aim3_cleaned <- shiri_aim3_cleaned %>% rename(PHU_BYAM_LOG_NAI = PHU_BYAM_NAI)
shiri_aim3_cleaned <- shiri_aim3_cleaned %>% rename(BRIS_BVIC_LOG_NAI = BRIS_BVIC_NAI)
#create log base-2 transformed titer variables
shiri_aim3_cleaned$PHU_BYAM_NAI <- 5*2^(shiri_aim3_cleaned$PHU_BYAM_LOG_NAI)
shiri_aim3_cleaned$BRIS_BVIC_NAI <- 5*2^(shiri_aim3_cleaned$BRIS_BVIC_LOG_NAI)
shiri_aim3_cleaned$PHU_BYAM_NAI_LOG2 <- log2(shiri_aim3_cleaned$PHU_BYAM_NAI)
shiri_aim3_cleaned$BRIS_BVIC_NAI_LOG2 <- log2(shiri_aim3_cleaned$BRIS_BVIC_NAI)

#HAI
#create log base-2 transformed titer variables
shiri_aim3_cleaned$PHU_BYAM_HAI_LOG2 <- log2(shiri_aim3_cleaned$PHU_BYAM_HAI)
shiri_aim3_cleaned$BRIS_BVIC_HAI_LOG2 <- log2(shiri_aim3_cleaned$BRIS_BVIC_HAI)

#log2 transformation for figure axis
log_breaks <- log2(c(5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120))
log_labels <- c("5", "10", "20", "40", "80", "160", "320", "640", "1280", "2560", "5120")

#Long to wide dataset
shiri_aim3_wide <- shiri_aim3_cleaned %>%
  dplyr::select(STUDY_ID, COLLECTION_DATE, 
                COLLECTION_POINT, ONSET_DATE, TYPE_FINAL, 
                SUBTYPE_FINAL, CASE_CONTROL, SPECIMEN_ID, 
                ALIQUOT_ID, AGE_ENROLL, SEX, STUDY_YEAR, 
                VACC_DATE_2017, VACC_STAT_2013, VACC_STAT_2014, 
                VACC_STAT_2015, VACC_STAT_2016, VACC_STAT_2017, 
                BRIS_BVIC_LOG_NAI:BRIS_BVIC_HAI_LOG2) %>%
  pivot_wider(
    names_from = COLLECTION_POINT,
    values_from = c(COLLECTION_DATE, SPECIMEN_ID, ALIQUOT_ID, BRIS_BVIC_LOG_NAI:BRIS_BVIC_HAI_LOG2),
    id_cols = c(STUDY_ID, STUDY_YEAR, ONSET_DATE, VACC_DATE_2017, TYPE_FINAL, SUBTYPE_FINAL, 
                VACC_STAT_2013:VACC_STAT_2017, SEX, AGE_ENROLL, CASE_CONTROL)) #262 observations, 51 variables

# filter study ID's of participants
# shiri_aim3_studyids <- shiri_aim3_wide %>%
#   dplyr::select(STUDY_ID)

#infection status variable
shiri_aim3_wide <- mutate(shiri_aim3_wide, infected = ifelse(TYPE_FINAL=="FluB", 1, 0))

#subset datasets by vaccination status
shiri_aim3_wide_vax <- shiri_aim3_wide %>%
  filter(VACC_STAT_2017 !=0)

# Descriptive statistics --------------------------------------------------
table(shiri_aim3_wide$TYPE_FINAL) #64 IBV cases, 113 IBV negative, 85 no symptomatic illness reported
table(shiri_aim3_wide$CASE_CONTROL) #64 cases, 198 controls
table(shiri_aim3_wide$VACC_STAT_2017) #187 influenza vaccinated HCP, 75 unvaccinated against influenza

#Descriptive stats by cases:
shiri_aim3_wide_case <- shiri_aim3_wide %>%
  filter(CASE_CONTROL== "Case") 
as.numeric(median(shiri_aim3_wide_case$AGE_ENROLL, na.rm = TRUE)) #median age= 43
min(shiri_aim3_wide_case$AGE_ENROLL) #min age: 26
max(shiri_aim3_wide_case$AGE_ENROLL) #max age: 67
table(shiri_aim3_wide_case$SEX) #51 female, 13 male
prop.table(table(shiri_aim3_wide_case$SEX)) #79.7% female, 20.3% male
table(shiri_aim3_wide_case$VACC_STAT_2017) #37 unvax, 27 vax
prop.table(table(shiri_aim3_wide_case$VACC_STAT_2017)) #57.8% unvax, 42.2% vax

#Descriptive stats by controls:
shiri_aim3_wide_controls <- shiri_aim3_wide %>%
  filter(CASE_CONTROL != "Case")
as.numeric(median(shiri_aim3_wide_controls$AGE_ENROLL, na.rm = TRUE)) #median age= 43
min(shiri_aim3_wide_controls$AGE_ENROLL) #min age: 23
max(shiri_aim3_wide_controls$AGE_ENROLL) #max age: 70
table(shiri_aim3_wide_controls$SEX) #139 female, 59 male
prop.table(table(shiri_aim3_wide_controls$SEX)) #70.2% female, 29.8% male
table(shiri_aim3_wide_controls$VACC_STAT_2017) #38 unvax, 160 vax
prop.table(table(shiri_aim3_wide_controls$VACC_STAT_2017)) #19.2% unvax, 80.8% vax

#prior vax history
#Previous season vaccinations
table(shiri_aim3_wide$VACC_STAT_2016) #186 vaccinated in 2016/17 season
table(shiri_aim3_wide$VACC_STAT_2016, shiri_aim3_wide$VACC_STAT_2017) #163 HCP vaccinated in 2016/17 and 2017/18 seasons
sum(shiri_aim3_wide_vax$VACC_STAT_2016==1 &shiri_aim3_wide_vax$VACC_STAT_2017==1) /sum(shiri_aim3_wide$VACC_STAT_2016==1) #87.6% of HCP vaccinated in 2016-17 season also vaccinated in 2017/18 season

#how many participants that were vaccinated in the 2017-18 season have recorded prior vaccination in prior 4 years
prior_vax_count <-as.numeric(sum(apply(shiri_aim3_wide_vax[, c("VACC_STAT_2013", "VACC_STAT_2014", "VACC_STAT_2015", "VACC_STAT_2016")], 1, function(row) any(row == 1, na.rm = TRUE)))) # 183 vaccinated participants received prior vaccination
prior_vax_count /sum(shiri_aim3_wide_vax$VACC_STAT_2017==1) #97.9% of participants vaccinated in 2017-18 season had recieved at least one additional vaccine in the prior 4 years

#frequent vax participants (recieved at least 3/4 prior season vaccines)
shiri_aim3_wide_vax$num_vax <- rowSums(shiri_aim3_wide_vax[, c("VACC_STAT_2013", "VACC_STAT_2014", "VACC_STAT_2015", "VACC_STAT_2016")], na.rm = TRUE)
shiri_aim3_wide_vax <- shiri_aim3_wide_vax %>%
  mutate(freq_vax = ifelse(num_vax>=3, 1,0))
table(shiri_aim3_wide_vax$freq_vax) #108 vaccinated participants received 3/4 prior season vaccines
sum(shiri_aim3_wide_vax$freq_vax)/sum(shiri_aim3_wide_vax$VACC_STAT_2017) #57.8% of HCP vaccinated in 2017-18 season received influenza vaccine in 3/4 of previous seasons


#check the time between vaccination and S2 collection in vaccinated HCP (should be collected ~30 days post-vax)
shiri_aim3_wide_vax <- shiri_aim3_wide_vax %>%
  mutate(vax_time = COLLECTION_DATE_S2 - VACC_DATE_2017)
as.numeric(median(shiri_aim3_wide_vax$vax_time, na.rm = TRUE)) #median= 28 days
quantile(shiri_aim3_wide_vax$vax_time, probs = c(0.25, 0.75), na.rm= TRUE) #IQR:(24-34) days

#did any vaccinated participants receive vaccination after infection
ifelse(is.na(shiri_aim3_wide$ONSET_DATE) | is.na(shiri_aim3_wide$VACC_DATE_2017), 0, ifelse(shiri_aim3_wide$VACC_DATE_2017 >= shiri_aim3_wide$ONSET_DATE, 1, 0)) #no

#Create descriptive statistics table (Table 1)
#re-order IBV infection status
shiri_aim3_wide <- shiri_aim3_wide %>%
  mutate(flub_type = ifelse(is.na(TYPE_FINAL) | TYPE_FINAL == "N/A", "No symptomatic illness", as.character(TYPE_FINAL))) %>%
  mutate(flub_type = ifelse(flub_type == "FluB", "Influenza B", as.character(flub_type))) %>%
  mutate(flub_type = ifelse(flub_type == "Negative", "Influenza B Negative", as.character(flub_type)))
shiri_aim3_wide$flub_type <- factor(shiri_aim3_wide$flub_type, levels = c("Influenza B", "Influenza B Negative", "No symptomatic illness"))
#vaccination status variable labels 
shiri_aim3_wide <-  shiri_aim3_wide %>%
  mutate(vacc_stat_2017_label = recode(VACC_STAT_2017, 
                                       "0" = "Unvaccinated",
                                       "1" = "Vaccinated"))
table1 <- shiri_aim3_wide %>%
  dplyr::select(SEX, AGE_ENROLL, vacc_stat_2017_label, CASE_CONTROL,flub_type) %>%
  tbl_summary(by = CASE_CONTROL, digits = list(SEX~c(0,1), AGE_ENROLL ~c(0,1,1), vacc_stat_2017_label~c(0,1), flub_type~c(0,1)),
              label = list(SEX ~ "Sex",
                           AGE_ENROLL ~ "Age",
                           vacc_stat_2017_label ~ "Influenza Vaccination Status",
                           flub_type ~ "Infection Status"))%>%
  bold_labels() %>%
  modify_footnote(everything() ~ NA)%>%
  modify_footnote_header(
    footnote = "Median (Q1, Q3); n(%) reported with column percentages.",
    columns = all_stat_cols(),
    replace = FALSE) %>%
  modify_footnote_body(
    footnote = "Healthcare personnel with “no symptomatic illness” did not report illness during study season.",
    columns = "label",
    rows = label  == "No symptomatic illness"
  )%>%
  modify_footnote_body(
    footnote = "Healthcare personnel with “Influenza B Negative” reported illness during study season but tested negative for influenza B infection.",
    columns = "label",
    rows = label  == "Influenza B Negative"
  )

table1 <- table1 %>%
  modify_table_body(
    function(table1) {
      table1 <- table1 %>%
        dplyr::add_row(
          variable = "flub_type",  
          label = "– B/non-repeatable",
          row_type = "level",     
          stat_1 = "1 (1.6%)", stat_2 = "0 (0.0%)",
          .after = which(table1$label == "Influenza B")
        ) %>%
        dplyr::add_row(
          variable = "flub_type",  
          label = "– B/not able to lineage",
          row_type = "level",     
          stat_1 = "6 (9.4%)",stat_2 = "0 (0.0%)",
          .after = which(table1$label == "Influenza B")
        ) %>%
        dplyr::add_row(
          variable = "flub_type",
          label = "– B/Yamagata",
          row_type = "level",
          stat_1 = "57 (89.1%)", stat_2 = "0 (0.0%)", 
          .after = which(table1$label == "Influenza B")
        )
      table1 })

gt_table1 <- as_gt(table1) %>%
  cols_width(
    everything() ~ px(250))
gt_table1
gtsave(gt_table1, filename = file.path(secure_data, "Tables + Figures/table1_case_stat.docx"))

# Baseline Titer Analysis  ----------------------------------------------------------------
##Q1: Is S1 titer values different between influenza vaccinated and unvaccinated HCP?

#remove participants with missing S1 collection from analysis
shiri_aim3_wide_q1 <- shiri_aim3_wide %>%
  drop_na(BRIS_BVIC_HAI_LOG2_S1, PHU_BYAM_NAI_LOG2_S1, PHU_BYAM_HAI_LOG2_S1) #one participant missing only B/Brisbane S1 NAI titers, dont remove in analysis subset for inclusion of participant in B/Phuket S1 analysis
#260 participants in B/Phuket and 259 B/Brisbane analysis 

#compare baseline NAI by vaccination status
#B/Phuket NAI
t.test(shiri_aim3_wide_q1$PHU_BYAM_NAI_LOG2_S1 ~ shiri_aim3_wide_q1$VACC_STAT_2017) #p-value = 1.169e-05
#wilcox.test(shiri_aim3_wide_q1$PHU_BYAM_NAI_LOG2_S1 ~ shiri_aim3_wide_q1$VACC_STAT_2017) #p-value = 2.412e-05

phuket_nai_vax <- ggplot(shiri_aim3_wide_q1, aes(x=factor(VACC_STAT_2017), y=PHU_BYAM_NAI_LOG2_S1)) +
  geom_violin(fill = "#f28e69") +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.75, color = "black", linetype = "dashed", linewidth= 0.3) +
  geom_jitter(shape=16, alpha=0.25,position=position_jitter(0.2), size=1)+
  #  ggtitle("NAI B/Phuket/3073/2013") +
  xlab("Influenza Vaccination Status") +
  scale_x_discrete(labels = c("0" = "Unvaccinated", "1" = "Vaccinated")) +
  # ylab(expression(log[2]~'NAI Titer')) + #log2 transformed label y axis
  ylab(expression('NAI Titer')) + #non-transformed titer label y axis
  #scale_y_continuous(breaks = seq(0, 13, 1), limits = c(0, 13)) +   #log2 transformed scale y axis
  scale_y_continuous(breaks = log_breaks, labels = log_labels, limits=c(1.75,12.75)) + #non-transformed titer scale y axis
  theme_minimal() +
  theme(legend.position = "none")
phuket_nai_vax
#check for normality
hist(shiri_aim3_wide_q1$PHU_BYAM_NAI_LOG2_S1)
qqnorm(shiri_aim3_wide_q1$PHU_BYAM_NAI_LOG2_S1)
qqline(shiri_aim3_wide_q1$PHU_BYAM_NAI_LOG2_S1)

#B/Brisbane NAI
t.test(shiri_aim3_wide_q1$BRIS_BVIC_NAI_LOG2_S1 ~ shiri_aim3_wide_q1$VACC_STAT_2017) #p-value = 0.05243
#wilcox.test(shiri_aim3_wide_q1$BRIS_BVIC_NAI_LOG2_S1 ~ shiri_aim3_wide_q1$VACC_STAT_2017) #p-value = 0.03903

brisbane_nai_vax <-ggplot(shiri_aim3_wide_q1, aes(x=factor(VACC_STAT_2017), y=BRIS_BVIC_NAI_LOG2_S1)) +
  geom_violin(fill = "#1E76C3")+
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.75, color = "black", linetype = "dashed", linewidth= 0.3) +
  geom_jitter(shape=16, alpha=0.25,position=position_jitter(0.2), size=1)+
  # ggtitle("NAI B/Brisbane/60/2008") +
  xlab("Influenza Vaccination Status") +
  scale_x_discrete(labels = c("0" = "Unvaccinated", "1" = "Vaccinated")) +
  # ylab(expression(log[2]~'NAI Titer')) + #log2 transformed label y axis
  ylab(expression('NAI Titer')) + #non-transformed titer label y axis
  #scale_y_continuous(breaks = seq(0, 13, 1), limits = c(0, 13)) +   #log2 transformed scale y axis
  scale_y_continuous(breaks = log_breaks, labels = log_labels, limits=c(1.75,12.75)) + #non-transformed titer scale y axis
  theme_minimal() +
  theme(legend.position = "none")
brisbane_nai_vax
#check for normality
hist(shiri_aim3_wide_q1$BRIS_BVIC_NAI_LOG2_S1)
qqnorm(shiri_aim3_wide_q1$BRIS_BVIC_NAI_LOG2_S1)
qqline(shiri_aim3_wide_q1$BRIS_BVIC_NAI_LOG2_S1)

#compare baseline HAI by vaccination status
#B/Phuket HAI
t.test(shiri_aim3_wide_q1$PHU_BYAM_HAI_LOG2_S1~ shiri_aim3_wide_q1$VACC_STAT_2017) #p-value = 0.02051
#wilcox.test(shiri_aim3_wide_q1$PHU_BYAM_HAI_LOG2_S1 ~ shiri_aim3_wide_q1$VACC_STAT_2017) #p-value = 0.02236

#check for normality
hist(shiri_aim3_wide_q1$PHU_BYAM_HAI_LOG2_S1, breaks=9)
qqnorm(shiri_aim3_wide_q1$PHU_BYAM_HAI_LOG2_S1)
qqline(shiri_aim3_wide_q1$PHU_BYAM_HAI_LOG2_S1)

#B/Brisbane HAI
t.test(shiri_aim3_wide_q1$BRIS_BVIC_HAI_LOG2_S1 ~ shiri_aim3_wide_q1$VACC_STAT_2017) #p-value = 0.1916
#wilcox.test(shiri_aim3_wide_q1$BRIS_BVIC_HAI_LOG2_S1 ~ shiri_aim3_wide_q1$VACC_STAT_2017) #p-value = 0.08192

#check for normality 
hist(shiri_aim3_wide_q1$BRIS_BVIC_HAI_LOG2_S1)
qqnorm(shiri_aim3_wide_q1$BRIS_BVIC_HAI_LOG2_S1)
qqline(shiri_aim3_wide_q1$BRIS_BVIC_HAI_LOG2_S1)

#format figures using patchwork
nai_vax_boxplots <-phuket_nai_vax + brisbane_nai_vax
nai_vax_boxplots
ggsave(plot = nai_vax_boxplots, width = 10, height = 5, dpi = 300, filename = file.path(secure_data, "Tables + Figures/nai_vax_boxplots.png"))

#GMT at enrollment by vaccination status
#stratify Q1 datasets by 2017-18 vaccination status
shiri_aim3_wide_q1_unvax <- shiri_aim3_wide_q1 %>%
  filter(VACC_STAT_2017 !=1)
shiri_aim3_wide_q1_vax <- shiri_aim3_wide_q1 %>%
  filter(VACC_STAT_2017 !=0)

#GMT and geometric mean fold ratio calculations
##formula for unpaired samples = 2^((gm_mean_vax-gm_mean_unvax) +/- 1.96*sqrt((se_vax)^2 + (se_unvax)^2))

#NAI titers
#B/Phuket 
phuk_vax_nai_gmt <- geometric_summary_log2(shiri_aim3_wide_q1_vax$PHU_BYAM_NAI_LOG2_S1) #vax GMT (95% CI): 190.6 (157.7, 230.5) n=186
phuk_unvax_nai_gmt <- geometric_summary_log2(shiri_aim3_wide_q1_unvax$PHU_BYAM_NAI_LOG2_S1) #unvax GMT (95% CI): 77.1 (54.9,108.1) n=74
#B/Brisbane
bris_vax_nai_gmt <- geometric_summary_log2(shiri_aim3_wide_q1_vax$BRIS_BVIC_NAI_LOG2_S1) #vax GMT (95% CI): 367.6 (298.2, 453.1) n=185 
bris_unvax_nai_gmt <- geometric_summary_log2(shiri_aim3_wide_q1_unvax$BRIS_BVIC_NAI_LOG2_S1) #unvax GMT (95% CI): 250.8 (182.0, 345.6) n=74

# NAI geometric mean titer ratios 
#B/Phuket 
geometric_mean_fold_ratio(phuk_vax_nai_gmt, phuk_unvax_nai_gmt) #2.5 times higher NAI titers in vaccinated HCP, 95% CI (1.7, 3.6)
#B/Brisbane
geometric_mean_fold_ratio(bris_vax_nai_gmt, bris_unvax_nai_gmt) #1.5 times higher NAI titers in vaccinated HCP, 95% CI (1.0, 2.1)


#HAI titers
#B/Phuket 
phuk_vax_hai_gmt <- geometric_summary_log2(shiri_aim3_wide_q1_vax$PHU_BYAM_HAI_LOG2_S1) #vax GMT (95% CI): 40.4 (33.3, 49.1) n=186
phuk_unvax_hai_gmt <-geometric_summary_log2(shiri_aim3_wide_q1_unvax$PHU_BYAM_HAI_LOG2_S1) #unvax GMT (95% CI): 26.5 (19.7,35.6) n=74
#B/Brisbane
bris_vax_hai_gmt <- geometric_summary_log2(shiri_aim3_wide_q1_vax$BRIS_BVIC_HAI_LOG2_S1) #vax GMT (95% CI): 79.7 (69.0, 92.1) n=186 
bris_unvax_hai_gmt <-geometric_summary_log2(shiri_aim3_wide_q1_unvax$BRIS_BVIC_HAI_LOG2_S1) #unvax GMT (95% CI): 65.1 (49.9, 84.9) n=74

# HAI geometric mean titer ratio 
#B/Phuket 
geometric_mean_fold_ratio(phuk_vax_hai_gmt, phuk_unvax_hai_gmt) #1.5 times higher HAI titers in vaccinated HCP, 95%CI (1.1,2.2)
#B/Brisbane
geometric_mean_fold_ratio(bris_vax_hai_gmt, bris_unvax_hai_gmt) #1.2 times higher HAI titers in vaccinated HCP, 95%CI (0.9, 1.7)

#Baseline GMT table by vaccination status
#NAI
table_gmt_q1 <- data.frame(
  Variables = c("B/Phuket",
                "B/Brisbane"),
  Column_1 = c("77.1 (54.9, 108.1) n=74","250.8 (182.0, 345.6) n=74"),
  Column_2 = c("190.6 (157.7, 230.5) n=186","367.6 (298.2, 453.1) n=185"),
  Column_3 = c("2.5 (1.7, 3.6)", "1.5 (1.0, 2.1)"),
  Column_4 = c("<0.001***", "0.052"))
table_gmt_q1 <- table_gmt_q1 %>%
  gt() %>%
  cols_align(
    align = "left",   
    columns = everything())%>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())) %>%
  cols_label(
    Variables = "NAI Titer Measurements",
    Column_1 = "S1 GMT Influenza Unvaccinated",
    Column_2 = "S1 GMT Influenza Vaccinated",
    Column_3 = "GMT Ratio",
    Column_4 = "p-value") %>%
  tab_options(
    heading.align = "left")%>%
  tab_footnote(
    footnote = md("p-value for test of baseline GMT by ifluenza vaccination status for each antigen using independent t-tests."),
    locations = cells_column_labels(columns = Column_4))%>%
  tab_footnote(
    footnote = md("*p* < 0.05 ( * ), *p* < 0.01 ( ** ), *p* < 0.001 ( *** )"),
    locations = cells_column_labels(columns = Column_4))%>%
  tab_footnote(
    footnote = md("GMT (95% CI)"),
    locations = cells_column_labels(columns = c(Column_1,Column_2))) %>%
  tab_footnote(
    footnote = md("GMT Ratio (95% CI)"),
    locations = cells_column_labels(columns = Column_3))
table_gmt_q1
gtsave(table_gmt_q1, filename = file.path(secure_data, "Tables + Figures/q1_gmt_table.pdf"))

#HAI
table_gmt_q1_hai <- data.frame(
  Variables = c("B/Phuket",
                "B/Brisbane"),
  Column_1 = c("26.5 (19.7, 35.6) n=74","65.1 (49.9, 84.9) n=74"),
  Column_2 = c("40.4 (33.3, 49.1) n=186","79.7 (69.0, 92.1) n=186"),
  Column_3 = c("1.5 (1.1, 2.2)", "1.2 (0.9, 1.7)"),
  Column_4 = c("0.021*", "0.192"))
table_gmt_q1_hai <- table_gmt_q1_hai %>%
  gt() %>%
  cols_align(
    align = "left",   
    columns = everything())%>%
#  tab_style(
#    style = cell_text(weight = "bold"),
#    locations = cells_column_labels(columns = everything())) %>%
  cols_label(
    Variables = md("**HAI Titer Measurements**"),
    Column_1 = md("**S1 GMT Influenza Unvaccinated**"),
    Column_2 = md("**S1 GMT Influenza Vaccinated**"),
    Column_3 = md("**GMT Ratio**"),
    Column_4 = md("**p-value**")) %>%
  tab_options(
    heading.align = "left")%>%
  tab_footnote(
    footnote = md("p-value for test of baseline GMT by influenza vaccination status for each antigen using independent t-tests."),
    locations = cells_column_labels(columns = Column_4))%>%
  tab_footnote(
    footnote = md("p < 0.05 ( * ), p < 0.01 ( ** ), p < 0.001 ( *** )"),
    locations = cells_column_labels(columns = Column_4))%>%
  tab_footnote(
    footnote = md("GMT (95% CI)"),
    locations = cells_column_labels(columns = c(Column_1,Column_2))) %>%
  tab_footnote(
    footnote = md("GMT Ratio (95% CI)"),
    locations = cells_column_labels(columns = Column_3))
table_gmt_q1_hai
gtsave(table_gmt_q1_hai, filename = file.path(secure_data, "Tables + Figures/q1_gmt_table_hai.docx"))


# Association between NAI titers and IBV infection --------------------

## Time varying left-truncated cox proportional-hazards model with robust sandwich estimators ----------------------------------------------
#Data cleaning
#vaccinated participants have 2 rows in analysis (unvaccinated time vs vaccinated time)
#if vacc_stat_2017 = 1
# AND collection_point = S1, then time in study is COLLECTION_DATE to VACC_DATE_2017 (would then be considered unvaccinated time)
# AND collection_point= S2, then time in study is COLLECTION_DATE to ONSET if CASE or 4-1-2018 if CONTROL (would then be considered vaccinated time)
#if vacc_stat_2017 = 0, then time in study is from COLLECTION_DATE to ONSET if CASE or 4-1-2018 if CONTROL

study_end_date <- as.Date("2018-04-01")
shiri_aim3_cleaned$ONSET_DATE <- as.Date(shiri_aim3_cleaned$ONSET_DATE)
shiri_aim3_cleaned$VACC_DATE_2017 <- as.Date(shiri_aim3_cleaned$VACC_DATE_2017)
shiri_aim3_cleaned$COLLECTION_DATE <- as.Date(shiri_aim3_cleaned$COLLECTION_DATE)

shiri_aim3_cleaned <- shiri_aim3_cleaned %>%
  mutate(infected = ifelse(is.na(ONSET_DATE), 0,
                           ifelse(ONSET_DATE>1, 1, 0)))

#create infection and vaccination indicator variables for specific time intervals
shiri_aim3_cleaned <- shiri_aim3_cleaned %>%
  mutate(study_time = ifelse(VACC_STAT_2017==1 & COLLECTION_POINT=="S1", VACC_DATE_2017-COLLECTION_DATE, 
                             ifelse(CASE_CONTROL=="Case", ONSET_DATE-COLLECTION_DATE, study_end_date-COLLECTION_DATE))) 
#if vaccinated, the unvax study time is time from S1 to vax date. If infected (case), study time is from collection to onset, everyone else is collection to end of follow up 

shiri_aim3_cleaned <- shiri_aim3_cleaned %>%
  mutate(vacc_stat_17_long = ifelse(VACC_STAT_2017==0, 0,
                                    ifelse(VACC_STAT_2017 ==1 & COLLECTION_POINT=="S1",0,1)))
#0 for unvaccinated participants and S1 vax participants, 1 for vaccinated S2 time points

shiri_aim3_cleaned <- shiri_aim3_cleaned %>%
  mutate(infected_long = ifelse(is.na(ONSET_DATE), 0,
                                ifelse(VACC_STAT_2017 ==1 & COLLECTION_POINT=="S1"& ONSET_DATE>COLLECTION_DATE & ONSET_DATE<VACC_DATE_2017,1,
                                       ifelse(VACC_STAT_2017 ==1 & COLLECTION_POINT=="S2"& ONSET_DATE>COLLECTION_DATE & ONSET_DATE<study_end_date,1,
                                              ifelse(VACC_STAT_2017==0 & COLLECTION_POINT=="S1" & ONSET_DATE>COLLECTION_DATE & ONSET_DATE<study_end_date,1,0)))))

#time-varying survival time
study_start_date <- as.Date("2017-06-01") # before first collection 
study_end_date <- as.Date("2018-04-01") # after last reported onset 

shiri_aim3_cleaned <- shiri_aim3_cleaned %>%
  arrange(STUDY_ID, ONSET_DATE) %>%  
  mutate(
    # Calculate numeric time differences, define risk start as: serum collection date - study start date, define risk end as: infection date/ end of follow up - study start date
    time_start_numeric = as.numeric(COLLECTION_DATE-study_start_date),
    time_stop_numeric = ifelse(infected_long==1, as.numeric(ONSET_DATE-study_start_date), 
                               ifelse(COLLECTION_POINT=="S1" &VACC_STAT_2017==0, as.numeric(study_end_date-study_start_date), 
                                      ifelse(COLLECTION_POINT=="S1" & VACC_STAT_2017==1, as.numeric(VACC_DATE_2017-study_start_date), as.numeric(study_end_date-study_start_date)))))

# If original calculated values are NA, fill them in with logical defaults
shiri_aim3_cleaned<- shiri_aim3_cleaned %>%
  mutate(
    time_start_numeric = ifelse(is.na(time_start_numeric), 
                                as.numeric(COLLECTION_DATE - study_start_date), # if no serum collected, the risk start = risk end (avoid missing value)
                                time_start_numeric),
    time_stop_numeric = ifelse(is.na(time_stop_numeric), 
                               as.numeric(study_end_date - study_start_date), 
                               time_stop_numeric)) #no original calculated values are NA in dataset

#separate dataset by antigen for analysis
#remove participants previously excluded from baseline titer analysis
shiri_tvcoxph_phuk <- shiri_aim3_cleaned %>%
  drop_na(PHU_BYAM_NAI_LOG2, PHU_BYAM_HAI_LOG2)
shiri_tvcoxph_bris <- shiri_aim3_cleaned %>%
  drop_na(BRIS_BVIC_NAI_LOG2, BRIS_BVIC_HAI_LOG2)

#if participants are vaccinated and dont have both S1 and S2 collection data, remove them
shiri_tvcoxph_phuk <- shiri_tvcoxph_phuk %>%
  group_by(STUDY_ID) %>%
  filter(!(VACC_STAT_2017 == 1 & n() < 2)) %>%
  ungroup()
shiri_tvcoxph_bris <- shiri_tvcoxph_bris %>%
  group_by(STUDY_ID) %>%
  filter(!(VACC_STAT_2017 == 1 & n() < 2)) %>%
  ungroup()

# Select relevant columns
shiri_tvcoxph_phuk <- shiri_tvcoxph_phuk %>%
  filter(
    !is.na(time_start_numeric) & !is.na(time_stop_numeric) & 
      time_stop_numeric > time_start_numeric
  )  %>%
  dplyr::select(STUDY_ID,SPECIMEN_ID, COLLECTION_POINT,time_start_numeric, time_stop_numeric, BRIS_BVIC_NAI_LOG2, BRIS_BVIC_HAI_LOG2, PHU_BYAM_NAI_LOG2, PHU_BYAM_HAI_LOG2,AGE_ENROLL,
                VACC_STAT_2017, ONSET_DATE, COLLECTION_DATE, VACC_DATE_2017, infected, infected_long, study_time, vacc_stat_17_long) %>% 
  distinct()
#check which participants have missing data for analysis
anti_join(shiri_aim3_cleaned, shiri_tvcoxph_phuk, by = "STUDY_ID") #check excluded participants, 1 additional person removed in B/Phuket analysis, 3 total
anti_join(shiri_aim3_cleaned, shiri_tvcoxph_phuk) #unvaccinated time of 6 removed for S1 occuring after vaccination
length(unique(shiri_tvcoxph_phuk$STUDY_ID)) #259 participants in B/Phuket cox-ph analysis

shiri_tvcoxph_bris <- shiri_tvcoxph_bris %>%
  filter(
    !is.na(time_start_numeric) & !is.na(time_stop_numeric) & 
      time_stop_numeric > time_start_numeric
  )  %>%
  dplyr::select(STUDY_ID,SPECIMEN_ID, COLLECTION_POINT,time_start_numeric, time_stop_numeric, BRIS_BVIC_NAI_LOG2, BRIS_BVIC_HAI_LOG2,AGE_ENROLL,
                VACC_STAT_2017, ONSET_DATE, COLLECTION_DATE, VACC_DATE_2017, infected, infected_long, study_time, vacc_stat_17_long) %>% 
  distinct() 
#check which participants have missing data for analysis
anti_join(shiri_aim3_cleaned, shiri_tvcoxph_bris, by = "STUDY_ID") #check excluded participants, 1 additional person removed in B/Brisbane analysis, 4 total
check <-anti_join(shiri_aim3_cleaned, shiri_tvcoxph_bris) #unvaccinated time of 5 removed for S1 occurring after vaccination
length(unique(shiri_tvcoxph_bris$STUDY_ID)) #258 participants

#save dataset for use in percent hazard reduction code
write_csv(shiri_tvcoxph_phuk, file.path(secure_data, "Datasets/shiri_aim3_coxph_phuk.csv"))
write_csv(shiri_tvcoxph_bris, file.path(secure_data, "Datasets/shiri_aim3_coxph_bris.csv"))

#left truncated time-varying cox-ph model (with time_start and time_stop variables)
#crude model (NAI-only predictior)
coxphtv_inf_bris_left_crude_nai <- as.formula("Surv(time_start_numeric, time_stop_numeric,infected_long)~ BRIS_BVIC_NAI_LOG2 + cluster(STUDY_ID)")
coxphtv_inf_phuk_left_crude_nai <- as.formula("Surv(time_start_numeric, time_stop_numeric,infected_long)~ PHU_BYAM_NAI_LOG2 + cluster(STUDY_ID)")
fit.coxphtv_inf_bris_left_crude_nai <- coxph(coxphtv_inf_bris_left_crude_nai, data = shiri_tvcoxph_bris)
fit.coxphtv_inf_phuk_left_crude_nai <- coxph(coxphtv_inf_phuk_left_crude_nai, data = shiri_tvcoxph_phuk)
summary(fit.coxphtv_inf_bris_left_crude_nai)
summary(fit.coxphtv_inf_phuk_left_crude_nai)
#crude model (HAI-only predictior)
coxphtv_inf_bris_left_crude_hai <- as.formula("Surv(time_start_numeric, time_stop_numeric,infected_long)~  BRIS_BVIC_HAI_LOG2   + cluster(STUDY_ID)")
coxphtv_inf_phuk_left_crude_hai <- as.formula("Surv(time_start_numeric, time_stop_numeric,infected_long)~ PHU_BYAM_HAI_LOG2   + cluster(STUDY_ID)")
fit.coxphtv_inf_bris_left_crude_hai <- coxph(coxphtv_inf_bris_left_crude_hai, data = shiri_tvcoxph_bris)
fit.coxphtv_inf_phuk_left_crude_hai <- coxph(coxphtv_inf_phuk_left_crude_hai, data = shiri_tvcoxph_phuk)
summary(fit.coxphtv_inf_bris_left_crude_hai)
summary(fit.coxphtv_inf_phuk_left_crude_hai) 

#adjusted models (including NAI and HAI modeled by each antigen)
coxphtv_inf_bris_left_adj <- as.formula("Surv(time_start_numeric, time_stop_numeric,infected_long)~ BRIS_BVIC_NAI_LOG2 + BRIS_BVIC_HAI_LOG2   + cluster(STUDY_ID)")
coxphtv_inf_phuk_left_adj <- as.formula("Surv(time_start_numeric, time_stop_numeric,infected_long)~ PHU_BYAM_NAI_LOG2 +PHU_BYAM_HAI_LOG2   + cluster(STUDY_ID)")
fit.coxphtv_inf_bris_left_adj <- coxph(coxphtv_inf_bris_left_adj, data = shiri_tvcoxph_bris)
fit.coxphtv_inf_phuk_left_adj <- coxph(coxphtv_inf_phuk_left_adj, data = shiri_tvcoxph_phuk)
summary(fit.coxphtv_inf_bris_left_adj) 
summary(fit.coxphtv_inf_phuk_left_adj) 

## Sensitivity analysis -------------------------------------
# risk start on specificity date, assume no waning before
# data cleaning 
shiri_tvcoxph_phuk$start_circulate <- as.Date("2017-12-01") #define a finer wave within the study period, start date = before the first IBV case detected in the study period
shiri_tvcoxph_bris$start_circulate <- as.Date("2017-12-01") #define a finer wave within the study period, start date = before the first IBV case detected in the study period

shiri_tvcoxph_sens_phuk <- shiri_tvcoxph_phuk %>%
  arrange(SPECIMEN_ID, COLLECTION_DATE) %>%  # Ensure the data is sorted by ID and date
  group_by(STUDY_ID) %>%
  filter(COLLECTION_DATE < start_circulate) %>% #make sure only use serum collected before the wave start
  ungroup() 

shiri_tvcoxph_sens_bris <- shiri_tvcoxph_bris %>%
  arrange(SPECIMEN_ID, COLLECTION_DATE) %>%  # Ensure the data is sorted by ID and date
  group_by(STUDY_ID) %>%
  filter(COLLECTION_DATE < start_circulate) %>% #make sure only use serum collected before the wave start
  ungroup() 

#participants missing in sensitivity analysis:
#phuket: 
anti_join(shiri_tvcoxph_phuk, shiri_tvcoxph_sens_phuk, by = "STUDY_ID") #1 additional HCP removed in sensitivity analysis for outside of collection period
check_phuk <- anti_join(shiri_aim3_cleaned, shiri_tvcoxph_sens_phuk) #unvax time of 5 and vax time of 38 removed 
length(unique(shiri_tvcoxph_sens_phuk$STUDY_ID)) #258 participants

#brisbane:
anti_join(shiri_tvcoxph_bris, shiri_tvcoxph_sens_bris, by = "STUDY_ID") #1 additional HCP removed in sensitivity analysis for outside of collection period
check_bris<- anti_join(shiri_aim3_cleaned, shiri_tvcoxph_sens_bris) #unvax time of 4 and vax time of 38 removed 
length(unique(shiri_tvcoxph_sens_bris$STUDY_ID)) #257 participants

#sensitivity analysis: cox proportional-hazards models
#crude NAI
fit.coxphtv_inf_bris_left_crude_nai_sens <- coxph(coxphtv_inf_bris_left_crude_nai, data = shiri_tvcoxph_sens_bris)
fit.coxphtv_inf_phuk_left_crude_nai_sens <- coxph(coxphtv_inf_phuk_left_crude_nai, data = shiri_tvcoxph_sens_phuk)
summary(fit.coxphtv_inf_bris_left_crude_nai_sens) 
summary(fit.coxphtv_inf_phuk_left_crude_nai_sens) 
#crude HAI
fit.coxphtv_inf_bris_left_crude_hai_sens <- coxph(coxphtv_inf_bris_left_crude_hai, data = shiri_tvcoxph_sens_bris)
fit.coxphtv_inf_phuk_left_crude_hai_sens <- coxph(coxphtv_inf_phuk_left_crude_hai, data = shiri_tvcoxph_sens_phuk)
summary(fit.coxphtv_inf_bris_left_crude_hai_sens) 
summary(fit.coxphtv_inf_phuk_left_crude_hai_sens) 

#adjusted models (including NAI and HAI modeled by each antigen)
fit.coxphtv_inf_bris_left_adj_sens <- coxph(coxphtv_inf_bris_left_adj, data = shiri_tvcoxph_sens_bris)
fit.coxphtv_inf_phuk_left_adj_sens <- coxph(coxphtv_inf_phuk_left_adj, data = shiri_tvcoxph_sens_phuk)
summary(fit.coxphtv_inf_bris_left_adj_sens) 
summary(fit.coxphtv_inf_phuk_left_adj_sens) 

### Cox proportional-hazards model tables ------------------------------------------------------------------
#Time-varying left truncated cox proportional hazards model
table_tvcoxph <- data.frame(
  Variables = c("Phuket NAI",
                "Phuket HAI",
                "Brisbane NAI",
                "Brisbane HAI"),
  Column_1 = c("0.58 (0.52,0.65)", "0.89 (0.78,1.01)", "0.75 (0.66,0.85)", "0.82 (0.70,0.97)"),
  Column_2 = c( "<0.001***",
                "0.068",
                "<0.001***",
                "0.022*"),
  Column_3 = c("0.58 (0.52,0.64)", "0.99 (0.87,1.13)", "0.77 (0.67,0.87)", "0.88 (0.74,1.05)"),
  Column_4 = c( "<0.001***",
                "0.925",
                "<0.001***",
                "0.155"),
  Block = c("B/Phuket (n=259)", "B/Phuket (n=259)", "B/Brisbane (n=258)", "B/Brisbane (n=258)"))

#add space between columns
table_tvcoxph$Spacer1 <- ""
table_tvcoxph$Spacer2 <- ""

#generate table
table_tvcoxph <- table_tvcoxph %>%
  gt(groupname_col = "Block") %>%
  cols_label(
    Variables = "Antigen",    
    Column_1 = "HR (95% CI)" ,
    Spacer1 = "",
    Column_2 = "p-value",
    Column_3 = "HR (95% CI)",
    Spacer2 = "",
    Column_4 = "p-value")%>%
  tab_spanner(
    label = "Crude",
    columns = c(Column_1, Spacer1, Column_2), replace = TRUE
  ) %>%
  tab_spanner(
    label = "Adjusted",
    columns = c(Column_3,Spacer2, Column_4), replace = TRUE
  ) %>%
  tab_style(
    style = cell_text(size = px(14), weight = "bold"),
    locations = cells_column_spanners(spanners = everything())) %>%
  tab_style(
    style = cell_text(size = px(14), weight = "bold"),
    locations = cells_column_labels(columns = everything())) %>%
  fmt_number(columns = c("Column_2", "Column_4"), decimals = 3) %>%
  tab_options(
    table.font.size = px(12),
    table.width = pct(80),
    row_group.font.weight = "bold")%>%
  tab_footnote(
    footnote = md("p < 0.05 ( * ), p < 0.01 ( ** ), p < 0.001 ( *** )"),
    locations = cells_column_labels(columns = c(Column_2,Column_4)))%>%
  tab_footnote(
    footnote = md("Adjusted models contained both NAI and HAI titer values of matched antigens as covariates."),
    locations = cells_column_spanners(spanners = "Adjusted"))%>%
  tab_footnote(
    footnote = "The unvaccinated time period of 6 HCP were removed for influenza vaccination received before/ on the date of S1 collection.",
    locations = cells_row_groups(groups = "B/Phuket (n=259)"))%>%
  tab_footnote(
    footnote = "The unvaccinated time period of 5 HCP were removed for influenza vaccination received before/ on the date of S1 collection.",
    locations = cells_row_groups(groups = "B/Brisbane (n=258)")
  )
table_tvcoxph
gtsave(table_tvcoxph, filename = file.path(secure_data, "Tables + Figures/table_tvcoxph.docx"))

# sensitivity analysis table
table_coxph_sens <- data.frame(
  Variables = c("Phuket NAI",
                "Phuket HAI",
                "Brisbane NAI",
                "Brisbane HAI"),
  Column_1 = c("0.59 (0.53,0.66)", "0.90 (0.79,1.02)", "0.76 (0.68,0.87)", "0.86 (0.73,1.01)"),
  Column_2 = c( "<0.001***",
                "0.108",
                "<0.001***",
                "0.073"),
  Column_3 = c("0.59 (0.53,0.66)", "0.99 (0.87,1.14)", "0.78 (0.68,0.88)", "0.91 (0.77,1.08)"),
  Column_4 = c( "<0.001***",
                "0.911",
                "<0.001***",
                "0.302"),
  Block = c("B/Phuket (n=258)", "B/Phuket (n=258)", "B/Brisbane (n=257)", "B/Brisbane (n=257)"))

#add space between columns
table_coxph_sens$Spacer1 <- ""
table_coxph_sens$Spacer2 <- ""

# Generate sensitivity analysis table
table_coxph_sens <- table_coxph_sens %>%
  gt(groupname_col = "Block") %>%
  cols_label(
    Variables = "Antigen",    
    Column_1 = "HR (95% CI)" ,
    Spacer1 = "",
    Column_2 = "p-value",
    Column_3 = "HR (95% CI)",
    Spacer2 = "",
    Column_4 = "p-value")%>%
  tab_spanner(
    label = "Crude",
    columns = c(Column_1, Spacer1, Column_2)
  ) %>%
  tab_spanner(
    label = "Adjusted",
    columns = c(Column_3,Spacer2, Column_4)
  ) %>%
  tab_style(
    style = cell_text(size = px(14), weight = "bold"),
    locations = cells_column_spanners(spanners = everything())) %>%
  tab_style(
    style = cell_text(size = px(14), weight = "bold"),
    locations = cells_column_labels(columns = everything())) %>%
  fmt_number(columns = c("Column_2", "Column_4"), decimals = 3) %>%
  tab_options(
    table.font.size = px(12),
    table.width = pct(80),
    row_group.font.weight = "bold")%>%
  tab_footnote(
    footnote = md("p < 0.05 ( * ), p < 0.01 ( ** ), p < 0.001 ( *** )"),
    locations = cells_column_labels(columns = c(Column_2,Column_4)))%>%
  tab_footnote(
    footnote = md("Adjusted models contained both NAI and HAI titer values of matched antigens as covariates."),
    locations = cells_column_spanners(spanners = "Adjusted"))%>%
  tab_footnote(
    footnote = "The unvaccinated time period of 5 HCP were removed for influenza vaccination received before/ on the date of S1 collection.",
    locations = cells_row_groups(groups = "B/Phuket (n=258)"))%>%
  tab_footnote(
    footnote = "The unvaccinated time period of 4 HCP were removed for influenza vaccination received before/ on the date of S1 collection.",
    locations = cells_row_groups(groups = "B/Brisbane (n=257)"))%>%
  tab_footnote(
    footnote = "The vaccinated period of 38 HCP were removed due to S2 collection date occurring after the start of observed influenza B circulation in study.",
    locations = cells_row_groups(groups = c("B/Phuket (n=258)", "B/Brisbane (n=257)")))
table_coxph_sens
gtsave(table_coxph_sens, filename = file.path(secure_data, "Tables + Figures/table_sensitivity.docx"))

# Vaccine Immunogenicity Analysis -------------------------------------------------------------
#Did vaccination change antibody levels?
#187 vaccinated participants

### Create NAI titer indicator variables 
#* 2-fold increase (difference of 1) 
#* 4-fold increase (difference of 2)

#calculate difference in NAI titers at each timepoint
shiri_aim3_wide_vax$ph_dif <- as.numeric(round(shiri_aim3_wide_vax$PHU_BYAM_NAI_LOG2_S2-shiri_aim3_wide_vax$PHU_BYAM_NAI_LOG2_S1))
shiri_aim3_wide_vax$bris_dif <- as.numeric(round(shiri_aim3_wide_vax$BRIS_BVIC_NAI_LOG2_S2-shiri_aim3_wide_vax$BRIS_BVIC_NAI_LOG2_S1))
shiri_aim3_wide_vax$ph_dif_hai <- shiri_aim3_wide_vax$PHU_BYAM_HAI_LOG2_S2-shiri_aim3_wide_vax$PHU_BYAM_HAI_LOG2_S1
shiri_aim3_wide_vax$bris_dif_hai <- shiri_aim3_wide_vax$BRIS_BVIC_HAI_LOG2_S2-shiri_aim3_wide_vax$BRIS_BVIC_HAI_LOG2_S1

#create variables of fold-rise/fall >= 2 and >=4
#B/Phuket
shiri_aim3_wide_vax <- fold_change_ind(shiri_aim3_wide_vax, PHU_BYAM_LOG_NAI_S1, PHU_BYAM_LOG_NAI_S2, ph_dif, ph_ge2, ph_ge4, ph_ge2_decrease, ph_ge4_decrease)
shiri_aim3_wide_vax <- fold_change_ind(shiri_aim3_wide_vax, PHU_BYAM_LOG_HAI_S1, PHU_BYAM_LOG_HAI_S2, ph_dif_hai, ph_ge2_hai, ph_ge4_hai, ph_ge2_hai_decrease, ph_ge4_hai_decrease)
#B/Brisbane
shiri_aim3_wide_vax <- fold_change_ind(shiri_aim3_wide_vax, BRIS_BVIC_LOG_NAI_S1, BRIS_BVIC_LOG_NAI_S2, bris_dif, bris_ge2, bris_ge4, bris_ge2_decrease, bris_ge4_decrease)
shiri_aim3_wide_vax <- fold_change_ind(shiri_aim3_wide_vax, BRIS_BVIC_LOG_HAI_S1, BRIS_BVIC_LOG_HAI_S2, bris_dif_hai, bris_ge2_hai, bris_ge4_hai, bris_ge2_hai_decrease, bris_ge4_hai_decrease)

#remove participants with missing S1 and S2 collection from analysis
shiri_aim3_wide_vax_bris <- shiri_aim3_wide_vax %>%
  drop_na(BRIS_BVIC_LOG_NAI_S1, BRIS_BVIC_LOG_NAI_S2, BRIS_BVIC_LOG_HAI_S1, BRIS_BVIC_LOG_HAI_S2)
shiri_aim3_wide_vax_phuk <- shiri_aim3_wide_vax %>%
  drop_na(PHU_BYAM_LOG_NAI_S1, PHU_BYAM_LOG_NAI_S2, PHU_BYAM_LOG_HAI_S1, PHU_BYAM_LOG_HAI_S2) 

#confirm only vaccinated participants
table(shiri_aim3_wide_vax_bris$VACC_STAT_2017) #184 participants in B/Brisbane analysis
table(shiri_aim3_wide_vax_phuk$VACC_STAT_2017) #185 participants in B/Phuket analysis

#number of HCP with 2/4-fold increasing/decreasing NAI titers 
#increasing titers
table(shiri_aim3_wide_vax_phuk$ph_ge2) 
sum(shiri_aim3_wide_vax_phuk$ph_ge2==1)/ nrow(shiri_aim3_wide_vax_phuk) #76 (41.1%)
table(shiri_aim3_wide_vax_phuk$ph_ge4) 
sum(shiri_aim3_wide_vax_phuk$ph_ge4==1)/ nrow(shiri_aim3_wide_vax_phuk) #41 (22.2%)
table(shiri_aim3_wide_vax_bris$bris_ge2) 
sum(shiri_aim3_wide_vax_bris$bris_ge2==1)/ nrow(shiri_aim3_wide_vax_bris) #78 (42.4%)
table(shiri_aim3_wide_vax_bris$bris_ge4) 
sum(shiri_aim3_wide_vax_bris$bris_ge4==1)/ nrow(shiri_aim3_wide_vax_bris) #40 (21.7%)
#decreasing titers
table(shiri_aim3_wide_vax_phuk$ph_ge2_decrease)
sum(shiri_aim3_wide_vax_phuk$ph_ge2_decrease==1)/ nrow(shiri_aim3_wide_vax_phuk) #55 (29.7%)
table(shiri_aim3_wide_vax_phuk$ph_ge4_decrease)
sum(shiri_aim3_wide_vax_phuk$ph_ge4_decrease==1)/ nrow(shiri_aim3_wide_vax_phuk) #30 (16.2%)
table(shiri_aim3_wide_vax_bris$bris_ge2_decrease)
sum(shiri_aim3_wide_vax_bris$bris_ge2_decrease==1)/ nrow(shiri_aim3_wide_vax_bris) #61 (33.2%)
table(shiri_aim3_wide_vax_bris$bris_ge4_decrease)
sum(shiri_aim3_wide_vax_bris$bris_ge4_decrease==1)/ nrow(shiri_aim3_wide_vax_bris) #36 (19.6%)

#number of HCP with 2/4-fold increasing/decreasing HAI titers 
#increasing titers
table(shiri_aim3_wide_vax_phuk$ph_ge2_hai) 
sum(shiri_aim3_wide_vax_phuk$ph_ge2_hai==1)/ nrow(shiri_aim3_wide_vax_phuk) #45 (24.3%)
table(shiri_aim3_wide_vax_phuk$ph_ge4_hai) 
sum(shiri_aim3_wide_vax_phuk$ph_ge4_hai==1)/ nrow(shiri_aim3_wide_vax_phuk) #13 (7.0%)
table(shiri_aim3_wide_vax_bris$bris_ge2_hai) 
sum(shiri_aim3_wide_vax_bris$bris_ge2_hai==1)/ nrow(shiri_aim3_wide_vax_bris) #64 (34.8%)
table(shiri_aim3_wide_vax_bris$bris_ge4_hai) 
sum(shiri_aim3_wide_vax_bris$bris_ge4_hai==1)/ nrow(shiri_aim3_wide_vax_bris) #16 (8.7%)
#decreasing titers
table(shiri_aim3_wide_vax_phuk$ph_ge2_hai_decrease)
sum(shiri_aim3_wide_vax_phuk$ph_ge2_hai_decrease==1)/ nrow(shiri_aim3_wide_vax_phuk) #30 (16.2%)
table(shiri_aim3_wide_vax_phuk$ph_ge4_hai_decrease)
sum(shiri_aim3_wide_vax_phuk$ph_ge4_hai_decrease==1)/ nrow(shiri_aim3_wide_vax_phuk) #6 (3.2%)
table(shiri_aim3_wide_vax_bris$bris_ge2_hai_decrease)
sum(shiri_aim3_wide_vax_bris$bris_ge2_hai_decrease==1)/ nrow(shiri_aim3_wide_vax_bris) #13 (7.1%)
table(shiri_aim3_wide_vax_bris$bris_ge4_hai_decrease)
sum(shiri_aim3_wide_vax_bris$bris_ge4_hai_decrease==1)/ nrow(shiri_aim3_wide_vax_bris) #2 (1.1%)

#Table of HCP titer fold changes 
#NAI
table_fold_changes <- data.frame(
  Variables = c("4-fold Titer Decrease",
                "2-fold Titer Decrease",
                "2-fold Titer Increase",
                "4-fold Titer Increase"),
  Column_1 = c("30 (16.2%)", "55 (29.7%)", "76 (41.1%)" , "41 (22.2%)"),
  Column_2 = c("36 (19.6%)", "61 (33.2%)", "78 (42.4%)", "40 (21.7%)")
)
table_fold_changes <- table_fold_changes %>%
  gt() %>%  
  cols_align(
    align = "left",   
    columns = everything())%>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())) %>%
  cols_label(
    Variables = "NAI Titer Fold Changes from Influenza Vaccination",
    Column_1 = "B/Phuket n=185",
    Column_2 = "B/Brisbane n=184") %>%
  tab_options( 
    heading.align = "left") %>%
  tab_footnote(
    footnote = md("Frequencies and column percentages reported as n(%)."),
    locations = cells_column_labels(columns = c(Column_1,Column_2)))
table_fold_changes

#HAI
table_fold_changes_hai <- data.frame(
  Variables = c("4-fold Titer Decrease",
                "2-fold Titer Decrease",
                "2-fold Titer Increase",
                "4-fold Titer Increase"),
  Column_1 = c("6 (3.2%)", "30 (16.2%)", "45 (24.3%)" , "13 (7.0%)"),
  Column_2 = c("2 (1.1%)", "13 (7.1%)", "64 (34.8%)", "16 (8.7%)")
)
table_fold_changes_hai <- table_fold_changes_hai %>%
  gt() %>%  
  cols_align(
    align = "left",   
    columns = everything())%>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())) %>%
  cols_label(
    Variables = "HAI Titer Fold Changes from Influenza Vaccination",
    Column_1 = "B/Phuket n=185",
    Column_2 = "B/Brisbane n=184") %>%
  tab_options( 
    heading.align = "left") %>%
  tab_footnote(
    footnote = md("Frequencies and column percentages reported as n(%)."),
    locations = cells_column_labels(columns = c(Column_1,Column_2)))
table_fold_changes_hai

#save tables
gtsave(table_fold_changes, filename = file.path(secure_data, "Tables + Figures/table_fold_changes.pdf"))
gtsave(table_fold_changes_hai, filename = file.path(secure_data, "Tables + Figures/table_fold_changes_hai.pdf"))


#plot changes in titer values from pre and post vaccination for each antigen
#NAI
#B/Phuket
phuket_fold_NAI_change <- ggplot(shiri_aim3_wide_vax_phuk, aes(x = reorder(STUDY_ID,ph_dif), y = ph_dif)) +
  geom_bar(stat = "identity", fill = "#f28e69") +
  #  ggtitle("B/Phuket NAI titer fold changes") +
  ylab("NAI Titer Fold Change") + 
  xlab("Participant") +
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(breaks = seq(-4, 6, 1), limits = c(-4, 6))
table(shiri_aim3_wide_vax_phuk$ph_dif)

#B/Brisbane
bris_fold_NAI_change <- ggplot(shiri_aim3_wide_vax_bris, aes(x = reorder(STUDY_ID,bris_dif), y = bris_dif)) +
  geom_bar(stat = "identity", fill= "#1E76C3") +
  #  ggtitle("B/Brisbane NAI titer fold changes") +
  ylab("NAI Titer Fold Change") + 
  xlab("Participant") +
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(breaks = seq(-4, 6, 1), limits = c(-4, 6)) 
table(shiri_aim3_wide_vax_bris$bris_dif)

#HAI
#B/Phuket
phuket_fold_HAI_change <- ggplot(shiri_aim3_wide_vax_phuk, aes(x = reorder(STUDY_ID,ph_dif_hai), y = ph_dif_hai)) +
  geom_bar(stat = "identity", fill = "#f28e69") +
  #  ggtitle("B/Phuket HAI titer fold changes") +
  ylab("HAI Titer Fold Change") + 
  xlab("Participant") +
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(breaks = seq(-7, 7, 1), limits = c(-7, 7))
table(shiri_aim3_wide_vax_phuk$ph_dif_hai)

#B/Brisbane
bris_fold_HAI_change <- ggplot(shiri_aim3_wide_vax_bris, aes(x = reorder(STUDY_ID,bris_dif_hai), y = bris_dif_hai)) +
  geom_bar(stat = "identity", fill= "#1E76C3") +
  #  ggtitle("B/Brisbane HAI titer fold changes") +
  ylab("HAI Titer Fold Change") + 
  xlab("Participant") +
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(breaks = seq(-7, 7, 1), limits = c(-7, 7)) 
table(shiri_aim3_wide_vax_bris$bris_dif_hai)

foldchange_patchwork_NAI<- phuket_fold_NAI_change + bris_fold_NAI_change
foldchange_patchwork_NAI
ggsave(plot = foldchange_patchwork_NAI, width = 12, height = 5, dpi = 300, filename = file.path(secure_data, "Tables + Figures/foldchange_nai.png"))
foldchange_patchwork_HAI <- phuket_fold_HAI_change + bris_fold_HAI_change
foldchange_patchwork_HAI
ggsave(plot = foldchange_patchwork_HAI, width = 12, height = 5, dpi = 300, filename = file.path(secure_data, "Tables + Figures/foldchange_hai.png"))

#compare paired S1 vs S2 NAI
#B/Phuket
t.test(shiri_aim3_wide_vax_phuk$PHU_BYAM_NAI_LOG2_S1, shiri_aim3_wide_vax_phuk$PHU_BYAM_NAI_LOG2_S2, paired = TRUE) #p-value = 0.08462
#wilcox.test(shiri_aim3_wide_vax_phuk$PHU_BYAM_NAI_LOG2_S1, shiri_aim3_wide_vax_phuk$PHU_BYAM_NAI_LOG2_S2, paired = TRUE) #p-value = 0.08062
#B/Brisbane
t.test(shiri_aim3_wide_vax_bris$BRIS_BVIC_NAI_LOG2_S1, shiri_aim3_wide_vax_bris$BRIS_BVIC_NAI_LOG2_S2, paired = TRUE) #p-value = 0.258
#wilcox.test(shiri_aim3_wide_vax_bris$BRIS_BVIC_NAI_LOG2_S1, shiri_aim3_wide_vax_bris$BRIS_BVIC_NAI_LOG2_S2, paired = TRUE) #p-value = 0.3582

#compare paired S1 and S2 HAI
#B/Phuket
t.test(shiri_aim3_wide_vax_phuk$PHU_BYAM_HAI_LOG2_S1,shiri_aim3_wide_vax_phuk$PHU_BYAM_HAI_LOG2_S2,paired = TRUE) #p-value= 0.0681
#wilcox.test(shiri_aim3_wide_vax_phuk$PHU_BYAM_HAI_LOG2_S1,shiri_aim3_wide_vax_phuk$PHU_BYAM_HAI_LOG2_S2,paired = TRUE) #p-value: 0.04792
#B/Brisbane
t.test(shiri_aim3_wide_vax_bris$BRIS_BVIC_HAI_LOG2_S1,shiri_aim3_wide_vax_bris$BRIS_BVIC_HAI_LOG2_S2,paired = TRUE) #p-value= 2.183e-08
#wilcox.test(shiri_aim3_wide_vax_bris$BRIS_BVIC_HAI_LOG2_S1,shiri_aim3_wide_vax_bris$BRIS_BVIC_HAI_LOG2_S2,paired = TRUE) #p-value: 8.413e-09

#GMT NAI in vaccinated HCP S1 vs. S2 
#B/Phuket S1
geometric_summary_log2(shiri_aim3_wide_vax_phuk$PHU_BYAM_NAI_LOG2_S1)# GMT(95% CI): 191.5 (158.3, 231.7) n=185, se = 0.1402476
#B/Phuket S2 
geometric_summary_log2(shiri_aim3_wide_vax_phuk$PHU_BYAM_NAI_LOG2_S2)#GMT(95% CI): 222.5 (185.0, 267.6) n=185, se = 0.1358419

#B/Brisbane S1
geometric_summary_log2(shiri_aim3_wide_vax_bris$BRIS_BVIC_NAI_LOG2_S1) #GMT(95% CI): 367.9 (298.1, 454.0) n=184, se= 0.1548649, geom_mean=8.5230151
#B/Brisbane S2 
geometric_summary_log2(shiri_aim3_wide_vax_bris$BRIS_BVIC_NAI_LOG2_S2) #GMT(95% CI): 407.2 (338.8, 489.5) n=184, se = 1.837287, geom_mean=8.6697542

#GMFR from vaccination
#NAI
gmfr(shiri_aim3_wide_vax_phuk$PHU_BYAM_NAI_LOG2_S1, shiri_aim3_wide_vax_phuk$PHU_BYAM_NAI_LOG2_S2) #B/Phuket NAI GMFR (95% CI) code: 1.2 (1.0, 1.4)
gmfr(shiri_aim3_wide_vax_bris$BRIS_BVIC_NAI_LOG2_S1, shiri_aim3_wide_vax_bris$BRIS_BVIC_NAI_LOG2_S2)#B/Brisbane NAI GMFR (95% CI) code: 1.1 (0.9, 1.3)


#GMT HAI in vaccinated HCP S1 vs. S2 
#B/Phuket S1
geometric_summary_log2(shiri_aim3_wide_vax_phuk$PHU_BYAM_HAI_LOG2_S1)# GMT(95% CI): 40.2 (33.0, 48.8) n=185, se = 0.1437081
#B/Phuket S2 
geometric_summary_log2(shiri_aim3_wide_vax_phuk$PHU_BYAM_HAI_LOG2_S2) #GMT(95% CI): 45.1 (37.1, 54.8) n=185, se = 0.1440624

#B/Brisbane S1
geometric_summary_log2(shiri_aim3_wide_vax_bris$BRIS_BVIC_HAI_LOG2_S1) #GMT(95% CI): 80.0 (69.4, 92.2) n=184, se = 0.1042572, geom_mean = 6.3219281
#B/Brisbane S2 
geometric_summary_log2(shiri_aim3_wide_vax_bris$BRIS_BVIC_HAI_LOG2_S2) #GMT(95% CI): 114.4 (100.1, 130.7) n=184, se = 0.0980856, geom_mean = 6.8382324

#GMFR from vaccination
gmfr(shiri_aim3_wide_vax_phuk$PHU_BYAM_HAI_LOG2_S1,shiri_aim3_wide_vax_phuk$PHU_BYAM_HAI_LOG2_S2) #B/Phuket HAI GMFR (95% CI) code: 1.1 (1.0, 1.3)
gmfr(shiri_aim3_wide_vax_bris$BRIS_BVIC_HAI_LOG2_S1, shiri_aim3_wide_vax_bris$BRIS_BVIC_HAI_LOG2_S2) #B/Brisbane HAI GMFR (95% CI) code: 1.4 (1.3, 1.6)

#GMT table S1 and S2 collection points
#NAI
table_gmt_q2 <- data.frame(
  Variables = c("B/Phuket, n=185",
                "B/Brisbane, n=184"),
  Column_1 = c("191.5 (158.3, 231.7)","367.9 (298.1, 454.0)" ),
  Column_2 = c("222.5 (185.0, 267.6)","407.2 (338.8, 489.5)"),
  Column_3 = c("1.2 (1.0, 1.4)", "1.1 (0.9, 1.3)"),
  Column_4 = c("0.085","0.258"))
table_gmt_q2 <- table_gmt_q2 %>%
  gt() %>%
  cols_align(
    align = "left",   
    columns = everything())%>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())) %>%
  cols_label(
    Variables = "NAI Titer Measurements",
    Column_1 = "S1 GMT Influenza Vaccinated",
    Column_2 = "S2 GMT Influenza Vaccinated",
    Column_3 = "GMFR",
    Column_4 = "p-value") %>%
  tab_options(
    heading.align = "left")%>%
  tab_footnote(
    footnote = md("p-value for test of baseline GMT by influenza vaccination status for each antigen using paired t-tests."),
    locations = cells_column_labels(columns = Column_4))%>%
  tab_footnote(
    footnote = md("*p* < 0.05 ( * ), *p* < 0.01 ( ** ), *p* < 0.001 ( *** )"),
    locations = cells_column_labels(columns = Column_4))%>%
  tab_footnote(
    footnote = md("GMT (95% CI)"),
    locations = cells_column_labels(columns = c(Column_1,Column_2))) %>%
  tab_footnote(
    footnote = md("GMFR (95% CI)"),
    locations = cells_column_labels(columns = Column_3))
table_gmt_q2
gtsave(table_gmt_q2, filename = file.path(secure_data, "Tables + Figures/q2_gmt_s1s2_table.pdf"))

table_gmt_q2_hai <- data.frame(
  Variables = c("B/Phuket, n=185",
                "B/Brisbane, n=184"),
  Column_1 = c("40.2 (33.0, 48.8)","45.1 (37.1, 54.8)" ),
  Column_2 = c("80.0 (69.4, 92.2)","114.4 (100.1, 130.7)"),
  Column_3 = c("1.1 (1.0, 1.3)", "1.4 (1.3, 1.6)"),
  Column_4 = c("0.068","<0.001***"))
table_gmt_q2_hai <- table_gmt_q2_hai %>%
  gt() %>%
  cols_align(
    align = "left",   
    columns = everything()
  )%>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())) %>%
  cols_label(
    Variables = "HAI Titer Measurements",
    Column_1 = "S1 GMT Influenza Vaccinated",
    Column_2 = "S2 GMT Influenza Vaccinated",
    Column_3 = "GMFR",
    Column_4 = "p-value") %>%
  tab_options(
    heading.align = "left")%>%
  tab_footnote(
    footnote = md("p-value for test of baseline GMT by influenza vaccination status for each antigen using paired t-tests."),
    locations = cells_column_labels(columns = Column_4))%>%
  tab_footnote(
    footnote = md("*p* < 0.05 ( * ), *p* < 0.01 ( ** ), *p* < 0.001 ( *** )"),
    locations = cells_column_labels(columns = Column_4))%>%
  tab_footnote(
    footnote = md("GMT (95% CI)"),
    locations = cells_column_labels(columns = c(Column_1,Column_2))) %>%
  tab_footnote(
    footnote = md("GMFR (95% CI)"),
    locations = cells_column_labels(columns = Column_3))
table_gmt_q2_hai
gtsave(table_gmt_q2_hai, filename = file.path(secure_data, "Tables + Figures/q2_gmt_s1s2_table_hai.pdf"))

## Supplement immunogenicity analysis --------------------------------------
#filter NA values (same population as immunogenicity, need to remove 1 from Phuket analysis for missing Brisbane data, n=184 for both antigens)
shiri_aim3_wide_vax_cor <- shiri_aim3_wide_vax %>%
  drop_na(PHU_BYAM_LOG_NAI_S1, PHU_BYAM_LOG_NAI_S2, PHU_BYAM_LOG_HAI_S1, PHU_BYAM_LOG_HAI_S2,BRIS_BVIC_LOG_NAI_S1, BRIS_BVIC_LOG_NAI_S2, BRIS_BVIC_LOG_HAI_S1, BRIS_BVIC_LOG_HAI_S2) 
#keep variables for correlation analysis
shiri_wide_cor <- shiri_aim3_wide_vax_cor %>%
  select(PHU_BYAM_NAI_LOG2_S1, BRIS_BVIC_NAI_LOG2_S1, PHU_BYAM_HAI_LOG2_S1, BRIS_BVIC_HAI_LOG2_S1,ph_dif, ph_dif_hai, bris_dif, bris_dif_hai)
cor_table <-cor(shiri_wide_cor)

#long format tables for merging 
p_calc <- rcorr(as.matrix(shiri_wide_cor))
p_val<-p_calc$P
p_val_table <- p_val %>%
  as_tibble(rownames = 'col_a')%>%
  pivot_longer(
    -col_a, 
    names_to = "col_b",
    values_to = "pval"
  )
cor_table <- cor_table %>%
  as_tibble(rownames = 'col_a')%>%
  pivot_longer(
    -col_a, 
    names_to = "col_b",
    values_to = "correlation"
  )
cor_table_comb <- left_join(cor_table, p_val_table)
cor_table_comb <- cor_table_comb %>%
  mutate(sig = ifelse(pval < .001, "***", 
                      ifelse(pval < .01, "**", 
                             ifelse(pval < .05, "*", ""))))%>%
  mutate(sig = replace_na(as.character(sig), ""))%>%
  mutate(corr_rounded = format(round(correlation, 2)))
#merge correlation and p-value significance cells
cor_table_comb$cor_sig <- paste(cor_table_comb$corr_rounded, cor_table_comb$sig)
#factor variables
variable_names <- unique(cor_table_comb$col_a)
cor_table_factored <- cor_table_comb %>%
  mutate(col_a= factor(col_a, levels = variable_names),
          col_b = factor(col_b, levels = rev(variable_names)))
#relabel variable names for heatmap
cor_table_relabeled <- cor_table_factored %>%
  mutate(col_a = fct_relabel(col_a,\(x) recode_values(x,
                                                      'PHU_BYAM_NAI_LOG2_S1' ~'S1 B/Phuket NAI',
                                                      'BRIS_BVIC_NAI_LOG2_S1' ~ 'S1 B/Brisbane NAI',
                                                      'PHU_BYAM_HAI_LOG2_S1' ~ 'S1 B/Phuket HAI',
                                                      'BRIS_BVIC_HAI_LOG2_S1' ~ 'S1 B/Brisbane HAI',
                                                      'ph_dif' ~ 'B/Phuket NAI titer fold-change',
                                                      'ph_dif_hai' ~ 'B/Phuket HAI titer fold-change',
                                                      'bris_dif' ~ 'B/Brisbane NAI titer fold-change',
                                                      'bris_dif_hai' ~ 'B/Brisbane HAI titer fold-change'
  )),
  col_b = fct_relabel(col_b,\(x) recode_values(x,
                                               'PHU_BYAM_NAI_LOG2_S1' ~'S1 B/Phuket NAI',
                                               'BRIS_BVIC_NAI_LOG2_S1' ~ 'S1 B/Brisbane NAI',
                                               'PHU_BYAM_HAI_LOG2_S1' ~ 'S1 B/Phuket HAI',
                                               'BRIS_BVIC_HAI_LOG2_S1' ~ 'S1 B/Brisbane HAI',
                                               'ph_dif' ~ 'B/Phuket NAI titer fold-change',
                                               'ph_dif_hai' ~ 'B/Phuket HAI titer fold-change',
                                               'bris_dif' ~ 'B/Brisbane NAI titer fold-change',
                                               'bris_dif_hai' ~ 'B/Brisbane HAI titer fold-change'
  )))

cor_table_relabeled_leveled <- cor_table_relabeled %>%
  mutate(sec_a = as.numeric(col_a),
         sec_b = as.numeric(col_b %>% fct_rev()),
         correlation = ifelse(sec_a<=sec_b, correlation, NA )) %>%
  mutate(cor_sig = if_else(is.na(correlation), NA, cor_sig))

cor_table_heatmap <- cor_table_relabeled_leveled%>%
  ggplot(aes(col_a, col_b))+
  geom_tile(aes(fill = correlation)) +
  geom_text(
    aes(label = cor_sig), 
    color = ifelse(abs(cor_table$correlation) > 0.8,
                   'white', 'black'))+
  theme_minimal(base_size = 16) +
  labs(fill = expression("Correlation (" ~ italic(r) ~")"),
       caption = "p < 0.05 ( * ), p < 0.01 ( ** ), p < 0.001 ( *** )")+
  scale_fill_gradient2(high = 'blue',
                       low= 'orange1',
                       na.value = 'white',
                       limits = c(-1,1))+
  coord_cartesian(expand = FALSE) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text.y = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y =  element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.9, 0.75),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 13, barheight = 2.5,
                               title.position = "top", title.hjust = 0.5))
cor_table_heatmap
ggsave(plot = cor_table_heatmap, width = 12, height = 9, dpi = 600, filename = file.path(secure_data, "Tables + Figures/rogers_supplement_figure2.tiff"))

#save datasets
#write_csv(shiri_aim3_studyids, file.path(secure_data, "Datasets/shiri_aim3_studyids.csv"))
#write_csv(shiri_aim3_wide, file.path(secure_data, "Datasets/SHIRI_AIM3_WIDE.csv"))
#write_csv(shiri_aim3_wide_vax, file.path(secure_data, "Datasets/SHIRI_AIM3_vaccinated.csv"))
#write_csv(shiri_aim3_cleaned, file.path(secure_data, "Datasets/SHIRI_AIM3_CLEANED_UPDATED.csv"))
#write_csv(shiri_aim3_wide_vax_phuk, file.path(secure_data, "Datasets/shiri_aim3_wide_vax_phuk.csv"))
#write_csv(shiri_aim3_wide_vax_bris, file.path(secure_data, "Datasets/shiri_aim3_wide_vax_bris.csv"))
#write_csv(shiri_aim3_wide_vax, file.path(secure_data, "Datasets/shiri_aim3_wide_vax.csv"))
