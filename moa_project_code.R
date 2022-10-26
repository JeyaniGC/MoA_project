rm(list=ls())

library(readr)
library(caret)
library(rpart.plot)
library(pROC)
library(knitr)
library(corrplot)
library(FactoMineR)
library(tidyr)
library(mixOmics)
library(dplyr)
library(ggplot2)
library(viridis)

################################################################################
#                                   Load Data                                  #
################################################################################
train_features <- read_csv("data/train_features.csv")
train_score <- read_csv("data/train_targets_scored.csv")
train_drug <- read_csv("data/train_drug.csv")

test_features <- read_csv("data/test_features.csv")


################################################################################
#                                   Clean data                                 #
################################################################################
# MoA
summary(train_score)

summary(train_features)
# features with 'g-' are gene expression data
# features with 'c-' are cell viability data
spec(train_features)
# columns with 'cp-' are information on the drugs used or not

which(is.na(train_features)==T)
# No NA in the dataset

# Selection of variables
train_ft <- train_features %>% 
    dplyr::select(-cp_type, -cp_dose, -cp_time)

# correlation matrix
train_cor <- cor(train_ft)
#corrplot(train_cor, tl.cex = 0.5)  # plot doesn't work
var <- findCorrelation(train_cor, cutoff = 0.95, names = T)
# pas de VA correlees

# Colinearity
train_colin <- caret::findLinearCombos(train_ft)
train_colin$remove
# pas de colin

################################################################################
#                               Prepare Data                                   #
################################################################################

# separe train_score in placebo and treatment
score_placebo <- train_score[train_features$cp_type == "ctl_vehicle",]
score_treatment <- train_score[train_features$cp_type == "trt_cp",]

# take out the samples without MoA 
score_treatment_noMoA <- score_treatment %>% 
  filter(rowSums(across(where(is.numeric))) != 0)

# Add a column noMoa for the placebo
score_placebo$noMoa <- rep(1, times=dim(score_placebo)[1])
score_treatment_noMoA$noMoa <- rep(0, times=dim(score_treatment_noMoA)[1])

# Add the placebo
score_clean <- rbind(score_treatment_noMoA, score_placebo)

# reduction y
# change values to numeric
score_t <- score_clean %>% 
  t()

target_sums <- score_clean %>% 
  dplyr::select(-sig_id) %>% 
  # compte l'occurence de chaque MoA equivalent d'un colSums
  summarise(across(everything(), sum)) %>% 
  # Passe les rows en col et inversement avec 2 colonnes target et sum
  pivot_longer(everything(), names_to = "target", values_to = "sum")

type <- target_sums %>%  
  separate(target, into = c("a", "b", "c", "d", "e", "type"), fill = "left") %>% 
  select(type) 

# transposée
colnames(score_t) <- score_t[1, ]
score_t <- score_t[-1, ]
dim(score_t)
# change values to numeric
score_t <- as.data.frame(apply(score_t, 2, as.numeric))
# target and type
target_type <- data.frame(target_sums$target, type)
# add target and their type to scores
type_score <- data.frame(target_type, score_t)
# tot score by id to verify
tot_score <- apply(type_score[-c(1, 2)], 2, sum)

inhibitor <- type_score[which(type_score$type == 'inhibitor'),]
inhib_score <- apply(inhibitor[-c(1, 2)], 2, sum)

antagonist <- type_score[which(type_score$type == 'antagonist'),]
anta_score <- apply(antagonist[-c(1, 2)], 2, sum)

agonist <- type_score[which(type_score$type == 'agonist'),]
ago_score <- apply(agonist[-c(1, 2)], 2, sum)

placebo <- type_score[which(type_score$type == 'noMoa'),]
pla_score <- apply(placebo[-c(1, 2)], 2, sum)

other <- type_score[!(type_score$type %in% c('inhibitor', 'antagonist', 'agonist', 'noMoa')),]
other_score <- apply(other[-c(1, 2)], 2, sum)


type <- c('inhibitor', 'antagonist', 'agonist', 'placebo', 'other')
score_resum <- as.data.frame(rbind(inhib_score, anta_score, ago_score, pla_score, other_score))
# totaux des scores resumés pr verif
tot_resum <- apply(score_resum, 2, sum)

# transpose pr avoir dtf comme au debut
score_resum.T <- as.data.frame(t(score_resum))
names(score_resum.T) <- type
# 0 ou 1
score_resum.T$inhibitor <- ifelse(score_resum.T$inhibitor>=1, 1, 0)
score_resum.T$antagonist <- ifelse(score_resum.T$antagonist>=1, 1, 0)
score_resum.T$agonist <- ifelse(score_resum.T$agonist>=1, 1, 0)
score_resum.T$placebo <- ifelse(score_resum.T$placebo>=1, 1, 0)
score_resum.T$other <- ifelse(score_resum.T$other>=1, 1, 0)



# clean train_features in same way
train_features_clean <- train_features %>%
  filter(train_features$sig_id %in% score_clean$sig_id)%>%
    dplyr::select(-cp_type, -cp_dose, -cp_time)
# normalise
train_features_clean <- cbind(train_features_clean$sig_id, 
                              data.frame(scale(train_features_clean[, -1])))
colnames(train_features_clean)[1] <- "sig_id"

for (i in 2:dim(train_features_clean)[2]){
  train_features_clean[,i] <- ifelse(train_features_clean[, i]<(-4), 0, train_features_clean[, i])
  train_features_clean[,i] <- ifelse(train_features_clean[, i]>4, 0, train_features_clean[, i])
}


# binarise
train_features_bool <- train_features_clean
for (i in 2:dim(train_features_bool)[2]){
  train_features_bool[,i] <- ifelse(train_features_bool[, i]>0, 1, 0)
}


# clean train_drug
train_drugs_clean <- train_drug %>%
  filter(train_drug$sig_id %in% score_clean$sig_id)

################################################################################
#                                   plot                                       #
################################################################################
colour <- c("Indianred", "royalblue", "Gold", "Olivedrab4")

# Proportion of placebo and treatment
# in rought data
train_features %>%
  ggplot(aes(x = cp_type, fill=cp_type)) +
  geom_bar(stat="count")
# a lot more treatement data than control


# Proportion of sample for the 3 time stamps
# in rought data
train_features %>%
  ggplot(aes(x = as.factor(cp_time)),) +
  geom_bar(stat="count", fill=c("#999999", "#E69F00", "#56B4E9"))
# 3 duration of treatment


# Proportion of the two dosage
# in rought data
train_features %>%
  ggplot(aes(x = cp_dose, fill=cp_dose)) +
  geom_bar(stat="count")
# high & low dose


# Histogram of the first gene expression
# in rought data
hist(train_features$'g-0', col=colour[1], 
     main="Distribution de l'expression du gène 'g-0' des molécules thérapeutiques")
# distribution of g-0
hist(train_features_clean$g.0, col=colour[1], 
     main="Distribution de l'expression du gène 'g-0' des molécules thérapeutiques")

barplot(summary(as.factor(train_features_clean$g.0)), col=colour,
        main="Distribution de l'expression du gène 'g-0'\n des molécules thérapeutiques après pre-processing")
# distribution of g-0


#  Histogram of the first cell viability
# in rought data
hist(train_features$'c-0', col=colour[1],
     main="Distribution de la viabilité cellulaire de 'c-0'\n des molécules thérapeutiques")
# distribution of c-0
# in clean data
hist(train_features_clean$c.0, col=colour[1],
     main="Distribution de la viabilité cellulaire de 'c-0'\n des molécules thérapeutiques")

barplot(summary(as.factor(train_features_clean$c.0)), col=colour,
        main="Distribution de la viabilité cellulaire de 'c-0'\n des molécules thérapeutiques après pre-processing")
# distribution of c-0


# Distribution of MoA attributed
# In rought score data
prop_moa <- rep(NA, size = length(train_score)[1])
for(i in 1:dim(train_score)[1]){
  prop_moa[i] = sum(train_score[i, -1])
}

barplot(summary(as.factor(prop_moa))*100/dim(train_score)[1], 
        ylim=c(0, 60), col=colour, main="Compte du nombre de MoA par échantillon avant nettoyage")
# We have almost 40% of our data that has 0 MoA
# Almost 50% have 1 MoA
# In clean score data
prop_moa_clean <- rep(NA, size = length(score_clean)[1])
for(i in 1:dim(score_clean)[1]){
  prop_moa_clean[i] = sum(score_clean[i, -1])
}

barplot(summary(as.factor(prop_moa_clean))*100/dim(score_clean)[1], 
        ylim=c(0, 100), col=colour, main="Compte du nombre de MoA par échantillon après nettoyage")



prop_moa_bool <- rep(NA, size = length(score_resum.T)[1])
for(i in 1:dim(score_resum.T)[1]){
  prop_moa_bool[i] = sum(score_resum.T[i, -1])
}

barplot(summary(as.factor(prop_moa_bool))*100/dim(score_clean)[1], 
        ylim=c(0, 100), col=colour, main="Compte du nombre de MoA par échantillon après nettoyage")



# Summary of the different MoA
target_sums <- train_score %>% 
  dplyr::select(-sig_id) %>% 
  # compte l'occurence de chaque MoA equivalent d'un colSums
  summarise(across(everything(), sum)) %>% 
  # Passe les rows en col et inversement avec 2 colonnes target et sum
  pivot_longer(everything(), names_to = "target", values_to = "sum")

target_sums %>% 
  separate(target, into = c("a", "b", "c", "d", "e", "type"), fill = "left") %>% 
  count(type) %>% 
  add_tally(n, name = "total") %>% 
  mutate(perc = n/total) %>% 
  filter(n > 1) %>% 
  ggplot(aes(reorder(type, n, FUN = min), n, fill = n)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.2f%%", perc*100)), nudge_y = 6) +
  theme_minimal()

target_sums_clean <- score_clean %>% 
  dplyr::select(-sig_id) %>% 
  # compte l'occurence de chaque MoA equivalent d'un colSums
  summarise(across(everything(), sum)) %>% 
  # Passe les rows en col et inversement avec 2 colonnes target et sum
  pivot_longer(everything(), names_to = "target", values_to = "sum")


################################################################################
#                                    ACP                                       #
################################################################################

train_pca_data <- train_features %>% 
  dplyr::select(-cp_time, -cp_type, -cp_dose)

color_cptype <- train_features %>% 
  filter(train_features$sig_id %in% train_features_clean$sig_id) %>% 
  dplyr::select(cp_type)

color_cptime <- train_features %>% 
  filter(train_features$sig_id %in% train_features_clean$sig_id) %>% 
  dplyr::select(cp_time)

color_cpdose <- train_features %>% 
  filter(train_features$sig_id %in% train_features_clean$sig_id) %>% 
  dplyr::select(cp_dose)

train_pca <- pca(train_features_clean[,-1])
plotIndiv(train_pca, var.names = F, group = color_cptype$cp_type)
plotIndiv(train_pca, var.names = F, group = color_cpdose$cp_dose)
plotIndiv(train_pca, var.names = F, group = color_cptime$cp_time)

plotVar(train_pca, var.names = F)
# on montre que la premiere dim est essentielle
plot(train_pca, ylim=c(0, 0.4))

# PLSDA


################################################################################
#                               Upload Donnees                                 #
################################################################################

# upload train_features_clean
write.csv(train_features_clean, "train_features_clean.csv")
# upload train_features_bool (binarised)
write.csv(train_features_bool, 'train_features_boolean.csv')

# upload score_clean (normalised and with a column no MoA)
write.csv(score_clean, "train_score_clean.csv")
# upload score_reum.T (binarised)
write.csv(score_resum.T, 'train_score_clean_boolean.csv')

# upload train_drug_clean
write.csv(train_drugs_clean, "train_drug_clean.csv")




