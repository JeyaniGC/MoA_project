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

score_sums <- score_clean %>% 
  dplyr::select(-sig_id) %>% 
  # compte l'occurence de chaque MoA equivalent d'un colSums
  summarise(across(everything(), sum)) %>% 
  pivot_longer(everything(), names_to = "target", values_to = "sum")

target <- head(sort(as.numeric(score_sums$sum), decreasing = T), 20)
score_sums_clean <- score_sums[which(score_sums$sum %in% target),]
sum(head(sort(as.numeric(score_sums$sum), decreasing = T), 20))

other_sum <- sum(as.numeric(score_sums$sum[!(score_sums$sum %in% target)]))
other <- c('other', other_sum)
score_sums_clean <- rbind(score_sums_clean, other)

names <- colnames(score_clean)
other_target <- score_clean[,!(names %in% score_sums_clean$target)]
sum_otar <- rowSums(other_target[-1])
other <- data.frame(other = sum_otar)

train_score_clean <- cbind(sig_id = score_clean$sig_id, 
                           score_clean[,names %in% score_sums_clean$target],
                           other)
train_score_clean$other <- ifelse(train_score_clean$other>=1, 1, 0)



# clean train_features in same way
train_features_clean <- train_features %>%
  filter(train_features$sig_id %in% score_clean$sig_id)%>%
    dplyr::select(-cp_type, -cp_dose, -cp_time)

# cut out outlayers
for (i in 2:dim(train_features_clean)[2]){
  train_features_clean[,i] <- ifelse(train_features_clean[, i]<(-4), 0, train_features_clean[, i])
  train_features_clean[,i] <- ifelse(train_features_clean[, i]>4, 0, train_features_clean[, i])
}

# normalise
train_features_clean <- cbind(train_features_clean$sig_id, 
                              data.frame(scale(train_features_clean[, -1])))
colnames(train_features_clean)[1] <- "sig_id"


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
     main="Distribution de l'expression du g??ne 'g-0' des mol??cules th??rapeutiques")
# distribution of g-0
hist(train_features_clean$g.0, col=colour[1], 
     xlim=c(-5, 5), 
     main="Distribution de l'expression du g??ne 'g-0' des mol??cules th??rapeutiques")


#  Histogram of the first cell viability
# in rought data
hist(train_features$'c-0', col=colour[1],
     main="Distribution de la viabilit?? cellulaire de 'c-0'\n des mol??cules th??rapeutiques")
# distribution of c-0
# in clean data
hist(train_features_clean$c.0, col=colour[1], 
     xlim=c(-5, 5),
     main="Distribution de la viabilit?? cellulaire de 'c-0'\n des mol??cules th??rapeutiques")


# Distribution of MoA attributed
# In rought score data
prop_moa <- rep(NA, size = length(train_score)[1])
for(i in 1:dim(train_score)[1]){
  prop_moa[i] = sum(train_score[i, -1])
}

barplot(summary(as.factor(prop_moa))*100/dim(train_score)[1], 
        ylim=c(0, 60), col=colour, main="Compte du nombre de MoA par ??chantillon avant nettoyage")
# We have almost 40% of our data that has 0 MoA
# Almost 50% have 1 MoA
# In clean score data
prop_moa_clean <- rep(NA, size = length(score_clean)[1])
for(i in 1:dim(score_clean)[1]){
  prop_moa_clean[i] = sum(score_clean[i, -1])
}

barplot(summary(as.factor(prop_moa_clean))*100/dim(score_clean)[1], 
        ylim=c(0, 100), col=colour, main="Compte du nombre de MoA par ??chantillon apr??s nettoyage")



prop_moa_bool <- rep(NA, size = length(score_resum.T)[1])
for(i in 1:dim(score_resum.T)[1]){
  prop_moa_bool[i] = sum(score_resum.T[i, -1])
}

barplot(summary(as.factor(prop_moa_bool))*100/dim(score_clean)[1], 
        ylim=c(0, 100), col=colour, main="Compte du nombre de MoA par ??chantillon apr??s nettoyage")



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

# upload score_clean (normalised and with a column no MoA)
write.csv(score_clean, "train_score_clean.csv")
# upload score_reum.T (binarised)
write.csv(train_score_clean, 'train_score_clean_best.csv')

# upload train_drug_clean
write.csv(train_drugs_clean, "train_drug_clean.csv")





