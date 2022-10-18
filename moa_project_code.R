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
#                                   plot                                       #
################################################################################

barplot(summary(as.factor(train_features$cp_type)))
# a lot more treatement data than control

barplot(summary(as.factor(train_features$cp_time)))
# 3 duration of treatment

barplot(summary(as.factor(train_features$cp_dose)))
# high & low dose

hist(train_features$'g-0')
# distribution of g-0

hist(train_features$`c-0`)
# distribution of c-0

########## MoA
summary(train_score)

prop_moa <- rep(NA, size = length(train_score)[1])
for(i in 1:dim(train_score)[1]){
    prop_moa[i] = sum(train_score[i, -1])
}

barplot(summary(as.factor(prop_moa))*100/dim(train_score)[1], ylim=c(0, 60))
# We have almost 40% of our data that has 0 MoA
# Almost 50% have 1 MoA
# cleaning data en virant les médicmaent qui sont associées à aucun MoA
placebo <- train_features %>% 
  filter(cp_type == "ctl_vehicle")

# vire les placebo de train_score
score_no_vehicule <- train_score[train_features$cp_type == "trt_cp",]

# vire les MoA qui sont assignés à aucune catégorie 
score_no_vehicule <- score_no_vehicule %>% 
  filter(rowSums(across(where(is.numeric))) != 0)
  
# train features avec pas de placebo et pas de MoA
train_ft_noplacebo_nomoa <- score_no_vehicule %>% 
  rbind(train_score[which(train_score$sig_id %in% placebo$sig_id),])


train_clean <- train_ft %>% 
  filter(train_ft$sig_id %in% train_ft_noplacebo_nomoa$sig_id)

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


train_pca_data <- train_features %>% 
  dplyr::select(-cp_time, -cp_type, -cp_dose)

# reduction des dimension et ACP
# visualisation
color_cptype <- train_features %>% 
  filter(train_features$sig_id %in% train_clean$sig_id) %>% 
  dplyr::select(cp_type)

color_cptime <- train_features %>% 
  filter(train_features$sig_id %in% train_clean$sig_id) %>% 
  dplyr::select(cp_time)

color_cpdose <- train_features %>% 
  filter(train_features$sig_id %in% train_clean$sig_id) %>% 
  dplyr::select(cp_dose)

train_pca <- pca(train_clean[,-1])
plotIndiv(train_pca, var.names = F, group = color_cptype$cp_type)
plotIndiv(train_pca, var.names = F, group = color_cpdose$cp_dose)
plotIndiv(train_pca, var.names = F, group = color_cptime$cp_time)

plotVar(train_pca, var.names = F)
# on montre que la première dim est essentielle
plot(train_pca)

# PLS

# upload data
write.csv(train_clean, "train_clean.csv")
write.csv(train_ft, "train_ft.csv")
write.csv(train_features_clean, "train_feature_clean.csv")