rm(list=ls()) 

library(readr)
library(caret)
library(rpart.plot)
library(pROC)
library(knitr)
library(corrplot)
library(FactoMineR)
library(tidyr)

################################################################################
#                                   Load Data                                  #
################################################################################
train_features <- read_csv("data/train_features.csv")
train_score <- read_csv("data/train_targets_scored.csv")
train_drug <- read_csv("data/train_drug.csv")

test_features <- read_csv("data/test_features.csv")


################################################################################
#                                   Exploration                                #
################################################################################
summary(train_features)
spec(train_features)

# Sélectionn des variables 
# matrice de corrélation
train_ft <- train_features %>% 
  dplyr::select(-sig_id, -cp_type, -cp_dose, -cp_time)

train_cor <- cor(train_ft)

# plot qui ne marche point
corrplot(train_cor, tl.cex = 0.5)

var <- findCorrelation(train_cor, cutoff = 0.95, names = T)
# pas de VA corrélés

train_colin <- caret::findLinearCombos(train_ft)
train_colin$remove
# pas de colin 

# clean dataset

# réduction des dimension et ACP
# visuslisation
train_pca <- PCA(train_ft)


          