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
    dplyr::select(-sig_id, -cp_type, -cp_dose, -cp_time)

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

# reduction des dimension et ACP
# visualisation
train_pca <- PCA(train_ft)


          