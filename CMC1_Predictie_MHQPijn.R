# Development, training, and validation of prediction model for MHQ Pain Score
# after surgery for CMC1 OA
#
##########################################################################################
# Code written by N.L. Loos & L. Hoogendam for the Hand Wrist Study Group, 
# Erasmus Medical Center, Rotterdam
##########################################################################################
# 
# This script prepares data for analysis. It splits the data into a training, validation and test
# data set. And fits several models to this data: Generalized Linear Model (GLM), Gradient Boosting 
# Machine (GBM) and Random Forest (RF) using different resampling methods. The models are trained 
# and validated on separate data sets. 
# Hopefully this script might provide those who are interested in developing prediction models with 
# some tips on how to perform this on their own data. 
#
# Note: this script is specific to our data. For this script to be used on different data set some 
# adaptation will be needed.

rm(list = ls())

library(caret)
library(stringdist)
library(gbm)
library(tidyverse)
library(here)
library(RANN)
library(ggpubr)
library(ggplot2)
library(pROC)

select <- dplyr::select
rename <- dplyr::rename

# load data set for developing prediction models. Note: this script is not an ETL script. 
# Some pre-processing will probably be required. 

##########################################################################################
# Pre-processing 
##########################################################################################

# Removing redundant outcome variables and checking whether variables are in the data type.
# Remove all variables you won't be using for the development of your models. We also created
# a new variable whether the dominant hand of the patient was treated. 
CMC_OK_dataset <- CMC_OK_dataset %>%
  select(-c(Patient.traject.ID, Respondent.ID, delta_MHQPijn, Score_Pijn_MHQ12m, behandeling, 
            ipqscore, contains('rounddescription'), contains('926'), Vingers, type_behandeling, 
            organization, trackname, Survey.ID, location, Invuldatum, mrsaAsielzoeker, drugs, contains('12m'), contains('uitslag'))) %>%
  droplevels() %>% 
  mutate(jamar_punt2_links = ifelse(nchar(jamar_punt2_links)>2, 
                                    as.numeric(str_sub(jamar_punt2_links, 
                                                       start = 2)), as.numeric(jamar_punt2_links))) %>%
  mutate(jamar_punt2_rechts = ifelse(nchar(jamar_punt2_rechts)>2, 
                                     as.numeric(str_sub(jamar_punt2_rechts, 
                                                        start = 2)), as.numeric(jamar_punt2_rechts))) %>%
  mutate(pinch_sleutel_rechts = ifelse(nchar(pinch_sleutel_rechts)>2, 
                                       as.numeric(str_sub(pinch_sleutel_rechts, 
                                                          start = 2)), as.numeric(pinch_sleutel_rechts))) %>%
  mutate(pinch_punt3_links = ifelse(nchar(pinch_punt3_links)>2, 
                                    as.numeric(str_sub(pinch_punt3_links, 
                                                       start = 2)), as.numeric(pinch_punt3_links))) %>%
  mutate(pinch_punt3_rechts = ifelse(nchar(pinch_punt3_rechts)>2, 
                                     as.numeric(str_sub(pinch_punt3_rechts, 
                                                        start = 2)), as.numeric(pinch_punt3_rechts))) %>%
  mutate(pinch_punt2_rechts = ifelse(nchar(pinch_punt2_rechts)>2, 
                                     as.numeric(str_sub(pinch_punt2_rechts,
                                                        start = 2)), as.numeric(pinch_punt2_rechts))) %>% 
  mutate(Dominant_behandeld = ifelse(zijde_EQ5D == "Rechter hand" & dominant == "Rechts", "Ja", 
                                     ifelse(zijde_EQ5D == "Linker hand" & dominant == "Links", "Ja", "Nee"))) %>% 
  mutate(hoeLangKlacht = replace(hoeLangKlacht, hoeLangKlacht>200, NA))


CMC_OK_dataset <- rename(CMC_OK_dataset, c("zijde" = "zijde_EQ5D"))
CMC_OK_dataset$Geslacht <- as.factor(CMC_OK_dataset$Geslacht)
CMC_OK_dataset$side <- as.factor(CMC_OK_dataset$side)
CMC_OK_dataset$Dominant_behandeld <- factor(CMC_OK_dataset$Dominant_behandeld, levels = c("Nee", "Ja"))

##########################################################################################

# Removing impossible hands strength measurements
CMC_OK_dataset <- CMC_OK_dataset %>%
  mutate(jamar_punt2_links  = ifelse(jamar_punt2_links > 990, NA, jamar_punt2_links)) %>%
  mutate(jamar_punt2_rechts  = ifelse(jamar_punt2_rechts > 990, NA, jamar_punt2_rechts)) %>% 
  mutate(pinch_sleutel_links = ifelse(pinch_sleutel_links > 90, NA, pinch_sleutel_links)) %>%
  mutate(pinch_sleutel_rechts = ifelse(pinch_sleutel_rechts > 90, NA, pinch_sleutel_rechts)) 

# Creating new variables with strength measurements for the affected side
# Jamar punt2
CMC_OK_dataset <- CMC_OK_dataset %>%
  dplyr::mutate(jamar_punt2_aangedane_hand = 
                  ifelse(CMC_OK_dataset$zijde == "Linker hand", jamar_punt2_links, 
                         ifelse(CMC_OK_dataset$zijde == "Rechter hand", jamar_punt2_rechts, NA))) 
CMC_OK_dataset <- CMC_OK_dataset %>%
  dplyr::mutate(jamar_punt2_niet_aangedane_hand = 
                  ifelse(CMC_OK_dataset$zijde == "Linker hand", jamar_punt2_rechts, 
                         ifelse(CMC_OK_dataset$zijde == "Rechter hand", jamar_punt2_links, NA))) 

# Pinch_sleutel
CMC_OK_dataset <- mutate(CMC_OK_dataset, pinch_sleutel_aangedane_hand = 
                           ifelse(CMC_OK_dataset$zijde == "Linker hand", pinch_sleutel_links, 
                                  ifelse(CMC_OK_dataset$zijde == "Rechter hand", pinch_sleutel_rechts, NA)))

CMC_OK_dataset <- mutate(CMC_OK_dataset, pinch_sleutel_niet_aangedane_hand = 
                           ifelse(CMC_OK_dataset$zijde == "Linker hand", pinch_sleutel_rechts, 
                                  ifelse(CMC_OK_dataset$zijde == "Rechter hand", pinch_sleutel_links, NA)))

# Pinch_punt3
CMC_OK_dataset <- mutate(CMC_OK_dataset, pinch_punt3_aangedane_hand = 
                           ifelse(CMC_OK_dataset$zijde == "Linker hand", pinch_punt3_links, 
                                  ifelse(CMC_OK_dataset$zijde == "Rechter hand", pinch_punt3_rechts, NA)))

CMC_OK_dataset <- mutate(CMC_OK_dataset, pinch_punt3_niet_aangedane_hand = 
                           ifelse(CMC_OK_dataset$zijde == "Linker hand", pinch_punt3_rechts, 
                                  ifelse(CMC_OK_dataset$zijde == "Rechter hand", pinch_punt3_links, NA)))

# Pinch_punt2
CMC_OK_dataset <- mutate(CMC_OK_dataset, pinch_punt2_aangedane_hand = 
                           ifelse(CMC_OK_dataset$zijde == "Linker hand", pinch_punt2_links, 
                                  ifelse(CMC_OK_dataset$zijde == "Rechter hand", pinch_punt2_rechts, NA)))

CMC_OK_dataset <- mutate(CMC_OK_dataset, pinch_punt2_niet_aangedane_hand = 
                           ifelse(CMC_OK_dataset$zijde == "Linker hand", pinch_punt2_rechts, 
                                  ifelse(CMC_OK_dataset$zijde == "Rechter hand", pinch_punt2_links, NA)))

# Removing strength measurements from the not affected side, and old strength measurements
CMC_OK_dataset <- CMC_OK_dataset %>%
  select(-c(jamar_punt2_links:pinch_punt2_rechts, contains('niet'))) %>%
  droplevels()

##########################################################################################
# Random Train-Test-Validation Split

# The data (100%) is split into:
# training (~60%)
# validation (~20%)
# testing (~20%)
##########################################################################################

set.seed(123)
inTrain <- as.data.frame(createDataPartition(CMC_OK_dataset$MCID_MHQPijn, 
                                             p = 0.60, 
                                             list = FALSE, 
                                             times = 1))

train <- CMC_OK_dataset[inTrain$Resample1, ]
test_plus_val <- CMC_OK_dataset[-inTrain$Resample1, ]

set.seed(123)
in_test_plus_val <- as.data.frame(createDataPartition(test_plus_val$MCID_MHQPijn,
                                                      p = 0.5,
                                                      list = FALSE,
                                                      times = 1))

test <- test_plus_val[-in_test_plus_val$Resample1, ]
val <- test_plus_val[in_test_plus_val$Resample1, ]

# Making Dummy variables
fit <- MCID_MHQPijn ~.
options(na.action='na.pass')
dummies <- model.matrix(fit, data = train) %>%  
  as.data.frame()
y <- train$MCID_MHQPijn
dummies <- cbind(dummies, y)
names(dummies)[ncol(dummies)] <- "MCID_MHQPijn"

# Centering, scaling and imputation
preProcValues <- preProcess(as.data.frame(dummies), method = c("center", "scale", "nzv", "knnImpute")) # Geeft error als df nog een tbl_df is

trainTransformed <- predict(preProcValues, dummies) %>% 
  select(MCID_MHQPijn, everything()) %>% 
  select_if(~ !any(is.na(.)))
x_frame <- trainTransformed %>% 
  select(-MCID_MHQPijn)
y <- trainTransformed$MCID_MHQPijn
trainTransformed$MCID_MHQpijn <- factor(trainTransformed$MCID_MHQPijn)
sapply(trainTransformed, function(x) class(x))

##########################################################################################
# Model Development 
    # Recursive Feature Elimination (RFE) 
    # Without sampling, with downsampling and with upsampling
##########################################################################################

# RFE for GLM model without sampling
set.seed(123)
trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummary)

new_lr <- lrFuncs
new_lr$summary <- twoClassSummary

ctrl<-rfeControl(functions=new_lr, method = "repeatedcv", repeats = 5)
subsets <- c(1:10, 15, 30)
set.seed(123)

glmProfile <-rfe(x = data.frame(model.matrix(~., data = x_frame)[,-1]),
                                  y,
                                  sizes =subsets,
                                  rfeControl = ctrl,
                                  trControl =trainctrl,
                                  metric = "ROC")
trellis.par.set(caretTheme())
plot(glmProfile, type = c("g", "o"),  main = "RFE without Sampling for GLM")
predictors(glmProfile)

# RFE for GLM model with downsampling
lr_fit <- function (x, y, first, last, ...)
{tmp_2 <- caret::downSample(x, y)
tmp <- if (is.data.frame(tmp_2))
  tmp_2
else as.data.frame(tmp_2, stringsAsFactors = TRUE)
a <- glm(Class ~ ., data = tmp, family = "binomial")
}

new_lr_down <- lrFuncs
new_lr_down$summary <- twoClassSummary
new_lr_down$fit <- lr_fit

ctrl<-rfeControl(functions=new_lr_down, method = "repeatedcv", repeats = 5)
set.seed(123)

glmProfile_down <-rfe(x = data.frame(model.matrix(~., data = x_frame)[,-1]),
                              y,
                              sizes =subsets,
                              rfeControl = ctrl,
                              trControl =trainctrl,
                              metric = "ROC")
trellis.par.set(caretTheme())
plot(glmProfile_down, type = c("g", "o"),  main = "RFE with Downsampling for GLM")
predictors(glmProfile_down)

# RFE for GLM model with upsampling
lr_fit <- function (x, y, first, last, ...)
{tmp_2 <- caret::upSample(x, y)
tmp <- if (is.data.frame(tmp_2))
  tmp_2
else as.data.frame(tmp_2, stringsAsFactors = TRUE)
a <- glm(Class ~ ., data = tmp, family = "binomial")
}

new_lr_up <- lrFuncs
new_lr_up$summary <- twoClassSummary
new_lr_up$fit <- lr_fit

ctrl<-rfeControl(functions=new_lr_up, method = "repeatedcv", repeats = 5)
set.seed(123)

glmProfile_up <-rfe(x = data.frame(model.matrix(~., data = x_frame)[,-1]),
                            y,
                            sizes =subsets,
                            rfeControl = ctrl,
                            trControl =trainctrl,
                            metric = "ROC")
trellis.par.set(caretTheme())
plot(glmProfile_up, type = c("g", "o"),  main = "RFE with Upsampling for GLM")
predictors(glmProfile_up)

##########################################################################################

# RFE for GBM model without sampling
trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummary)
caretFuncs$summary <- twoClassSummary
ctrl<-rfeControl(functions=caretFuncs, method = "repeatedcv", repeats=5)

set.seed(123)
gbmProfile<-rfe(data.frame(model.matrix(~., data = x_frame)[,-1]),
                y,
                sizes =subsets, 
                rfeControl = ctrl,
                trControl =trainctrl,
                method = "gbm",
                metric = "ROC")

trellis.par.set(caretTheme())
plot(gbmProfile, type = c("g", "o"), main = "RFE without sampling for GBM")
predictors(gbmProfile)

# RFE for GBM model with downsampling
gbm_down <- caretFuncs
gbm_down$summary <- twoClassSummary
gbm_down$fit <- function(x, y, first, last, ...){
  
  df_down <- caret::downSample(x, y)
  
  train(select(df_down, -Class),
        df_down$Class,
        ...)
}
ctrl<-rfeControl(functions=gbm_down, method = "repeatedcv", repeats=5)

set.seed(123)
gbmProfile_down<-rfe(data.frame(model.matrix(~., data = x_frame)[,-1]),
                     y,
                     sizes =subsets,
                     rfeControl = ctrl,
                     trControl =trainctrl,
                     method = "gbm",
                     metric = "ROC",
                     verbose = F)

trellis.par.set(caretTheme())
plot(gbmProfile_down, type = c("g", "o"), main = "RFE with downsampling for GBM")
predictors(gbmProfile_down)

# RFE for GBM model with upsampling
gbm_up <- caretFuncs
gbm_up$summary <- twoClassSummary
gbm_up$fit <- function(x, y, first, last, ...){
  
  df_up <- caret::upSample(x, y)
  
  train(select(df_up, -Class),
        df_up$Class,
        ...)
}
ctrl<-rfeControl(functions=gbm_up, method = "repeatedcv", repeats=5)

set.seed(123)
gbmProfile_up<-rfe(data.frame(model.matrix(~., data = x_frame)[,-1]),
                   y,
                   sizes =subsets,
                   rfeControl = ctrl,
                   trControl =trainctrl,
                   method = "gbm",
                   metric = "ROC",
                   verbose = F)

trellis.par.set(caretTheme())
plot(gbmProfile_up, type = c("g", "o"), main = "RFE with upsampling for GBM")
predictors(gbmProfile_up)

##########################################################################################

# RFE for RF model without sampling
rfFuncs$summary <- twoClassSummary
ctrl<-rfeControl(functions=rfFuncs,method = "repeatedcv", repeats = 5)
subsets <- c(1:10, 15, 30)
set.seed(123)
rfProfile<-rfe(data.frame(model.matrix(~., data = x_frame)[,-1]),
               y,
               sizes =subsets, 
               rfeControl = ctrl,
               trControl =trainctrl,
               metric = "ROC")

trellis.par.set(caretTheme())
plot(rfProfile, type = c("g", "o"), main  = "RFE without sampling for RF")
predictors(rfProfile)

# RFE for RF model with downsampling
library(randomForest)
set.seed(123)
trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummary)
rfFuncs_down <- rfFuncs
rfFuncs_down_fit <- function(x, y, first, last, ...){
  loadNamespace("randomForest")
  
  df_down <- caret::downSample(x, y)
  
  randomForest::randomForest(
    select(df_down, -Class),
    df_down$Class,
    importance = (first | last),
    ...)
}

rfFuncs_down$fit <- rfFuncs_down_fit
rfFuncs_down$summary <- twoClassSummary
ctrl<-rfeControl(functions=rfFuncs_down,method = "repeatedcv", repeats =5)
subsets <- c(1:10, 15, 30)
set.seed(123)

rfProfile_down<-rfe(data.frame(model.matrix(~., data = x_frame)[,-1]),
                    y,
                    sizes =subsets, 
                    rfeControl = ctrl,
                    trControl =trainctrl,
                    metric = "ROC")

trellis.par.set(caretTheme())
plot(rfProfile_down, type = c("g", "o"), main = "RFE with downsampling for RF")

predictors(rfProfile_down)

# RFE for RF model with upsampling
set.seed(123)
trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummary)
rfFuncs_up <- rfFuncs
rfFuncs_up_fit <- function(x, y, first, last, ...){
  loadNamespace("randomForest")
  
  df_up <- caret::upSample(x, y)
  
  randomForest::randomForest(
    select(df_up, -Class),
    df_up$Class,
    importance = (first | last),
    ...)
}

rfFuncs_up$fit <- rfFuncs_up_fit
rfFuncs_up$summary <- twoClassSummary

ctrl<-rfeControl(functions=rfFuncs_up, method = "repeatedcv", repeats =5)
subsets <- c(1:10, 15, 30)
set.seed(123)
rfProfile_up<-rfe(data.frame(model.matrix(~., data = x_frame)[,-1]),
                  y,
                  sizes =subsets, 
                  rfeControl = ctrl,
                  trControl =trainctrl,
                  metric = "ROC")

trellis.par.set(caretTheme())
plot(rfProfile_up, type = c("g", "o"), main = "RFE with upsampling for RF")

predictors(rfProfile_up)

##########################################################################################

# Store results of RFE
rfe_res <- list(
  "glmProfile_up" = glmProfile_up,
  "glmProfile_down" = glmProfile_down,
  "glmProfile" = glmProfile,
  "rfProfile_up" = rfProfile_up,
  "rfProfile_down" = rfProfile_down,
  "rfProfile" = rfProfile,
  "gbmProfile_up" = gbmProfile_up,
  "gbmProfile_down" = gbmProfile_down,
  "gbmProfile" = gbmProfile
)

save(rfe_res, file = "Results_rfe_MCID_MHQ_CMC_231220.Rdata")

# Checking performance of RFE models and choose the best performing sampling method for each model.
# This model will be used for further analysis
h <- lapply(rfe_res, FUN = rfe_results, train = train)

##########################################################################################
# Model Development
    # Model Training and Validation 
##########################################################################################

# GLM MODEL
# From Dummy variables to original names 
vars <- predictors(rfe_res$glmProfile)
check <- names(train)
pred_glmProf <- getOriginalNames(vars = vars, check = check)

# Check whether variables names are correct. If not, manually correct them 
pred_glmProf[5] <- 'zwaarteBeroep'

# Training GLM model
fm <- as.formula(paste("MCID_MHQPijn", paste(pred_glmProf, collapse = " + "), sep = " ~ ")) 
set.seed(123)

na_rfe_control <- trainControl(method = 'repeatedcv',repeats = 5, number = 10,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE)

glm_no_sampling <- train(fm, 
                  train, 
                  method = 'glm', 
                  trControl = na_rfe_control, 
                  preProcess=c("scale","center","nzv", "knnImpute"), 
                  na.action = na.pass,
                  metric = "ROC")

# Evaluation of GLM model on  validation data set 
predicted.probability <- as.vector(predict.train(glm_no_sampling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_glm <- pROC::roc(MCID_MHQPijn ~ prob, data = val, AUC=TRUE)
plot(g_glm)

# GBM MODEL
# From Dummy variables to original names 
vars <- predictors(rfe_res$gbmProfile_up)
check <- names(train)
pred_gbmProf <- getOriginalNames(vars = vars, check = check)

# Check whether variables names are correct. If not, manually correct them
pred_gbmProf[24] <- 'zwaarteBeroep'

# Training GBM model
fm <- as.formula(paste("MCID_MHQPijn", paste(pred_gbmProf, collapse = " + "), sep = " ~ ")) 
set.seed(123)

na_rfe_control <- trainControl(method = 'repeatedcv',repeats = 5, number = 10,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE)

gbm_up_sampling <- train(fm, 
                          train, 
                          method = 'gbm', 
                          trControl = na_rfe_control, 
                          preProcess=c("scale","center","nzv", "knnImpute"), 
                          na.action = na.pass,
                          metric = "ROC")

# Evaluation of GBM model onvalidation data set 
predicted.probability <- as.vector(predict.train(gbm_up_sampling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_gbm <- pROC::roc(MCID_MHQPijn ~ prob, data = val, AUC=TRUE)
plot(g_gbm)

# RF MODEL
# From Dummy variables to original names 
vars <- predictors(rfe_res$rfProfile_down)
check <- names(train)
pred_rfProf <- getOriginalNames(vars = vars, check = check)

# Check whether variables names are correct. If not, manually correct them 
pred_rfProf [1] <- 'zwaarteBeroep'

# Training RF model
fm <- as.formula(paste("MCID_MHQPijn", paste(pred_rfProf, collapse = " + "), sep = " ~ ")) 
set.seed(123)

na_rfe_control <- trainControl(method = 'repeatedcv',repeats = 5, number = 10,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE)

rf_down_sampling <- train(fm, 
                          train, 
                          method = 'rf', 
                          trControl = na_rfe_control, 
                          preProcess=c("scale","center","nzv", "knnImpute"), 
                          na.action = na.pass,
                          metric = "ROC")

# Evaluation of RF model onvalidation data set 
predicted.probability <- as.vector(predict.train(rf_down_sampling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_rf <- pROC::roc(MCID_MHQPijn ~ prob, data = val, AUC=TRUE)
plot(g_rf)

# Store results
save(g_rf, file = "Results_rf_MCID_MHQ_Pijn_CMC1_070121.Rdata")
save(g_gbm, file = "Results_gbm_MCID_MHQ_Pijn_CMC1_070121.Rdata")
save(g_glm, file = "Results_glm_MCID_MHQ_Pijn_CMC1_070121.Rdata")

# Save the best performing model as the final model of the analysis
save(rf_down_sampling, file = "Final_model_pain_CMC1.Rdata")


##########################################################################################
# Model Development
    # RF Model Validation on Test set  
##########################################################################################

# Determining NPV and PPV for the RF model 
predicted.probability <- as.vector(predict.train(rf_down_sampling, test, type = "prob", na.action = na.pass)[,2])
test$prob <- predicted.probability
test$pred <- as.factor(ifelse(test$prob >= 0.719, "MCID_Gehaald", "MCID_NietGehaald"))

confusionMatrix(data = test$pred, reference = test$MCID_MHQPijn, positive = "MCID_Gehaald")

# Evaluation of GBM model on test set
predicted.probability <- as.vector(predict.train(rf_down_sampling, test, type = "prob", na.action = na.pass)[,2])
test$prob=predicted.probability
par(pty = 's')
g_rf_test <- pROC::roc(MCID_MHQPijn ~ prob, data = test, AUC=TRUE, plot = TRUE, legacy.axes = TRUE)
save(g_rf_test, file = "Results_rf_final_MCID_MHQ_Pijn_CMC1_070121.Rdata")
coords(g_rf_test, 'best')

##########################################################################################
# Model Development
    # Plotting calibration belts for visual inspection of calibration
##########################################################################################

val_2 <- val
test_2 <- test

# GLM model validation set -----

# Calibration belt for GLM model on VALIDATION data set
predicted.probability <-  as.vector(predict.train(glm_no_sampling, val, type = "prob", na.action = na.pass)[,2])
val_2$prob=predicted.probability
val_2$outcome <- ifelse(val_2$MCID_MHQPijn == "MCID_Gehaald", 1, 0)

tiff("GLM_Val.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = val_2$outcome, e = val_2$prob, devel = "external")
plot(cb, main = "Calibration GLM model without sampling",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

# Calibration belt for GBM model on VALIDATION data set
predicted.probability <-  as.vector(predict.train(gbm_up_sampling, val, type = "prob", na.action = na.pass)[,2])
val_2$prob=predicted.probability
val_2$outcome <- ifelse(val_2$MCID_MHQPijn == "MCID_Gehaald", 1, 0)

tiff("GBM_Val.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = val_2$outcome, e = val_2$prob, devel = "external")
plot(cb, main = "Calibration upsampled GBM model",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

# Calibration belt for RF model on VALIDATION data set
predicted.probability <-  as.vector(predict.train(rf_down_sampling, val, type = "prob", na.action = na.pass)[,2])
predicted.probability[predicted.probability == 1.000] <- 0.9999

val_2$prob=predicted.probability
val_2$outcome <- ifelse(val_2$MCID_MHQPijn == "MCID_Gehaald", 1, 0)

tiff("RF_Val.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = val_2$outcome, e = val_2$prob, devel = "external")
plot(cb, main = "Calibration RF model",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

# Calibration belt for RF model (our final model) on TEST data set
predicted.probability <-  as.vector(predict.train(rf_down_sampling, test, type = "prob", na.action = na.pass)[,2])
test_2$prob=predicted.probability
test_2$outcome <- ifelse(test_2$MCID_MHQPijn == "MCID_Gehaald", 1, 0)

tiff("RF_Test.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = test_2$outcome, e = test_2$prob, devel = "external")
plot(cb, main = "Calibration GLM model - Testset",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

