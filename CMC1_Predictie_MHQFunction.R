# Development, training, and validation of prediction model for MHQ Function Score
# after surgery for CMC1 OA
#
##########################################################################################
# Code written by N.L. Loos & L. Hoogendam for the Hand Wrist Study Group, 
# Erasmus Medical Center, Rotterdam
##########################################################################################
# 
# This script prepares data for analysis. It splits the data into a training, validation and test
# data set. And fits several models to this data: Generalized Linear Model (GLM), Gradient Boosting 
# Machine (GBM) and Random Forest (RF). The models are trained and validated on separate data sets. 
# Hopefully this script might provide those who are interested in developing prediction models with 
# some tips on how to perform this on their own data. 
#
# Note: this script is specific to our data. For this script to be used on different data set some 
# adaptation will be needed.

rm(list = ls())

# Load packages 
library(caret)
library(stringdist)
library(gbm)
library(tidyverse)
library(here)
library(RANN)
library(ggpubr)
library(ggplot2)
library(givitiR)
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
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  select(-c(Patient.traject.ID, Respondent.ID, delta_MHQFunctie, score_handfunctie_aangedane_hand_MHQ12m, behandeling, 
            ipqscore, contains('rounddescription'), contains('926'), Vingers, type_behandeling, 
            organization, trackname, Survey.ID, location, Invuldatum, mrsaAsielzoeker, drugs, contains('12m'), contains('uitslag'), contains('Gonio'))) %>%
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


CMC_OK_Functie_ds <- rename(CMC_OK_Functie_ds, c("zijde" = "zijde_EQ5D"))
CMC_OK_Functie_ds$Geslacht <- as.factor(CMC_OK_Functie_ds$Geslacht)
CMC_OK_Functie_ds$side <- as.factor(CMC_OK_Functie_ds$side)
CMC_OK_Functie_ds$Dominant_behandeld <- factor(CMC_OK_Functie_ds$Dominant_behandeld, levels = c("Nee", "Ja"))

##########################################################################################

# Removing impossible hands strength measurements
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  mutate(jamar_punt2_links  = ifelse(jamar_punt2_links > 990, NA, jamar_punt2_links)) %>%
  mutate(jamar_punt2_rechts  = ifelse(jamar_punt2_rechts > 990, NA, jamar_punt2_rechts)) %>% 
  mutate(pinch_sleutel_links = ifelse(pinch_sleutel_links > 90, NA, pinch_sleutel_links)) %>%
  mutate(pinch_sleutel_rechts = ifelse(pinch_sleutel_rechts > 90, NA, pinch_sleutel_rechts)) 

# Creating new variables with strength measurements for the affected side
# Jamar punt2
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  dplyr::mutate(jamar_punt2_aangedane_hand = 
                  ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", jamar_punt2_links, 
                         ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", jamar_punt2_rechts, NA))) 
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  dplyr::mutate(jamar_punt2_niet_aangedane_hand = 
                  ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", jamar_punt2_rechts, 
                         ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", jamar_punt2_links, NA))) 

# Pinch_sleutel
CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_sleutel_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_sleutel_links, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_sleutel_rechts, NA)))

CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_sleutel_niet_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_sleutel_rechts, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_sleutel_links, NA)))

# Pinch_punt3
CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_punt3_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_punt3_links, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_punt3_rechts, NA)))

CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_punt3_niet_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_punt3_rechts, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_punt3_links, NA)))

# Pinch_punt2
CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_punt2_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_punt2_links, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_punt2_rechts, NA)))

CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_punt2_niet_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_punt2_rechts, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_punt2_links, NA)))

# Removing strength measurements from the not affected side, and old strength measurements
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
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
inTrain <- as.data.frame(createDataPartition(CMC_OK_Functie_ds$MCID_MHQFunctie, 
                                             p = 0.60, 
                                             list = FALSE, 
                                             times = 1))

train <- CMC_OK_Functie_ds[inTrain$Resample1, ]
test_plus_val <- CMC_OK_Functie_ds[-inTrain$Resample1, ]

set.seed(123)
in_test_plus_val <- as.data.frame(createDataPartition(test_plus_val$MCID_MHQFunctie,
                                                      p = 0.5,
                                                      list = FALSE,
                                                      times = 1))

test <- test_plus_val[-in_test_plus_val$Resample1, ]
val <- test_plus_val[in_test_plus_val$Resample1, ]

# Making Dummy variables
fit <- MCID_MHQFunctie ~.
options(na.action='na.pass')
dummies <- model.matrix(fit, data = train) %>%  
  as.data.frame()
y <- train$MCID_MHQFunctie
dummies <- cbind(dummies, y)
names(dummies)[ncol(dummies)] <- "MCID_MHQFunctie"

# Centering, scaling and imputation
preProcValues <- preProcess(as.data.frame(dummies), method = c("center", "scale", "nzv", "knnImpute")) # Geeft error als df nog een tbl_df is

trainTransformed <- predict(preProcValues, dummies) %>% 
  select(MCID_MHQFunctie, everything()) %>% 
  select_if(~ !any(is.na(.)))
x_frame <- trainTransformed %>% 
  select(-MCID_MHQFunctie)
y <- trainTransformed$MCID_MHQFunctie
trainTransformed$MCID_MHQFunctie <- factor(trainTransformed$MCID_MHQFunctie)
sapply(trainTransformed, function(x) class(x))

##########################################################################################
# Model Development
    # Recursive Feature Elimination (RFE)
##########################################################################################

# RFE for GLM model
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
plot(glmProfile, type = c("g", "o"),  main = "RFE for GLM model")
predictors(glmProfile)

##########################################################################################

# RFE for GBM model
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
plot(gbmProfile, type = c("g", "o"), main = "RFE for GBM model")
predictors(gbmProfile)

##########################################################################################

# RFE for RF model
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
plot(rfProfile, type = c("g", "o"), main  = "RFE for RF model")
predictors(rfProfile)

##########################################################################################

# Store results of RFE
rfe_res_func <- list(
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

save(rfe_res_func, file = "Results_rfe_func_MCID_MHQ_CMC_231220.Rdata")

# Checking performance of RFE models
h <- lapply(rfe_res_func, FUN = rfe_results, train = train)

##########################################################################################
# Model Development
    # Model Training and Validation 
##########################################################################################

# GLM MODEL
# From Dummy variables to original names 
vars <- predictors(rfe_res_func$glmProfile)
check <- names(train)
pred_glmProf <- getOriginalNames(vars = vars, check = check)

# Check whether variables names are correct. If not, manually correct them 

# Training GLM model
fm <- as.formula(paste("MCID_MHQFunctie", paste(pred_glmProf, collapse = " + "), sep = " ~ "))
set.seed(123)

na_rfe_control <- trainControl(method = 'repeatedcv',repeats = 5, number = 10,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE)

glm_no_sampiling <- train(fm,
                          train,
                          method = 'glm',
                          trControl = na_rfe_control,
                          preProcess=c("scale","center","nzv", "knnImpute"),
                          na.action = na.pass,
                          metric = "ROC")

# Evaluation of GLM model on  validation data set 
predicted.probability <- as.vector(predict.train(glm_no_sampiling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_glm <- pROC::roc(MCID_MHQFunctie ~ prob, data = val, AUC=TRUE)
plot(g_glm)

# RF MODEL
# From Dummy variables to original names 
vars <- predictors(rfe_res_func$rfProfile)
check <- names(train)
pred_rfProf <- getOriginalNames(vars = vars, check = check)

# Check whether variables names are correct. If not, manually correct them 

# Training RF model
fm <- as.formula(paste("MCID_MHQFunctie", paste(pred_rfProf, collapse = " + "), sep = " ~ "))
set.seed(123)

na_rfe_control <- trainControl(method = 'repeatedcv',repeats = 5, number = 10,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE)

rf_no_sampling <- train(fm,
                  train,
                  method = 'rf',
                  trControl = na_rfe_control,
                  preProcess=c("scale","center","nzv", "knnImpute"),
                  na.action = na.pass,
                  metric = "ROC")

# Evaluation of RF model on  validation data set 
predicted.probability <- as.vector(predict.train(rf_no_sampling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_rf <- pROC::roc(MCID_MHQFunctie ~ prob, data = val, AUC=TRUE)
plot(g_rf)

# GBM MODEL
# From Dummy variables to original names 
vars <- predictors(rfe_res_func$gbmProfile)
check <- names(train)
pred_gbmProf <- getOriginalNames(vars = vars, check = check)

# Check whether variables names are correct. If not, manually correct them 

# Training GBM model
fm <- as.formula(paste("MCID_MHQFunctie", paste(pred_gbmProf, collapse = " + "), sep = " ~ "))
set.seed(123)

na_rfe_control <- trainControl(method = 'repeatedcv',repeats = 5, number = 10,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE)

gbm_no_sampling <- train(fm,
                        train,
                        method = 'gbm',
                        trControl = na_rfe_control,
                        preProcess=c("scale","center","nzv", "knnImpute"),
                        na.action = na.pass,
                        metric = "ROC")

# Evaluation of GBM model on validation data set 
predicted.probability <- as.vector(predict.train(gbm_no_sampling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_gbm <- pROC::roc(MCID_MHQFunctie ~ prob, data = val, AUC=TRUE)
plot(g_gbm)

# Store results
save(g_rf, file = "Results_rf_MCID_MHQ_Functie_CMC1_070121.Rdata")
save(g_gbm, file = "Results_gbm_MCID_MHQ_Functie_CMC1_070121.Rdata")
save(g_glm, file = "Results_glm_MCID_MHQ_Functie_CMC1_070121.Rdata")

# Save the best performing model as the final model of the analysis
save(gbm_no_sampling, file = "Final_model_function_CMC1.Rdata")

##########################################################################################
# Model Development
    # GBM Model Validation on Test set  
##########################################################################################

# Determining NPV and PPV for the GBM model 
predicted.probability <- as.vector(predict.train(gbm_no_sampling, test, type = "prob", na.action = na.pass)[,2])
test$prob <- predicted.probability
test$pred <- as.factor(ifelse(test$prob >= 0.622842, "MCID_Gehaald", "MCID_NietGehaald"))

confusionMatrix(data = test$pred, reference = test$MCID_MHQFunctie, positive = "MCID_Gehaald")

# Evaluation of GBM model on test set
predicted.probability <- as.vector(predict.train(gbm_no_sampling, test, type = "prob", na.action = na.pass)[,2])
test$prob=predicted.probability
par(pty = 's')
g_gbm_test <- pROC::roc(MCID_MHQFunctie ~ prob, data = test, AUC=TRUE, plot = TRUE, legacy.axes = TRUE)
save(g_gbm_test, file = "Results_gbm_final_MCID_MHQ_Functie_CMC1_070121.Rdata")
coords(g_gbm_test, 'best')

##########################################################################################
# Model Development
    # Plotting calibration belts for visual inspection of calibration
##########################################################################################

val_2 <- val
test_2 <- test

# Calibration belt for GLM model on VALIDATION data set
predicted.probability <-  as.vector(predict.train(glm_no_sampiling, val, type = "prob", na.action = na.pass)[,2])
val_2$prob=predicted.probability
val_2$outcome <- ifelse(val_2$MCID_MHQFunctie == "MCID_Gehaald", 1, 0)

tiff("GLM_Function_Val.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = val_2$outcome, e = val_2$prob, devel = "external")
plot(cb, main = "Calibration GLM model without sampling",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

# Calibration belt for GBM model on VALIDATION data set
predicted.probability <-  as.vector(predict.train(gbm_no_sampling, val, type = "prob", na.action = na.pass)[,2])
val_2$prob=predicted.probability
val_2$outcome <- ifelse(val_2$MCID_MHQFunctie == "MCID_Gehaald", 1, 0)

tiff("GBM_Function_Val.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = val_2$outcome, e = val_2$prob, devel = "external")
plot(cb, main = "Calibration GBM model without sampling",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

# Calibration belt for RF model on VALIDATION data set
predicted.probability <-  as.vector(predict.train(rf_no_sampling, val, type = "prob", na.action = na.pass)[,2])
predicted.probability[predicted.probability == 1.000] <- 0.9999
predicted.probability[predicted.probability == 0.0000] <- 0.0001
val_2$prob=predicted.probability
val_2$outcome <- ifelse(val_2$MCID_MHQFunctie == "MCID_Gehaald", 1, 0)


tiff("RF_Function_Val.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = val_2$outcome, e =val_2$prob, devel = "external")
plot(cb, main = "Calibration RF model without sampling",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

# Calibration belt for GBM model (our final model) on TEST data set
predicted.probability <-  as.vector(predict.train(gbm_no_sampling, test, type = "prob", na.action = na.pass)[,2])
test_2$prob=predicted.probability
test_2$outcome <- ifelse(test_2$MCID_MHQFunctie == "MCID_Gehaald", 1, 0)

tiff("GBM_Function_test.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = test_2$outcome, e = test_2$prob, devel = "external")
plot(cb, main = "Calibration GBM model without sampling",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()
