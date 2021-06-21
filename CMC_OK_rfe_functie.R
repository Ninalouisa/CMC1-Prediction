rm(list = ls())

library(caret)
library(stringdist)
library(gbm)
library(tidyverse)
setwd("/Volumes/Passport NL")
library(here)
library(RANN)
library(ggpubr)
library(ggplot2)
library(givitiR)

select <- dplyr::select
rename <- dplyr::rename

load("/Volumes/Passport NL/Datasets/CMC_OK_Functie_ds.Rdata")

# Load the specific functions written by Lisa Hoogendam
source(here("R scripts", "Functies_Lisa.R"))


# Preprocessing ----

# Verwijderen van uitkomst data
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  select(-c(Patient.traject.ID, Respondent.ID, delta_MHQFunctie, score_handfunctie_aangedane_hand_MHQ12m, behandeling, 
            ipqscore, contains('rounddescription'), contains('926'), Vingers, type_behandeling, 
            organization, trackname, Survey.ID, location, Invuldatum, mrsaAsielzoeker, drugs, contains('12m'), contains('uitslag'), contains('Gonio'))) %>%
  droplevels() %>% 
  #mutate(uitslag_MCP_AROMfl_aangedane_hand = ifelse(nchar(uitslag_MCP_AROMfl_aangedane_hand)>2, 
  #as.numeric(str_sub(uitslag_MCP_AROMfl_aangedane_hand, 
  #start = 2)), as.numeric(uitslag_MCP_AROMfl_aangedane_hand))) %>% 
  #mutate(uitslag_IP_AROMex_aangedane_hand = ifelse(nchar(uitslag_IP_AROMex_aangedane_hand)>2, 
  #   as.numeric(str_sub(uitslag_IP_AROMex_aangedane_hand, 
  #                      start = 2)), as.numeric(uitslag_IP_AROMex_aangedane_hand))) %>% 
  #mutate(uitslag_MIP_AROMfl_aangedane_hand = ifelse(nchar(uitslag_MIP_AROMfl_aangedane_hand)>2, 
  #   as.numeric(str_sub(uitslag_MIP_AROMfl_aangedane_hand, 
  #                      start = 2)), as.numeric(uitslag_MIP_AROMfl_aangedane_hand))) %>% 
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

# Kracht aangedane hand
## jamar punt2
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  mutate(jamar_punt2_links  = ifelse(jamar_punt2_links > 990, NA, jamar_punt2_links)) %>%
  mutate(jamar_punt2_rechts  = ifelse(jamar_punt2_rechts > 990, NA, jamar_punt2_rechts)) %>% 
  mutate(pinch_sleutel_links = ifelse(pinch_sleutel_links > 90, NA, pinch_sleutel_links)) %>%
  mutate(pinch_sleutel_rechts = ifelse(pinch_sleutel_rechts > 90, NA, pinch_sleutel_rechts)) 

CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  dplyr::mutate(jamar_punt2_aangedane_hand = 
                  ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", jamar_punt2_links, 
                         ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", jamar_punt2_rechts, NA))) 
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  dplyr::mutate(jamar_punt2_niet_aangedane_hand = 
                  ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", jamar_punt2_rechts, 
                         ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", jamar_punt2_links, NA))) 

## pinch_sleutel
CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_sleutel_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_sleutel_links, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_sleutel_rechts, NA)))

CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_sleutel_niet_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_sleutel_rechts, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_sleutel_links, NA)))

## pinch_punt3
CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_punt3_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_punt3_links, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_punt3_rechts, NA)))

CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_punt3_niet_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_punt3_rechts, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_punt3_links, NA)))

## pinch_punt2
CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_punt2_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_punt2_links, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_punt2_rechts, NA)))

CMC_OK_Functie_ds <- mutate(CMC_OK_Functie_ds, pinch_punt2_niet_aangedane_hand = 
                           ifelse(CMC_OK_Functie_ds$zijde == "Linker hand", pinch_punt2_rechts, 
                                  ifelse(CMC_OK_Functie_ds$zijde == "Rechter hand", pinch_punt2_links, NA)))

# Kracht niet aangedane hand + dubbele variabelen verwijderen
CMC_OK_Functie_ds <- CMC_OK_Functie_ds %>%
  select(-c(jamar_punt2_links:pinch_punt2_rechts, contains('niet'))) %>%
  droplevels()

# Random Train-Test-Validation Split
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

# The data (100%) is split into:
# training (~60%)
# validation (~20%)
# testing (~20%)

# Dummy variables aanmaken----
fit <- MCID_MHQFunctie ~.
options(na.action='na.pass')
dummies <- model.matrix(fit, data = train) %>%  
  as.data.frame()
y <- train$MCID_MHQFunctie
dummies <- cbind(dummies, y)
names(dummies)[ncol(dummies)] <- "MCID_MHQFunctie"

# Centering, scaling en imputation-----
preProcValues <- preProcess(as.data.frame(dummies), method = c("center", "scale", "nzv", "knnImpute")) # Geeft error als df nog een tbl_df is

trainTransformed <- predict(preProcValues, dummies) %>% 
  select(MCID_MHQFunctie, everything()) %>% 
  select_if(~ !any(is.na(.)))
x_frame <- trainTransformed %>% 
  select(-MCID_MHQFunctie)
y <- trainTransformed$MCID_MHQFunctie
trainTransformed$MCID_MHQFunctie <- factor(trainTransformed$MCID_MHQFunctie)
sapply(trainTransformed, function(x) class(x))

#########################################################################

# RFE model zonder sampling ----
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
plot(glmProfile, type = c("g", "o"),  main = "RFE without Sampling")
predictors(glmProfile)

df <- test  %>%
  count(MCID_MHQFunctie) %>%
  mutate(Percent = n / sum(n)*100)
df %>%
  ggplot(aes(x = MCID_MHQFunctie, y = Percent, fill = MCID_MHQFunctie))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste(round(Percent),"%"), y = Percent),
            position = position_stack(vjust = 0.5)) +
  scale_x_discrete(labels=c("MCID not reached", "MCID Reached")) +
  labs(x = "", y = "Percentage", title = "Percentage MCID not reached vs. percentage MCID reached", subtitle = "12 months post-surgery")+
  theme(plot.title = element_text(hjust = 0.5, size = 26), plot.subtitle = element_text(hjust = 0.5, size = 18), axis.title.y = element_text(size = 50) ) +
  theme_minimal() +
  scale_fill_brewer(palette="Blues")

# RFE model met downsampling ----
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
plot(glmProfile_down, type = c("g", "o"),  main = "RFE with Downsampling")
predictors(glmProfile_down)

# RFE model met upsampling ----
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
plot(glmProfile_up, type = c("g", "o"),  main = "RFE with Upsampling")
predictors(glmProfile_up)

#########################################################################

# GBM met Downsampling ----
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

# GBM met Upsampling ----
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


# GBM zonder Sampling ----
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

##################################################################################

# RF met Downsampling
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

# RF met Upsampling 
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

# RF zonder Sampling
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

# Store results -------
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
load("/Volumes/Passport NL/Results_rfe_func_MCID_MHQ_CMC_231220.Rdata")
##############################################################################

# Checking performance of RFE models ----
h <- lapply(rfe_res_func, FUN = rfe_results, train = train)

# Van dummy vars naar originele namen ----
vars <- predictors(rfe_res_func$glmProfile)
check <- names(train)
pred_glmProf <- getOriginalNames(vars = vars, check = check)
vars # Checken of de namen goed zijn omgezet en indien nodig handmatig wijzigen
pred_glmProf  <- pred_glmProf[2:79]

# Formule voor in train functie maken GLM no sampling----
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

#evaluatie op validatieset
predicted.probability <- as.vector(predict.train(glm_no_sampiling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_glm <- pROC::roc(MCID_MHQFunctie ~ prob, data = val, AUC=TRUE)
plot(g_glm)

# Van dummy vars naar originele namen ----
vars <- predictors(rfe_res_func$rfProfile)
check <- names(train)
pred_rfProf <- getOriginalNames(vars = vars, check = check)
vars # Checken of de namen goed zijn omgezet en indien nodig handmatig wijzigen
pred_rfProf

# Formule voor in train functie maken RF no sampling----
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

#evaluatie op validatieset
predicted.probability <- as.vector(predict.train(rf_no_sampling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_rf <- pROC::roc(MCID_MHQFunctie ~ prob, data = val, AUC=TRUE)
plot(g_rf)

# Van dummy vars naar originele namen ----
vars <- predictors(rfe_res_func$gbmProfile)
check <- names(train)
pred_gbmProf <- getOriginalNames(vars = vars, check = check)
vars # Checken of de namen goed zijn omgezet en indien nodig handmatig wijzigen
pred_gbmProf

# Formule voor in train functie maken GBM no sampling----
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

#evaluatie op validatieset
predicted.probability <- as.vector(predict.train(gbm_no_sampling, val, type = "prob", na.action = na.pass)[,2])
val$prob=predicted.probability
g_gbm <- pROC::roc(MCID_MHQFunctie ~ prob, data = val, AUC=TRUE)
plot(g_gbm)


save(g_rf, file = "Results_rf_MCID_MHQ_Functie_CMC1_070121.Rdata")
save(g_gbm, file = "Results_gbm_MCID_MHQ_Functie_CMC1_070121.Rdata")
save(g_glm, file = "Results_glm_MCID_MHQ_Functie_CMC1_070121.Rdata")

save(gbm_no_sampling, file = "Final_model_function_CMC1_100321.Rdata")
load("Final_model_function_CMC1_100321.Rdata")

#NPV/PPV GBM Model Test set
predicted.probability <- as.vector(predict.train(gbm_no_sampling, test, type = "prob", na.action = na.pass)[,2])
test$prob <- predicted.probability
test$pred <- as.factor(ifelse(test$prob >= 0.622842, "MCID_Gehaald", "MCID_NietGehaald"))

confusionMatrix(data = test$pred, reference = test$MCID_MHQFunctie, positive = "MCID_Gehaald")

# Evaluatie op testset GBM
predicted.probability <- as.vector(predict.train(gbm_no_sampling, test, type = "prob", na.action = na.pass)[,2])
test$prob=predicted.probability
par(pty = 's')
g_gbm_test <- pROC::roc(MCID_MHQFunctie ~ prob, data = test, AUC=TRUE, plot = TRUE, legacy.axes = TRUE)
save(g_gbm_test, file = "Results_gbm_final_MCID_MHQ_Functie_CMC1_070121.Rdata")
library(pROC)
coords(g_gbm_test, 'best')

# Calibration belts and application of function
library(givitiR)
# GLM model validationset -----
val_2 <- val

predicted.probability <-  as.vector(predict.train(glm_no_sampiling, val, type = "prob", na.action = na.pass)[,2])
val_2$prob=predicted.probability
val_2$outcome <- ifelse(val_2$MCID_MHQFunctie == "MCID_Gehaald", 1, 0)

tiff("GLM_Function_Val.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = val_2$outcome, e = val_2$prob, devel = "external")
plot(cb, main = "Calibration GLM model without sampling",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

# GBM model validation set
val_2 <- val
predicted.probability <-  as.vector(predict.train(gbm_no_sampling, val, type = "prob", na.action = na.pass)[,2])
val_2$prob=predicted.probability
val_2$outcome <- ifelse(val_2$MCID_MHQFunctie == "MCID_Gehaald", 1, 0)

tiff("GBM_Function_Val.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = val_2$outcome, e = val_2$prob, devel = "external")
plot(cb, main = "Calibration GBM model without sampling",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()

# RF model validation set
val_2 <- val
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

# GBM model test set
test_2 <- test
predicted.probability <-  as.vector(predict.train(gbm_no_sampling, test, type = "prob", na.action = na.pass)[,2])
test_2$prob=predicted.probability
test_2$outcome <- ifelse(test_2$MCID_MHQFunctie == "MCID_Gehaald", 1, 0)

tiff("GBM_Function_test.tiff", units = "in", width = 5, height = 5, res = 300)
cb <- givitiCalibrationBelt(o = test_2$outcome, e = test_2$prob, devel = "external")
plot(cb, main = "Calibration GBM model without sampling",
     xlab =  "Predicted probability",
     ylab = "Observed probability of reaching MCID")
dev.off()
