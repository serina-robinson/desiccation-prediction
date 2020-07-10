# Install packages
pacman::p_load("tidyverse", "caret", "rsample", "ranger", "e1071")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/desiccation-prediction/")

# Read in the dataset
rawdat <- read_csv("data/desiccationgenomes.csv")
colnames(rawdat)[1] <- "nams"

# Remove variables with nonzero variance (don't expect any - just a sanity check)
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] # none! 

# Check for duplicates (don't expect any - just a sanity check)
dat <- rawdat[!duplicated(rawdat),]

# Actually the random forest code throws an error with 0/1s so convert to a factor
dat$binary <- factor(ifelse(dat$binary == 1, "Tol", "Sen"))
table(dat$binary)

# Set random seed 
set.seed(1234)

# Split into test and training data 
# 80% training
# 20% test
dat_split <- rsample::initial_split(dat, strata = "binary", prop = 0.75)
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)
nrow(dat_test)
nrow(dat_train)/nrow(dat) # 82.3% of the data is in dat_train

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("nams", "binary")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nams", "binary")]

# Dependent variable
y_train <- dat_train$binary
y_test <- dat_test$binary
y_test

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nams)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nams)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$nams)

# Optional tuning different mtrys
mtrys <- c(round(log2(ncol(df_train)), 0), round(sqrt(ncol(df_train)), 0), round(ncol(df_train)/2, 0))
mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = c("gini", "extratrees"),
                       min.node.size = 1)

# Train a machine learning model
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, # this is how many folds
                           repeats = 1, # increase this to 3 when you run the code 
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

# Training set accuracy
getTrainPerf(rf) # Training set accuracy 
rf$finalModel$prediction.error # out-of-bag error

#Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                   surrogates = FALSE, 
                   competes = FALSE)
rf_imp
ggplot(rf_imp, top = 20) + 
  xlab("") +
  theme_classic()

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred

cm_rf <- confusionMatrix(rf_pred, y_test)
cm_rf

# ROC curve
rf$pred$obs
rf$pred$pred


rf_roc <- pROC::roc(response = ifelse(rf$pred$obs == "Sen", 0, 1),
                   predictor = ifelse(rf$pred$pred == "Sen", 0, 1),
                   plot = TRUE)
rf_roc

# AUC
plot(rf_roc, type = "s", 
     col = "#529DCF", xaxs = "i", yaxs="i",
     print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.6)


# Write model to file
saveRDS(rf, "output/desiccation_model_20200708.rds")
# Read model back in
my_mod <- readRDS("output/desiccation_model_20200708.rds")
sort(desc(my_mod$finalModel$variable.importance))[1:10]

## Now try tuning a model using feature selection!
feats_keep <- varImp(rf, 30, scale = FALSE, 
                         surrogates = FALSE, 
                         competes = FALSE)



ggplot(feats_keep, top = 30) + 
  xlab("") +
  theme_classic()

# Select top 30 most important variables
which_keep <- feats_keep$importance %>%
  arrange(desc(.)) %>% 
  dplyr::slice(1:30) %>%
  rownames()
which_keep # These are the 30 variables to keep

# Independent variables
x_train <- dat_train[,colnames(dat_train) %in% which_keep] 
x_test <- dat_test[,colnames(dat_test) %in% which_keep]
dim(x_test) # I'm not sure why this is 25 instead of 30?? # Maia, any ideas

# Dependent variable
y_train <- dat_train$binary
y_test <- dat_test$binary
y_test

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nams)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nams)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$nams)


# Train a machine learning model
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10, # this is how many folds
                           repeats = 1, # increase this to 3 when you run the code 
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

# Training set accuracy
getTrainPerf(rf) # Training set accuracy 
rf$finalModel$prediction.error # out-of-bag error     

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred

cm_rf <- confusionMatrix(rf_pred, y_test)
cm_rf # Accuracy looks way lower...want to do this in a loop to test
