# This function does feature selection using the caret and mbench packages
# https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
library(mlbench)
library(caret)
library(rmcfs)
# load the data

# function to tune hyperparameter and calculate accuracy for different feature selections. stolen from here:
# https://www.kaggle.com/reisel/how-to-handle-correlated-features

#this first function helps us to find the optimal number of features
calculateAccFeatureSelect <- function(featureSelectName, predictorName, nFeature, hp){
  # featureSelectName is the type of feature estimation
  # predictorName is either Random forest or support vector machine
  
  # get order of features
  if (featureSelectName == "GE"){
    featureOrder <- featureGE
    useFeatureOrder <- TRUE
  } else if (featureSelectName == "LG"){
    featureOrder <- featureLG
    useFeatureOrder <- TRUE
  } else if (featureSelectName == "RFE"){
    featureOrder <- featureRFE
    useFeatureOrder <- TRUE
  } else if (featureSelectName == "PCA"){
    useFeatureOrder <- FALSE
  } else {
    stop("unknown feature selection")
  }
  
  # number of hyperparameter
  nHp <- length(hp)
  
  nNFeature <- length(nFeature)
  accuracyNFeature <- data.frame(nFeature = nFeature)
  accuracyNFeature$predictor <- predictorName
  accuracyNFeature$featureSelect <- featureSelectName
  accuracyNFeature$hyperparameter <- NA
  accuracyNFeature$accVal <- NA
  accuracyNFeature$accTest <- NA
  
  # iterate through number of features
  for (iFeature in 1:nNFeature){
    if (featureSelectName %in% c("GE", "LG", "RFE")){
      # names of selected features
      featureSelect <- featureOrder[1:nFeature[iFeature]]
    }
    
    # average accuracy of hyperparameter
    accHp <- rep(NA, nHp)    
    
    # iterate through hyper parameter
    for (iHp in 1:nHp){
      # skip this iteration if random forest is chosen and number of features exceeds hyperparameter
      if ((predictorName == "RF") && (nFeature[iFeature] < hp[iHp])) {
        accHp[iHp] <- -1
        next
      }
      
      # accuracy of set
      accFold <- rep(NA, nFold)
      
      # iterate through folds
      for (iFold in 1:nFold){
        XTrainSelect <- XTrainFold[[iFold]]
        XValSelect <- XValFold[[iFold]]                
        yTrainSelect <- yTrainFold[[iFold]]
        yValSelect <- yValFold[[iFold]]
        if (useFeatureOrder){
          # prepare training and validation data
          XTrainSelect <- XTrainSelect[, colnames(XTrainSelect) %in% featureSelect]
          XValSelect <- XValSelect[, colnames(XValSelect) %in% featureSelect]
        } else {
          # transform data with PCA
          pcaModel <- prcomp(XTrainSelect, scale. = TRUE)
          XTrainSelect <- pcaModel$x[, 1:nFeature[iFeature]]
          XValFull <- predict(pcaModel, XValSelect)
          XValSelect <- XValFull[, 1:nFeature[iFeature]]
        }
        
        # fit model
        if (predictorName == "SVM"){
          model <- svm(XTrainSelect, as.factor(yTrainSelect), scale = TRUE, kernel = "radial", gamma = hp[iHp])
        } else if (predictorName == "RF") {
          model <- randomForest(XTrainSelect, as.factor(yTrainSelect), ntree = 100, mtry = hp[iHp])
        } else {
          stop("Unknown predictor")
        }
        # predict survival
        predVal <- as.numeric(predict(model, XValSelect)) - 1
        # calculate accuracy
        accFold[iFold] <- mean(yValSelect == predVal)
      }
      
      # calcualte average accuracy of hyperparameter
      accHp[iHp] <- mean(accFold)
    }
    
    # select hyperparameter with highest accuracy
    iHpBest <- which.max(accHp)
    
    # get best hyperparameter
    hpBest <- hp[iHpBest]
    
    # insert validation accuracy to data frame
    accuracyNFeature$accVal[iFeature] <- accHp[iHpBest]
    
    # combine training and validation data and select features
    XTrainVal <- rbind(XTrainFold[[1]], XValFold[[1]])
    yTrainVal <- c(yTrainFold[[1]], yValFold[[1]])
    
    if (useFeatureOrder){
      XTrainVal <- XTrainVal[, colnames(XTrainVal) %in% featureSelect]
      XTestSelect <- XTest[, colnames(XTest) %in% featureSelect]
    } else {
      # transform data with PCA
      pcaModel <- prcomp(XTrainVal, scale. = TRUE)
      XTrainVal <- pcaModel$x[, 1:nFeature[iFeature]]
      XTestFull <- predict(pcaModel, XTest)
      XTestSelect <- XTestFull[, 1:nFeature[iFeature]]
    }
    
    # train model with best hyperparameter on train and validation data
    if (predictorName == "SVM"){
      model <- svm(XTrainVal, as.factor(yTrainVal), scale = TRUE, kernel = "radial", gamma = hpBest)
    } else if (predictorName == "RF"){
      # check if hyper parameter is valid
      if (hpBest > nFeature[iFeature]) {
        # set hyperparameter to maximum value
        hpBest <- nFeature[iFeature]
      }
      model <- randomForest(XTrainVal, as.factor(yTrainVal), ntree = 100, mtry = hpBest)
    }
    
    # insert best hyperparameter
    accuracyNFeature$hyperparameter[iFeature] <- hpBest
    
    # predict survival on test set
    predTest <- as.numeric(predict(model, XTestSelect)) - 1
    
    # calculate accuracy
    accTest <- mean(yTest == predTest)
    
    # insert accuracy to data frame
    accuracyNFeature$accTest[iFeature] <- accTest
  }
  accuracyNFeature
}




feature_selection <- function(data_frame_predictors, resp_vector, corr.cutoff = 0.75){
  #data_frame_predictors is a dataframe in which columns are the variables and rows are observations
  #corr.cutoff is the cutoff beyond which should include
  set.seed(4)
  # define the control using a random forest selection function
  control <- rfeControl(functions=rfFuncs, method="cv", number=10)
  # run the RFE algorithm
  results <- rfe(data_frame_predictors, response_vector, sizes=c(1:dim(data_frame_predictors)[1]), rfeControl=control)
  # summarize the results
  # list the chosen features
  #predictors(results)
  # plot the results
  plot(results, type=c("g", "o"))
  return(results)
}

#Notes: rfe is a simple backwards selection (recursive feature elimination algorith)

# problems with feature selection: https://stats.stackexchange.com/questions/27750/feature-selection-and-cross-validation