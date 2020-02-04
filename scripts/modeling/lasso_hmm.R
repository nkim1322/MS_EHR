library(expm)

expit <- function(v) {
  1/(1+exp(-v))
}

trainTransitionMatrix <- function(trainset){
  counts <- matrix(rep(0,4),2,2)
  rownames(counts) <- c('Start0','Start1')
  colnames(counts) <- c('End0','End1')
  
  currentPat <- trainset$PatientNum[1]
  currentPer <- trainset$Period[1]
  currentCC <- trainset$CC[1]
  for (i in 2:nrow(trainset)){
    if (trainset$PatientNum[i] != currentPat){
      counts[currentCC+1,currentCC+1] <- counts[currentCC+1,currentCC+1] + trainset$Period[i-1] - currentPer
      currentPat <- trainset$PatientNum[i]
      currentPer <- trainset$Period[i]
      currentCC <- trainset$CC[i]
    }
    else if (trainset$CC[i] != currentCC){
      counts[currentCC+1,currentCC+1] <- counts[currentCC+1,currentCC+1] + trainset$Period[i-1] - currentPer
      counts[currentCC+1,trainset$CC[i]+1] <- counts[currentCC+1,trainset$CC[i]+1] + 1
      currentPer <- trainset$Period[i]
      currentCC <- trainset$CC[i]
    }
  }
  
  transition <- matrix(nrow=2,ncol=2)
  rownames(transition) <- c('Start0','Start1')
  colnames(transition) <- c('End0','End1')
  transition['Start0','End0'] <- counts['Start0','End0']/(counts['Start0','End0']+counts['Start0','End1'])
  transition['Start0','End1'] <- counts['Start0','End1']/(counts['Start0','End0']+counts['Start0','End1'])
  transition['Start1','End0'] <- counts['Start1','End0']/(counts['Start1','End0']+counts['Start1','End1'])
  transition['Start1','End1'] <- counts['Start1','End1']/(counts['Start1','End0']+counts['Start1','End1'])
  
  return(list('counts'=counts,'tmat'=transition))
}

trainEmissionParams <- function(Xtrain){
  model <- lm(X~Y,Xtrain)
  return(list('mu0'=coef(model)[1], 'mu1'=coef(model)[1]+coef(model)[2], 'sigma'=summary(model)$sigma))
}

predictY <- function(Xtest,transition,emission,startProb){
  unlist(sapply(unique(Xtest$PatientNum), function(patient){
    keep <- which(Xtest$PatientNum==patient)
    Xi <- Xtest$X[keep]
    Ti <- Xtest$Period[keep]
    
    fwd <- matrix(0,2,length(Xi))
    f_prev <- startProb
    t_prev <- Ti[1]-1
    for (i in 1:length(Xi)){
      fwd[,i] <- (t(transition)%^%(Ti[i]-t_prev) %*% f_prev) * dnorm(Xi[i],c(emission$mu0,emission$mu1),emission$sigma)
      f_prev <- fwd[,i]
      t_prev <- Ti[i]
    }
    
    bkw <- matrix(0,2,length(Xi))
    b_next <- bkw[,length(Xi)] <- startProb
    t_next <- Ti[length(Xi)]
    if (length(Xi) > 1){
      for (i in (length(Xi)-1):1){
        bkw[,i] <- (transition%^%(t_next-Ti[i]) %*% b_next) * dnorm(Xi[i+1],c(emission$mu0,emission$mu1),emission$sigma)
        b_next <- bkw[,i]
        t_next <- Ti[i]
      }
    }
    
    post <- sapply(1:length(Xi), function(i){
      fwd[2,i]*bkw[2,i]/(fwd[,i]%*%bkw[,i])
    })

    post
  }))
}


hmm <- function(trainset,fitted_train,fitted_test){
  transition_matrix <- trainTransitionMatrix(trainset)$tmat
  emission_params <- trainEmissionParams(fitted_train)
  start_probs <- c(1-mean(trainset$CC),mean(trainset$CC))
  test_predictions <- predictY(fitted_test,transition_matrix,emission_params,start_probs)
  AUC <- auc(fitted_test$Y,test_predictions)
  return(list('transition_matrix'=transition_matrix, 'emission_params'=emission_params,
              'test_predictions'=test_predictions, 'AUC'=AUC))
}


runAll <- function(train,test){
  train <- train[order(train$Period),]; train <- train[order(train$PatientNum),]
  test <- test[order(test$Period),]; test <- test[order(test$PatientNum),]
  logreg <- glm(train$CC ~ as.matrix(train[,-c(1:9,206)]), family='binomial')
  coefs <- logreg$coefficients
  coefs[is.na(coefs)] <- 0
  logreg_prediction <- c(expit(as.matrix(cbind(1,test[,-c(1:9,206)]))%*%coefs))
  logreg_auc <- auc(test$CC,c(logreg_prediction))
  
  print('Logreg done')
  
  lasso <- cv.glmnet(as.matrix(train[,-c(1:9,206)]), train$CC, family='binomial', type.measure='auc')
  lasso_prediction <- predict(lasso,newx=as.matrix(test[,-c(1:9,206)]), s='lambda.1se', type='response')[,1]
  lasso_auc <- auc(test$CC,c(lasso_prediction))
  
  print('LASSO done')
  
  logreg_fitted_train <- predict(lasso,newx=as.matrix(train[,-c(1:9,206)]))[,1]
  logreg_fitted_test <- predict(lasso,newx=as.matrix(test[,-c(1:9,206)]))[,1]
  logreg_fitted_train_noGP <- as.data.frame(cbind(train$PatientNum,logreg_fitted_train,train$CC,train$Period))
  colnames(logreg_fitted_train_noGP) <- c('PatientNum','X','Y','Period')
  logreg_fitted_test_noGP <- as.data.frame(cbind(test$PatientNum,logreg_fitted_test,test$CC,test$Period))
  colnames(logreg_fitted_test_noGP) <- c('PatientNum','X','Y','Period')
  
  hmm12_result <- hmm(train,logreg_fitted_train_noGP,logreg_fitted_test_noGP)
  hmm12_auc <- hmm12_result$AUC
  hmm12_prediction <- hmm12_result$test_predictions
  
  print('HMM done')
  
  compiled_aucs <- c(logreg_auc,lasso_auc,hmm12_auc)
  names(compiled_aucs) <- c('LogReg','LASSO','HMM')
  print(compiled_aucs)
  output <- cbind(logreg_prediction,lasso_prediction,hmm12_prediction)
  colnames(output) <- c('LogReg', 'LASSO', 'HMM')
  
  return(list(output, compiled_aucs))
}


# Example use
# predictions <- runAll(CC_comb,CC_comb_val)



# Nicole's modification for Causal project --------------------------------

hmmSS <- function(trainset,fitted_train,fitted_test){
  transition_matrix <- trainTransitionMatrix(trainset)$tmat
  emission_params <- trainEmissionParams(fitted_train)
  start_probs <- c(1-mean(trainset$CC),mean(trainset$CC))
  test_predictions <- predictY(fitted_test,transition_matrix,emission_params,start_probs)
  return(list('transition_matrix'=transition_matrix, 'emission_params'=emission_params,
              'test_predictions'=test_predictions))
}

runAllSS <- function(train,test){
  train <- train[order(train$Period),]; train <- train[order(train$PatientNum),]
  test <- test[order(test$Period),]; test <- test[order(test$PatientNum),]
  # Model predictors, excludes ID/dates columns
  predictors <- readRDS("intermediate_data/model_predictors.rds")
  
  lasso <- cv.glmnet(x=as.matrix(train[,predictors]), y=train$CC, family='binomial', type.measure='auc')
  print('LASSO done')
  
  logreg_fitted_train <- predict(lasso,newx=as.matrix(train[,predictors]))[,1]
  logreg_fitted_test <- predict(lasso,newx=as.matrix(test[,predictors]))[,1]
  logreg_fitted_train_noGP <- as.data.frame(cbind(train$PatientNum,logreg_fitted_train,train$CC,train$Period))
  colnames(logreg_fitted_train_noGP) <- c('PatientNum','X','Y','Period')
  logreg_fitted_test_noGP <- as.data.frame(cbind(test$PatientNum,logreg_fitted_test,test$Period))
  colnames(logreg_fitted_test_noGP) <- c('PatientNum','X','Period')
  
  hmm12_result <- hmmSS(train,logreg_fitted_train_noGP,logreg_fitted_test_noGP)
  hmm12_prediction <- hmm12_result$test_predictions
  print('HMM done')
  
  return(hmm12_prediction)
}

