# This script defines functions for use with the SDM. The majority of the 
# functions are adapted from previously existing code to more specifically fit 
# our purposes

# Modified from https://github.com/bobmuscarella/ENMeval/blob/master/R/calc.aicc.R
calc.AICc <- function(nparam, occ, model, predictors, predictive.maps = NULL, multiple_models = F) {
  if(is.null(predictive.maps) & !is.null(model) & !is.null(predictors)){
    predictive.maps = dismo::predict(predictors, model = model, type = "exponential")
  }
  AIC.valid <- nparam < nrow(occ)
  if (nlayers(predictive.maps) == 0) {
    res <- data.frame(cbind(AICc=NA, delta.AICc=NA, w.AIC=NA, parameters=nparam))
    warning("Cannot calculate AICc when rasterPreds = FALSE... returning NA's.")
  } else {
    vals <- extract(predictive.maps, occ)
    probsum <- cellStats(predictive.maps, sum)
    # The log-likelihood was incorrectly calculated (see next line) in ENMeval v.1.0.0 when working with >1 model at once.
    #   LL <- colSums(log(vals/probsum), na.rm=T)
    # The corrected calculation (since v.0.1.1) is:
    LL <- colSums(log(t(t(vals)/probsum)), na.rm=T)
    AICc <- (2*nparam - 2*LL) + (2*(nparam)*(nparam+1)/(nrow(occ)-nparam-1))
    AICc[AIC.valid==FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if(sum(is.na(AICc))==length(AICc)){
      warning("AICc not valid... returning NA's.")
      res <- data.frame(cbind(AICc, delta.AICc=NA, w.AIC=NA, parameters=nparam))
    } else {
      delta.AICc <- (AICc - min(AICc, na.rm=TRUE))
      w.AIC <- (exp(-0.5*delta.AICc))/(sum(exp(-0.5*delta.AICc), na.rm=TRUE))
      res <- data.frame(AICc, delta.AICc, w.AIC, parameters=nparam)
      rownames(res) <- NULL
    }    
  }
  rownames(res) <- NULL
  return(res)
}
  
  
# Code from https://rdrr.io/github/adamlilith/enmSdm/src/R/trainMaxNet.r
trainMaxNet <- function(
  data,
  resp=names(data)[1],
  preds=names(data)[2:ncol(data)],
  regMult=c(seq(0.5, 5, by=0.5), 6:10, 12.5, 15, 17.5, 20),
  classes='default',
  testClasses=TRUE,
  out='model',
  anyway=TRUE,
  verbose=FALSE,
  ...
) {
  
  ###########
  ## setup ##
  ###########
  
  # response and predictors
  if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
  if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]
  
  # get response and predictors
  presentBg <- data[ , resp]
  data <- data[ , preds, drop=FALSE]
  
  ### get combinations of features to test for each regularization multiplier
  
  if (classes == 'default') {
    classesToTest <- if (ncol(data) > 1) {
      c('l', 'p', 'q', 'h')
    } else {
      c('l', 'q', 'h')
    }
  } else {
    classesToTest <- rep(NA, nchar(classes))
    for (i in 1:nchar(classes)) classesToTest[i] <- substr(classes, i, i)
  }
  
  # create df of 1/0 to indicate each combination of classes to test
  if (testClasses) {
    classGrid <- expand.grid(rep(list(c(1, 0)), length(classesToTest)))
    classGrid <- classGrid[-which(rowSums(classGrid) == 0), ]
    classGrid <- as.data.frame(classGrid)
  } else {
    classGrid <- data.frame(matrix(rep(1, length(classesToTest)), nrow=1))
  }
  
  names(classGrid) <- classesToTest
  if (any(classGrid$l == 0)) classGrid <- classGrid[-which(classGrid$l == 0), ]
  
  ### collate presences and BG sites
  presences <- data[which(presentBg == 1), ]
  if (class(presences) != 'data.frame') presences <- as.data.frame(presences)
  names(presences) <- names(data) # names of data to which to predict
  
  bg <- data[which(presentBg == 0), ]
  if (class(bg) != 'data.frame') bg <- as.data.frame(bg)
  names(bg) <- names(data)
  
  ##########
  ## MAIN ##
  ##########
  
  tuning <- data.frame()
  
  # for each regularization multiplier
  for (thisRegMult in regMult) {
    
    if (verbose) print(paste0('Calculating AICc for multipler ', thisRegMult, ' with features:'))
    
    # for each combination of class features
    for (countCombo in 1:nrow(classGrid)) {
      
      if (verbose) print(paste0(classesToTest[c(classGrid[countCombo, ]) == 1]))
      
      theseClasses <- paste(classesToTest[as.logical(unlist(classGrid[countCombo, ]))], collapse='')
      
      # add dummy column if doing univariate model to avoid error in maxnet.default.regularizationMOD
      if (ncol(data) == 1 & theseClasses == 'l') {
        thisData <- data
        thisPresences <- presences
        thisBg <- bg
        
        thisData$DUMMY <- rep(1, nrow(thisData))
        thisPresences$DUMMY <- rep(1, nrow(presences))
        thisBg$DUMMY <- rep(1, nrow(bg))
        
      } else {
        thisData <- data
        thisPresences <- presences
        thisBg <- bg
        
      }
      
      # train model
      model <- maxnet::maxnet(
        p=as.vector(presentBg),
        data=thisData,
        f=maxnet::maxnet.formula(p=as.vector(presentBg), data=thisData, classes=theseClasses),
        regfun=maxnet::maxnet.default.regularization,
        regmult=thisRegMult,
        ...
      )
      
      # predict presences
      predPres <- stats::predict(
        object=model,
        newdata=presences,
        type='exponential',
        ...
      )
      
      # predict to background
      predBg <- stats::predict(
        object=model,
        newdata=bg,
        type='exponential',
        ...
      )
      
      rawSum <- sum(c(predPres, predBg), na.rm=TRUE)
      
      ## calculate log likelihood
      ll <- sum(log(predPres / rawSum), na.rm=TRUE)
      
      ## number of parameters
      K <- length(model$betas)
      
      # AICc
      AICc <- -2 * ll + 2 * K + (2 * K * (K + 1)) / (sum(presentBg) - K - 1)
      
      # remember
      thisAicFrame <- data.frame(
        regMult=thisRegMult,
        n=sum(presentBg),
        classes=theseClasses,
        logLik=ll,
        K=K,
        AICc=AICc
      )
      
      tuning <- rbind(tuning, thisAicFrame)
      
    } # next combination of class features
    
    if (verbose) print(paste0(''))
    
  } # next reg mult
  
  # remove models with more parameters than data points that have more than 0 parameters
  tuning <- tuning[which(tuning$n >= tuning$K & tuning$K > 0), ]
  
  # re-order frame so one with lowest AICc, number of parameters, and reg mult are used (in that order, used to break ties)
  if (nrow(tuning) > 0) {
    
    tuning <- tuning[order(tuning$regMult, decreasing=TRUE), ]
    tuning <- tuning[order(tuning$AICc, tuning$K, tuning$regMult), ]
    
    tuning$deltaAICc <- tuning$AICc - min(tuning$AICc)
    tuning$relLike <- exp(-0.5 * tuning$deltaAICc)
    tuning$aicWeight <- tuning$relLike / sum(tuning$relLike)
    
  }
  
  if (verbose) {
    
    print(paste0(''))
    print(tuning)
    print(paste0(''))
    
  }
  
  # if user wants best model returned
  if ('model' %in% out) {
    
    # train model
    if (nrow(tuning) > 0) {
      
      if (anyway) {
        
        # add dummy column if doing univariate model to avoid error in maxnet.default.regularizationMOD
        if (ncol(data) == 1 & theseClasses == 'l') {
          
          thisData <- data
          thisPresences <- presences
          thisBg <- bg
          
          thisData$DUMMY <- rep(1, nrow(thisData))
          thisPresences$DUMMY <- rep(1, nrow(presences))
          thisBg$DUMMY <- rep(1, nrow(bg))
          
        } else {
          
          thisData <- data
          thisPresences <- presences
          thisBg <- bg
          
        }
        
        model <- maxnet::maxnet(
          p=as.vector(presentBg),
          data=thisData,
          f=maxnet::maxnet.formula(p=as.vector(presentBg), data=thisData, classes=if (nrow(tuning) == 0) { 1 } else { tuning$classes[1] }),
          regfun=maxnet::maxnet.default.regularization,
          regmult=if (nrow(tuning) == 0) { 1 } else { tuning$regMult[1] },
          ...
        )
        
        if (nrow(tuning) == 0) warning('No models had fewer coefficients than predictors. No model returned.', immediate.=TRUE)
        
      } else {
        
        warning('No models had fewer coefficients than predictors. No model returned.', immediate.=TRUE)
        model <- 'No MAXENT model had number of parameters < number of training presences.'
        
      }
      
    }
    
  }
  
  # return stuff
  if ('model' %in% out & !('tuning' %in% out)) {
    model
  } else if (!('model' %in% out) & 'tuning' %in% out) {
    tuning
  } else {
    list(tuning=tuning, model=model)
  }
  
}


gain = function(maxnet_model, p, data, occ.data = data[p == "1",,drop = F], 
                bg.data = data[p == "0",, drop = F],
                gaintype = c("glmnet", "maxent")){
  penalty <- maxnet_model$lambda[200]*abs(maxnet_model$beta[,200])%*%maxnet_model$penalty.factor
  # negative binomial ll from https://web.stanford.edu/~hastie/Papers/Glmnet_Vignette.pdf
  if(gaintype == "glmnet") {
    sum1 = p%*%predAll/length(p) - mean(log(1+exp(predAll)))
    # predict to all
    predAll = stats::predict(
      object=maxnet_model,
      newdata=data,
      type='link'
    )
    return(sum1 - penalty)
  }
  
  # Alternatively, the gain function from Merow 2013
  if(gaintype == "maxent"){
    predPres = stats::predict(
      object=maxnet_model,
      newdata=occ.data,
      type='link'
    )
    predBg = stats::predict(
      object=maxnet_model,
      newdata=bg.data,
      type='exponential'
    )
    return(mean(predPres) - log(mean(predBg)) - penalty)
  }
}

# Attempted to define a variable jackknife function to calculate variable 
# importance as is calculated by MaxEnt usually, but it didn't quite work
# variable_jackknife = function(p, data, gaintype = c("glmnet", "maxent")){
#   gains = matrix(NA, nrow = ncol(data), ncol = 2)
#   for (i in 1:ncol(data)){
#     mod_i <- maxnet(p = p, data = data.frame(data[,i, drop = F], DUMMY = rep(1, nrow(data))))
#     mod_noi <- maxnet(p = p, data = data[,-i, drop = F])
#     gains[i,] = c(gain(mod_i, p, data[,i, drop = F], 
#                               gaintype = gaintype),
#                          gain(mod_noi, p, data[,-i, drop = F], gaintype = gaintype ))
#   }
#   
#   gains <- cbind(colnames(data), gains) %>% as.data.frame 
#   colnames(gains) = c("variable","only", "without")
#   rownames(gains) = seq(1:nrow(gains))
#   
#   gains %<>%reshape(direction = "long", varying = 2:3, v.names = "gain", 
#                             timevar = "model_type", times = c("only", "without")) %>%
#     mutate_if(is.factor, as.character)
#   
#   gains$variable = as.character(gains$variable)
#   mod_all <- maxnet(p = p, data = data)
#   gains %<>% rbind(c(variable = "full_model", model_type = "full", 
#                      gain = gain(mod_all, p, data, gaintype = gaintype), 
#                      id = 0))
#   return(gains)
# }



