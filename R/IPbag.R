
#' Predicting outcomes using non-parametric Inverse-Probability bagging
#' @author Norbert Krautenbacher, Kevin Strauss, Maximilian Mandl, Christiane Fuchs
#' @description This method trains classifiers by correcting them for sample selection bias via non-parametric
#' inverse-probability bagging. This method fits classifiers from different resampled data whose observations
#' are increased per stratum to correct for the bias in the original sample. The so attained ensemble of predictors
#' is aggregated by averaging.
#' @references Krautenbacher, N., Theis, F. J., & Fuchs, C. (2017). Correcting Classifiers for Sample Selection Bias in Two-Phase Case-Control Studies. Computational and mathematical methods in medicine, 2017.
#' @param ... see the parameter ipOversampling() of package sambia.
#' @param learner a character indicating which classifier is used to train. Note: set learner='rangerTree'
#' if random forest should be applied as in Krautenbacher et al. (2017), i.e. the correction step is incorporated in the inherent random forest resampling procedure.
#' @param list.train.learner a list of parameters specific to the classifier that will be trained. Note that the
#' parameter 'data' need not to be provided in this list since the training data which the model will learn on is already attained by new sampled data produced by the method genSample().
#' @param list.predict.learner a list of parameters specifiying how to predict new data given the learned model.
#' (This learned model is uniquely determined by parameters 'learner' and 'list.train.learner').
#' @param n.bs number of bootstramp samples.
#' @examples
#' ## simulate data for a population
#' require(pROC)
#'
#' set.seed(1342334)
#' N = 100000
#' x1 <- rnorm(N, mean=0, sd=1) 
#' x2 <- rt(N, df=25)
#' x3 <- x1 + rnorm(N, mean=0, sd=.6)
#' x4 <- x2 + rnorm(N, mean=0, sd=1.3)
#' x5 <- rbinom(N, 1, prob=.6)
#' x6 <- rnorm(N, 0, sd = 1) # noise not known as variable
#' x7 <- x1*x5 # interaction
#' x <- cbind(x1, x2, x3, x4, x5, x6, x7)
#'
#' ## stratum variable (covariate)
#' xs <- c(rep(1,0.1*N), rep(0,(1-0.1)*N))
#'
#' ## effects
#' beta <- c(-1, 0.2, 0.4, 0.4, 0.5, 0.5, 0.6)
#' beta0 <- -2
#'
#' ## generate binary outcome
#' linpred.slopes <-  log(0.5)*xs + c(x %*% beta)
#' eta <-  beta0 + linpred.slopes
#'
#' p <- 1/(1+exp(-eta)) # this is the probability P(Y=1|X), we want the binary outcome however:
#' y<-rbinom(n=N, size=1, prob=p) #
#'
#' population <- data.frame(y,xs,x)
#'
#' #### draw "given" data set for training
#' sel.prob <- rep(1,N)
#' sel.prob[population$xs == 1] <- 9
#' sel.prob[population$y == 1] <- 8
#' sel.prob[population$y == 1 & population$xs == 1] <- 150
#' ind <- sample(1:N, 200, prob = sel.prob)
#'
#' data = population[ind, ]
#'
#' ## calculate weights from original numbers for xs and y
#' w.matrix <- table(population$y, population$xs)/table(data$y, data$xs)
#' w <- rep(NA, nrow(data))
#' w[data$y==0 & data$xs ==0] <- w.matrix[1,1]
#' w[data$y==1 & data$xs ==0] <- w.matrix[2,1]
#' w[data$y==0 & data$xs ==1] <- w.matrix[1,2]
#' w[data$y==1 & data$xs ==1] <- w.matrix[2,2]
#'
#' ### draw a test data set
#' newdata = population[sample(1:N, size=200 ), ]
#'
#' n.bs = 5
#' pred_nb <- sambia::IPbag(data = data, weights = w,
#'            learner='naiveBayes', list.train.learner = list(formula=formula(y~.)),
#'            list.predict.learner = list(newdata=newdata, type="raw"),
#'            n.bs = n.bs)
#' roc(newdata$y, pred_nb, direction="<")
#' @export
#' @import e1071 ranger


IPbag <- function(...,learner, list.train.learner, list.predict.learner, n.bs){
  # get all data types of list.predict.learner
  cats <- sapply(list.predict.learner,class)
  # get name of data frame
  #which(class.list == 'factor')
  data_frame_name <-  names(cats[match('data.frame',cats)])
  # get newdata of list
  newdata <- list.predict.learner[[data_frame_name]]
  # create an empty matrix containing n.bs boostrap predictions
  pred <- matrix(NA,nrow(newdata),n.bs)
  for(i in 1:n.bs){
    sample <- sambia::ipOversampling(...)
    ind.bs <- sample(1:nrow(sample), replace = TRUE, size = nrow(sample))
    sampled_data <- sample[ind.bs, ]
    list.train.learner$data <- sampled_data
    if(learner == 'naiveBayes'){
      fit <- do.call(learner, list.train.learner)
      list.predict.learner$object=fit
      pred[,i] <- do.call(predict, list.predict.learner)[,2]
    }else if(learner == 'rangerTree'){

      # modify list.train.learner
      learner1 = 'ranger'
      list.train.def.par <- c('write.forest','probability','num.trees','replace','sample.fraction')
      list.train.def.val <- c(TRUE,TRUE,1,FALSE,1)
      for(ind in 1:length(list.train.def.par)){
        list.train.learner[[list.train.def.par[ind]]] <- list.train.def.val[ind]
      }
      # apply model
      #list.train.learner$data[,'y'] <- as.factor(list.train.learner$data[,'y'])
      fit <- do.call(learner1, list.train.learner)
      list.predict.learner$object=fit
      pred[,i] <- do.call(predict, list.predict.learner)$predictions[,2]
    }else{
      fit <- do.call(learner, list.train.learner)
      list.predict.learner$object=fit
      pred[,i] <- do.call(predict, list.predict.learner)
    }
  }
  pred <- rowMeans(pred)
  return(pred)
}

