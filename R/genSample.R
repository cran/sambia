#' Generate synthetic observations using inverse-probability weights
#' @author Norbert Krautenbacher, Kevin Strauss, Maximilian Mandl, Christiane Fuchs
#' @description
#' This method corrects a given data set for sample selection bias by
#' generating new observations via Stochastic inverse-probability
#' oversampling or parametric inverse-probability sampling using
#' inverse-probability weights and information on covariance structure of the given strata (Krautenbacher et al, 2017).
#' @references Krautenbacher, N., Theis, F. J., & Fuchs, C. (2017). Correcting Classifiers for Sample Selection Bias in Two-Phase Case-Control Studies. Computational and mathematical methods in medicine, 2017.
#' @param data a data frame containing the observations rowwise, along with their corresponding categorical strata feature.
#' @param strata.variables a character vector of the names determined by the categorical stratum features.
#' @param stratum a numerical vector of the length of the number of rows of the data specifying 
#' the stratum ID. Either 'strata.variables' or 'stratum' has to be provided.
#' This vector will not be included as a column in the resulting data set.
#' @param weights a numerical vector whose length must coincide with the number of the rows of data. The i-th value contains the inverse-probability e.g. determines how often the i-th observation of data shall be replicated.
#' @param distr character object that describes the distribution
#' @param type character which decides which method is used to correct a given data set for sample selection bias. Stochastic Inverse-Probabiltiy oversampling is applied if type = 'stochIP' or Parametric Inverse-Probability Bagging if type = 'parIP'.
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
#' #### draw "given" data set 
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
#' ## parametric IP bootstrap sample
#' sample1 <- sambia::genSample(data=data, strata.variables = c('y', 'xs'),
#'                           weights = w, type='parIP')
#' ## stochastic IP oversampling; treating 'y' and 'xs' as usual input variable
#' ## but using strata info unambiguously defined by the weights w                        
#' sample2 <- sambia::genSample(data=data,
#'                             weights = w, type='stochIP', stratum= round(w))
#' @return $data data frame containing synthetic data which is corrected for sample selection bias by generating new observations via Stochastic inverse-probability oversampling or parametric inverse-probability oversampling.
#' @return $orig.data original data frame which shall to corrected
#' @return $stratum vector containing the stratum of each observation
#' @return $method a character indicating which method was used. If method = 'stochIP' then Stochastic Inverse-Probabiltiy oversampling was used, 
#' if method = 'parIP' the Parametric Inverse-Probability sampling was used.
#' @return $strata.tbl a data frame containing all variables and their feature occurences
#' @return $N number of rows in data
#' @return $n number of rows in original data
#' @export
#' @import dplyr mvtnorm


genSample <- function (data, strata.variables = NULL, stratum = NULL, weights = rep(1,nrow(data)), distr = "mvnorm", type = c("parIP", "stochIP"))
{
  if((is.null(strata.variables) & is.null(stratum) ) | ((!is.null(strata.variables) & !is.null(stratum) ))) {
    stop("Either argument 'strata.variables' or 'stratum' has to be given!")
  }
  if (any(sapply(weights, is.na))) {
    stop("weights sholud not contain any NA")
  }
  if (nrow(data) != length(weights)) {
    stop("The number of observations must be equal to the number of weights")
  }
  if (!(is.null(stratum))) {
    if (nrow(data) != length(stratum)) {
      stop("Size of stratum and data must coincide!")
    }
  }
  if (!(type %in% c("parIP", "stochIP"))) {
    stop("type parameter must be eithter set to 'parIP' or 'stochIP'")
  }

  #### nk: following function not needed
  # strataTbl <- function(strata.features) {
  #   colnames <- colnames(strata.features)
  #   l <- list()
  #   for (name in colnames) {
  #     l <- c(l, unique(strata.features[name]))
  #   }
  #   return(expand.grid(l))
  # }
  strataIndex <- function(strata.tbl, strata.features, i) {
    colnames <- colnames(strata.features)
    if (length(colnames) == 1) {
      return(which(strata.tbl[, colnames[1]] == strata.features[i,colnames[1]]))
    }
    sub.strata.tbl <- strata.tbl
    for (name in colnames) {
      sub.strata.tbl <- sub.strata.tbl[sub.strata.tbl[,name] == strata.features[i, name], ]
    }
    if (nrow(sub.strata.tbl) != 1) {
      stop("There could not be found a unique stratum.")
    }
    return(as.numeric(rownames(sub.strata.tbl)))
  }

  gen.stratum <- function(strata.tbl, strata.features) {
    n <- nrow(strata.features)
    stratum <- rep(NA, n)
    for (i in 1:n) {
      stratum[i] <- strataIndex(strata.tbl = strata.tbl,
                                strata.features = strata.features, i = i)
    }
    return(stratum)
  }
  weights.min <- weights/min(weights)
  #if (is.null(strata.variables)) {
  #  stop("Vector strats indicating the names of the strats need to be provided!")
  #}
  if (!all(strata.variables %in% colnames(data))) {
    stop("Not all strata variables could be detected!") # deleted one 'not'
  }
  strata.features <- data[strata.variables]

  # transfers strata variables and their exposures into a dummy coded representation
  transferFactor <- function(df) {
    factor.names <- c()
    class.list <- sapply(df, class)
    newdf <- df
    if ("factor" %in% class.list) {
      warning('Data contains factor variables: converted to dummy coded variables ')
      factor.names <- names(class.list[which(class.list ==
                                               "factor")])
      newdf <- df[setdiff(colnames(df), factor.names)]
      for (name in factor.names){
        dum.coded <- model.matrix(formula(paste0('~',name)),df)
        dum.coded <- as.data.frame(dum.coded)
        colnames <- colnames(dum.coded)[-1]
        dum.coded <- as.data.frame(dum.coded[,-1])
        colnames(dum.coded) <- colnames
        newdf <- cbind(newdf,dum.coded)
      }
    }
    return(list(factor.names = factor.names, df.nofactors = newdf))
  }

  tr.fac <- transferFactor(strata.features)
  strata.features <- tr.fac$df.nofactors
  #factor.names <- tr.fac$factor.names
  strata.tbl <- unique(strata.features) # nk edit: simply uniqueness of the COMBINATIONS; not expand.grid!!
  if (is.null(stratum)) {
    stratum <- gen.stratum(strata.tbl, strata.features)
  }

  data.new <- data.frame(matrix(NA, sum(round(weights.min)),
                                ncol(data)))
  ind <- 0
  subs <- data[setdiff(colnames(data), strata.variables)]

  # transforms all factors of data frame df into numeric
  factorsToNumeric <- function(df){
    indx <- sapply(df, is.factor)
    if(any(indx)){
      df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
      factornames <- names(indx[indx==TRUE])
      warning(paste0('The Columns [',paste(factornames,collapse = ' '),'] have been transformed from factor to numeric values'))
    }
    return(df)
  }
  # note that subs should merely contain numerics (as subs contains is data less the
  # strata variables), if not so all factors will be corrected for numerics
  subs <- factorsToNumeric(subs)

  # this function accounts for changes in the output variable. The output variable
  # may not fit to the data. As thereby new combinations of the strata.variables evolved
  # new strata were produced.
  # this functions returns a stratum vector correted for these new evolved strata.
  # adjStrata <- function(stratum){
  #   # here we consider the case where strata.variable just holds one varaible
  #   if(length(strata.variables)!=1){
  #     return(stratum)
  #   }
  #   # strata.tbl <- unique(data.frame(strata.features,stratum))
  #   # #to.stratify <- strata.tbl[which(duplicated(strata.tbl$stratum,fromLast = TRUE)),]
  #   # 
  #   # max.strat <- max(strata.tbl$stratum)
  #   # strata.tbl <- strata.tbl %>% group_by(stratum) %>% filter(n()>1)
  #   # strata.tbl <- as.data.frame(strata.tbl)
  #   # dupl <- duplicated(strata.tbl[,'stratum'])
  #   # if(length(dupl)==0){
  #   #   return(stratum) # no duplicates could be detected
  #   # }
  #   # warning("further combinations from \'stratum\' and strata.varibales evolved: new strata were produced!")
  #   # 
  #   # new.strata <- seq(max(stratum)+1,length.out = sum(dupl))
  #   # strata.tbl[dupl,'stratum'] <- new.strata
  #   # 
  #   # adapt <- function(stratum, new.strata.tbl){
  #   #   for(i in 2:nrow(new.strata.tbl)){
  #   #     stratum[strata.features[,1]==new.strata.tbl[i,1] & new.strata.tbl[1,2] == stratum] = new.strata.tbl[i,2]
  #   #   }
  #   #   return(stratum)
  #   # }
  #   # 
  #   # orig.strata.idx <- which(dupl==FALSE)
  #   # #for(i in 1:length(orig.strata.idx)-1){
  #   #   # TODO
  #   # #}
  #   # new.strata.idx <- seq(orig.strata.idx[length(orig.strata.idx)],length(dupl))
  #   # new.strata.tbl <- strata.tbl[new.strata.idx,]
  #   # stratum <- adapt(stratum,new.strata.tbl)
  #   # return(stratum)
  # }
  # 
  # 
  # stratum <- adjStrata(stratum)

  for (h in unique(stratum)) {
    ind.weights.h <- which(stratum == h)
    subs.h <- subset(subs, stratum == h)
    weights.h <- round(weights.min[ind.weights.h])
    N.h <- sum(weights.h)
    n.h <- nrow(subs.h)

    dfEmpty <- function(df){
      return(ncol(df)==0)
    }

    if(!dfEmpty(strata.features)){
      strata.features.h <- subset(strata.features, stratum ==h)
    }else{
      strata.features.h <- data.frame(a=rep(0,sum(stratum==h)))[c()]
    }

    if (distr == "mvnorm") {
      pos.var.0 <- which(apply(subs.h[setdiff(colnames(subs.h),
                                              strata.variables)], 2, var) == 0)
      if (length(pos.var.0) != 0)
        subs.h[, pos.var.0] <- subs.h[, pos.var.0] +
          rnorm(nrow(subs.h), 1e-04)
      if (type == "stochIP") {
        cor.fac <- (weights.h[1] - 1)/(weights.h[1] *
                                         n.h - 1)
        term.add.h <- rmvnorm(N.h, mean = rep(0, ncol(subs.h)),sigma = cor.fac * cov(subs.h))
        new.obs.h <- cbind(ipOversampling(strata.features.h,weights.min[ind.weights.h]),
                           ipOversampling(subs.h,weights.min[ind.weights.h]) + term.add.h)
      }
      else {
        if(dim(subs.h)[1]==1){
          stop(paste0('stratum ', h, ' just has 1 observation.  Each stratum should have at least 2 observations.')) # there is just one observation
        }
        new.obs.h <- cbind(ipOversampling(strata.features.h, weights.min[ind.weights.h]),
                           data.frame(rmvnorm(N.h,mean = colMeans(subs.h), sigma = cov(subs.h))))
      }
    }
    data.new[(ind + 1):(ind + N.h), ] <- new.obs.h
    ind <- ind + N.h
    colnames(data.new) <- colnames(new.obs.h)
  }

  return(list(data = data.new, orig.data = data, stratum=stratum, method = type,
              strata.tbl = strata.tbl, N = nrow(data.new), n = nrow(data)))
}
