#' Rejection Sampling is a method used in sambia's function 'costing' (Krautenbacher et al, 2017).
#' @author Norbert Krautenbacher, Kevin Strauss, Maximilian Mandl, Christiane Fuchs
#' @description Rejection Sampling is a method used in sambias costing function. It is sampling scheme that allows us
#' to draw examples independently from a distribution X, given examples drawn independently from distribution Y.
#' @references Krautenbacher, N., Theis, F. J., & Fuchs, C. (2017). Correcting Classifiers for Sample Selection Bias in Two-Phase Case-Control Studies. Computational and mathematical methods in medicine, 2017.
#' @param data a data frame containing the observations rowwise, along with their corresponding categorical strata feature
#' @param weights a numerical vector whose length must coincide with the number of the rows of data. The i-th value contains the inverse-probability e.g. determines how often the i-th observation of data shall be replicated.
#' @examples
#' library(smotefamily)
#' library(sambia)
#' data.example <- sample_generator(100,ratio = 0.80)
#' result <- gsub('n','0',data.example[,'result'])
#' result <- gsub('p','1',result)
#' data.example[,'result'] <- as.numeric(result)
#' weights <- data.example[,'result']
#' weights <- ifelse(weights==1,1,4)
#' rej.sample <- sambia:::rejSamp(data=data.example, weights = weights)



rejSamp <- function(data, weights){
  index.if = rep(FALSE,nrow(data))
  for(i in 1:nrow(data)){
    prob = weights[i]/max(weights)
    if(rbinom(size = 1, prob = prob, n = 1) == 1){
      index.if[i] = TRUE
    }
  }
  return(data[index.if, , drop=FALSE])
}

