#' Plain replication of each observation by inverse-probability weights
#' @author Norbert Krautenbacher, Kevin Strauss, Maximilian Mandl, Christiane Fuchs
#' @description This method corrects for the sample selection bias by the plain replication of each observation in the sample according to its IP weight,
#' i.e. in a stratified random sample one replicates an observation of stratum h by the factor w_h.
#' @param data a data frame containing the observations rowwise, along with their corresponding categorical strata feature(s).
#' @param weights a numerical vector whose length must coincide with the number of the rows of data. The i-th value contains the inverse-probability e.g. determines how often the i-th observation of data shall be replicated.
#' @param normalize If weight vector should be normalized, i.e. the smallest entry of the vector will be set to 1.
#' @details If the numeric vector contains numbers which are not natural but real, they will be rounded.
#' @examples
#' library(smotefamily)
#' library(sambia)
#' data.example <- sample_generator(100,ratio = 0.80)
#' result <- gsub('n','0',data.example[,'result'])
#' result <- gsub('p','1',result)
#' data.example[,'result'] <- as.numeric(result)
#' weights <- data.example[,'result']
#' weights <- ifelse(weights==1,1,4)
#' sample <- sambia::ipOversampling(data.example,weights)
#' @export


ipOversampling <- function(data, weights, normalize = FALSE){

  if(nrow(data)!=length(weights)){
    stop('Length of weights and observations in data must coincide.')
  }
  if(!is.numeric(weights)){
    stop('Weights are required to be numeric.')
  }
  if(min(weights)<=0){
    stop('All weights are required to be greater than 0.')
  }

  if(normalize){
    weights <- weights/min(weights) # let minimum weight equal one --> as less observations as possible (as much as necessary)
  }
  pos <- rep(1:nrow(data), round(weights))
  data <- data.frame(data[pos, , drop=FALSE])
  rownames(data) <- 1:nrow(data)
  return(data)
}
