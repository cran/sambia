#' smoteNew is a necessary function that modifies the SMOTE algorithm.
#' @author Norbert Krautenbacher, Kevin Strauss, Maximilian Mandl, Christiane Fuchs
#' @description smoteNewis a necessary function that modifies the SMOTE algorithm in the following ways: (1) correct bug in original
#' smotefamily::SMOTE() function and (2) lets the user specifiy which class to be oversampled.
#' @param data.x A data frame or matrix of numeric-attributed dataset
#' @param data.y A vector of a target class attribute corresponding to a dataset X
#' @param K The number of nearest neighbors during sampling process
#' @param dup_size The number or vector representing the desired times of synthetic minority instances over the original number of majority instances
#' @param class.to.oversample Class to be oversampled
#' @examples
#' library(smotefamily)
#' library(sambia)
#' data.example = sample_generator(10000,ratio = 0.80)
#' genData = sambia:::smoteNew(data.example[,-3],data.example[,3],K = 5,class.to.oversample = 'p')


smoteNew <- function (data.x, data.y, K, dup_size = 0, class.to.oversample){
  ncD = ncol(data.x)
  n_data.y = table(data.y)
  classP = as.character(class.to.oversample)
  P_set = subset(data.x, data.y == as.character(class.to.oversample) )[sample(n_data.y[ classP]    ),
    ]
  N_set = subset(data.x, data.y != classP)
  P_class = rep(classP, nrow(P_set))
  N_class = data.y[data.y != classP]
  sizeP = nrow(P_set)
  sizeN = nrow(N_set)
  knear = smotefamily::knearest(P_set, P_set, K)
  sum_dup = smotefamily::n_dup_max(sizeP + sizeN, sizeP, sizeN, dup_size = dup_size)# corrected line
  syn_dat = NULL
  for (i in 1:sizeP) {
    if (is.matrix(knear)) {
      pair_idx = knear[i, ceiling(runif(sum_dup) * K)]
    }
    else {
      pair_idx = rep(knear[i], sum_dup)
    }
    g = runif(sum_dup)
    P_i = matrix(unlist(P_set[i, ]), sum_dup, ncD, byrow = TRUE)
    Q_i = as.matrix(P_set[pair_idx, ])
    syn_i = P_i + g * (Q_i - P_i)
    syn_dat = rbind(syn_dat, syn_i)
  }
  P_set[, ncD + 1] = P_class
  colnames(P_set) = c(colnames(data.x), "stratum")
  N_set[, ncD + 1] = N_class
  colnames(N_set) = c(colnames(data.x), "stratum")
  rownames(syn_dat) = NULL
  syn_dat = data.frame(syn_dat)
  syn_dat[, ncD + 1] = rep(classP, nrow(syn_dat))
  colnames(syn_dat) = c(colnames(data.x), "stratum")
  NewD = rbind(P_set, syn_dat, N_set)
  rownames(NewD) = NULL
  D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set,
    orig_P = P_set, K = K, K_all = NULL, dup_size = sum_dup,
    outcast = NULL, eps = NULL, method = "SMOTE")
  class(D_result) = "gen_data"
  return(D_result)
}

