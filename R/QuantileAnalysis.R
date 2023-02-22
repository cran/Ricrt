#' @title Quantile Effect Analysis
#' @description
#' A function for quantile analysis that provides a p-value and a point estimation.
#'
#' @param S A numeric column vector with individuals' stratum number
#' @param C A numeric column vector with individuals' cluster number
#' @param Z A numeric column vector with individuals' treatment assignment
#' (binary)
#' @param R A numeric column vector with individuals' outcome
#' @param quantile A numeric value for intended quantile in the input data for
#' running the quantile test, the default of this value is 0.5.
#' @param delta_hyp A numeric value that is the hypothesized treatment effect
#' given the quantile. The default of this value is 0.
#' @importFrom randomForest randomForest
#' @importFrom stats lm
#' @importFrom stats ecdf
#' @importFrom stats median
#' @importFrom stats predict
#' @import tidyverse
#' @import dplyr
#' @return A list of the outputs
#' @examples
#' # First we need to obtain the vectors for the inputs.
#' S = example1$S
#' C = example1$C
#' Z = example1$Z
#' R = example1$R
#' \donttest{QuantileAnalysis(S, C, Z, R)}
#' @export
QuantileAnalysis <- function(S, C, Z, R, quantile = 0.5, delta_hyp = 0) {
  # check whether the inputs are valid:
  check = length(S)
  if (check != length(S) || check != length(C) || check != length(Z) || check != length(R)) {
    stop("The input vectors have different length")
  }

  df = cbind(S, C, Z, R)
  df = as.data.frame(df)

  # getting basic information
  col = ncol(df)
  S = max(df$S) # number of stratum
  # getting the number of clusters in one stratum and the number of treatment group in each stratum
  k_s = rep(0, S)
  m_s = rep(0, S)
  for (i in 1:S) {
    temp = df %>% select(S, C, Z) %>% group_by(S, C) %>% summarise(Z = max(Z), .groups = 'drop') %>% filter(S == i)
    if (max(temp$Z) != 1 || min(temp$Z) != 0) {
      stop("Please input valid treatment assignment. ")
    }
    k_s[i] = max(temp$C)
    m_s[i] = sum(temp$Z)
  }
  K = max(k_s)


  n_s = array(0, S)# get n_s
  for (s in 1:S) {
    n_s[s] = nrow(df[df$S == s,])
  }
  df_lm = df
  cov_col = 0
  indicator = TRUE


  # compute the quantile effect
  get_e_quantile <- function(df, df_lm, delta_hyp, cov_col, indicator, quantile) {
    df_lm$A = df$R - (1 - df$Z) * delta_hyp
    if (indicator) {
      df$A = df_lm$A
    } else {
      model = lm(A ~ ., df_lm)
      df$A = df_lm$A - predict(model, df_lm[,1:cov_col])
    }
    A_quantile = quantile(df$A, quantile)[[1]]
    df$e = 0
    for (i in 1:nrow(df)) {
      if (df[i, ]$A >= A_quantile) {
        df[i, ]$e = 1
      } else {
        df[i, ]$e = 0
      }
    }
    return(df)
  }

  get_H <- function(df, S, df_lm, delta_hyp, cov_col, indicator, quantile) {

    dff = df
    # calculate q, q_bar, and weight
    q = array(0, dim = c(S, K))
    q_bar = array(0, dim = S)
    weight = array(0, dim = S)
    for (i in 1:S) {
      n = array(0, k_s[i])
      for (j in 1:k_s[i]) {
        temp = dff %>% filter(dff$S == i, dff$C == j)
        n_sk = nrow(temp)
        q[i, j] = (1 / n_sk) * sum(temp$e)
        n[j] = nrow(temp)
      }
      q_bar[i] = 0.5 * (q[i, 1] + q[i, 2])
      weight[i] = (2 * n[1] * n[2]) / (n[1] + n[2])
    }

    # calculating H-statistic
    H = 0
    for (i in 1:S) {
      for (j in 1:k_s[i]) {
        temp = dff %>% filter(dff$S == i, dff$C == j)
        H = H + temp$Z[1] * weight[i] * (q[i, j] - q_bar[i])
      }
    }
    return(H)
  }

  get_H_P <- function(df, S, df_lm, delta_hyp, cov_col, indicator, quantile, n_s, n_permute = 500) {

    dff = get_e_quantile(df, df_lm, delta_hyp, cov_col, indicator, quantile)
    H_actual = get_H(dff, S, df_lm, delta_hyp, cov_col, indicator, quantile)
    df_temp = dff
    dd = dff %>% group_by(S, C) %>% summarize(Z = max(Z), .groups = 'drop')

    # permutate treatments
    H_sim = rep(0, n_permute)
    for (iter in 1:n_permute){
      for (s in 1:S) {
        ind = unique(sample(rep(1:k_s[s]), m_s[s], replace = FALSE))
        for (c in 1:k_s[s]) {
          if (c %in% ind) {
            dd[dd$S == s, ]$Z[c] = 1
          } else {dd[dd$S == s, ]$Z[c] = 0}
        }
      }
      df_temp = df_temp %>% select(-Z)
      df_temp = left_join(df_temp, dd, by = c("S", "C"))
      H_sim[iter] = get_H(df_temp, S, df_lm, delta_hyp, cov_col, indicator, quantile)
    }

    if (H_actual > median(H_sim)) {
      return((1 - ecdf(H_sim)(H_actual)) * 2)
    } else if (H_actual < median(H_sim)) {
      return(ecdf(H_sim)(H_actual) * 2)
    } else {return(1)}
  }

  H_estimation <- function(df, S, df_lm, delta_hyp, cov_col, indicator, quantile) {
    ## A function that calculate the expected ranksum value
    total = length(df$S)
    average = ((1 + total) * total * 0.5) / total
    Sumexp = 0

    delta_hat = delta_hyp
    distance = 0.1
    if (get_H(get_e_quantile(df, df_lm, delta_hat, cov_col, indicator, quantile), S, df_lm, delta_hat, cov_col, indicator, quantile) < Sumexp) {
      while (1 > 0) {
        if (get_H(get_e_quantile(df, df_lm, delta_hat + distance, cov_col, indicator, quantile), S, df_lm, delta_hat, cov_col, indicator, quantile) < Sumexp) {
          delta_hat = delta_hat + distance
          distance = distance * 2
        } else {distance = distance / 2}
        if (distance < 0.0001) {break}
      }
    }
    if (get_H(get_e_quantile(df, df_lm, delta_hat, cov_col, indicator, quantile), S, df_lm, delta_hat, cov_col, indicator, quantile) > Sumexp) {
      while (1 > 0) {
        if (get_H(get_e_quantile(df, df_lm, delta_hat - distance, cov_col, indicator, quantile), S, df_lm, delta_hat, cov_col, indicator, quantile) > Sumexp) {
          delta_hat = delta_hat - distance
          distance = distance * 2
        } else {distance = distance / 2}
        if (distance < 0.0001) {break}
      }
    }
    return(delta_hat)
  }

  dff = get_e_quantile(df, df_lm, delta_hyp, cov_col, indicator, quantile)
  H_stats = get_H(dff, S, df_lm, delta_hyp, cov_col, indicator, quantile)
  Quantile_P = get_H_P(df, S, df_lm, delta_hyp, cov_col, indicator, quantile, n_s, n_permute = 100)
  estimation = H_estimation(df, S, df_lm, delta_hyp, cov_col, indicator, quantile)
  result = list(delta_hyp = delta_hyp,
                quantile = quantile,
                H_statistics = H_stats,
                p_value = Quantile_P,
                estimation = estimation)
}
