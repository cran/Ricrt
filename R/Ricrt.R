#' @title Ricrt
#' @description
#' This package can use Mann-Whitney-Wilcoxon or signed-rank test to perform
#' randomization inference. The statistics, p-value, point estimation, and a
#' two-sided 95% confidence interval will be returned.
#'
#' @param S A numeric column vector with individuals' stratum number
#' @param C A numeric column vector with individuals' cluster number
#' @param Z A numeric column vector with individuals' treatment assignment
#' (binary)
#' @param R A numeric column vector with individuals' outcome
#' @param X A numeric matrix with each column being a covariate
#' @param tau_hyp A numeric value for hypothesized treatment effect,
#' the default for this value is 0.
#' @param method A string being either "W" or "sr", indicating either
#' weighted sum of S Mann–Whitney–Wilcoxon statistics will be used or
#' signed-rank test will be used
#' @param reg A string being either "lm" or "rf," indicating either
#' linear model or random forest model being used for fitting the data with
#' covariates. The default is "lm."
#' @param permutation A numeric value indicating the number of permutation
#' inside the function when using permutation tests for p-values,
#' the default is 50.
#' @import dplyr
#' @importFrom randomForest randomForest
#' @importFrom stats lm
#' @importFrom stats ecdf
#' @importFrom stats median
#' @importFrom stats predict
#' @importFrom rlang is_null
#' @import tidyverse
#' @import SuperLearner
#' @import glmnet
#' @importFrom Rdpack reprompt
#' @return A list of the outputs
#' @examples
#' # First we need to obtain the vectors for the inputs.
#' S = example1$S
#' C = example1$C
#' Z = example1$Z
#' R = example1$R
#' X = cbind(example1$X1, example1$X2, example1$X3, example1$X4, example1$X5)
#'
#' # Let's see the first example with method = W and reg = lm.
#' set.seed(123)
#' \donttest{Ricrt(S, C, Z, R, X, tau_hyp = 10, method = "W", reg = "lm", permutation = 5)}
#'
#' # Let's see the second example with method = W and reg = rf
#' \donttest{Ricrt(S, C, Z, R, X, tau_hyp = 10, method = "W", reg = "rf", permutation = 5)}
#' @export
Ricrt <- function(S, C, Z, R, X = NULL, tau_hyp = 0, method = "W", reg = "lm", permutation = 100){

  # check whether the inputs are valid:
  check = length(S)
  if (check != length(S) || check != length(C) || check != length(Z) || check != length(R)) {
    stop("The input vectors have different length.")
  }
  if (check != nrow(X) && !is_null(X)) {
    stop("The input matrix has wrong row number.")
  }
  if (reg != "lm" && reg != "rf") {
    stop("Please input valid regression method.")
  }
  # combine vectors to dataframe
  if (is.null(X)) {
    df = cbind(S, C, Z, R)
    indicator = TRUE
  } else {
    df = cbind(S, C, Z, R, X)
    indicator = FALSE
  }
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
      return("Please input valid treatment assignment. ")
    }
    k_s[i] = max(temp$C)
    m_s[i] = sum(temp$Z)
  }
  n_s = array(0, S)# get n_s
  for (s in 1:S) {
    n_s[s] = nrow(df[df$S == s,])
  }
  # separate the covariates:
  if (indicator == FALSE) {
    df_lm = df[,5:col]
    cov_col = col - 4 # the number of covariance
  } else {
    df_lm = df
    cov_col = 0
  }
  df = df[, 1:4]

  get_e <- function(df, df_lm, tau_hyp, cov_col, indicator) {
    df_lm$e = df$R - df$Z * tau_hyp
    e = as.matrix(df_lm$e)
    cov = df_lm %>% select(-e)
    if (indicator) {
      df$e = df_lm$e
    } else {
      if (reg == "lm") {
        model = SuperLearner(e, cov, SL.library = c("SL.lm"))
        y_pred = predict(model, onlySL = TRUE, newdata = df_lm[, 1:cov_col])
        df$e = df_lm$e - y_pred$pred
      } else if (reg == "rf") {
        rf.fit = SuperLearner(e, cov, SL.library = c("SL.mean", "SL.lm", "SL.randomForest"))
        y_pred = predict(rf.fit, onlySL = TRUE, newdata = df_lm[, 1:cov_col])
        df$e = df_lm$e - y_pred$pred
      }
    }
    return(df)
  }

  if (method == "W") {
    # a function the return df with updated e
    df = get_e(df, df_lm, tau_hyp, cov_col, indicator)
    # calculate the W-statistic:
    get_W <- function(df, n_s, S) {
      u = 0
      for (s in 1:S) {
        ds = df[df$S == s,]
        ds1 = ds[ds$Z == 1, ]
        temp = ds1$e
        for (i in 1:nrow(ds1)) {
          u = u + sum(-sign(ds$e - temp[i])) * (1 / (1 + n_s[s]))
        }
      }
      return(u)
    }

    W_stats = get_W(df, n_s, S)
    # calculate the variance of W-statistic:
    get_W_variance <- function(df, n_s, S) {
      u = 0
      for (s in 1:S) {
        u1 = 0
        u2 = 0
        ds = df[df$S == s,]
        ds1 = ds[ds$Z == 1, ]
        ds2 = ds[ds$Z == 0, ]
        temp = ds1$e
        temp2 = ds2$e
        for (i in 1:nrow(ds1)) {
          u1 = u1 + sum(-sign(ds$e - temp[i])) * (1 / (1 + n_s[s]))
        }
        for (i in 1:nrow(ds2)) {
          u2 = u2 + sum(-sign(ds$e - temp2[i])) * (1 / (1 + n_s[s]))
        }
        u = u + (u1^2 + u2^2) * 0.5
      }
      return(u)
    }
    Variance = get_W_variance(df, n_s, S)
    # Calculate the p-value
    get_W_P <- function(df, df_lm, tau_hyp, n_s, cov_col, S, n_permute = permutation) {
      df = get_e(df, df_lm, tau_hyp, cov_col, indicator)
      W_actual = get_W(df, n_s, S)
      df_temp = df
      dd = df %>% group_by(S, C) %>% summarize(Z = max(Z), .groups = 'drop')

      # permutate treatments
      W_sim = rep(0, n_permute)
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
        W_sim[iter] = get_W(df_temp, n_s, S)
      }
      if (W_actual > median(W_sim)) {
        return((1 - ecdf(W_sim)(W_actual)) * 2)
      } else {return(ecdf(W_sim)(W_actual) * 2)}
    }

    P_value = get_W_P(df, df_lm, tau_hyp, n_s, cov_col, S)
    get_W_estimation <- function(df, df_lm, tau_hyp, n_s, cov_col, S) {
      tau_hat = tau_hyp
      distance = 0.01
      if (get_W(get_e(df, df_lm, tau_hat - distance, cov_col, indicator), n_s, S) < 0) {
        while (1 > 0) {
          if (get_W(get_e(df, df_lm, tau_hat - distance, cov_col, indicator), n_s, S) < 0) {
            tau_hat = tau_hat - distance
            distance = distance * 2
          }else {distance = distance / 2}
          if (distance < 0.00001) {break}
        }
      }
      if (get_W(get_e(df, df_lm, tau_hat + distance, cov_col, indicator), n_s, S) > 0) {
        while (1 > 0) {
          if (get_W(get_e(df, df_lm, tau_hat + distance, cov_col, indicator), n_s, S) > 0) {
            tau_hat = tau_hat + distance
            distance = distance * 2
          }else {distance = distance / 2}
          if (distance < 0.00001) {break}
        }
      }
      return(tau_hat)
    }
    estimation = get_W_estimation(df, df_lm, tau_hyp, n_s, cov_col, S)
    lower_bound = estimation
    upper_bound = estimation
    distance = 0.1
    if (get_W_P(df, df_lm, lower_bound, n_s, cov_col, S, n_permute = 50) < 0.05) {
      warning("Estimation has a p-value of less than 0.05, unable to obtain feasible a confidence interval")
    }
    while(1 > 0) {
      if (get_W_P(df, df_lm, lower_bound - distance, n_s, cov_col, S, n_permute = permutation) >= 0.05) {
        lower_bound = lower_bound - distance
        distance = distance * 2
      }else {distance = distance / 2}
      if (distance < 0.01) {break}
    }

    distance = 0.1
    while(1 > 0) {
      if (get_W_P(df, df_lm, upper_bound + distance, n_s, cov_col, S, n_permute = permutation) >= 0.05) {
        upper_bound = upper_bound + distance
        distance = distance * 2
      }else {distance = distance / 2}
      if (distance < 0.01) {break}
    }
    result <- list(tau_hyp = tau_hyp,
                   W_statistic = W_stats,
                   W_statistic_Variance = Variance,
                   p_value = P_value,
                   estimation = estimation,
                   lower_bound_confidence_interval = lower_bound,
                   upper_bound_confidence_interval = upper_bound)
    return(result)
  } else if (method == "sr") {

    # a function the return df with updated e
    df = get_e(df, df_lm, tau_hyp, cov_col, indicator)

    ## A function that calculates the Hodges-Lemman statistic given the hypthesis
    get_HL <- function(df, df_lm, cov_col, indicator) {
      stratamean = rep(0, S)
      for (x in 1:S) {
        temp = df %>% filter(df$S == x)
        stratamean[x] = mean(temp$e)
      }
      df = df[order(df$S),]
      ## subtract response from mean
      df$Anew = df$e
      ind = 0
      for (x in 1:S) {
        for (y in 1:n_s[x]) {
          df$Anew[ind + y] = df$Anew[ind + y] - stratamean[x]
        }
        ind = ind + n_s[x]
      }
      df$q = rank(df$Anew)
      ## calculate Hodges-Lehmann statistics
      return(sum(df$q * df$Z))
    }
    HL_statistics = get_HL(df, df_lm, cov_col, indicator)

    ## A function that outputs the p-value using Hodges-Lehmman given the hypothesized tau value
    twosided_HodgesLehmann_p <- function(df, df_lm, tau_hyp, cov_col, indicator, n_permute = permutation) {
      df = get_e(df, df_lm, tau_hyp, cov_col, indicator)
      HL_actual = get_HL(df, df_lm, cov_col, indicator)
      df_temp = df
      dd = df %>% group_by(S, C) %>% summarize(Z = max(Z), .groups = 'drop')

      # permutate treatments
      HL_sim = rep(0, n_permute)
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
        HL_sim[iter] = get_HL(df_temp, df_lm, cov_col, indicator)
      }

      if (HL_actual > median(HL_sim)) {
        return((1 - ecdf(HL_sim)(HL_actual)) * 2)
      } else {return(ecdf(HL_sim)(HL_actual) * 2)}
    }
    HL_P_value = twosided_HodgesLehmann_p(df, df_lm, tau_hyp, cov_col, indicator, n_permute = permutation)

    ## A function that outputs the estimation of the treatment effect
    Hodges_lehmman_estimation <- function(df, df_lm, tau_hyp, cov_col, indicator) {
      ## A function that calculate the expected ranksum value
      total = length(df$S)
      average = ((1 + total) * total * 0.5) / total
      Sumexp = average * sum(df$Z)

      tau_hat = tau_hyp
      distance = 0.01
      if (get_HL(get_e(df, df_lm, tau_hat - distance, cov_col, indicator), df_lm, cov_col, indicator) < Sumexp) {
        while (1 > 0) {
          if (get_HL(get_e(df, df_lm, tau_hat - distance, cov_col, indicator), df_lm, cov_col, indicator) < Sumexp) {
            tau_hat = tau_hat - distance
            distance = distance * 2
          }else {distance = distance / 2}
          if (distance < 0.00001) {break}
        }
      }
      if (get_HL(get_e(df, df_lm, tau_hat + distance, cov_col, indicator), df_lm, cov_col, indicator) > Sumexp) {
        while (1 > 0) {
          if (get_HL(get_e(df, df_lm, tau_hat + distance, cov_col, indicator), df_lm, cov_col, indicator) > Sumexp) {
            tau_hat = tau_hat + distance
            distance = distance * 2
          } else {distance = distance / 2}
          if (distance < 0.00001) {break}
        }
      }
      return(tau_hat)
    }
    HL_estimation = Hodges_lehmman_estimation(df, df_lm, tau_hyp, cov_col, indicator)
    ### Getting the 95 percent confidence interval:
    lower_bound = HL_estimation
    upper_bound = HL_estimation
    distance = 0.01
    if (twosided_HodgesLehmann_p(df, df_lm, lower_bound, cov_col, indicator, n_permute = 100) < 0.05) {
      warning("Estimation has a p-value of less than 0.05, unable to obtain feasible a confidence interval")
    }

    while(1 > 0) {
      if (twosided_HodgesLehmann_p(df, df_lm, lower_bound - distance, cov_col, indicator, n_permute = 50) >= 0.05) {
        lower_bound = lower_bound - distance
        distance = distance * 2
      } else {
        distance = distance / 2
      }
      if (distance < 0.0001) {
        break
      }
    }

    distance = 0.01
    while(1 > 0) {
      if (twosided_HodgesLehmann_p(df, df_lm, upper_bound + distance, cov_col, indicator, n_permute = 50) >= 0.05) {
        upper_bound = upper_bound + distance
        distance = distance * 2
      } else {
        distance = distance / 2
      }
      if (distance < 0.0001) {
        break
      }
    }
    result <- list(tau_hyp = tau_hyp,
                   signed_rank_stats = HL_statistics,
                   p_value = HL_P_value,
                   estimation = HL_estimation,
                   lower_bound_confidence_interval = lower_bound,
                   upper_bound_confidence_interval = upper_bound)
    return(result)
  } else {
    stop("Please input valid method: either W or sr")
  }
}
