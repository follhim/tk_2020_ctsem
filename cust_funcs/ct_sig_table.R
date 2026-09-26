#' Significance table for a ctsem fit (95% intervals + stars)
#'
#' Pulls the population means, TIpred effects, random-effect SDs and
#' random-effect correlations out of `summary(fit)` and flags every parameter
#' whose 95% interval excludes 0.
#'
#' @param fit      a ctsem fit (julia or stan backend)
#' @param blocks   which summary blocks to include
#' @param only_sig TRUE = keep only the starred rows
#' @param digits   decimals shown in the gt table
#' @param gt       TRUE = return a gt table, FALSE = a data frame
ct_sig_table <- function(fit,
                         blocks = c("popmeans", "tipreds", "popsd", "rawpopcorr"),
                         only_sig = FALSE, digits = 3, gt = TRUE) {
  s <- summary(fit)
  labs <- c(popmeans   = "Population means (fixed effects)",
            tipreds    = "TIpred effects (moderation)",
            popsd      = "Random-effect SDs",
            rawpopcorr = "Random-effect correlations")
  pull <- function(b) {
    x <- s[[b]]
    if (is.null(x) || !(is.data.frame(x) || is.matrix(x)) || nrow(x) == 0) return(NULL)
    x <- as.data.frame(x, check.names = FALSE)
    num <- function(col) {
      hit <- names(x)[tolower(names(x)) == tolower(col)]
      if (length(hit) == 0) return(rep(NA_real_, nrow(x)))
      suppressWarnings(as.numeric(x[[hit[1]]]))
    }
    data.frame(block = unname(labs[b]), param = rownames(x),
               mean = num("mean"), sd = num("sd"),
               lower = num("2.5%"), upper = num("97.5%"),
               stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, lapply(intersect(blocks, names(s)), pull))
  if (is.null(out) || nrow(out) == 0) stop("No summary blocks with intervals found.")
  # SDs are bounded at 0, so "excludes 0" is not a meaningful test for them
  testable <- out$block != labs[["popsd"]]
  out$sig <- ifelse(testable & !is.na(out$lower) & !is.na(out$upper) &
                      (out$lower > 0 | out$upper < 0), "*", "")
  if (only_sig) out <- out[out$sig == "*", , drop = FALSE]
  rownames(out) <- NULL
  if (!gt) return(out)
  out |>
    gt::gt(groupname_col = "block") |>
    gt::cols_label(param = "Parameter", mean = "Estimate", sd = "SD",
                   lower = "2.5%", upper = "97.5%", sig = "") |>
    gt::fmt_number(columns = c(mean, sd, lower, upper), decimals = digits) |>
    gt::tab_style(style = gt::cell_text(weight = "bold"),
                  locations = gt::cells_body(rows = sig == "*")) |>
    gt::tab_footnote("* = 95% interval excludes 0 (not tested for random-effect SDs).")
}

#' Any results table -> gt with a significance-star column and bold sig rows
#'
#' @param df     data frame (e.g. TIpred effects, cross-effect peaks, half-lives)
#' @param lower  column holding the lower interval bound
#' @param upper  column holding the upper interval bound
#' @param digits decimals for numeric columns
ct_sig_gt <- function(df, lower = "2.5%", upper = "97.5%", digits = 3) {
  df <- as.data.frame(df, check.names = FALSE)
  lo <- suppressWarnings(as.numeric(df[[lower]]))
  hi <- suppressWarnings(as.numeric(df[[upper]]))
  df$sig <- ifelse(!is.na(lo) & !is.na(hi) & (lo > 0 | hi < 0), "*", "")
  num <- names(df)[vapply(df, is.numeric, logical(1))]
  g <- gt::gt(df) |>
    gt::fmt_number(columns = dplyr::all_of(num), decimals = digits) |>
    gt::cols_label(sig = "")
  sig_rows <- which(df$sig == "*")
  if (length(sig_rows))
    g <- gt::tab_style(g, style = gt::cell_text(weight = "bold"),
                       locations = gt::cells_body(rows = sig_rows))
  g |> gt::tab_footnote("* (bold) = 95% interval excludes 0.")
}
