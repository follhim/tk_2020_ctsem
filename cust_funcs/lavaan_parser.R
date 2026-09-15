ct_lavaan <- function(spec, latentNames = NULL, manifestNames = NULL, data = NULL) {

  raw_lines <- sub("#.*$", "", strsplit(spec, "\n")[[1]])
  raw_lines <- trimws(raw_lines)
  raw_lines <- raw_lines[nzchar(raw_lines)]

  meas_lines <- list(); struct_lines <- list(); cov_lines <- list()
  t0cov_lines <- list()
  cint_lines <- list(); t0means_lines <- list(); mm_lines <- list()
  intercept_lines_raw <- list()

  # ---- manifest type keywords: BINARY:/ORDINAL:/COUNT:/CENSOR:/CONTINUOUS: ----
  # manifesttype codes per ctsem (juliaFit branch): 0 continuous, 1 binary,
  # 3 Poisson/count, 4 censored -- all four confirmed working in this project.
  # 2 = ordinal is ctsem's documented code but has NOT been exercised in this
  # project yet -- verify against `?ctModel`/your installed ctsem version
  # before trusting ORDINAL: for a real fit.
  type_keywords <- list(
    BINARY     = 1L,
    ORDINAL    = 2L,
    COUNT      = 3L,
    CENSOR     = 4L,
    CONTINUOUS = 0L
  )
  type_tokens  <- setNames(vector("list", length(type_keywords)), names(type_keywords))
  for (kw in names(type_tokens)) type_tokens[[kw]] <- character(0)
  censor_bounds <- list()  # named by manifest: list(min=, max=)

  seen_latent <- character(0); seen_manifest <- character(0)
  add_latent   <- function(v) if (!v %in% seen_latent)   seen_latent   <<- c(seen_latent, v)
  add_manifest <- function(v) if (!v %in% seen_manifest) seen_manifest <<- c(seen_manifest, v)

  parse_intercept_term <- function(rhs, default_label) {
    rhs <- trimws(rhs)
    if (grepl("\\*", rhs)) {
      pieces <- trimws(strsplit(rhs, "\\*")[[1]])
      label <- pieces[1]; val <- pieces[2]
    } else { label <- NA; val <- rhs }

    if (val == "1")      list(free = TRUE,  label = if (is.na(label)) default_label else label)
    else if (val == "0") list(free = FALSE, label = "0")
    else stop("Intercept-style line must end in 0 or 1 (optionally 'label*1'): got '", rhs, "'")
  }

  parse_rhs <- function(rhs_string, default_fn) {
    terms <- trimws(strsplit(rhs_string, "\\+")[[1]])
    lapply(seq_along(terms), function(i) {
      term <- terms[i]
      if (grepl("\\*", term)) {
        pieces <- trimws(strsplit(term, "\\*")[[1]])
        list(var = pieces[2], label = pieces[1])
      } else list(var = term, label = default_fn(term, i))
    })
  }

  # Shared range-expansion for any TYPE: keyword's token list, e.g.
  # "esm19-esm22" -> c("esm19","esm20","esm21","esm22"). Used by
  # BINARY/ORDINAL/COUNT/CENSOR/CONTINUOUS alike.
  expand_range_token <- function(token) {
    token <- trimws(token)
    if (!grepl("-", token, fixed = TRUE)) return(token)

    parts <- trimws(strsplit(token, "-", fixed = TRUE)[[1]])
    if (length(parts) != 2) return(token)

    m1 <- regmatches(parts[1], regexec("^(.*?)(\\d+)$", parts[1]))[[1]]
    m2 <- regmatches(parts[2], regexec("^(.*?)(\\d+)$", parts[2]))[[1]]

    if (length(m1) != 3 || length(m2) != 3 || m1[2] != m2[2]) return(token)

    prefix <- m1[2]
    lo <- as.integer(m1[3]); hi <- as.integer(m2[3])
    if (lo > hi) stop("Invalid range (start > end): ", token)

    width <- nchar(m1[3])
    paste0(prefix, formatC(lo:hi, width = width, flag = "0"))
  }

  # ---- PASS 1: classify every line ----
  for (line in raw_lines) {

    kw_match <- regmatches(line, regexec("^([A-Za-z]+)\\s*:(.*)$", line))[[1]]
    if (length(kw_match) == 3 && toupper(kw_match[2]) %in% names(type_keywords)) {
      kw  <- toupper(kw_match[2])
      rhs <- kw_match[3]
      raw_toks <- trimws(strsplit(rhs, "[+,]")[[1]])
      raw_toks <- raw_toks[nzchar(raw_toks)]

      for (raw_tok in raw_toks) {
        tok <- raw_tok
        bound_str <- NULL

        # optional "[min=...,max=...]" suffix -- CENSOR: only
        bmatch <- regmatches(raw_tok, regexec("^(.*)\\[(.*)\\]\\s*$", raw_tok))[[1]]
        if (length(bmatch) == 3) {
          tok <- trimws(bmatch[1])
          bound_str <- bmatch[2]
        }

        expanded <- expand_range_token(tok)

        if (!is.null(bound_str)) {
          if (kw != "CENSOR")
            stop("'[min=/max=]' bounds are only valid on CENSOR: tokens -- got on ",
                 kw, ": '", raw_tok, "'")
          if (length(expanded) > 1)
            stop("'[min=/max=]' bounds cannot be combined with a range-expansion ",
                 "token: '", raw_tok, "'")

          b <- list(min = NA_real_, max = NA_real_)
          for (p in trimws(strsplit(bound_str, ",")[[1]])) {
            kvm <- regmatches(p, regexec("^(min|max)\\s*=\\s*(-?Inf|-?[0-9.]+)$", p, ignore.case = TRUE))[[1]]
            if (length(kvm) != 3) stop("Could not parse CENSOR bound '", p, "' in '", raw_tok, "'")
            b[[tolower(kvm[2])]] <- as.numeric(kvm[3])
          }
          censor_bounds[[tok]] <- b
        }

        type_tokens[[kw]] <- c(type_tokens[[kw]], expanded)
      }
      next
    }

    # ---- T0VAR covariance operator: "%~%" ----
    # Deliberately not built from "~" alone -- "~~t0~~" shares a substring
    # with "~t0~" (T0MEANS), which made parsing order-dependent. "%~%"
    # shares no substring with any other operator in this DSL, so this
    # check can sit anywhere in the chain with no ordering risk.
    if (grepl("%~%", line, fixed = TRUE)) {
      parts <- strsplit(line, "%~%", fixed = TRUE)[[1]]
      lhs <- trimws(parts[1]); add_latent(lhs)
      terms <- parse_rhs(parts[2], function(v, i)
        if (v == lhs) paste0("T0var_", lhs) else paste0("T0var_", lhs, "_", v))
      for (t in terms) add_latent(t$var)
      t0cov_lines[[length(t0cov_lines)+1]] <- list(lhs = lhs, terms = terms)
      next
    }

    if (grepl("~t0~", line, fixed = TRUE)) {
      parts <- strsplit(line, "~t0~", fixed = TRUE)[[1]]
      lhs <- trimws(parts[1])
      intercept_lines_raw[[length(intercept_lines_raw)+1]] <-
        list(op = "t0", lhs = lhs, rhs = parts[2])
      next
    }

    if (grepl("=~", line, fixed = TRUE)) {
      parts <- strsplit(line, "=~", fixed = TRUE)[[1]]
      lhs <- trimws(parts[1]); add_latent(lhs)
      terms <- parse_rhs(parts[2], function(v, i) if (i == 1) "1" else paste0("lam_", v))
      for (t in terms) add_manifest(t$var)
      meas_lines[[length(meas_lines)+1]] <- list(lhs = lhs, terms = terms)
      next
    }

    if (grepl("~~", line, fixed = TRUE)) {
      parts <- strsplit(line, "~~", fixed = TRUE)[[1]]
      lhs <- trimws(parts[1]); add_latent(lhs)
      terms <- parse_rhs(parts[2], function(v, i)
        if (v == lhs) paste0("diff_", lhs) else paste0("diff_", lhs, "_", v))
      for (t in terms) add_latent(t$var)
      cov_lines[[length(cov_lines)+1]] <- list(lhs = lhs, terms = terms)
      next
    }

    parts <- strsplit(line, "~", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("Could not parse line (unexpected '~' count): ", line)
    lhs <- trimws(parts[1]); rhs <- trimws(parts[2])
    is_intercept_rhs <- grepl("^(\\S+\\*)?[01]$", rhs)

    if (is_intercept_rhs) {
      intercept_lines_raw[[length(intercept_lines_raw)+1]] <-
        list(op = "cint_or_mm", lhs = lhs, rhs = rhs)
    } else {
      add_latent(lhs)
      terms <- parse_rhs(rhs, function(v, i)
        if (v == lhs) paste0("auto_", lhs) else paste0("dr_", lhs, "_", v))
      for (t in terms) add_latent(t$var)
      struct_lines[[length(struct_lines)+1]] <- list(lhs = lhs, terms = terms)
    }
  }

  if (is.null(latentNames))   latentNames   <- seen_latent
  if (is.null(manifestNames)) manifestNames <- seen_manifest

  # ---- PASS 2: resolve intercept lines now that latent/manifest sets are known ----
  for (l in intercept_lines_raw) {
    if (l$op == "t0") {
      if (!l$lhs %in% latentNames)
        stop("'", l$lhs, " ~t0~ ...' (T0MEANS) requires a latent name. ",
             "Define it via a ~ or =~ line first, or pass latentNames= explicitly.")
      parsed <- parse_intercept_term(l$rhs, paste0("T0m_", l$lhs))
      t0means_lines[[length(t0means_lines)+1]] <-
        list(var = l$lhs, label = if (parsed$free) parsed$label else "0")

    } else {
      if (l$lhs %in% latentNames) {
        parsed <- parse_intercept_term(l$rhs, paste0("cint_", l$lhs))
        cint_lines[[length(cint_lines)+1]] <-
          list(var = l$lhs, label = if (parsed$free) parsed$label else "0")
      } else if (l$lhs %in% manifestNames) {
        parsed <- parse_intercept_term(l$rhs, paste0("mm_", l$lhs))
        mm_lines[[length(mm_lines)+1]] <-
          list(var = l$lhs, label = if (parsed$free) parsed$label else "0")
      } else {
        stop("'", l$lhs, " ~ ", l$rhs, "': unknown variable -- not seen as a latent ",
             "or manifest elsewhere in spec. Define it first, or pass ",
             "latentNames=/manifestNames= explicitly.")
      }
    }
  }

  # ---- validate TYPE: tokens now that manifestNames is known ----
  for (kw in names(type_tokens)) {
    type_tokens[[kw]] <- unique(type_tokens[[kw]])
    unknown <- setdiff(type_tokens[[kw]], manifestNames)
    if (length(unknown) > 0) {
      stop(kw, ": references manifest(s) not seen elsewhere in spec: ",
           paste(unknown, collapse = ", "),
           " -- define them via a =~ line first, or pass manifestNames= explicitly.")
    }
  }

  # a manifest typed by more than one keyword is ambiguous -- refuse rather
  # than silently letting one keyword win (see the zpos/zpos_effect saga in
  # 03_ple_model.qmd for why silent precedence rules are worth avoiding here)
  all_typed <- unlist(type_tokens, use.names = FALSE)
  dup <- unique(all_typed[duplicated(all_typed)])
  if (length(dup) > 0) {
    stop("Manifest(s) assigned more than one type keyword (",
         paste(names(type_keywords), collapse = "/"), "): ",
         paste(dup, collapse = ", "))
  }

  LAMBDA <- matrix("0", length(manifestNames), length(latentNames),
                   dimnames = list(manifestNames, latentNames))
  for (b in meas_lines) for (t in b$terms) LAMBDA[t$var, b$lhs] <- t$label

  DRIFT <- matrix("0", length(latentNames), length(latentNames),
                  dimnames = list(latentNames, latentNames))
  for (b in struct_lines) for (t in b$terms) DRIFT[b$lhs, t$var] <- t$label

  DIFFUSION <- matrix("0", length(latentNames), length(latentNames),
                      dimnames = list(latentNames, latentNames))
  for (b in cov_lines) for (t in b$terms) {
    i <- match(b$lhs, latentNames)
    j <- match(t$var, latentNames)
    row <- max(i, j); col <- min(i, j)
    DIFFUSION[row, col] <- t$label
  }

  # ---- build T0VAR (same lower-triangle routing as DIFFUSION above).
  #      Any cell not mentioned via "%~%" defaults to fixed "0" -- same
  #      convention as DIFFUSION. Diagonal (variance) cells are NOT free
  #      by default here, unlike ctsem's own ctModel() default -- if you
  #      want a latent's T0 variance estimated, say so explicitly,
  #      e.g. "cog %~% cog".
  T0VAR <- matrix("0", length(latentNames), length(latentNames),
                  dimnames = list(latentNames, latentNames))
  for (b in t0cov_lines) for (t in b$terms) {
    i <- match(b$lhs, latentNames)
    j <- match(t$var, latentNames)
    row <- max(i, j); col <- min(i, j)
    T0VAR[row, col] <- t$label
  }

  CINT <- setNames(rep("0", length(latentNames)), latentNames)
  for (l in cint_lines) CINT[l$var] <- l$label

  T0MEANS <- setNames(rep("0", length(latentNames)), latentNames)
  for (l in t0means_lines) T0MEANS[l$var] <- l$label

  MANIFESTMEANS <- setNames(rep("0", length(manifestNames)), manifestNames)
  for (l in mm_lines) MANIFESTMEANS[l$var] <- l$label

  manifesttype <- setNames(rep(0L, length(manifestNames)), manifestNames)
  for (kw in names(type_keywords)) {
    manifesttype[type_tokens[[kw]]] <- type_keywords[[kw]]
  }

  # ---- CENSOR: bounds -- default -Inf/Inf (no censoring) for every
  #      manifest, floor (and optionally ceiling) filled in for CENSOR:
  #      tokens. A bound can be given explicitly in the DSL
  #      ("CENSOR: var[min=2]"), or left to be computed from `data` (the
  #      observed minimum, na.rm=TRUE -- this project's floor-effect
  #      convention throughout: censormax stays Inf unless [max=...] is
  #      given explicitly).
  censormin <- setNames(rep(-Inf, length(manifestNames)), manifestNames)
  censormax <- setNames(rep(Inf,  length(manifestNames)), manifestNames)

  for (v in type_tokens$CENSOR) {
    b <- censor_bounds[[v]]

    if (!is.null(b) && !is.na(b$min)) {
      censormin[v] <- b$min
    } else if (!is.null(data)) {
      censormin[v] <- min(data[[v]], na.rm = TRUE)
    } else {
      stop("CENSOR: '", v, "' has no explicit '[min=...]' bound and no `data=` ",
           "was passed to ct_lavaan() to auto-compute its floor from the observed ",
           "minimum. Either pass data = <your data frame>, or write ",
           "'CENSOR: ", v, "[min=<value>]' explicitly.")
    }

    if (!is.null(b) && !is.na(b$max)) censormax[v] <- b$max
    # else stays Inf -- floor-only censoring, matching every use in this project
  }

  list(
    LAMBDA        = LAMBDA,
    DRIFT         = DRIFT,
    DIFFUSION     = DIFFUSION,
    T0VAR         = T0VAR,
    CINT          = unname(CINT),
    T0MEANS       = unname(T0MEANS),
    MANIFESTMEANS = unname(MANIFESTMEANS),
    manifesttype  = unname(manifesttype),
    censormin     = unname(censormin),
    censormax     = unname(censormax),
    latentNames   = latentNames,
    manifestNames = manifestNames
  )
}
