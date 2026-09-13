ct_lavaan <- function(spec, latentNames = NULL, manifestNames = NULL) {
  
  raw_lines <- sub("#.*$", "", strsplit(spec, "\n")[[1]])
  raw_lines <- trimws(raw_lines)
  raw_lines <- raw_lines[nzchar(raw_lines)]
  
  meas_lines <- list(); struct_lines <- list(); cov_lines <- list()
  t0cov_lines <- list()
  cint_lines <- list(); t0means_lines <- list(); mm_lines <- list()
  binary_tokens <- character(0)
  intercept_lines_raw <- list()
  
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
  
  expand_binary_token <- function(token) {
    token <- trimws(token)
    if (!grepl("-", token, fixed = TRUE)) return(token)
    
    parts <- trimws(strsplit(token, "-", fixed = TRUE)[[1]])
    if (length(parts) != 2) return(token)
    
    m1 <- regmatches(parts[1], regexec("^(.*?)(\\d+)$", parts[1]))[[1]]
    m2 <- regmatches(parts[2], regexec("^(.*?)(\\d+)$", parts[2]))[[1]]
    
    if (length(m1) != 3 || length(m2) != 3 || m1[2] != m2[2]) return(token)
    
    prefix <- m1[2]
    lo <- as.integer(m1[3]); hi <- as.integer(m2[3])
    if (lo > hi) stop("Invalid BINARY range (start > end): ", token)
    
    width <- nchar(m1[3])
    paste0(prefix, formatC(lo:hi, width = width, flag = "0"))
  }
  
  # ---- PASS 1: classify every line ----
  for (line in raw_lines) {
    
    if (grepl("^BINARY\\s*:", line, ignore.case = TRUE)) {
      rhs <- sub("^BINARY\\s*:\\s*", "", line, ignore.case = TRUE)
      toks <- trimws(strsplit(rhs, "[+,]")[[1]])
      toks <- toks[nzchar(toks)]
      binary_tokens <- c(binary_tokens, toks)
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
  
  binary_manifests <- unique(unlist(lapply(binary_tokens, expand_binary_token)))
  unknown_binary <- setdiff(binary_manifests, manifestNames)
  if (length(unknown_binary) > 0) {
    stop("BINARY references manifest(s) not seen elsewhere in spec: ",
         paste(unknown_binary, collapse = ", "),
         " -- define them via a =~ line first, or pass manifestNames= explicitly.")
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
  manifesttype[binary_manifests] <- 1L
  
  list(
    LAMBDA        = LAMBDA,
    DRIFT         = DRIFT,
    DIFFUSION     = DIFFUSION,
    T0VAR         = T0VAR,
    CINT          = unname(CINT),
    T0MEANS       = unname(T0MEANS),
    MANIFESTMEANS = unname(MANIFESTMEANS),
    manifesttype  = unname(manifesttype),
    latentNames   = latentNames,
    manifestNames = manifestNames
  )
}