# cust_funcs/ct_latex_html.R
#
# Renders a ctsem model/fit's LaTeX equation (via ctsem::ctModelLatex()) as a
# genuine vector SVG that shows up both when you run the chunk by hand (in
# the RStudio Viewer) and when the .qmd is rendered to HTML (embedded
# directly as inline SVG markup) -- crisp at any zoom, no PNG rasterizing.
# HTML has no native LaTeX rendering, and the equations ctModelLatex()
# builds (parbox, flalign*, custom \vect{} macros, underbrace stacks) are
# not things MathJax/KaTeX can reliably reproduce, so this leans on
# ctModelLatex()'s own working pdflatex/tinytex compilation rather than
# hand-translating the LaTeX.
#
# A real .pdf is compiled first (pdflatex only ever produces a PDF, and
# plain `latex` -> .dvi was tried and rejected: it uses a different font
# encoding than pdflatex and produced garbled glyphs for this document).
# The PDF is then converted to SVG and immediately deleted -- it's a
# means-to-an-end intermediate, never left on disk or shown to the user.
# The PDF -> SVG conversion tries, in order:
#   1. pdftocairo -svg  (part of poppler-utils; very likely present already,
#      since the pdftools R package links against the same poppler library)
#   2. dvisvgm --pdf    (bundled with TinyTeX, but its PDF backend depends
#      on being linked against a PDF library and can be broken/missing on
#      some TinyTeX builds -- kept as a fallback, not the primary path)
# Every intermediate file lives in a tempfile() folder deleted before the
# function returns.
#
# Usage inside a .qmd chunk:
#   source("cust_funcs/ct_latex_html.R")
#   ctModelLatexHTML(psx_twopart_mod)
# or, for a fitted model (population-level equation with TIpred effects):
#   ctModelLatexHTML(psx_twopart_fit)
# with extra top margin (default 1cm):
#   ctModelLatexHTML(psx_twopart_fit, margin_cm = 1.5)

# ctsem's own LaTeX template uses the "preview"/tightpage mechanism, which
# crops the compiled output (and therefore the SVG made from it) exactly to
# the ink bounding box -- there's no built-in margin to configure upstream.
# This pads the SVG after the fact: bumps its height/viewBox and shifts all
# existing content down inside a wrapping translate() group.
.ct_svg_add_top_margin <- function(svg_text, margin_cm = 1) {
  margin_pt <- margin_cm * 72 / 2.54

  parts <- regmatches(svg_text,
    regexec("(?s)(<svg\\b[^>]*>)(.*)(</svg>\\s*)$", svg_text, perl = TRUE))[[1]]
  if (length(parts) != 4) return(svg_text)  # unrecognized structure -- leave untouched
  open_tag <- parts[2]; body <- parts[3]; close_tag <- parts[4]

  hmatch <- regmatches(open_tag, regexec('height="([0-9.]+)([a-z%]*)"', open_tag))[[1]]
  if (length(hmatch) == 3) {
    new_h <- as.numeric(hmatch[2]) + margin_pt
    open_tag <- sub('height="[0-9.]+[a-z%]*"',
      paste0('height="', format(new_h, trim = TRUE), hmatch[3], '"'), open_tag)
  }

  vmatch <- regmatches(open_tag, regexec(
    'viewBox="([0-9.eE+-]+) ([0-9.eE+-]+) ([0-9.eE+-]+) ([0-9.eE+-]+)"', open_tag))[[1]]
  if (length(vmatch) == 5) {
    new_vh <- as.numeric(vmatch[5]) + margin_pt
    open_tag <- sub('viewBox="[0-9.eE+-]+ [0-9.eE+-]+ [0-9.eE+-]+ [0-9.eE+-]+"',
      paste0('viewBox="', vmatch[2], " ", vmatch[3], " ", vmatch[4], " ",
             format(new_vh, trim = TRUE), '"'), open_tag)
  }

  body <- paste0('<g transform="translate(0,', format(margin_pt, trim = TRUE), ')">',
                  body, "</g>")

  paste0(open_tag, body, close_tag)
}

# Namespaces every id="..." in a compiled SVG and rewrites every reference
# to it (href="#id", xlink:href="#id", url(#id)) so that multiple
# independently-compiled SVGs -- e.g. two ctModelLatexHTML() calls on the
# same rendered .qmd page -- can be embedded inline in the same HTML
# document without their ids colliding.
#
# pdftocairo/dvisvgm number font glyphs (id="glyph-0-0", "glyph-0-1", ...)
# fresh on every single compile, starting from whatever character happens
# to appear first *in that one document*. SVG/HTML ids are global to the
# whole page, not scoped to their own <svg> -- so if two such SVGs land on
# the same page, a later <use xlink:href="#glyph-0-0"> silently resolves to
# the *first* SVG's glyph-0-0 (same position, wrong glyph shape) instead of
# its own. That's what turned real parameter values into a mix of correct
# digits and stray letters once the model and fit equations were both on
# 03_ple_model.qmd -- a lone SVG (e.g. the standalone popup this function
# opens when run outside knitr) never collides with anything, which is why
# it always looked fine locally.
.ct_svg_namespace_ids <- function(svg_text, prefix) {
  ids <- unique(regmatches(svg_text,
    gregexpr('(?<=id=")[^"]+', svg_text, perl = TRUE))[[1]])
  for (id in ids) {
    esc <- gsub("([.*+?^${}()|\\[\\]\\\\])", "\\\\\\1", id, perl = TRUE)
    repl <- paste0(prefix, "-", id)
    svg_text <- gsub(paste0('id="', esc, '"'),
      paste0('id="', repl, '"'), svg_text, perl = TRUE)
    svg_text <- gsub(paste0('(xlink:href|href)="#', esc, '"'),
      paste0('\\1="#', repl, '"'), svg_text, perl = TRUE)
    svg_text <- gsub(paste0('url\\(#', esc, '\\)'),
      paste0('url(#', repl, ')'), svg_text, perl = TRUE)
  }
  svg_text
}

# Resolves a bundled binary's absolute path, trying Sys.which() first and
# falling back to a direct search of TinyTeX's own bin directory (dvisvgm
# only -- pdftocairo isn't a TeX Live tool, so no such fallback for it).
.ct_find_tex_bin <- function(name, tinytex_fallback = FALSE) {
  bin <- unname(Sys.which(name))
  if (!identical(bin, "")) return(bin)
  if (tinytex_fallback && requireNamespace("tinytex", quietly = TRUE)) {
    root <- tryCatch(tinytex::tinytex_root(), error = function(e) NULL)
    if (!is.null(root)) {
      candidates <- Sys.glob(file.path(root, "bin", "*", name))
      if (length(candidates) > 0) return(candidates[1])
    }
  }
  ""
}

ctModelLatexHTML <- function(x, margin_cm = 1, ...) {

  if (!requireNamespace("ctsem", quietly = TRUE)) {
    stop("ctsem package required.", call. = FALSE)
  }
  has_pdflatex <- !identical(Sys.which("pdflatex"), "")
  has_tinytex  <- requireNamespace("tinytex", quietly = TRUE)
  if (!has_pdflatex && !has_tinytex) {
    stop("No LaTeX compiler found (pdflatex not on PATH, tinytex not installed). ",
         "Install one, e.g.: tinytex::install_tinytex()", call. = FALSE)
  }

  pdftocairo_bin <- .ct_find_tex_bin("pdftocairo")
  dvisvgm_bin    <- .ct_find_tex_bin("dvisvgm", tinytex_fallback = TRUE)
  if (identical(pdftocairo_bin, "") && identical(dvisvgm_bin, "")) {
    stop("No PDF-to-SVG converter found (tried pdftocairo and dvisvgm). ",
         "Install poppler (e.g. `brew install poppler` on macOS) or ",
         "run tinytex::tlmgr_install('dvisvgm').", call. = FALSE)
  }

  folder <- tempfile("ctsemEq_")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  filename <- "eq"

  # compile to PDF only -- ctModelLatex()'s own popup logic (PNG-based,
  # gated on interactive()) is skipped entirely; we handle both the popup
  # and the render case ourselves below, from the SVG. ctModelLatex() (and
  # the tinytex/pdflatex compile it triggers) both print a fair bit of
  # noise straight to the console/output -- a "Computing quantities..."
  # message() plus the full pdflatex compile log -- none of which is
  # useful here, so both are swallowed.
  invisible(suppressMessages(utils::capture.output(
    ctsem::ctModelLatex(
      x,
      folder   = folder,
      filename = filename,
      tex      = TRUE,
      compile  = TRUE,
      open     = FALSE,
      savepng  = FALSE,
      ...
    )
  )))

  pdf_path <- file.path(folder, paste0(filename, ".pdf"))
  if (!file.exists(pdf_path)) {
    stop("LaTeX did not compile -- check pdflatex/tinytex output above.",
         call. = FALSE)
  }

  svg_path <- file.path(folder, paste0(filename, ".svg"))
  conv_errs <- character(0)

  if (!identical(pdftocairo_bin, "")) {
    status <- tryCatch(
      system2(pdftocairo_bin,
        c("-svg", shQuote(pdf_path), shQuote(svg_path)),
        stdout = TRUE, stderr = TRUE),
      error = function(e) { conv_errs[[1]] <<- conditionMessage(e); NULL })
    if (!file.exists(svg_path)) conv_errs <- c(conv_errs, paste("pdftocairo:", paste(status, collapse = " ")))
  }

  if (!file.exists(svg_path) && !identical(dvisvgm_bin, "")) {
    status <- tryCatch(
      system2(dvisvgm_bin,
        c("--pdf", shQuote(pdf_path), "-o", shQuote(svg_path)),
        stdout = TRUE, stderr = TRUE),
      error = function(e) { conv_errs[[length(conv_errs) + 1]] <<- conditionMessage(e); NULL })
    if (!file.exists(svg_path)) conv_errs <- c(conv_errs, paste("dvisvgm:", paste(status, collapse = " ")))
  }

  # the PDF was only ever a means to an end -- gone as soon as we have (or
  # have failed to get) an SVG from it.
  unlink(pdf_path)

  if (!file.exists(svg_path)) {
    stop("Could not convert the equation PDF to SVG:\n", paste(conv_errs, collapse = "\n"),
         call. = FALSE)
  }

  svg_text <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  if (margin_cm > 0) svg_text <- .ct_svg_add_top_margin(svg_text, margin_cm)
  svg_text <- .ct_svg_namespace_ids(svg_text, basename(folder))

  if (isTRUE(getOption("knitr.in.progress"))) {
    # inside a knitr/Quarto render -- emit the raw SVG markup directly so
    # the browser renders it natively (crisp at any zoom, no chunk-figure
    # PNG capture involved). asis_output() means results: asis isn't
    # required on the chunk, but it doesn't hurt to set it explicitly.
    return(knitr::asis_output(svg_text))
  }

  # run by hand -- pop the SVG up in the RStudio Viewer (or the system
  # browser if the Viewer isn't available), via a tiny wrapper HTML page.
  tmp_html <- tempfile(fileext = ".html")
  writeLines(c('<html><body style="margin:0">', svg_text, '</body></html>'), tmp_html)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    rstudioapi::viewer(tmp_html)
  } else {
    utils::browseURL(tmp_html)
  }

  invisible(NULL)
}
