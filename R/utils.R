
#' Not in operator
#' @param x left-hand side argument
#' @param y right-hand side argument
'%not in%' <- function(x,y)!('%in%'(x,y))


factorize <- function(tbl, col, levels = unique(tbl[["col"]])) {
  tbl[[col]] <- parse_factor(tbl[[col]], levels = levels)
  tbl
}


get_verbosity <- function() {
  v <- verbose::verbose("cevoverse")
  if (is.null(v)) {
    0
  } else {
    v
  }
}


join_aes <- function(aes_default, aes_2) {
  aes <- c(as.list(aes_default[names(aes_default) %not in% names(aes_2)]), aes_2)
  class(aes) <- 'uneval'
  aes
}


msg <- function(...,
                collapse = "",
                col = "steelblue3",
                new_line = TRUE,
                verbose = get_verbosity()) {
  msg <- str_c(list(...), collapse = collapse)
  if (verbose && new_line) {
    cli::cat_line(msg, col = col)
  } else if (verbose) {
   cat(crayon::blue(msg))
  }
}


#' @export
print.cevo_snvs <- function(x, ...) {
  msg("<cevo_snvs> tibble")
  NextMethod()
}


#' @describeIn quick_save Quick load of ~/.cevomod/object.Rds
#' @export
quick_load <- function() {
  read_rds("~/.cevodata/object.Rds")
}


#' Quick save to ~/.cevomod directory
#' @param object object to save
#' @export
quick_save <- function(object) {
  dir.create("~/.cevodata", showWarnings = FALSE)
  write_rds(object, "~/.cevodata/object.Rds")
}


require_packages <- function(...) {
  pkgs <- list(...)
  missing <- !map_lgl(pkgs, requireNamespace, quietly = TRUE)
  if (any(missing)) {
    stop(
      paste0("Package '", pkgs[missing], "' must be installed to use this function.\n"),
      call. = FALSE
    )
  }
}


require_columns <- function(tbl, ...) {
  cols <- list(...) |>
    unlist()
  tbl_name <- deparse(substitute(tbl))
  missing <- cols %not in% colnames(tbl)
  if (sum(missing) > 0) {
    stop("The following columns are missing in the ", tbl_name, ": ", str_c(cols[missing], collapse = ", "))
  }
}


segment <- function(vec) {
  x <- vec != lag(vec)
  x[1] <- 0
  cumsum(x)
}


verbose_down <- function(verbose) {
  if (isTRUE(verbose) || isFALSE(verbose) || verbose == 0) {
    FALSE
  } else if (is.numeric(verbose) && verbose > 0) {
    verbose - 1
  } else {
    stop("Verbose should be logical or positive")
  }
}
