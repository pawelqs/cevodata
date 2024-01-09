
#' Not in operator
#' @param x left-hand side argument
#' @param y right-hand side argument
'%not in%' <- function(x,y)!('%in%'(x,y))


join_aes <- function(aes_default, aes_2) {
  aes <- c(as.list(aes_default[names(aes_default) %not in% names(aes_2)]), aes_2)
  class(aes) <- 'uneval'
  aes
}


drop_na_columns <- function(.data) {
  .data |>
    keep(~all(!is.na(.x)))
}


get_f_range <- function(snvs, pct_left = 0.05, pct_right = 0.95) {
  bounds <- snvs |>
    filter(.data$f > 0.00001, !is.na(.data$f)) |>
    group_by(.data$sample_id) |>
    summarise(
      lower_bound = stats::quantile(.data$f, pct_left),
      higher_bound = stats::quantile(.data$f, pct_right)
    )
  bounds
}


msg <- function(...,
                collapse = "",
                col = "steelblue3",
                new_line = TRUE,
                verbose = get_cevomod_verbosity()) {
  msg <- str_c(list(...), collapse = collapse)
  if (verbose && new_line) {
    cli::cat_line(msg, col = col)
  } else if (verbose) {
   cat(crayon::blue(msg))
  }
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


#' @export
print.cevo_snvs <- function(x, ...) {
  msg("<cevo_snvs> tibble")
  NextMethod()
}


#' Fill na values in the object
#' @param object object
#' @param val value to fill the NAs
#' @export
fill_na <- function(object, val) {
  object[is.na(object)] <- val
  object
}


segment <- function(vec) {
  x <- vec != lag(vec)
  x[1] <- 0
  cumsum(x)
}


#' Quick save to ~/.cevomod directory
#' @param object object to save
#' @export
quick_save <- function(object) {
  dir.create("~/.cevodata", showWarnings = FALSE)
  write_rds(object, "~/.cevodata/object.Rds")
}


#' @describeIn quick_save Quick load of ~/.cevomod/object.Rds
#' @export
quick_load <- function() {
  read_rds("~/.cevodata/object.Rds")
}