
# ------------------------- cv_subitem class -----------------------------------

new_cv_subitem <- function(x) {
  classes <- c("cv_subitem", class(x)) |>
    unique()
  structure(x, class = classes)
}


as_cv_subitem <- function(x) {
  validate_cv_subitem(new_cv_subitem(x))
}


validate_cv_subitem <- function(cv_subitem) {
  stopifnot(is.list(cv_subitem))
  cv_subitem
}


# --------------------------- core methods  ------------------------------------

#' @export
semi_join.cv_subitem <- function(x, y, by = NULL, copy = FALSE, ...) {
  keep <- c("settings", "info")
  tibbles_to_filter <- setdiff(names(x), keep)

  for (i in tibbles_to_filter) {
    join_by_vars <- c("patient_id", "sample_id") |>
      intersect(names(y)) |>
      intersect(names(x[[i]]))
    x[[i]] <- semi_join(x[[i]], y, by = join_by_vars)
  }
  x
}


#' @export
merge.cv_subitem <- function(x, y, ...) {
  const <- intersect(c("settings", "info"), names(x))
  tibbles_to_merge <- names(x) |>
    union(names(y)) |>
    setdiff(const)

  res <- vector("list", length(tibbles_to_merge))
  names(res) <- tibbles_to_merge
  for (i in tibbles_to_merge) {
    res[[i]] <- bind_rows(x[[i]], y[[i]])
  }
  res[const] <- x[const]

  as_cv_subitem(res)
}

