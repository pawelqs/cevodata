
# ----------------------------- Filter -----------------------------------------

#' @export
dplyr::filter


#' Filter/subset cevodata object
#'
#' This is a wrapper around dplyr::filter function which can be used to subset
#' cevodata object. Works like dplyr::filter, performs the filtering on metadata,
#' then filters SNVs, CNAs, clones and models keeping samples kept in metadata
#'
#' @param .data cevodata object
#' @param ... expression passed to dplyr::filter(...)
#' @inheritParams dplyr::filter
#' @return cevodata object
#' @export
filter.cevodata <- function(.data, ..., .preserve = FALSE) {
  cd <- .data
  cd$metadata <- get_metadata(cd) |>
    filter(...)
  id_cols <- c("patient_id", "sample_id")
  keep <- cd$metadata |>
    select(any_of(id_cols))

  # Lists of tibbles
  cd$SNVs  <- map(cd$SNVs,  \(x) semi_join(x, keep, by = "sample_id"))
  cd$CNAs  <- map(cd$CNAs,  \(x) semi_join(x, keep, by = "sample_id"))
  cd$stats <- map(cd$stats, \(x) semi_join(x, keep, by = intersect(names(x), id_cols)))

  # Lists of cv_subitems
  cd$models <- map(cd$models, \(x) semi_join(x, keep))
  cd$misc   <- map(cd$misc,   \(x) semi_join(x, keep))

  cd
}


# ------------------------------ Split -----------------------------------------

#' Split object
#' @param object object to split
#' @param ... other arguments
#' @export
split_by <- function(object, ...) {
  UseMethod("split_by")
}


#' @describeIn split_by Split cevodata object
#'
#' This function is based on `filter.cevodata()`, thus might be inefficient
#' if there are too many levels in the var column.
#'
#' @param object cevodata object
#' @param var name of column in metadata
#' @export
split_by.cevodata <- function(object, var, ...) {
  split_names <- object$metadata |>
    pull({{var}}) |>
    unique()
  splits <- split_names |>
    set_names(split_names) |>
    map(~filter(object, {{var}} == .x))
  class(splits) <- c("cevo_splits", "list")
  splits
}


# ----------------------------- Merge -----------------------------------------

#' Merge two cevodata objects
#' @inheritParams base::merge
#' @param name Name of the merged object
#' @param verbose Show messages?
#' @param .id datasets names will be saved to this metadata column, if provided
#' @export
merge.cevodata <- function(x, y,
                           ...,
                           name = "Merged datasets",
                           verbose = get_verbosity(),
                           .id = NULL) {
  if (!is.null(.id)) {
    x$metadata[[.id]] <- x$name
    y$metadata[[.id]] <- y$name
  }
  meta <- bind_rows(x$metadata, y$metadata)
  cd <- init_cevodata(name = name) |>
    add_metadata(meta)

  # Lists of tibbles
  cd$SNVs  <- merge_items(x$SNVs, y$SNVs)
  cd$CNAs  <- merge_items(x$CNAs, y$CNAs)
  cd$stats <- merge_items(x$stats, y$stats)

  # Lists of cv_subitems
  cd$models <- merge_items(x$models, y$models)
  cd$misc   <- merge_items(x$misc, y$misc)

  msg("Default SNVs: ", default_SNVs(x), verbose = verbose)
  msg("Default CNAs: ", default_CNAs(x), verbose = verbose)
  cd$settings <- x$settings
  cd
}


merge_items <- function(x, y) {
  if (length(x) == 0 && length(y) == 0) {
    return(list())
  } else if (length(x) == 0) {
    return(y)
  } else if (length(y) == 0) {
    return(x)
  }

  if ("tbl_df" %in% class(x[[1]])) {
    fun <- bind_rows
  } else if ("cv_subitem" %in% class(x[[1]])) {
    fun <- merge
  } else {
    stop("Unknown item type")
  }

  item_names <- union(names(x), names(y))
  res <- vector("list", length(item_names))
  names(res) <- item_names
  for (i in item_names) {
    if (i %in% names(x) && i %in% names(y)) {
      res[[i]] <- fun(x[[i]], y[[i]])
    } else if (i %in% names(y)) {
      res[[i]] <- y[[i]]
    } else {
      res[[i]] <- x[[i]]
    }
  }
  res
}


# ----------------------------- Update -----------------------------------------

#' Update cevodata object with values from another object
#' @param object object to update
#' @param object2 object to use
#' @param ... other args, unused now
#' @export
update.cevodata <- function(object, object2, ...) {
  sample_ids <- union(object$metadata$sample_id, object2$metadata$sample_id)
  object <- object |>
    filter(.data$sample_id %not in% object2$metadata$sample_id) |>
    merge(object2)
  object$metadata <- object$metadata |>
    mutate(sample_id = parse_factor(.data$sample_id, levels = sample_ids)) |>
    arrange(.data$sample_id) |>
    mutate(across("sample_id", as.character))
  object
}
