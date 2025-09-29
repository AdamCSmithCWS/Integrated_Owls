
check_range <- function(arg, range) {
  if(!all(arg >= range[1] & arg <= range[2])) {
    stop("`", deparse(substitute(arg)),"` ",
         "must range between ", range[1], " and ", range[2],
         call. = FALSE)
  }
}


check_logical <- function(..., allow_null = FALSE) {
  check_type(msg = "logical (TRUE/FALSE)", check_fun = is.logical, allow_null,
             ...)
}


check_numeric <- function(..., allow_null = FALSE) {
  check_type(msg = "a number", check_fun = is.numeric, allow_null, ...)
}

check_in <- function(arg, opts) {
  if(!arg %in% opts) {
    if(is.character(opts)) sep <- "'" else if(is.numeric(opts)) sep <- ""
    stop("`", deparse(substitute(arg)),"` ",
         "must be one of ", sep, paste0(opts, collapse = paste0(sep, ", ", sep)),
         sep, ".",
         call. = FALSE)
  }
}



check_type <- function(msg, check_fun, allow_null, ...) {
  args <- list(...)
  if(is.null(names(args))) {
    names(args) <- vapply(substitute(list(...))[-1], deparse, FUN.VALUE = "a")
  }
  ck <- vapply(args, check_fun, FUN.VALUE = TRUE)
  if(allow_null) ck <- ck | vapply(args, is.null, FUN.VALUE = TRUE)
  if(!all(ck)) {
    stop("`", paste0(names(ck[!ck]), collapse = "`, `"),
         "` must be ", msg, call. = FALSE)
  }
}
