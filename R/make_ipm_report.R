#' @title Generate an RMarkdown file with IPM metadata
#' @rdname ipm_report
#'
#' @description Generates a \code{.rmd} file containing a mathematical description
#' of the \code{proto_ipm} object.
#'
#' @param proto_ipm A proto_ipm object
#' @param rmd_dest  The folder to save the Rmd file at. The default is
#' \code{getwd()}. Alternatively, can be a complete file path that specifies
#' the location and title of the document with the extension \code{".rmd"}. in
#' this case, the current date will be appended to the title.
#' @param title The title to include in the document. This is not necessarily
#' the same as \code{rmd_dest}, as this appears at the top of the generated
#' report, and is not included in the file path!
#' @param output_format The format to include in the YAML header for the created
#' \code{.rmd} document.
#' @param render_output A logical indicating whether to call
#' \code{rmarkdown::render} on the generated \code{.rmd} file. Often times, the
#' \code{.rmd} file will need further editing before it's useful, so the default
#' is \code{FALSE}.
#'
#' @details \code{make_ipm_report_body} only translates the iteration
#'   expressions and vital rate expressions into Markdown with LaTeX, and does
#'   not produce any headers needed to knit the file. This function is exported
#'   mostly for re-usage in \code{\link[Rpadrino]{pdb_report}}, and isn't really
#'   intended for use by \code{ipmr} users.
#'
#' @return For \code{make_ipm_report}, the filepath to the \code{.rmd} file. The
#' default name is \code{ "ipmr_report_<current_date>.rmd"}. For
#' \code{make_ipm_report_body}, a character vector with Markdown and LaTeX
#' suitable for rendering, but without a header.
#'
#'
#'
#' @export


make_ipm_report <- function(proto_ipm,
                            rmd_dest      = getwd(),
                            title         = "",
                            output_format = "html",
                            render_output = FALSE) {

  if(!requireNamespace("rmarkdown", quietly = TRUE) && render_output) {
    stop("The 'rmarkdown' package is required for 'render_output = TRUE'.\n",
         "Please install it - install.packages('rmarkdown').",
         call. = FALSE)
  }

  # Create output file path from chosen destination folder, then yaml header,
  # then the actual report body. body contains translated equations and
  # parameter value lists. This doesn't try to insert species names or anything,
  # but rather punts to the user on that front. Rpadrino will handle automatic
  # insertions separately.

  rmd_dest     <- .make_ipm_report_fp(rmd_dest, keep_rmd = TRUE)
  header       <- .make_ipm_report_header(output_format, title)
  body         <- make_ipm_report_body(proto_ipm)

  rmd_contents <- c(header, body)

  fp <- file.create(rmd_dest, showWarnings = FALSE)

  writeLines(rmd_contents, con = rmd_dest)

  if(render_output) rmarkdown::render(rmd_dest)

  return(rmd_dest)

}

#' @rdname ipm_report
#' @export

make_ipm_report_body <- function(proto_ipm) {

  iter_exprs <- .make_ipm_report_iter_exprs(proto_ipm)
  vr_exprs   <- .make_ipm_report_vr_exprs(proto_ipm)
  all_params <- .make_ipm_report_params(proto_ipm)
  impl_args  <- .make_ipm_report_impl_args(proto_ipm)
  glossary   <- .make_ipm_report_glossary(proto_ipm)

  out <- c(iter_exprs, vr_exprs, all_params, impl_args, glossary)

  return(out)

}

.make_ipm_report_header <- function(output_format) {

  out <- paste("---",
               paste0("title: '", title, "'"),
               paste0("output: ", output_format, "_document\n"),
               paste0("date: '`r Sys.Date()`'"),
               paste0("urlcolor: blue"),
               "---\n",
               sep = "\n")

}

#' @noRd
# Carbon copy of Rpadrino:::pdb_rmd_dest, except it checks for 'tools'
# installation so that tools can remain in Suggests.

.make_ipm_report_fp <- function(rmd_dest, keep_rmd) {

  if(!requireNamespace("tools", quietly = TRUE)) {
    stop("The 'tools' package is required for 'make_ipm_report()'.\n",
         "Please install it - install.packages('tools').",
         call. = FALSE)
  }

  date <- gsub("-", "", Sys.Date())

  # Non-specified output directory gets a tempdir()
  if((is.null(rmd_dest) || is.na(rmd_dest) || rmd_dest == "") && keep_rmd) {

    rmd_dest <- tempfile(pattern = paste0("Rpadrino_report_", date),
                         fileext = ".Rmd")

    message("'keep_rmd = TRUE' and 'rmd_dest' is not specified! ",
            "Saving to a temporary file: \n", rmd_dest)

  } else if(tools::file_ext(rmd_dest) == "") {

    rmd_dest <- paste0(rmd_dest, "/ipmr_report_", date, ".Rmd")

    file.create(rmd_dest, showWarnings = FALSE)

  } else if(tools::file_ext(tolower(rmd_dest)) == "rmd") {

    rmd_dest <- gsub("\\.rmd$", paste0("_", date, ".Rmd"),
                     rmd_dest,
                     ignore.case = TRUE)

    file.create(rmd_dest, showWarnings = FALSE)

  } else {

    stop("'rmd_dest' must either be 'NULL', the name of a folder, or a file with",
         " a .rmd file extension!")
  }

  rmd_dest
}
