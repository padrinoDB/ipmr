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
#' @param translate_greek A logical. \code{TRUE} means the function will attempt
#' to translate symbols in vital rate expressions to Greek equivalents, which
#' can make reports prettier. \code{FALSE} will leave the expressions largely
#' as is before translating them to LaTeX equivalents.
#'
#' @details \code{make_ipm_report_body} only translates the iteration
#'   expressions and vital rate expressions into Markdown with LaTeX, and does
#'   not produce any headers needed to knit the file. This function is exported
#'   mostly for re-usage in \code{\link[Rpadrino]{pdb_report}}, and isn't really
#'   intended for use by \code{ipmr} users.
#'
#'   @section \strong{Translations}
#'
#'   For iteration expressions, vital rate expressions, and parameter names,
#'   \code{make_ipm_report} first attempts to translate functional forms into a
#'   common set of names. For example, \code{s = surv_int + surv_slope * z_1} is
#'   translated into \code{beta_0_s + beta_1_s * z_1}, and then is
#'   translated into LaTeX equations. This helps create useful Greek symbols
#'   more consistently. This may cause strange/incorrect results when parsing
#'   non-linear models.
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
                            rmd_dest        = getwd(),
                            title           = "",
                            output_format   = "html",
                            render_output   = FALSE,
                            translate_greek = TRUE) {

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
  body         <- make_ipm_report_body(proto_ipm, translate_greek)

  rmd_contents <- c(header, body)

  fp <- file.create(rmd_dest, showWarnings = FALSE)

  writeLines(rmd_contents, con = rmd_dest)

  if(render_output) rmarkdown::render(rmd_dest)

  return(rmd_dest)

}

#' @rdname ipm_report
#' @export

make_ipm_report_body <- function(proto_ipm, translate_greek) {

  iter_exprs <- .make_ipm_report_iter_exprs(proto_ipm, translate_greek)
  vr_exprs   <- .make_ipm_report_vr_exprs(proto_ipm,   translate_greek)
  all_params <- .make_ipm_report_params(proto_ipm,     translate_greek)
  impl_args  <- .make_ipm_report_impl_args(proto_ipm,  translate_greek)
  glossary   <- .make_ipm_report_glossary(proto_ipm,   translate_greek)

  out <- c(iter_exprs, vr_exprs, all_params, impl_args, glossary)

  return(out)

}

#' @noRd
# Translates:
# n_z_t_1 = P %*% n_z_t + F %*% n_z_t -> n(z', t + 1) = \int([P+F]n(z , t)dz)
# Use CC/DC to decide if \int applies!

.pdb_make_ipm_report_iter_exprs <- function(proto_ipm, translate_greek) {

  # Set up iteration expressions, and get sub-kernel families and start/end
  # state so we can append z/z' correctly.

  k_row      <- .init_iteration(proto_ipm, TRUE, "right")

  families   <- .sub_kernel_families(proto_ipm)
  start_end  <- .find_start_end(proto_ipm)

  iter_exprs <- k_row$params$vr_text

  latex_exprs <- .translate_iter_exprs(iter_exprs,
                                       families,
                                       start_end)

  if(translate_greek) {
    .update_make_ipm_report_proto_ipm(proto_ipm,
                                      translate_greek,
                                      latex_exprs)
  }

}

#' @noRd
# TODO need to work out how multi-state models fit into all this.

.translate_iter_exprs <- function(iter_exprs, families, start_end) {

  all_syms <- lapply(iter_exprs, .args_from_txt) %>%
    unlist(use.names = FALSE) %>%
    setNames(.) %>%
    as.list()


}

#' @noRd
# returns vector w/ family as value and kernel_id as name:
# c(P_yr = "CC", F_yr = "CC", goSB = "CD", ...)
.sub_kernel_families <- function(proto_ipm) {

  vapply(proto_ipm$params, function(x) x$family, character(1L)) %>%
    setNames(proto_ipm$kernel_id)

}

#' @noRd
# Translates:
# 1. if(translate_greek): linear functional forms into a common rep
#    - s_int + s_slope * z_1 -> beta_0_s + beta_1_s * z_1
# 2. translated forms into LaTex
#    - beta_0_s + beta_1_s * z_1 -> \beta_{0,s} + \beta_{1,s} * z
# Notes:
#     - if(translate_greek): attach dictionary mapping old names -> new names as
#       as an attribute to proto_ipm so subsequent functions have less work to
#       do.

.pdb_make_ipm_report_vr_exprs <- function(proto_ipm, translate_greek) {


  if(translate_greek) {
    .update_make_ipm_report_proto_ipm(proto_ipm,
                                      translate_greek,
                                      translated_values)
  }

  return(translated_values)

}

#' @noRd

.pdb_make_ipm_report_params <- function(proto_ipm, translate_greek) {

  # Define rest

  if(translate_greek) {
    .update_make_ipm_report_proto_ipm(proto_ipm,
                                      translate_greek,
                                      translated_values)
  }


}

#' @noRd

.pdb_make_ipm_report_impl_args <- function(proto_ipm, translate_greek) {


  if(translate_greek) {
    .update_make_ipm_report_proto_ipm(proto_ipm,
                                      translate_greek,
                                      translated_values)
  }
}

#' @noRd

.pdb_make_ipm_report_glossary <- function(proto_ipm, translate_greek) {



}

.update_make_ipm_report_proto_ipm <- function(proto_ipm,
                                              translate_greek,
                                              translated_values) {

    proto_ipm <- .append_greek_dictionary(proto_ipm, translated_values)

    # This will overwrite proto_ipm in make_report_body, as that is always
    # 2 frames up from the current call

    assign("proto_ipm", proto_ipm, envir = rlang::caller_env(n = 2))

}

#' @noRd
# to_append should be character vector or list w/ names of parameter values
# as values and translated values as names.

.append_greek_dictionary <- function(proto_ipm, to_append) {

  dict <- attr(proto_ipm, "translation_dictionary")

  if(is.null(dict)) {

    attr(proto_ipm, "translation_dictionary") <- to_append

  } else {

    to_append <- c(dict, to_append)

    attr(proto_ipm, "translation_dictionary") <- to_append

  }

  return(proto_ipm)
}

#' @noRd

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
