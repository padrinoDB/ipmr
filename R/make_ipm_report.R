#' @title Generate an RMarkdown file with IPM metadata
#' @rdname ipm_report
#'
#' @description Generates a \code{.rmd} file containing a mathematical
#'   description of the \code{proto_ipm} object.
#'
#' @param proto_ipm A proto_ipm object
#' @param rmd_dest  The folder to save the Rmd file at. The default is
#'   \code{getwd()}. Alternatively, can be a complete file path that specifies
#'   the location and title of the document with the extension \code{".rmd"}. in
#'   this case, the current date will be appended to the title.
#' @param title The title to include in the document. This is not necessarily
#'   the same as \code{rmd_dest}, as this appears at the top of the generated
#'   report, and is not included in the file path!
#' @param output_format The format to include in the YAML header for the created
#'   \code{.rmd} document.
#' @param render_output A logical indicating whether to call
#'   \code{rmarkdown::render} on the generated \code{.rmd} file. Often times,
#'   the \code{.rmd} file will need further editing before it's useful, so the
#'   default is \code{FALSE}.
#' @param block_eqs A logical. If \code{TRUE}, all equations will be inserted
#'   with blocks, and
#'
#' @details \code{make_ipm_report_body} only translates the iteration
#'   expressions and vital rate expressions into Markdown with LaTeX, and does
#'   not produce any headers needed to knit the file. This function is exported
#'   mostly for re-usage in \code{\link[Rpadrino]{pdb_report}}, and isn't really
#'   intended for use by \code{ipmr} users.
#'
#' @section \strong{Translations}
#'
#'   For iteration expressions, vital rate expressions, and parameter names,
#'   \code{make_ipm_report} first translates all values in the \code{data_list}
#'   to \code{beta_X}. For example, \code{s = surv_int + surv_slope * z_1} is
#'   translated into \code{beta_0 + beta_1 * z_1}, and then is translated into
#'   LaTeX equations. Since everything is call \code{beta_X}, a glossary is
#'   provided at the end of each report that matches \code{beta}s to their names
#'   in the \code{data_list}.
#'
#' @return For \code{make_ipm_report}, the filepath to the \code{.rmd} file. The
#'   default name is \code{ "ipmr_report_<current_date>.rmd"}. For
#'   \code{make_ipm_report_body}, a character vector with Markdown and LaTeX
#'   suitable for rendering, but without a header.
#'
#'
#' @export


make_ipm_report <- function(proto_ipm,
                            rmd_dest        = getwd(),
                            title           = "",
                            output_format   = "html",
                            render_output   = FALSE) {

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
  header       <- .make_ipm_report_header(output_format, title, block_eqs)
  body         <- make_ipm_report_body(proto_ipm, block_eqs)

  rmd_contents <- c(header, body)

  fp <- file.create(rmd_dest, showWarnings = FALSE)

  writeLines(rmd_contents, con = rmd_dest)

  if(render_output) rmarkdown::render(rmd_dest)

  return(rmd_dest)

}

#' @rdname ipm_report
#' @export

make_ipm_report_body <- function(proto_ipm, block_eqs) {

  # first, generate stack of translation environments. Each environment on the
  # stack uses the previous environment as the parent for symbol lookup and
  # translation.
  # 1. unk_env is first because it contains everything in the IPM. This guarantees
  # that every variable will at least get defined as "itself". Unknown functions
  # that aren't added to .math_env() will fall back to here.
  # 2. greek_env comes next, as this isn't so much unknown variables as it

  unk_env   <- .unk_env(proto_ipm)
  greek_env <- .greek_env(unk_env)
  math_env  <- .math_env(greek_env)
  par_env   <- .par_env(proto_ipm, math_env)
  pop_env   <- .pop_env(proto_ipm, par_env)


  iter_exprs <- .make_ipm_report_iter_exprs(proto_ipm, pop_env, block_eqs)
  vr_exprs   <- .make_ipm_report_vr_exprs(proto_ipm,   pop_env, block_eqs)
  all_params <- .make_ipm_report_params(proto_ipm,     pop_env, block_eqs)
  impl_args  <- .make_ipm_report_impl_args(proto_ipm,  pop_env, block_eqs)
  glossary   <- .make_ipm_report_glossary(proto_ipm,   pop_env, block_eqs)

  out <- c(iter_exprs, vr_exprs, all_params, impl_args, glossary)

  return(out)

}

.wrap_inline_eq <- function(eq, section, ind) {

  paste0("\n$",section, ".", ind, ": ", eq, "$\n")

}

.wrap_block_eq <- function(eq, section, ind) {

  paste0("$$\n\n", eq, "\\tag{", section, ".", ind, "}\n\n", "$$")


}

#' @noRd
# Translates:
# n_z_t_1 = P %*% n_z_t + F %*% n_z_t -> n(z', t + 1) = \int([P+F]n(z , t)dz)
# Use CC/DC to decide if \int applies!

.make_ipm_report_iter_exprs <- function(proto_ipm, pop_env, block_eqs) {

  # Set up iteration expressions, and get sub-kernel families and start/end
  # state so we can append z/z' correctly.

  k_row      <- .init_iteration(proto_ipm, TRUE, "right")

  families   <- .sub_kernel_families(proto_ipm)
  start_end  <- .find_start_end(proto_ipm)

  iter_exprs <- k_row$params[[1]]$formula %>%
    .add_right_mult_args(families, start_end, pop_env) %>%
    lapply(rlang::parse_expr)

  t_1_nms    <- names(iter_exprs)

  latex_exprs <- .translate_iter_exprs(iter_exprs,
                                       families,
                                       start_end,
                                       pop_env,
                                       block_eqs)

  c(.make_ipm_report_iter_exprs_header(), latex_exprs)

}

.add_right_mult_args <- function(iter_exprs, families, start_end, pop_env) {

  fam_env <- rlang::child_env(.parent = pop_env,
                              !!! families)

  se_env  <- rlang::child_env(.parent = pop_env,
                              !!! start_end)

  fam_env$right_mult <- function(kernel, vectr) {
    k2 <- rlang::enexpr(kernel)
    v2 <- rlang::enexpr(vectr)

    fam <- rlang::eval_bare(kernel, fam_env)

    rlang::expr_text(
      rlang::expr(
        right_mult(kernel = !! k2,
                   vectr  = !! v2,
                   family = !! fam)
      )
    )

  }

  se_env$right_mult <- function(kernel, vectr, family) {
    k2 <- rlang::enexpr(kernel)
    v2 <- rlang::enexpr(vectr)
    f2 <- rlang::enexpr(family)

    se <- rlang::eval_bare(kernel, se_env)[1]

    rlang::expr_text(
      rlang::expr(
        right_mult(kernel    = !! k2,
                   vectr     = !! v2,
                   family    = !! f2,
                   start_end = !! se)
      )
    )

  }

  temps <- lapply(iter_exprs, function(x, ev) {
    exs <- rlang::parse_expr(x)
    eval_bare(exs, ev)
  }, ev = fam_env)

  lapply(temps, function(x, ev) {
    x <- rlang::parse_expr(x)
    eval_bare(x, ev)
  }, ev = se_env)

}

#' @noRd
# TODO need to work out how multi-state models fit into all this.

.translate_iter_exprs <- function(full_exprs,
                                  families,
                                  start_end,
                                  pop_env,
                                  block_eqs) {

  out <- vector("list", length(full_exprs))

  for(i in seq_along(full_exprs)) {

    out[[i]] <- rlang::eval_bare(full_exprs[[i]], env = pop_env)

    nms[i] <- rlang::eval_bare(
      rlang::parse_expr(names(full_exprs)[i]),
      env = pop_env
    )

    out[[i]] <- paste0(nms[i], " = ", out[[i]])

    if(block_eqs) {
      out[[i]] <- .wrap_block_eq(out[[i]], 1, i)
    } else {
      out[[i]] <- .wrap_inline_eq(out[[i]], 1, i)
    }


  }

  names(out) <- names(full_exprs)

  return(out)

}

.add_integral <- function(eq, start_end) {

  dz <- paste0("d_", start_end)

  Lz <- paste0("L_{", start_end, "}")
  Uz <- paste0("U_{", start_end, "}")


  paste0(" \\int_{", Lz, "}^{", Uz, "}", eq, dz)

}

.unk_env <- function(proto_ipm) {

  all_forms <- c(vital_rate_exprs(proto_ipm),
                 kernel_formulae(proto_ipm))

  all_syms  <- lapply(all_forms,
                      function(x) {
    rlang::expr_text(x) %>%
      .args_from_txt()
  }) %>%
    unlist() %>%
    unique() %>%
  c(proto_ipm$kernel_id)

  vals <- gsub(pattern = "_", replacement = "\\_", x = all_syms, fixed = TRUE)

  rlang::as_environment(setNames(vals, all_syms))

}

.binary_group_op <- function(sep) {
  rlang::new_function(
    rlang::exprs(e1 = , e2 = ),
    rlang::expr(
      paste0(e1, !!sep, "{", e2, "}")
    ),
    rlang::caller_env()
  )
}

.math_env <- function(greek_env) {

  rlang::child_env(
    .parent = greek_env,
    `+` = .binary_op(" + "),
    `-` = .binary_op(" - "),
    `*` = .binary_op(" * "),
    `/` = function(a, b) {
      paste0("\\frac{", a, "}{", b, "}")
    },
    `^` = .binary_group_op("^"),
    `%*%` = .binary_op(""),

    paste = paste,
    exp = function(a) paste0("e^{", a, "}"),
    `[` = .unary_op("[", "]"),
    `(` = .unary_op("\\left(", "\\right)"),
    as.vector = .unary_op("as.vector(", ")"),
    predict   = function(a) paste0("\\mathrm{predict(", a, ")}"),
    right_mult = function(kernel, vectr, family, start_end) {

      out <- paste0(kernel, vectr)

      if(substr(family[1], 1, 1) == "C") {
        out <- .add_integral(out, start_end)
      }

      out
    },
    left_mult = function(kernel, vectr) paste0(vectr, kernel),

    # Other math functions
    sqrt = .unary_op("\\sqrt{", "}"),
    sin  = .unary_op("\\sin(", ")"),
    cos  = .unary_op("\\cos(", ")"),
    tan  = .unary_op("\\tan(", ")"),
    log  = .unary_op("\\log(", ")"),
    abs  = .unary_op("\\lvert|", "\rvert|")
  )

}

.greek_env <- function(unk_env) {
  greek <- c(
    "alpha", "theta", "tau", "beta", "vartheta", "pi", "upsilon",
    "gamma", "varpi", "phi", "delta", "kappa", "rho",
    "varphi", "epsilon", "lambda", "varrho", "chi", "varepsilon",
    "mu", "sigma", "psi", "zeta", "nu", "varsigma", "omega", "eta",
    "xi", "Gamma", "Lambda", "Sigma", "Psi", "Delta", "Xi",
    "Upsilon", "Omega", "Theta", "Pi", "Phi"
  )
  greek_list <- set_names(paste0("\\", greek), greek)
  greek_env <- rlang::as_environment(greek_list, parent = unk_env)

  return(greek_env)
}

.par_env <- function(proto_ipm, math_env) {

  params <- parameters(proto_ipm)

  # if(.has_par_sets(proto_ipm)) params <- .ipmr_report_par_set_nms(params,
  #                                                                 proto_ipm)

  beta_list <- vector("list", length(params))

  for(i in seq_along(beta_list)) {

    # This will capture all numbers, characters, etc. It won't capture
    # things like model objects. These will get replaced by their data_list
    # name.

    if(rlang::is_syntactic_literal(beta_list[[i]])) {

      beta_list[[i]] <- paste0("\\beta_{", i, "}")

    } else {

      beta_list[[i]] <- names(params)[i]

    }
  }

  beta_list <- setNames(beta_list, names(params))

  rlang::as_environment(beta_list, parent = math_env)

}

# Evaluates n_z_t_1 -> "n(z', t+1) = iter_expr * n_z_t"
# Evaluates n_z_t   -> "n_z_t      = n(z, t)"

.pop_env <- function(proto_ipm, par_env)  {

  dom_nms <- unique(unlist(proto_ipm$state_var))

  # ipmr names to be evaluated.
  t_nms   <- paste0("n_", dom_nms, "_t")
  t_1_nms <- paste0("n_", dom_nms, "_t_1")

  # Latex values output
  t_vals   <- paste0('n(', dom_nms, ", t)") %>%
    as.list() %>%
    set_names(t_nms)
  t_1_vals <- paste0('n(', dom_nms, "', t + 1)") %>%
    as.list() %>%
    set_names(t_1_nms)

  rlang::child_env(par_env,
                   !!! t_vals,
                   !!! t_1_vals)

}


# TODO: Devise way of grouping all parameters from an indexed symbol into
# one '\\beta_{x}^{ind_1, ind_2}'

.ipmr_report_par_set_nms <- function(params, proto_ipm) {



}

#' @noRd

.has_par_sets <- function(proto_ipm) {

  any(proto_ipm$uses_par_sets | proto_ipm$uses_age)

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


}

#' @noRd

.pdb_make_ipm_report_impl_args <- function(proto_ipm, translate_greek) {


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

.make_ipm_report_header <- function(output_format, title, block_eqs) {

  paste("---",
        paste0("title: '", title, "'"),
        paste0("output:\n  ",
               output_format,
               "_document:\n    toc: true\n    toc_depth: 3"),
        paste0("date: '`r Sys.Date()`'"),
        paste0("urlcolor: blue"),
        ifelse(block_eqs,
               paste0("header_includes:\n  - \\usepackage{amsmath}"),
               ""),
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
