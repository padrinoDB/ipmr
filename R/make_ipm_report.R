#' @title Generate an RMarkdown file with IPM metadata
#'
#' @rdname ipm_report
#'
#' @description Generates a \code{.rmd} file containing a mathematical
#'   description of the \code{proto_ipm} object.
#'
#' @param object A \code{proto_ipm} or output from \code{make_ipm()}.
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
#'   with blocks and numbered using \code{tag{}}. If \code{FALSE}, equations
#'   will be rendered as inline equations on a single line, and numbered as
#'   1.1, 1.2, 1.3 (iteration expressions), 2.1, 2.2 (vital rate expressions),
#'   etc.
#' @param long_eq_length For longer equations, \code{make_ipm_report} tries
#'   to wrap these into multiple lines using \code{\\\\}. This parameter controls
#'   the number of characters per line. Default is 65. Ignored when
#'   \code{block_eqs = FALSE}.
#'
#' @details \code{make_ipm_report_body} only translates the iteration
#'   expressions and vital rate expressions into Markdown with LaTeX, and does
#'   not produce any headers needed to knit the file. This function is exported
#'   mostly for re-usage in \code{\link[Rpadrino]{pdb_report}}, and isn't really
#'   intended for use by \code{ipmr} users.
#'
#' @section \strong{Translations}:
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

make_ipm_report <- function(object,
                            rmd_dest       = getwd(),
                            title          = "",
                            output_format  = "html",
                            render_output  = FALSE,
                            block_eqs      = TRUE,
                            long_eq_length = 65) {

  UseMethod("make_ipm_report")
}

#' @export
#' @rdname ipm_report

make_ipm_report.default <- function(object,
                                    rmd_dest      = getwd(),
                                    title         = "",
                                    output_format = "html",
                                    render_output = FALSE,
                                    block_eqs     = TRUE,
                                    long_eq_length = 65) {

  .make_ipm_report_impl(object, rmd_dest, title, output_format,
                        render_output, block_eqs, long_eq_length)
}

#' @export
#' @rdname ipm_report

make_ipm_report.ipmr_ipm <- function(object,
                                     rmd_dest      = getwd(),
                                     title         = "",
                                     output_format = "html",
                                     render_output = FALSE,
                                     block_eqs     = TRUE,
                                     long_eq_length = 65) {

  .make_ipm_report_impl(object$proto_ipm, rmd_dest, title, output_format,
                        render_output, block_eqs, long_eq_length)
}

#' @noRd

.make_ipm_report_impl <- function(proto_ipm,
                                  rmd_dest,
                                  title  ,
                                  output_format,
                                  render_output,
                                  block_eqs,
                                  long_eq_length) {

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
  body         <- make_ipm_report_body(proto_ipm, block_eqs,
                                       rmd_dest, long_eq_length)

  rmd_contents <- c(header, body)

  fp <- file.create(rmd_dest, showWarnings = FALSE)

  writeLines(rmd_contents, con = rmd_dest)

  if(render_output) rmarkdown::render(rmd_dest)

  return(rmd_dest)

}

#' @rdname ipm_report
#' @param proto_ipm A \code{proto_ipm} object. Only used for
#' \code{make_ipm_report_body}.
#' @export

make_ipm_report_body <- function(proto_ipm,
                                 block_eqs,
                                 rmd_dest,
                                 long_eq_length) {

  # first, generate stack of translation environments. Each environment on the
  # stack uses the previous environment as the parent for symbol lookup and
  # translation.
  # 1. unk_env is first because it contains everything in the IPM. This guarantees
  # that every variable will at least get defined as "itself". Unknown functions
  # that aren't added to .math_env() will fall back to here.

  unk_env   <- .unk_env(proto_ipm)
  greek_env <- .greek_env(unk_env)
  pop_env   <- .pop_env(proto_ipm, greek_env)
  math_env  <- .math_env(pop_env)
  par_env   <- .par_env(proto_ipm, math_env)


  iter_exprs <- .make_ipm_report_iter_exprs(proto_ipm, par_env,
                                            block_eqs, rmd_dest,
                                            long_eq_length)

  vr_exprs   <- .make_ipm_report_vr_exprs(proto_ipm,
                                          par_env,
                                          block_eqs,
                                          long_eq_length)
  all_params <- .make_ipm_report_params(proto_ipm,
                                        par_env)

  out <- c(iter_exprs, vr_exprs, all_params)

  return(out)

}

#' @noRd

.wrap_inline_eq <- function(eq, section, ind) {

  paste0("\n", section, ".", ind, ": ", "$", eq, "$\n")

}

#' @noRd

.wrap_block_eq <- function(eq, section, ind, long_eq_length) {

  if(nchar(eq) > long_eq_length) {

    eq <- .split_eq_lines(eq, long_eq_length)
  }

  paste0("\n$$", eq, "\\tag{", section, ".", ind, "}", "$$\n")

}

#' @noRd

# This is unlikely to be respected by output: pdf_document based on my
# experience, buy may as well try

.split_eq_lines <- function(eq, long_eq_length) {

  all_lines <- strwrap(eq, width = long_eq_length)

  paste(all_lines, collapse = " \\\\ ")

}

#' @noRd
# Translates:
# n_z_t_1 = P %*% n_z_t + F %*% n_z_t -> n(z', t + 1) = \int([P+F]n(z , t)dz)

.make_ipm_report_iter_exprs <- function(proto_ipm, pop_env,
                                        block_eqs, rmd_dest,
                                        long_eq_length) {

  # Set up iteration expressions, and get sub-kernel families and start/end
  # state so we can append z/z' correctly.

  k_row      <- .init_iteration(proto_ipm, TRUE, "right")

  families   <- .sub_kernel_families(proto_ipm)
  start_end  <- .find_start_end(proto_ipm)
  kern_forms <- kernel_formulae(proto_ipm)

  # This adds family and start_end to right_mult calls. These
  # are used by the function definition in 'math_env()' to decide
  # where a given sub-kernel/pop_state combo needs integration in the
  # latex expression.

  iter_exprs <- k_row$params[[1]]$formula %>%
    .add_right_mult_args(families, start_end, pop_env) %>%
    lapply(rlang::parse_expr)

  t_1_nms    <- names(iter_exprs)

  latex_exprs <- .translate_kern_exprs(iter_exprs,
                                       proto_ipm,
                                       families,
                                       start_end,
                                       pop_env,
                                       block_eqs,
                                       long_eq_length)

  unlist(c(.make_ipm_report_iter_exprs_header(rmd_dest),
           latex_exprs),
         use.names = FALSE)

}

#' @noRd

.make_ipm_report_iter_exprs_header <- function(rmd_dest) {

  paste0("\n## IPM Iteration Expressions\n\n",
         "These expressions iterate the IPM. Check translations from ",
         "R code to Latex for accuracy before distributing!",
         " If needed, edit the _Rmd_ file directly. It can be found here: ",
         paste0("`", rmd_dest, "`"),
         "\n")


}

#' @importFrom rlang eval_bare
#' @noRd

.add_right_mult_args <- function(iter_exprs, families, start_end, pop_env) {

  fam_env <- rlang::child_env(.parent = pop_env,
                              !!! families)

  se_env  <- rlang::child_env(.parent = pop_env,
                              !!! start_end)

  # Adds the kernel family (CC,CD,DC,DD) to the call to right_mult().
  # This determines whether or not a sub-kernel should get an integral.

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

  # Adds state_start to the call to right_mult. If the sub-kernel gets an
  # integral, this determines the domain that integration is performed over.

  se_env$right_mult <- function(kernel, vectr, family) {
    k2 <- rlang::enexpr(kernel)
    v2 <- rlang::enexpr(vectr)
    f2 <- rlang::enexpr(family)

    se <- rlang::eval_bare(kernel, se_env)

    rlang::expr_text(
      rlang::expr(
        right_mult(kernel    = !! k2,
                   vectr     = !! v2,
                   family    = !! f2,
                   start_end = !! se)
      )
    )

  }

  # This actually performs that evaluation.
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

.translate_kern_exprs <- function(full_exprs,
                                  proto_ipm,
                                  families,
                                  start_end,
                                  pop_env,
                                  block_eqs,
                                  long_eq_length) {

  kern_forms <- kernel_formulae(proto_ipm)
  names(kern_forms) <- .z_p_z_kernel_names(kern_forms, families, start_end)

  nms <- out <- vector("list", length(full_exprs) + length(kern_forms))

  # Handles the iteration expressions first. Next, we'll loop over kern_forms
  # to create the correct format for them.

  for(i in seq_along(full_exprs)) {

    # Translates right hand side of expression to latex

    out[[i]] <- rlang::eval_bare(full_exprs[[i]], env = pop_env)

    # Translate the left hand side to latex (i.e. n(z', t + 1))

    nms[i] <- rlang::eval_bare(
      rlang::parse_expr(names(full_exprs)[i]),
      env = pop_env
    )

    # And now put it all together and wrap w/ $$/$!
    out[[i]] <- paste0(nms[i], " = ", out[[i]])

    if(block_eqs) {
      out[[i]] <- .wrap_block_eq(out[[i]], 1, i, long_eq_length)
    } else {
      out[[i]] <- .wrap_inline_eq(out[[i]], 1, i)
    }

  }

  # Create iterator for (i + 1) from above. indexs
  for(i in (seq_along(kern_forms) + length(full_exprs))) {

    it <- i - length(full_exprs)

    out[[i]] <- eval_bare(kern_forms[[it]], pop_env) %>%
      .rm_dz(families[[it]], start_end[[it]])
    out[[i]] <- paste0(names(kern_forms)[it], " = ", out[[i]])

    if(block_eqs) {
      out[[i]] <- .wrap_block_eq(out[[i]], 1, i, long_eq_length)
    } else {
      out[[i]] <- .wrap_inline_eq(out[[i]], 1, i)
    }
  }

  return(out)

}

#' @noRd
# Removes d_z's that are already present in the iteration expressions.

.rm_dz <- function(txt, family, start_end) {

  if(family %in% c("DC", "DD")) return(txt)

  d_z <- paste0(c("\\* d\\\\_", "\\*d\\\\_"), start_end[1], "$")

  for(i in seq_along(d_z)) {
    txt <- gsub(d_z[i], "", txt)
  }

  return(txt)

}

#' @noRd
.z_p_z_kernel_names <- function(kern_forms, families, start_end) {

  stopifnot(
    all(names(families) %in% names(kern_forms)) &&
    length(families) == length(kern_forms)
  )

  stopifnot(
    all(names(start_end) %in% names(kern_forms)) &&
    length(start_end) == length(kern_forms)
  )

  # Make sure positional matching works
  families  <- families[names(kern_forms)]
  start_end <- start_end[names(kern_forms)]

  out <- character(length(kern_forms))
  for(i in seq_along(kern_forms)) {

    nm <- switch(
      families[[i]],
      "DD" = "",                       # Discrete to discrete get no modification
      "CC" = paste0("(",               # Create <(z',z)>
                    start_end[[i]][2],
                    "', ",
                    start_end[[i]][1],
                    ")"),
      "DC" = paste0("(", start_end[[i]][2], "')"), # Create <(z')>
      "CD" = paste0("(", start_end[[i]][1], ")")   # Create <(z)>
    )

    out[i] <- paste0(names(kern_forms)[i], nm)

  }

  # Escape underscored names
  out <- gsub("_", "\\_", out, fixed = TRUE)

  return(out)

}

#' @noRd

.add_integral <- function(eq, start_end) {

  dz <- paste0("d", start_end)

  Lz <- paste0("L_{", start_end, "}")
  Uz <- paste0("U_{", start_end, "}")


  paste0(" \\int_{", Lz, "}^{", Uz, "}", eq, dz)

}

#' @noRd

.unk_env <- function(proto_ipm) {

  # Hold symbols that we can't convert for some reason. Most of these
  # will get overridden by the .par_env object, but this houses all symbols
  # in the IPM so that at least they get a proper underscore

  all_forms <- c(vital_rate_exprs(proto_ipm),
                 kernel_formulae(proto_ipm))
  all_nms   <- gsub("_", "\\_", names(all_forms), fixed = TRUE)
  all_nms   <- setNames(as.list(all_nms), names(all_forms))

  all_syms  <- lapply(all_forms,
                      function(x) {
    rlang::expr_text(x) %>%
      .args_from_txt()
  }) %>%
    unlist() %>%
    unique()

  vals <- gsub(pattern = "_", replacement = "\\_", x = all_syms, fixed = TRUE)
  out <- c(setNames(vals, all_syms), all_nms)

  all_calls <- lapply(all_forms, function(x) {
    rlang::expr_text(x) %>%
    .calls_from_txt()
  }) %>%
    unlist() %>%
    unique()

  call_vals <- lapply(all_calls, .unknown_op) %>%
    setNames(all_calls)
  call_vals <- c(call_vals, list(c = base::c))

  out <- c(out, call_vals)

  rlang::as_environment(out[!duplicated(names(out))])

}

#' @noRd

.binary_group_op <- function(sep) {
  rlang::new_function(
    rlang::exprs(e1 = , e2 = ),
    rlang::expr(
      paste0(e1, !!sep, "{", e2, "}")
    ),
    rlang::caller_env()
  )
}

.binary_or_unary_op <- function(sep) {
  rlang::new_function(
    rlang::exprs(e1 = , e2 = ),
    rlang::expr(
      if(!missing(e2)){
        paste0(e1, !!sep, e2)
      } else {
        paste0(!!sep, e1)
      }

      ),
    rlang::caller_env()
  )
}

#' @noRd
# ... is not an essential argument of this function, but is required to silence
# cran NOTEs about 'possibly incorrect usage of ...' because it does not
# rlang::new_function() as something that creates new functions. Unfortunately,
# we need the quasi-quoting/unquoting version of function creation to correctly
# get 'op'/'fun' into the output function.

.unknown_op <- function(op, ...) {

  op <- gsub("_", "\\_", op, fixed = TRUE)

  rlang::new_function(
    exprs(... = ),
    expr({
      args <- paste(unlist(list(...)), collapse = ", ")
      paste0(!!paste0("\\mathrm{", op, "}("), args, ")")
    })
  )

}

#' @noRd
# ... is not an essential argument of this function, but is required to silence
# cran NOTEs about 'possibly incorrect usage of ...' because it does not
# rlang::new_function() as something that creates new functions. Unfortunately,
# we need the quasi-quoting/unquoting version of function creation to correctly
# get 'op'/'fun' into the output function.

.new_pdf <- function(fun, ...) {
  rlang::new_function(
    rlang::exprs(... = , .pop_ev = pop_env),
    rlang::expr(
      {
        args <- rlang::enexprs(...)

        for(i in seq_along(args)) args[[i]] <- rlang::eval_bare(args[[i]],
                                                                .pop_ev)


        paste0(!! fun,
               "(",
               args[[1]],
               " | ",
               paste(unlist(args[-1]), collapse = ", "),
               ")")

      }
    ),
    rlang::caller_env()
  )
}

#' @noRd

.math_env <- function(pop_env) {

  rlang::child_env(
    .parent = pop_env,

    # Usual operations
    `+` = .binary_op(" + "),
    `-` = .binary_or_unary_op(" - "),
    `*` = .binary_op(" * "),
    `/` = function(a, b) {
      paste0("\\frac{", a, "}{", b, "}")
    },
    `>` = function(a, b) paste0(a, "\\gt", b),
    `>=` = function(a, b) paste0(a, "\\gte", b),
    `<` = function(a, b) paste0(a, "\\lt", b),
    `<=` = function(a, b) paste0(a, "\\lte", b),

    `^` = .binary_group_op("^"),
    `%*%` = .binary_op(""),
    ifelse = function(a, tr, fa) {
      paste0("\\begin{cases} ", tr, " & \\text{if }", a, "\\\\",
             fa, "\\end{cases}")
    },
    t = function(a) paste0(a, "^T"),

    # Commonly used functions
    paste      = paste,
    exp        = function(a) paste0("e^{", a, "}"),
    `[`        = .unary_op("[", "]"),
    `(`        = .unary_op("\\left(", "\\right)"),
    predict    = .unknown_op("predict"),
    right_mult = function(kernel, vectr, family, start_end) {

      z_p_z <- switch(family,
                      "cc" = paste(rev(start_end), collapse = ", "),
                      "CD" = start_end[1],
                      "DC" = paste0(start_end[2], "'"),
                      "DD" = NA)

      out <- paste0("[", kernel,
                    ifelse(!is.na(z_p_z), paste0("(", z_p_z, ")"), ""),
                     "]", vectr)

      if(substr(family[1], 1, 1) == "C") {
        out <- .add_integral(out, start_end[1])
      }

      out
    },

    # Probability Density Functions
    dbeta     = .new_pdf("Beta"),
    dbinom    = .new_pdf("Binomial"),
    dcauchy   = .new_pdf("Cauchy"),
    dchisq    = .new_pdf("\\chi^2"),
    dexp      = .new_pdf("Exponential"),
    df        = .new_pdf("F"),
    dgeom     = .new_pdf("Geometric"),
    dhyper    = .new_pdf("Hypergeometric"),
    dlnorm    = .new_pdf("LogNormal"),
    dmultinom = .new_pdf("Multinomial"),
    dpois     = .new_pdf("Poisson"),
    dt        = .new_pdf("T"),
    dunif     = .new_pdf("Uniform"),
    dweibull  = .new_pdf("Weibull"),
    dmvnorm   = .new_pdf("MVNormal"),
    dnorm     = .new_pdf("Norm"),
    dgamma    = .new_pdf("Gamma"),
    dnbinom   = .new_pdf("NegBinomial"),


    # Other math functions
    sqrt = .unary_op("\\sqrt{", "}"),
    sin  = .unary_op("\\sin(", ")"),
    cos  = .unary_op("\\cos(", ")"),
    tan  = .unary_op("\\tan(", ")"),
    log  = .unary_op("\\log(", ")"),
    abs  = .unary_op("\\lvert|", "\rvert|")
  )

}

#' @noRd

.greek_env <- function(unk_env) {
  greek <- c(
    "alpha", "theta", "tau", "beta", "vartheta", "pi", "upsilon",
    "gamma", "varpi", "phi", "delta", "kappa", "rho",
    "varphi", "epsilon", "lambda", "varrho", "chi", "varepsilon",
    "mu", "sigma", "psi", "zeta", "nu", "varsigma", "omega", "eta",
    "xi", "Gamma", "Lambda", "Sigma", "Psi", "Delta", "Xi",
    "Upsilon", "Omega", "Theta", "Pi", "Phi"
  )
  greek_list <- stats::setNames(paste0("\\", greek), greek)
  greek_env <- rlang::as_environment(greek_list, parent = unk_env)

  return(greek_env)
}

#' @noRd
# TODO: Determine best way forward for parameter set indices.

.par_env <- function(proto_ipm, math_env) {

  params <- parameters(proto_ipm)
  vrs    <- vital_rate_exprs(proto_ipm)

  # These will get added afterwards to make sure they don't get beta'd
  doms   <- .domain_list_to_latex(proto_ipm)

  params <- c(params, vrs)

  # if(.has_par_sets(proto_ipm)) params <- .ipmr_report_par_set_nms(params,
  #                                                                 proto_ipm)

  beta_list <- vector("list", length(params))

  for(i in seq_along(beta_list)) {

    # This will capture all numbers, characters, etc. It won't capture
    # things like model objects. These will get replaced by their data_list
    # name.

    if(rlang::is_syntactic_literal(params)) {

      beta_list[[i]] <- paste0("\\beta_{", i, "}")

    } else {

      beta_list[[i]] <- gsub("_", "\\_",  names(params)[i], fixed = TRUE)

    }
  }

  beta_list <- setNames(beta_list, names(params)) %>%
    c(doms)


  rlang::as_environment(beta_list, parent = math_env)

}

# Evaluates n_z_t_1 -> "n(z', t+1) = iter_expr * n_z_t"
# Evaluates n_z_t   -> "n_z_t      = n(z, t)"

.pop_env <- function(proto_ipm, greek_env)  {

  sv_nms    <- unique(unlist(proto_ipm$state_var))
  cont_doms <- names(domains(proto_ipm))

  doms   <- .domain_list_to_latex(proto_ipm)

  out <- rlang::child_env(greek_env,
                          !!! doms)
  for(i in seq_along(sv_nms)) {

    t_nms   <- paste0("n_", sv_nms[i], "_t")
    t_1_nms <- paste0("n_", sv_nms[i], "_t_1")

    if(sv_nms[i] %in% cont_doms) {

      # Latex values output
      t_vals   <- paste0('n(', sv_nms[i], ", t)") %>%
        as.list() %>%
        setNames(t_nms)
      t_1_vals <- paste0('n(', sv_nms[i], "', t + 1)") %>%
        as.list() %>%
        setNames(t_1_nms)
    } else {
      t_vals <- paste0(paste0(sv_nms[i], "(t)")) %>%
        as.list() %>%
        setNames(t_nms)
      t_1_vals <- paste0(paste0(sv_nms[i], "(t + 1)")) %>%
        as.list() %>%
        setNames(t_1_nms)
    }

    rlang::env_bind(out, !!! t_vals, !!! t_1_vals)

  }

  out
}

.domain_list_to_latex <- function(proto_ipm) {

  doms <- domains(proto_ipm)

  out <- list()

  for(i in seq_along(doms)) {

    z_1 <- names(doms)[i]
    z_2 <- paste0(z_1, "'")

    z_1_nm <- paste0(z_1, "_1")
    z_2_nm <- paste0(z_1, "_2")

    out[[i]] <- rlang::list2(
      !! z_1_nm := z_1,
      !! z_2_nm := z_2
    )

  }

  .flatten_to_depth(out, 1L)
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
# TODO: add z',z if possible, to  vr_exprs

.make_ipm_report_vr_exprs <- function(proto_ipm, pop_env,
                                      block_eqs, long_eq_length) {

  vr_exprs <- vital_rate_exprs(proto_ipm)

  out <- nms <- vector("list", length(vr_exprs))

  for(i in seq_along(vr_exprs)) {

    out[[i]] <- rlang::eval_bare(vr_exprs[[i]], env = pop_env)

    # Translate the left hand side to latex (i.e. grow_mod = grow\\_mod)

    nms[i] <- rlang::eval_bare(
      rlang::parse_expr(names(vr_exprs)[i]),
      env = pop_env
    )

    # And now put it all together and wrap w/ $$/$!
    out[[i]] <- paste0(nms[i], " = ", out[[i]])

    if(block_eqs) {
      out[[i]] <- .wrap_block_eq(out[[i]], 2, i, long_eq_length)
    } else {
      out[[i]] <- .wrap_inline_eq(out[[i]], 2, i)
    }
  }

  unlist(c(.make_ipm_report_vr_exprs_header(), out), use.names = FALSE)

}

.make_ipm_report_vr_exprs_header <- function() {

  paste0("\n\n## IPM Vital Rate Expressions\n\n",
         "These expressions generate the vital rates the IPM. Check  ",
         "translations from R code to Latex for accuracy before distributing!\n")

}

#' @noRd
# TODO: add mesh points and int rule
.make_ipm_report_params <- function(proto_ipm, par_env) {

  pars <- parameters(proto_ipm)
  doms <- domains(proto_ipm)
  # Only one integration rule currently, might want to add it as an attribute
  # to domains() or something, so that it's easier to get text version directly
  # from the domain that is being integrated.
  # ints <- unique(proto_ipm$int_rule)

  par_trans <- rlang::env_get_list(par_env, names(pars))

  par_out <- character(length(pars))

  for(i in seq_along(pars)) {

    if(isTRUE(is.numeric(pars[[i]]))){
      par_out[i] <- paste(par_trans[[i]],
                          round(pars[[i]], 3),
                          sep = " = ")
    } else {

      par_out[i] <- paste(par_trans[[i]],
                          "is a user-defined object and has not been",
                          "translated. Please provide the details manually.")

    }

  }

  dom_out <- character(length(doms))

  for(i in seq_along(doms)) {

    d_nm <- names(doms)[i]
    l_nm <- paste0("L_{", d_nm, "}")
    u_nm <- paste0("U_{", d_nm, "}")

    dom_out[i] <- paste0("$", d_nm, " = ",
                         "[",
                         l_nm, " = ", round(doms[[i]][1], 3),
                         ", ",
                         u_nm, " = ", round(doms[[i]][2], 3),
                         "], ",
                         "n_{", d_nm, "} = ", doms[[i]][2],
                         "$\n\n",
                         " where $n_{x}$ denotes the number of meshpoints",
                         " for the midpoint rule for integration.")

  }

  par_out <- .pars_to_latex_list(par_out)
  dom_out <- .pars_to_latex_list(dom_out)

  par_hdr <- paste0(
    "\n\n## Implementation Details\n\n### Parameter values\n\n",
    "The following parameter values were used to implement this IPM: \n\n"
  )
  dom_hdr <- paste0(
    "\n\n### Domains and Integration Rules\n\n",
    "The following domains and integration rules were used to implement this IPM: \n\n")

  c(par_hdr, par_out, dom_hdr, dom_out)

}

.pars_to_latex_list <- function(pars) {

  pars <- paste("  -", pars)
  pars <- paste(pars, collapse = "\n\n")

  pars

}

#' @noRd

.make_ipm_report_header <- function(output_format, title, block_eqs) {

  if(block_eqs && ! output_format %in% c("pdf", "html")) {

    message("Block equation numbering may not work well in formats other than",
            " 'pdf' or 'html'!\nMake sure to inspect output.")
  }


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

    rmd_dest <- tempfile(pattern = paste0("ipmr_report_", date),
                         fileext = ".Rmd")

    message("'keep_rmd = TRUE' and 'rmd_dest' is not specified! ",
            "Saving to file: \n", rmd_dest)

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
