library(rlang)
library(ipmr)

intro <- tryCatch({
  source("tests/vignettes/ipmr-intro.R", local = TRUE, keep.source = TRUE)
},
error   = function(x) x,
warning = function(x) x)


general <- tryCatch({
  source("tests/vignettes/general-ipms.R", local = TRUE, keep.source = TRUE)
},
error   = function(x) x,
warning = function(x) x)


index <- tryCatch({
  source("tests/vignettes/index-notation.R", local = TRUE, keep.source = TRUE)
},
error   = function(x) x,
warning = function(x) x)

dens <- tryCatch({
  source("tests/vignettes/density-dependence.R", local = TRUE, keep.source = TRUE)
},
error   = function(x) x,
warning = function(x) x)

age <- tryCatch({
  source("tests/vignettes/age_x_size.R", local = TRUE, keep.source = TRUE)
},
error   = function(x) x,
warning = function(x) x)

probs <- list()
if(inherits(intro, c("simpleError", "simpleWarning"))) {
  probs$intro <- intro
}

if(inherits(general, c("simpleError", "simpleWarning"))) {
  probs$general <- general
}

if(inherits(index, c("simpleError", "simpleWarning"))) {
  probs$index <- index
}

if(inherits(dens, c("simpleError", "simpleWarning"))) {
  probs$dens <- dens
}

if(inherits(age, c("simpleError", "simpleWarning"))) {
  probs$age <- age
}

if(!rlang::is_empty(probs)) {

  msg <- c()
  for(i in seq_along(probs)) {

    temp <- paste(names(probs)[i],
                  "\n\tcall: ", probs[[i]]$call,
                  "\n\tmessage: ", probs[[i]]$message,
                  '\n\n')

    msg <- c(msg, temp)

  }


  cat("errors/warnings:\n",
      msg)

}

