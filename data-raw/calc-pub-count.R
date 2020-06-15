# Compute total IPM pubs

library(purrr)
library(dplyr)

pubs <- read.csv("data-raw/padr-pubs-feb-2020.csv",
                 stringsAsFactors = FALSE)

id <- pmap_chr(.l = data.frame(a = pubs$Authors, b = pubs$Journal, c = pubs$Year),
               .f = function(a, b, c) paste(a, b, c, sep = "_"))

length(unique(id))
