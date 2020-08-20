# subset padrino to aldo's model and only tables we need
library(dplyr)

padr <- readRDS("paper/Analyses/Data/padrino_for_synth.rds")

padr <- padr[-c(length(padr))]

id_ind <- which(padr[[1]]$corresponding_author == 'Compagnoni')

id <- padr[[1]]$ipm_id[id_ind]

use_padr      <- lapply(padr, function(x, id) filter(x, ipm_id == id), id = id)

# Add in model class for ipmr to work with. The next bit will be some heinous
# work to get the model into the right shape
use_padr[[1]] <- mutate(use_padr[[1]], model_class = 'general_di_stoch_kern')

keep_ind <- c(1, 5, 6, 8, 9, 10, 12)

padr_for_paper <- use_padr[keep_ind]

saveRDS(padr_for_paper, file = "paper/Analyses/Data/padrino_for_paper.rds")


