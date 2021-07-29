
library(ipmr)

data(iceplant_ex)

# growth model. 

grow_mod <- lm(log_size_next ~ log_size, data = iceplant_ex)
grow_sd  <- sd(resid(grow_mod))

# survival model

surv_mod <- glm(survival ~ log_size, data = iceplant_ex, family = binomial())

# Pr(flowering) model

repr_mod <- glm(repro ~ log_size, data = iceplant_ex, family = binomial())

# Number of flowers per plant model

flow_mod <- glm(flower_n ~ log_size, data = iceplant_ex, family = poisson())

# New recruits have no size(t), but do have size(t + 1)

recr_data <- subset(iceplant_ex, is.na(log_size))

recr_mu  <- mean(recr_data$log_size_next)
recr_sd  <- sd(recr_data$log_size_next)

# This data set doesn't include information on germination and establishment.
# Thus, we'll compute the realized recruitment parameter as the number
# of observed recruits divided by the number of flowers produced in the prior
# year.

recr_n   <- length(recr_data$log_size_next)

flow_n   <- sum(iceplant_ex$flower_n, na.rm = TRUE)

recr_pr  <- recr_n / flow_n


# Now, we put all parameters into a list. This case study shows how to use
# the mathematical notation, as well as how to use predict() methods

all_params <- list(
  surv_int = coef(surv_mod)[1],
  surv_slo = coef(surv_mod)[2],
  repr_int = coef(repr_mod)[1],
  grow_int = coef(grow_mod)[1],
  grow_slo = coef(grow_mod)[2],
  grow_sdv = grow_sd,
  repr_slo = coef(repr_mod)[2],
  flow_int = coef(flow_mod)[1],
  flow_slo = coef(flow_mod)[2],
  recr_n   = recr_n,
  flow_n   = flow_n,
  recr_mu  = recr_mu,
  recr_sd  = recr_sd,
  recr_pr  = recr_pr
)



L <- min(c(iceplant_ex$log_size,
           iceplant_ex$log_size_next),
         na.rm = TRUE) * 1.2

U <- max(c(iceplant_ex$log_size,
           iceplant_ex$log_size_next),
         na.rm = TRUE) * 1.2

n_mesh_p <- 100



carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")


## 
## carpobrotus_ipm <-  define_kernel(
##   proto_ipm = carpobrotus_ipm,
##   name      = "P",
##   formula   = s * G,
##   ...
## )
## 

## 
## carpobrotus_ipm <-  define_kernel(
##   proto_ipm = carpobrotus_ipm,
##   name      = "P",
##   formula   = s * G,
##   family    = "CC",
##   ...
## )
## 

## 
## carpobrotus_ipm <-  define_kernel(
##   proto_ipm = carpobrotus_ipm,
##   name      = "P",
##   formula   = s * G,
##   family    = "CC",
##   G         = dnorm(z_2, mu_g, grow_sdv),
##   mu_g      = grow_int + grow_slo * z_1,
##   s         = plogis(surv_int + surv_slo * z_1),
##   ...
## )
## 


carpobrotus_ipm <-  define_kernel(
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G,
  family    = "CC",
  G         = dnorm(z_2, mu_g, grow_sdv),
  mu_g      = grow_int + grow_slo * z_1,
  s         = plogis(surv_int + surv_slo * z_1),
  data_list = all_params,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun    = "norm", 
                                      target = "G")
)



carpobrotus_ipm <-  define_kernel(
  proto_ipm = carpobrotus_ipm,
  name      = "F",
  formula   = recr_pr * r_s * r_d * p_f,
  family    = "CC",
  r_s       = exp(flow_int + flow_slo * z_1),
  r_d       = dnorm(z_2, recr_mu, recr_sd),
  p_f       = plogis(repr_int + repr_slo * z_1),
  data_list = all_params,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun    = "norm", 
                                      target =  "r_d")
) 



carpobrotus_ipm <-  define_impl(
  proto_ipm = carpobrotus_ipm,
  make_impl_args_list(
    kernel_names = c("P", "F"),
    int_rule     = rep('midpoint', 2),
    state_start    = rep('z', 2),
    state_end      = rep('z', 2) 
  ) 
) 


carpobrotus_ipm <-  define_domains(
  proto_ipm = carpobrotus_ipm,
  z         = c(L, U, n_mesh_p)
)  



carpobrotus_ipm <-  define_pop_state(
  proto_ipm = carpobrotus_ipm,
    n_z     = rep(1/100, n_mesh_p)
) 



carpobrotus_ipm <-  make_ipm(
  proto_ipm       = carpobrotus_ipm,
  iterate         = TRUE,
  iterations      = 100,
  return_main_env = TRUE
)


asymp_grow_rate <- lambda(carpobrotus_ipm)
asymp_grow_rate



# Option 1: is_conv_to_asymptotic

is_conv_to_asymptotic(carpobrotus_ipm)

# Option 2: generate iteration kernel and compute eigenvalues

K <- make_iter_kernel(carpobrotus_ipm)

lam_eigen <- Re(eigen(K$mega_matrix)$values[1])

# If we've iterated our model enough, this should be approximately 0 (though 
# maybe a little off due to floating point errors).

asymp_grow_rate - lam_eigen



# Sub-kernels have their own print method to display the range of values 
# and some diagnotic information.

carpobrotus_ipm$sub_kernels

# Extract the time series of the population state (n_z), 
# and the n_t+1/n_t values (lambda)

pop_time_series    <- carpobrotus_ipm$pop_state$n_z
lambda_time_series <- carpobrotus_ipm$pop_state$lambda

# Next, we'll tweak the intercept of the p_f function and re-fit the model.

new_proto_ipm      <- carpobrotus_ipm$proto_ipm

# The parameters setter function takes a list. It can replace single values,
# create new values, or replace the entire parameter list, depending on how you
# set up the right hand side of the expression. 

parameters(new_proto_ipm) <- list(repr_int = -0.3)

new_carp_ipm <- make_ipm(new_proto_ipm, 
                         iterations = 100)

lambda(new_carp_ipm)



pred_par_list <- list(
  grow_mod = grow_mod,
  grow_sdv = grow_sd,
  surv_mod = surv_mod,
  repr_mod = repr_mod,
  flow_mod = flow_mod,
  recr_n   = recr_n,
  flow_n   = flow_n,
  recr_mu  = recr_mu,
  recr_sd  = recr_sd,
  recr_pr  = recr_pr
)

predict_method_carpobrotus <- init_ipm(sim_gen = "simple", 
                                       di_dd = "di",
                                       det_stoch = "det") %>%
  define_kernel(
    name      = "P",
    formula   = s * G,
    family    = "CC",
    G         = dnorm(z_2, mu_g, grow_sdv),
    mu_g      = predict(grow_mod, 
                        newdata = data.frame(log_size = z_1),
                        type = 'response'),
    s         = predict(surv_mod, 
                        newdata = data.frame(log_size = z_1),
                        type = "response"),
    data_list = pred_par_list,
    states    = list(c('z')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "G")
  ) %>%
  define_kernel(
    name      = "F",
    formula   = recr_pr * r_s * r_d * p_f,
    family    = "CC",
    r_s       = predict(flow_mod, 
                        newdata = data.frame(log_size = z_1),
                        type = "response"),
    r_d       = dnorm(z_2, recr_mu, recr_sd),
    p_f       = predict(repr_mod,
                        newdata = data.frame(log_size = z_1),
                        type = "response"),
    data_list = pred_par_list,
    states    = list(c("z")),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep('z', 2),
      state_end      = rep('z', 2) 
    ) 
  ) %>%
  define_domains(
    z = c(L, U, n_mesh_p)
  )  %>%
  define_pop_state(
    n_z = rep(1/100, n_mesh_p)
  ) %>%
  make_ipm(iterate    = TRUE,
           iterations = 100)



sens <- function(ipm_obj, d_z) {
  
  w <- right_ev(ipm_obj)[[1]]
  v <- left_ev(ipm_obj)[[1]]
  
  return(
    outer(v, w) / sum(v * w * d_z)
  )
  
}



elas <- function(ipm_obj, d_z) {
  
  K           <- make_iter_kernel(ipm_obj)$mega_matrix
  
  sensitivity <- sens(ipm_obj, d_z)
  
  lamb        <- lambda(ipm_obj)
  
  out         <- sensitivity * (K / d_z) / lamb
  
  return(out)
  
}



R_nought <- function(ipm_obj) {
  
  Pm <- ipm_obj$sub_kernels$P
  Fm <- ipm_obj$sub_kernels$F
  
  I  <- diag(dim(Pm)[1])
  
  N  <- solve(I - Pm)
  
  R  <- Fm %*% N
  
  return(
    Re(eigen(R)$values)[1]
  )
  
}

gen_time <- function(ipm_obj) {
  
  lamb     <- unname(lambda(ipm_obj))  
  
  r_nought <- R_nought(ipm_obj)
  
  return(log(r_nought) / log(lamb))
}




mesh_info <- int_mesh(carpobrotus_ipm)

sens_mat <- sens(carpobrotus_ipm, mesh_info$d_z)
elas_mat <- elas(carpobrotus_ipm, mesh_info$d_z)

R0    <- R_nought(carpobrotus_ipm)
gen_T <- gen_time(carpobrotus_ipm)

R0
gen_T



lab_seq  <- round(seq(L, U, length.out = 6), 2)
tick_seq <- c(1, 20, 40, 60, 80, 100)

par(mfrow = c(2, 2))

# Sub-kernels - ipmr contains plot methods for sub-kernels 

plot(carpobrotus_ipm$sub_kernels$P, 
     do_contour = TRUE,
     main       = "P",
     xlab       = "size (t)",
     ylab       = "size (t + 1)",
     yaxt       = "none", 
     xaxt       = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))

plot(carpobrotus_ipm$sub_kernels$F,
     do_contour = TRUE,
     main       = "F",
     xlab       = "size (t)",
     ylab       = "size (t + 1)",
     yaxt       = "none", 
     xaxt       = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))



# Sensitivity and elasticity

class(sens_mat) <- c("ipmr_matrix", class(sens_mat))
class(elas_mat) <- c("ipmr_matrix", class(elas_mat))

plot(sens_mat, 
     do_contour = TRUE,
     main       = "K Sensitivity",
     xlab       = "size (t)", 
     ylab       = "size (t + 1)",
     yaxt       = "none", 
     xaxt       = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))


plot(elas_mat, 
     do_contour = TRUE,
     main       = "K Elasticity",
     xlab       = "size (t)", 
     ylab       = "size (t + 1)",
     yaxt       = "none", 
     xaxt       = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))





par(mfrow = c(1, 1))
K <- make_iter_kernel(carpobrotus_ipm)

plot(K$mega_matrix,
     do_contour = TRUE,
     main       = "K",
     xlab       = "size (t)",
     ylab       = "size (t + 1)",
     yaxt       = "none", 
     xaxt       = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))






library(ggplot2)
library(gridExtra)

p_df    <- ipm_to_df(carpobrotus_ipm$sub_kernels$P)
f_df    <- ipm_to_df(carpobrotus_ipm$sub_kernels$F)
k_df    <- ipm_to_df(K$mega_matrix)

sens_df <- ipm_to_df(sens_mat)
elas_df <- ipm_to_df(elas_mat)

# Create a default theme for our plots

def_theme <- theme(
  panel.background  = element_blank(),
  axis.text         = element_text(size = 16),
  axis.ticks        = element_line(size = 1.5),
  axis.ticks.length = unit(0.08, "in"),
  axis.title.x      = element_text(
    size   = 20,
    margin = margin(
      t = 10,
      r = 0, 
      l = 0, 
      b = 2
    )
  ),
  axis.title.y = element_text(
    size   = 20,
    margin = margin(
      t = 0,
      r = 10,
      l = 2,
      b = 0
    )
  ),
  legend.text = element_text(size = 16)
)

p_plt <- ggplot(p_df) +
  geom_tile(aes(x    = t,
                y    = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("P kernel")

f_plt <- ggplot(f_df) +
  geom_tile(aes(x    = t,
                y    = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("F kernel")

k_plt <- ggplot(k_df) +
  geom_tile(aes(x    = t,
                y    = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("K kernel")

sens_plt <- ggplot(sens_df) +
  geom_tile(aes(x    = t,
                y    = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("K Sensitivity")

elas_plt <- ggplot(elas_df) +
  geom_tile(aes(x    = t,
                y    = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("K Elasticity")


grid.arrange(
  p_plt,    f_plt, k_plt,
  sens_plt, elas_plt,
             layout_matrix = matrix(c(1, 1, 2, 2,
                                      NA, 3, 3, NA, 
                                      4, 4, 5, 5),
                                    nrow = 3, 
                                    byrow = TRUE))



