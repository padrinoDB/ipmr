
devtools::load_all()
library(ggplot2)

# Figure top left

data(iceplant_ex)

mod <- lm(log_size_next ~ log_size, data = iceplant_ex)

xx <- seq(-5, 5, 0.1)

grow_pred <- data.frame(log_size = xx,
                        pred = predict(mod,
                                       data.frame(log_size = xx),
                                       type = 'response'))



theme_vr <- theme_bw() +
  theme(
    # Extras to theme_bw()
    axis.text.x       = element_text(size = 14),
    axis.text.y       = element_text(size = 14), # make y-axis text + title bigger
    axis.title.x      = element_text(size   = 20,
                                     margin = margin(
                                       t = 5,
                                       r = 0,
                                       l = 0,
                                       b = 7
                                     )
    ),
    axis.title.y     = element_text(size   = 20,
                                    margin = margin(
                                      t = 3,
                                      r = 10,
                                      l = 1,
                                      b = 0
                                    )
    ),
    strip.text        = element_text(size = 8), # Increases plot label size
    legend.background = element_rect(fill = NA,  # Box for the legend
                                     color = 'black'),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7)
  )

# Now, make the figure panel

grow_plot <- ggplot(iceplant_ex, aes(x = log_size,
                                     y = log_size_next)) +
  geom_point(color = 'black',
             size = 1.25) +
  geom_line(data = grow_pred,
            aes(x = log_size,
                y = pred),
            linetype = 'dashed',
            size = 1.25,
            color = 'red',
            show.legend = FALSE)  +
  theme_vr +
  scale_x_continuous('ln(Surface Area, t)',
                     limits = c(-6.1, 3.5)) +
  scale_y_continuous('ln(Surface Area, t + 1)',
                     limits = c(-6.1, 3.5))

png(filename = "paper/Figures/vr_panel.png",
    height = 4,
    width = 6,
    units = 'in',
    res = 300)

  print(grow_plot)

dev.off()


data(sim_di_det_ex)

K <- make_iter_kernel(sim_di_det_ex)

prot <- sim_di_det_ex$proto_ipm

new_ipm <- make_ipm(prot)

w <- right_ev(new_ipm)
v <- left_ev(new_ipm)

dom <- int_mesh(new_ipm) %>%
  lapply(unique)

dom_w <- rev(dom$dbh_1)


png(filename = "paper/Figures/sim_di_det_ex_K.png",
    height = 4,
    width = 6,
    units = 'in',
    res = 300)

  par(mar = c(5.1, 5.1, 1.1, 1.1))

  layout(mat = matrix(c(3, 0, 1, 2),
                      nrow = 2,
                      ncol = 2,
                      byrow = TRUE),
         heights = c(1, 2),
         widths = c(2, 1))


  plot.ipmr_matrix(A = K$mega_matrix ^ 0.15,
                   cex.lab = 2)

  par(mar = c(5.1, 1.1, 1.1, 5.1))

  plot(
    w$dbh_w,
    dom_w,
    type = 'l',
    ylab = "",
    xlab = expression(italic("w(z)")),
    yaxt = "n",
    cex.lab = 2
  )

  axis(side = 4, at = seq(0, 50, by = 10),
       labels = as.character(seq(100, 0, by = -20)))

  mtext(expression(paste("size ", italic(z))), side = 4, line = 3,
        cex = 2)


  par(mar = c(1.1, 5.1, 5.1, 1.1))

  plot(
    dom$dbh_1,
    v$dbh_v,
    type = "l",
    ylab = expression(italic("v(z)")),
    xlab = "",
    xaxt = 'n',
    cex.lab = 2
  )

  axis(side = 3, at = seq(0, 50, by = 10),
       labels = seq(0, 100, by = 20))

  mtext(expression(paste("size ", italic(z))), side = 3, line = 2,
        cex = 2)

dev.off()
