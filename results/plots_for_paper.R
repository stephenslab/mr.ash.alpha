


## function for boxplot
my.box <- function (dat, x, y,
                    values = c(1,2,0,3,4,5)) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           geom_jitter(aes(color = fit, shape = fit), width = .1) +
           scale_shape_manual(values = values) +
           geom_boxplot(alpha = 0.1, aes(color = fit), outlier.alpha = 0) +
           scale_alpha_manual(values = 0.1) +
           scale_fill_discrete(guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 18))
}

## function for error bar
my.errorbar <- function(dat, x, y,
                        values = c(1,2,0,3,4,5)) {
  
}


## function for violin plot
my.jitter <- function (dat, x, y,
                       values = c(1,2,0,3,4,5)) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           geom_jitter(aes(color = fit, shape = fit)) +
           scale_shape_manual(values = values) +
           geom_violin(alpha = 0.1, aes(color = fit, fill = fit), scale = "width") +
           scale_alpha_manual(0.1) +
           scale_fill_discrete(guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 18))
}

## function for color

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

##---------------------
## plot for figure 1
## initialization
##---------------------
setwd("~/git/caisar/result")
dat11 = read.table("dat11.txt", sep = ",")
dat12 = read.table("dat12.txt", sep = ",")

p11 = ggplot(dat11) + geom_line(aes(x = b, y = sb, color = method, linetype = method), size = .8) +
  theme_cowplot(font_size = 14) +
  labs(x  = "b", y = "S(b)",
       title = "MR-ASH shrinkage operator") +
  theme(axis.line    = element_blank(),
        legend.position = "none")
p12 = ggplot(dat12) + geom_line(aes(x = b, y = sb, color = method, linetype = method), size = .8) + 
  theme_cowplot(font_size = 14) +
  labs(x  = "b", y = "S(b)",
       title = "non-convex shrinkage operator") +
  theme(axis.line    = element_blank())
fig1 = plot_grid(p11,p12, nrow = 1, rel_widths = c(0.4,0.465))
ggsave("fig1.pdf", fig1, width = 12, height = 5)


##---------------------
## plot for figure 2
## Practical Consideration
##---------------------

setwd("~/git/caisar/result")
dat21 = read.table("dat21.txt", sep = ",", header = TRUE)[c(1:80,121:140),]
dat22 = read.table("dat22.txt", sep = ",", header = TRUE)[c(1:80,121:140),]
dat23 = read.table("dat23.txt", sep = ",", header = TRUE)[c(1:80,121:140),]
dat24 = read.table("dat24.txt", sep = ",", header = TRUE)[c(1:80,121:140),]
dat22$pred = dat22$pred / sqrt(574/2)
dat23$pred = dat23$pred / dat23$pred[1:20] * 5

p21 = my.box(dat21, "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(6)[4], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "EquiCorrGauss + SparseNormal",
       y     = "relative prediction error") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16), limits = c(1,13)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

p22 = my.box(dat22, "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(6)[4], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "RealGenotype + SparseNormal",
       y     = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16), limits = c(1,13)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.position = "none")

p23 = my.box(dat24, "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(6)[4], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "IndepLowdimGauss + SparseNormal",
       y     = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16), limits = c(1,13)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.position = "none")

p24 = my.box(dat23, "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(6)[4], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "Changepoint + FewSmallChange",
       y     = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16), limits = c(1,13)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.position = "none")

fig2 = plot_grid(p21,p22,p23,p24, nrow = 1, rel_widths = c(0.35,0.35))
ggsave("fig2.pdf", fig2, width = 18, height = 6)

##---------------------
## plot for figure 3
## CAISA versions
##---------------------

setwd("~/git/caisar/result")
dat31 = as.matrix(read.table("dat31.txt", sep = ","))
dat32 = as.matrix(read.table("dat32.txt", sep = ","))
dat33 = as.matrix(read.table("dat33.txt", sep = ","))
dat34 = read.table("dat34.txt", sep = ",")
dat34$y[2041:4000] = dat34$y[2040] # converged
dat34$y[4041:5000] = dat34$y[4040] # converged
dat34$y[7040:9000] = dat34$y[7039] # converged
dat34$y[5049:7000] = dat34$y[5048] # converged
dat35 = dat34[5001:10000,]
dat34 = dat34[1:5000,]

iter = 2000
p31  = ggplot()
time = matrix(1,20,3)
for (i in 1:20) {
  y  = dat31[i,1:(2.1 * iter)]
  y  = c(y, rep(min(y), 0.9 * iter))
  df = data.frame(x = c(1:iter, 1:iter, 1:iter),
                  y = (y - min(y)) / min(y) + 5e-4,
                  method = c(rep("CAISA-em",iter),rep("CAISA-g",iter),rep("CAISA-acc", iter)))
  p31 = p31 + geom_line(data = df, aes(x = x, y = y, color = method), alpha = 0.5)
}
p31 = p31 +
  theme_cowplot(font_size = 18) +
  theme(axis.line = element_blank(),
        legend.position = "none") +
  labs(x = "iteration", y = "relative objective value",
       title = "Orthogonal + SparseNormal") +
  scale_y_continuous(limits = c(5e-4,0.01), trans = "log10") +
  scale_x_continuous(limits = c(0,500))

iter = 2000
p32   = ggplot()
for (i in 1:20) {
  y  = dat32[i,1:(2.5 * iter)]
  df = data.frame(x = c(1:iter, 1:iter, 1:(iter/2)),
                  y = (y - min(y)) / min(y) + 5e-4,
                  method = c(rep("CAISA-em",iter),rep("CAISA-g",iter),rep("CAISA-acc", iter/2)))
  p32 = p32 + geom_line(data = df, aes(x = x, y = y, color = method), alpha = 0.5)
}
p32 = p32 +
  theme_cowplot(font_size = 18) +
  theme(axis.line = element_blank()) +
  labs(x = "iteration", y = "",
       title = "EquiCorrGauss + SparseNormal") +
  scale_y_continuous(limits = c(5e-4,10), trans = "log10") +
  scale_x_continuous(limits = c(0,900))

iter = 1000
p33   = ggplot()
for (i in 1:20) {
  y  = dat33[i,1:(2.1 * iter)]
  y  = c(y, rep(min(y), 0.9 * iter))
  df = data.frame(x = c(1:iter, 1:iter, 1:iter),
                  y = (y - min(y)) / min(y) + 5e-4,
                  method = c(rep("CAISA-em",iter),rep("CAISA-g",iter),rep("CAISA-acc", iter)))
  p33 = p33 + geom_line(data = df, aes(x = x, y = y, color = method), alpha = 0.5)
}
p33 = p33 +
  theme_cowplot(font_size = 18) +
  theme(axis.line = element_blank(),
        legend.position = "none") +
  labs(x = "iteration", y = "",
       title = "IndepLowdimGauss + ThreePointMass") +
  scale_y_continuous(limits = c(5e-4,2), trans = "log10") +
  scale_x_continuous(limits = c(0,120))

dat34$method = factor(dat34$method, levels = c("CAISA-em","CAISA-g","CAISA-acc"))
dat35$method = factor(dat35$method, levels = c("CAISA-em","CAISA-g","CAISA-acc"))
p34 = ggplot(dat34) + geom_line(aes(x = x, y = y, color = method), size = .7, alpha = 0.5) +
  theme_cowplot(font_size = 18) +
  theme(axis.line = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(1e-4,1), trans = "log10") +
  scale_x_continuous(limits = c(0,100)) +
  labs(x = "iteration",
       y = "relative error")
p35 = ggplot(dat35) + geom_line(aes(x = x, y = y, color = method), size = .7, alpha = 0.5) +
  theme_cowplot(font_size = 18) +
  theme(axis.line = element_blank()) +
  scale_y_continuous(limits = c(1e-4,1), trans = "log10") +
  scale_x_continuous(limits = c(0,100)) +
  labs(x = "iteration",
       y = "relative error")

fig3 = plot_grid(p31,p33,p32, nrow = 1, rel_widths = c(0.35,0.35,0.44))
fig31 = plot_grid(p34,p35, nrow = 1, rel_widths = c(0.35,0.45))
title = ggdraw() + draw_label("Two Cases showing different convergence", fontface = 'bold')
fig31 = plot_grid(title, fig31, ncol = 1, rel_heights = c(0.06,0.95))
ggsave("fig3.pdf", fig3, width = 18, height = 6)
ggsave("fig31.pdf", fig31, width = 12, height = 5)

##---------------------
## plot for figure 4
## nonconvex
##---------------------

setwd("~/git/caisar/result")
dat41 = read.table("dat41.txt", sep = ",")
dat42 = read.table("dat42.txt", sep = ",")
dat43 = read.table("dat43.txt", sep = ",")
dat44 = read.table("dat44.txt", sep = ",")
dat45 = read.table("dat45.txt", sep = ",")
dat46 = read.table("dat46.txt", sep = ",")

dat51 = read.table("dat51.txt", sep = ",")[1:505,]
dat52 = read.table("dat52.txt", sep = ",")[1:505,]
dat53 = read.table("dat53.txt", sep = ",")[1:505,]
dat54 = read.table("dat54.txt", sep = ",")[1:505,]
dat55 = read.table("dat55.txt", sep = ",")[1:505,]
dat56 = read.table("dat56.txt", sep = ",")[1:505,]

p41 = my.box(dat41, "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  theme_cowplot(font_size = 14) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(4)[1], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "EquiCorrGauss + SparseNormal",
       y     = "relative prediction error") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(1,16)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p42 = my.box(dat42,  "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  theme_cowplot(font_size = 14) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(4)[1], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "EquiCorrGauss + ThreePointMass", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,10,100), limits = c(1,400)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p43 = my.box(dat43,  "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  theme_cowplot(font_size = 14) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(4)[1], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "EquiCorrGauss + SparseHeavytail", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(1,20)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p44 = my.box(dat44,  "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  theme_cowplot(font_size = 14) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(4)[1], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "IndepLowdimGauss + SparseNormal", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(1,5)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p45 = my.box(dat45,  "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(4)[1], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = #"IndepLowdimGauss + 
       "ThreePointMass", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(1,16)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p46 = my.box(dat46,  "fit", "pred", c(3,2,1,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(4)[1], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = #"IndepLowdimGauss + 
         "SparseHeavytail", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(1,16)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

s = .8
p51 = ggplot(dat51) + geom_line(aes(x = b, y = sb, linetype = prior, color = prior), size = s) +
  theme_cowplot(font_size = 14) +
  #geom_point(aes(x = b, y = sb, color = prior), size = 1) + 
  labs(y = "S(b)") + 
  scale_linetype_manual(values=c("twodash", "dotted","dotdash","solid","dashed","longdash")) +
  theme(axis.line = element_blank(),
        legend.position = "none") +
  guides(color=guide_legend(title="method"), linetype=guide_legend(title="method"))
p52 = ggplot(dat52) + geom_line(aes(x = b, y = sb, linetype = prior, color = prior), size = s) +
  theme_cowplot(font_size = 14) +
  labs(y = "") + 
  scale_linetype_manual(values=c("twodash", "dotted","dotdash","solid","dashed","longdash")) +
  theme(axis.line = element_blank(),
        legend.position = "none") +
  guides(color=guide_legend(title="method"), linetype=guide_legend(title="method"))
p53 = ggplot(dat53) + geom_line(aes(x = b, y = sb, linetype = prior, color = prior), size = s) +
  theme_cowplot(font_size = 14) +
  #geom_point(aes(x = b, y = sb, color = prior), size = 1) + 
  labs(y = "") + 
  scale_linetype_manual(values=c("twodash", "dotted","dotdash","solid","dashed","longdash")) +
  theme(axis.line = element_blank(),
        legend.position = "none") +
  guides(color=guide_legend(title="method"), linetype=guide_legend(title="method"))
p54 = ggplot(dat54) + geom_line(aes(x = b, y = sb, linetype = prior, color = prior), size = s) +
  theme_cowplot(font_size = 14) +
  #geom_point(aes(x = b, y = sb, color = prior), size = 1) + 
  labs(y = "") + 
  scale_linetype_manual(values=c("twodash", "dotted","dotdash","solid","dashed","longdash")) +
  theme(axis.line = element_blank()) +
  guides(color=guide_legend(title="method"), linetype=guide_legend(title="method"))
p55 = ggplot(dat55) + geom_line(aes(x = b, y = sb, linetype = prior, color = prior), size = s) +
  #geom_point(aes(x = b, y = sb, color = prior), size = 1) + 
  labs(y = "") + 
  theme(axis.line = element_blank(),
        legend.position = "none") +
  scale_linetype_manual(values=c("twodash", "dotted","dotdash","solid","dashed","longdash")) +
  guides(color=guide_legend(title="method"), linetype=guide_legend(title="method"))
p56 = ggplot(dat56) + geom_line(aes(x = b, y = sb, linetype = prior, color = prior), size = s) +
  #geom_point(aes(x = b, y = sb, color = prior), size = 1) + 
  labs(y = "") +
  theme(axis.line = element_blank(),
        legend.position = "none") +
  scale_linetype_manual(values=c("twodash", "dotted","dotdash","solid","dashed","longdash"))# +
  guides(color=guide_legend(title="method"), linetype=guide_legend(title="method"))
pempt = ggplot() + geom_blank()
  
fig41 = plot_grid(p41,p42,p43,p44,pempt, nrow = 1, rel_widths = c(0.325,0.325,0.325,0.325,0.1))
fig42 = plot_grid(p51,p52,p53,p54, nrow = 1, rel_widths = c(0.325,0.325,0.325,0.42))
#fig42 = plot_grid(p44,p45,p46,p54,p55,p56, nrow = 2,
#                  rel_widths = c(0.325,0.325,0.325), rel_heights = c(0.5,0.45))
fig4 = plot_grid(fig41,fig42, nrow = 2, rel_heights = c(0.5,0.4))
ggsave("fig4.pdf", fig4, width = 18, height = 10)

##---------------------
## plot for figure 7
## Design 1
##---------------------

dat71 = read.table("dat71.txt", sep = ",")[c(1:20,41:200),]
dat72 = read.table("dat72.txt", sep = ",")[c(1:20,41:200),]
dat73 = read.table("dat73.txt", sep = ",")[c(1:20,41:200),]

p71 = my.box(dat71, "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "EquiCorrGauss + SparseNormal",
       y     = "relative variational objective") + 
  scale_y_continuous(trans = "log10", breaks = c(0.5,1,2,4,8,16), limits = c(0.9,20)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p72 = my.box(dat72,  "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "EquiCorrGauss + SingleEffect", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(0.5,1,2,4,8,16), limits = c(0.9,20)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p73 = my.box(dat73,  "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "EquiCorrGauss + Polygenic", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(0.5,1,2,4,8,16), limits = c(0.9,20)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

fig7 = plot_grid(p71,p72,p73, nrow = 1, rel_widths = c(0.35,0.325,0.325))
ggsave("fig7.pdf", fig7, width = 18, height = 6)

##---------------------
## plot for figure 8
## Design 2
##---------------------

dat81 = read.table("dat81.txt", sep = ",")[21:200,]
dat82 = read.table("dat82.txt", sep = ",")[21:200,]
dat83 = read.table("dat83.txt", sep = ",")[21:200,]

p81 = my.box(dat81, "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "RealGenotype + SparseNormal",
       y     = "relative prediction error") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(.9,26)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p82 = my.box(dat82,  "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "RealGenotype + SingleEffect", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(0.5,1,2,4,8,16), limits = c(.9,26)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p83 = my.box(dat83,  "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "RealGenotype + Polygenic", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(0.5,1,2,4,8,16), limits = c(.9,26)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

fig8 = plot_grid(p81,p82,p83, nrow = 1, rel_widths = c(0.35,0.325,0.325))
ggsave("fig8.pdf", fig8, width = 18, height = 6)


##---------------------
## plot for figure 9
## Design 3
##---------------------

setwd("~/git/caisar/result")
dat91 = read.table("dat91.txt", sep = ",")[21:260,]
dat92 = read.table("dat92.txt", sep = ",")[21:260,]
dat92$pred[dat92$pred < 0.21] = 0.21

p91 = my.box(dat91, "fit", "pred", c(1,2,5,4,3,6,0,8,9,10,7,12,13)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(7)[4], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "Changepoint + FewSmallChange",
       y     = "relative prediction error") + 
  scale_y_continuous(trans = "log10", breaks = c(0.5,1,2,4,8,16), limits = c(0.2,20)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p92 = my.box(dat92,  "fit", "pred", c(1,2,5,4,3,6,0,8,9,10,7,12,13)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(7)[4], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "Changepoint + BigSuddenChange", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(0.5,1,2,4,8,16), limits = c(0.2,20)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

fig9 = plot_grid(p91,p92, nrow = 1)
ggsave("fig9.pdf", fig9, width = 15, height = 5)

##---------------------
## plot for figure 9
## shape of
##---------------------
y1   = double(100); y1[1:10] = y1[1:10] + 1; y1[31:40]=y1[31:40]+4; y1[61:80]=y1[61:80]+2;
y2   = double(100); y2[50:51] = 8;
dat93 = data.frame(x = c(1:100,1:100),
                   y = c(y1, y2),
                   Scenario = rep(c("FewSmallChange", "BigSuddenChange"), each = 100))
p93 = ggplot(dat93) + geom_line(aes(x = x, y = y, color = Scenario)) +
  theme(axis.line    = element_blank()) +
  labs(x = "",
       title = "Shape of piecewise constant signals")
  

fig9 = plot_grid(p91,p92,p93, nrow = 1)
ggsave("fig9.pdf", fig9, width = 16, height = 6)

##---------------------
## plot for figure 10
## Design 4
##---------------------

setwd("~/git/caisar/result")
dat101 = read.table("dat101.txt", sep = ",")
dat102 = read.table("dat102.txt", sep = ",")
dat103 = read.table("dat103.txt", sep = ",")

p101 = my.box(dat101, "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "IndepLowdimGauss + SparseNormal",
       y     = "relative prediction error") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(0.9,41)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p102 = my.box(dat102,  "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "IndepLowdimGauss + SparseConstant", y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(0.9,41)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p103 = my.box(dat103,  "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "IndepLowdimGauss + SuBogdanCandes",
       y = "") + 
  scale_y_continuous(trans = "log10", breaks = c(1,2,4,8,16,32), limits = c(0.9,41)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

fig10 = plot_grid(p101,p102,p103, nrow = 1)
ggsave("fig10.pdf", fig10, width = 18, height = 6)

##---------------------
## plot for figure 11
## Bad case
##---------------------

setwd("~/git/caisar/result")
dat111 = read.table("dat111.txt", sep = ",")
dat112 = read.table("dat112.txt", sep = ",")
dat113 = read.table("dat113.txt", sep = ",")

p111 = my.box(dat111, "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "IndepLowdimGauss + NonsparseNormal",
       y     = "relative prediction error") + 
  scale_y_continuous(trans = "log10", breaks = 2^(0:7), limits = c(1,128)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p112 = my.box(dat112,  "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "EquiCorrGauss (rho = 0.95) + SparseConstant", y = "") + 
  scale_y_continuous(trans = "log10", breaks = 2^(0:7), limits = c(1,128)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p113 = my.box(dat113,  "fit", "pred", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "IndepLowdimGauss + SparseBimodal",
       y = "") + 
  scale_y_continuous(trans = "log10", breaks = 2^(0:7), limits = c(1,128)) +
  theme(axis.line    = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

fig11 = plot_grid(p111,p113,p112, nrow = 1)
ggsave("fig11.pdf", fig11, width = 18, height = 6)

##---------------------
## plot for figure 12
## Variable selection
##---------------------

setwd("~/git/caisar/result")
dat121 = read.table("dat121.txt", sep = ",")
dat122 = read.table("dat122.txt", sep = ",")
dat123 = read.table("dat123.txt", sep = ",")
dat124 = read.table("dat124.txt", sep = ",")
name   = dat121$fit[seq(1,160,20)]

a      = matrix(dat121$fdr1,20,8)
b      = matrix(dat121$power1,20,8)
fdr1.q = apply(a, 2, function(x) quantile(x, probs = c(0.1,0.5,0.9)))
pow1.q = apply(b, 2, function(x) quantile(x, probs = c(0.1,0.5,0.9)))

p121 = ggplot() + geom_point(aes(x = name, y = fdr1.q[2,], shape = "fdp", color = "fdp"), size = 2) +
  geom_errorbar(aes(x = name, ymin = fdr1.q[1,], ymax = fdr1.q[3,], color = "fdp"), width = 0.2) +
  geom_point(aes(x = name, y = pow1.q[2,], shape = "power", color = "power"), size = 2) +
  geom_errorbar(aes(x = name, ymin = pow1.q[1,], ymax = pow1.q[3,], color = "power"), width = 0.2) +
  theme_cowplot(font_size = 18) +
  labs(x = "", y = "FDP/POWER", title = "EquiCorrGauss + SparseNormal") +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

a      = matrix(dat121$fdr2,20,8)
b      = matrix(dat121$power2,20,8)
fdr2.q = apply(a, 2, function(x) quantile(x, probs = c(0.1,0.5,0.9)))
pow2.q = apply(b, 2, function(x) quantile(x, probs = c(0.1,0.5,0.9)))
p122 = ggplot() + geom_point(aes(x = name, y = fdr2.q[2,], shape = "fdp", color = "fdp"), size = 2) +
  geom_errorbar(aes(x = name, ymin = fdr2.q[1,], ymax = fdr2.q[3,], color = "fdp"), width = 0.2) +
  geom_point(aes(x = name, y = pow2.q[2,], shape = "power", color = "power"), size = 2) +
  geom_errorbar(aes(x = name, ymin = pow2.q[1,], ymax = pow2.q[3,], color = "power"), width = 0.2) +
  theme_cowplot(font_size = 18) +
  labs(x = "", y = "", title = "EquiCorrGauss + ThreePointMass") +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")

a      = matrix(dat121$fdr3,20,8)
b      = matrix(dat121$power3,20,8)
fdr3.q = apply(a, 2, function(x) quantile(x, probs = c(0.1,0.5,0.9)))
pow3.q = apply(b, 2, function(x) quantile(x, probs = c(0.1,0.5,0.9)))
p123 = ggplot() + geom_point(aes(x = name, y = fdr3.q[2,], shape = "fdp", color = "fdp"), size = 2) +
  geom_errorbar(aes(x = name, ymin = fdr3.q[1,], ymax = fdr3.q[3,], color = "fdp"), width = 0.2) +
  geom_point(aes(x = name, y = pow3.q[2,], shape = "power", color = "power"), size = 2) +
  geom_errorbar(aes(x = name, ymin = pow3.q[1,], ymax = pow3.q[3,], color = "power"), width = 0.2) +
  theme_cowplot(font_size = 18) +
  labs(x = "", y = "", title = "IndepLowdimGauss + SparseHeavytail") +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1)) +
  guides(color=guide_legend(title="stat"), shape=guide_legend(title="stat"))

p124 = ggplot(dat122) + geom_point(aes(x = x, y = y)) +
  theme_cowplot(font_size = 18) + 
  labs(x = "coefficient index",
       y = "local false discovery rate",
       title = "local false discovery rate") + 
  theme(axis.line    = element_blank())

p125 = ggplot(dat123) + geom_line(aes(x = x, y = y)) +
  theme_cowplot(font_size = 18) + 
  labs(x = "posterior quantile",
       y = "# of nonzero coefficients",
       title = "sparsity of posterior quantile") + 
  theme(axis.line    = element_blank())

p126 = ggplot(dat124) + geom_point(aes(x = x, y = 1-y)) +
  theme_cowplot(font_size = 18) + 
  labs(x = "coefficient index",
       y = "PIP",
       title = "posterior inclusion probability") + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(axis.line    = element_blank())

fig121 = plot_grid(p121,p122,p123, nrow = 1, rel_widths = c(.36,.36,.45))
fig122 = plot_grid(p126,p124,p125, nrow = 1)
ggsave("fig121.pdf", fig121, width = 15, height = 5)
ggsave("fig122.pdf", fig122, width = 15, height = 5)

##---------------------
## plot for figure 13
## Timing
##---------------------

setwd("~/git/caisar/result")
dat131 = read.table("dat131.txt", sep = ",")
dat133 = read.table("dat133.txt", sep = ",")
dat133$time = dat133$time / dat133$time[1:11]
dat131$n = factor(dat131$n, levels = c(250,500,1000,2000))
dat131$p = dat131$p/10

p133 = my.box(dat133,  "fit", "time", c(1,2,3,4,5,6,0,8,9)) +
  geom_abline(intercept = 0, slope = 0, color = gg_color_hue(8)[3], alpha = 0.8, linetype = 2, size = 0.2) +
  labs(title = "Average elapsed time (11 Scenarios)",
       y = "elapsed time (seconds)") +
  scale_y_continuous(trans = "log10", breaks = c(0.01,0.1,1,10,100), limits = c(0.04,103)) +
  theme(axis.line    = element_blank(),
        axis.text.x  = element_text(angle = 45,hjust = 1),
        legend.position = "none")
p131 = ggplot(dat131) + geom_line(aes(x = p, y = y, color = n)) +
  geom_point(aes(x = p, y = y, color = n, shape = n), size = 3) +
  theme_cowplot(font_size = 18) +
  scale_x_continuous(trans = "log10", breaks = .5 * 2^(1:10)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.line    = element_blank()) +
  labs(x = "p / 1000",
       y = "elapsed time (seconds)",
       title = "Average elasped time for varying (n,p)")
p132 = ggplot() + geom_blank()
fig13 = plot_grid(p133,p132,p131, nrow = 1, rel_widths = c(0.5,0.1,0.5))
ggsave("fig13.pdf", fig13, width = 16, height = 6)
