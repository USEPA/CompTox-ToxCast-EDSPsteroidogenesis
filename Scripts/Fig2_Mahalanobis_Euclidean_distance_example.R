#-----------------------------------------------------------------------------------#
# Haggard et al. "Using a high-throughput steroidogenesis assay to predict effects on 
# ...estrogen and androgen production or overall perturbation of the steroid biosynthesis pathway"

# Figure 2. Illustrating the difference between Mahalanabos and Eucvlidean distance

# This figure is to show the difference in Mahalanobis and Euclidean distance metrics

#Original 06/22/2017
#updated 08/14/2017

rm(list = ls())

library(ggplot2)
library(ggthemes)
library(ellipse)
library(grid)
library(gridExtra)

setwd("./Misc")

Cor <- matrix(c(1,.8,.8,1),nrow=2, ncol=2)
Cor0 <- matrix(c(1,0,0,1),nrow=2, ncol=2)
Scale <- c(0.08, 0.16)
Scale0 <- c(0.16,0.16)
Cov <- diag(Scale) %*% Cor %*% diag(Scale)
Cov0 <- diag(Scale) %*% Cor0 %*% diag(Scale)
iCov  <- solve(Cov)
MHD <- function(x1,x2,iCov) {
    (x2 - x1) %*% iCov %*% (x2 - x1)
}

ED <- function(x1,x2) {
    (x2 - x1) %*% (x2 - x1)
}

## x, y
z0 <- c(0,0)
z1 <- c(0.15, 0.25)
z2 <- c(0.2915476, 0.0) #z2 has same distance from z0 as z1

pts1 <- as.data.frame(rbind(z0,z1,z2))
pts1 <- sweep(pts1,2,c(-3,-3),"+")
colnames(pts1) <- c("x","y")
pts1T <- pts1
pts1T$labels=c("conc 1","conc 2","conc 3")
TL <- as.data.frame(ellipse(Cor, scale=Scale, centre=unlist(pts1[1,])))
p1 <- ggplot() +
  #geom_point(data=pts1T, aes(x=x,y=y, shape = factor(labels))) +
  geom_text(data=pts1T[2,], aes(x=x,y=y,label=labels), nudge_y = 0.02) +
  geom_text(data=pts1T[3,], aes(x=x,y=y,label=labels), nudge_x = 0.05) +
  geom_text(data=pts1T[1,], aes(x=x,y=y,label=labels), nudge_y = -0.02) +
  geom_path(data=TL, aes(x=x,y=y), linetype = "dashed") +
  geom_segment(aes(x=pts1$x[1], xend=pts1$x[2], y=pts1$y[1], yend=pts1$y[2]),
               arrow=arrow(length=unit(0.125,"inches")), size = 0.6) +
  geom_segment(aes(x=pts1$x[1], xend=pts1$x[3], y=pts1$y[1], yend=pts1$y[3]),
               arrow=arrow(length=unit(0.125,"inches")), size = 0.6) +
  scale_x_continuous("Measured Analyte A (log uM)", limits = c(-3.4, -2.6)) +
  scale_y_continuous("Measured Analyte B (log uM)", limits = c(-3.4, -2.6)) +
  theme_few() +
  #theme(aspect.ratio = 1)
  coord_fixed()


zz <- eigen(iCov)
tr0 <- t(zz$vectors %*% diag(sqrt(zz$values)))

pts2 <- as.data.frame(t(tr0 %*% t(data.matrix(pts1))))
colnames(pts2) <- c("x","y")
pts2T <- pts2
pts2T$labels=c("conc 1","conc 2","conc 3")
TL2 <- as.data.frame(t(tr0 %*% t(data.matrix(TL))))
colnames(TL2) <- c("x","y")
p2 <- ggplot() +
  #geom_point(data=pts2T, aes(x=x,y=y, shape = factor(labels))) +
  geom_text(data=pts2T[2,], aes(x=x,y=y,label=labels), nudge_y = -0.15) +
  geom_text(data=pts2T[3,], aes(x=x,y=y,label=labels), nudge_x = -0.6) +
  geom_text(data=pts2T[1,], aes(x=x,y=y,label=labels), nudge_y = 0.2) +
  geom_path(data=TL2, aes(x=x,y=y), linetype = "dashed") +
  geom_segment(aes(x=pts2$x[1], xend=pts2$x[2], y=pts2$y[1], yend=pts2$y[2]),
               arrow=arrow(length=unit(0.125,"inches")), size = 0.6) +
  geom_segment(aes(x=pts2$x[1], xend=pts2$x[3], y=pts2$y[1], yend=pts2$y[3]),
               arrow=arrow(length=unit(0.125,"inches")), size = 0.6) +
  scale_x_continuous("Rotated and Scaled Axis 1", limits = c(28, 38)) +
  scale_y_continuous("Rotated and Scaled Axis 2", limits = c(20, 30)) +
  theme_few() +
  #theme(aspect.ratio = 1)
  coord_fixed()


pdf("mahalanobis_example.pdf", width = 12, height =8)
multiplot(p1, p2, cols = 2)
dev.off()
