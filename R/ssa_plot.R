#########################################################################
## Project: Secondary Structure Analysis (SSA)
## Script purpose: Visualize Secondary Structure Analysis
## Date: Dec 19, 2020
## Author: Umut Gerlevik
#########################################################################
#' @title ggstride
#'
#' @description ssa_plot() returns a plot that contains percentage of
#' secondary structures observed along the trajectories per residue by
#' comparing the two ssa data frames obtained by ssa().
#'
#' @param ssa1,ssa2,name1,color_number1,name2,color_number2
#'
#' @return plot
#'
#' @examples
#' ssa_plot(ssa1 = ssa_WT, ssa2 = ssa_S50P,
#' name1 = "Wild-type", color_number1 = 1,
#' name2 = "Ser50Pro", color_number2 = 2)
#'
#' ggstride::ssa_plot(ssa1, ssa2, name1 = "Wild-type", color_number1 = 1, name2, color_number2)
#'
#' @export
#'
#' @import "ggplot2", "ggpubr", "ggsci", "reshape2", "stats"

# Requirements:
# Two ssa() returned data frames (ssa1 and ssa2)
# R packages: "ggplot2", "ggpubr", "ggsci", "reshape2"

# ssa_plot() will return a plot that contains
# percentage of secondary structures
# observed along the trajectories per residue
# by comparing the two ssa_df of the trajectories
# obtained by ssa()

# ssa1: the ssa data frame calculated for trajectory 1
# ssa2: the ssa data frame calculated for trajectory 1
# name1: name of your trajectory 1
# name2: name of your trajectory 2
# resid1: residue id of the first residue of the
	# region you want to visualize
# resid2: residue id of the last residue of the
	# region you want to visualize
# color_number1 (between 1 and 7):
	# select a color from jama palette of ggsci
	# for the ssa of trajectory 1
# color_number2 (between 1 and 7):
	# select a color from jama palette of ggsci
	# for the ssa of trajectory 2

# Example (do not run):
# plot_S50P <- ssa_plot(ssa1 = ssa_WT, ssa2 = ssa_S50P,
#						name1 = "Wild-type", color_number1 = 1,
#						name2 = "Ser50Pro", color_number2 = 2)
# plot_S50P

ssa_plot <- function(ssa1, ssa2,
                     name1 = "Wild-type", name2,
                     resid1 = as.numeric(gsub('\\D+','', colnames(ssa1)[1])),
                     resid2 = as.numeric(gsub('\\D+','', colnames(ssa1)[length(colnames(ssa1))])),
                     color_number1 = 1, color_number2 = 2) {
  palet <- ggsci::pal_jama()
  colorPalet <- palet(7)

  resids <- resid1:resid2

  out <- matrix(ncol = 7, nrow = length(resids))
  # H: AlphaHelix, E: ExtendedConformation(BetaSheet), B or b: Bridge,
  # T: Turn, C or " ": Coil, G: Helix310, I: PiHelix
  colnames(out) <- c("AlphaHelix", "BetaSheet", "Bridge", "Turn",
                     "Coil", "Helix310", "PiHelix")

  out1 <- as.data.frame(na.omit(out))
  for(i in resids) {
    perc <- (table(ssa1[i]) / sum(table(ssa1[i])))*100

    out1[i, "AlphaHelix"] <- perc["H"]
    out1[i, "BetaSheet"] <- perc["E"]
    out1[i, "Bridge"] <- sum(perc[c("B","b")])
    out1[i, "Turn"] <- perc["T"]
    out1[i, "Coil"] <- perc["C"]
    out1[i, "Helix310"] <- perc["G"]
    out1[i, "PiHelix"] <- perc["I"]
  }
  out2 <- as.data.frame(na.omit(out))
  for(i in resids) {
    perc <- (table(ssa2[i]) / sum(table(ssa2[i])))*100

    out2[i, "AlphaHelix"] <- perc["H"]
    out2[i, "BetaSheet"] <- perc["E"]
    out2[i, "Bridge"] <- sum(perc[c("B","b")])
    out2[i, "Turn"] <- perc["T"]
    out2[i, "Coil"] <- perc["C"]
    out2[i, "Helix310"] <- perc["G"]
    out2[i, "PiHelix"] <- perc["I"]
  }

  out1 <- out1[resid1:resid2, ]
  out1$resid <- resids
  out1 <- reshape2::melt(out1[1:8], id.var = "resid")
  out1$mut <- name1

  out2 <- out2[resid1:resid2, ]
  out2$resid <- resids
  out2 <- reshape2::melt(out2[1:8], id.var = "resid")
  out2$mut <- name2

  out <- rbind(out1, out2)

  require(ggplot2)
  require(ggpubr)
  p1 <- ggplot(out[out$variable == "AlphaHelix",], aes(x = resid, y = value, fill = mut,))
  p1 <- p1 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p1 <- p1 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p1 <- p1 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p1 <- p1 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p1 <- p1 + labs(x = "Residue Index", y = "Alpha Helix (%)")
  p1 <- p1 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p1 <- p1 + theme_minimal()
  p1 <- p1 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

  p2 <- ggplot(out[out$variable == "BetaSheet",], aes(x = resid, y = value, fill = mut,))
  p2 <- p2 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p2 <- p2 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p2 <- p2 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p2 <- p2 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p2 <- p2 + labs(x = "Residue Index", y = "Beta Sheet (%)")
  p2 <- p2 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p2 <- p2 + theme_minimal()
  p2 <- p2 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

  p3 <- ggplot(out[out$variable == "Bridge",], aes(x = resid, y = value, fill = mut,))
  p3 <- p3 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p3 <- p3 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p3 <- p3 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p3 <- p3 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p3 <- p3 + labs(x = "Residue Index", y = "Beta Bridge (%)")
  p3 <- p3 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p3 <- p3 + theme_minimal()
  p3 <- p3 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

  p4 <- ggplot(out[out$variable == "Turn",], aes(x = resid, y = value, fill = mut,))
  p4 <- p4 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p4 <- p4 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p4 <- p4 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p4 <- p4 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p4 <- p4 + labs(x = "Residue Index", y = "Turn (%)")
  p4 <- p4 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p4 <- p4 + theme_minimal()
  p4 <- p4 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

  p5 <- ggplot(out[out$variable == "Coil",], aes(x = resid, y = value, fill = mut,))
  p5 <- p5 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p5 <- p5 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p5 <- p5 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p5 <- p5 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p5 <- p5 + labs(x = "Residue Index", y = "Coil (%)")
  p5 <- p5 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p5 <- p5 + theme_minimal()
  p5 <- p5 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

  p6 <- ggplot(out[out$variable == "Helix310",], aes(x = resid, y = value, fill = mut,))
  p6 <- p6 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p6 <- p6 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p6 <- p6 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p6 <- p6 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p6 <- p6 + labs(x = "Residue Index", y = "3(10) Helix (%)")
  p6 <- p6 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p6 <- p6 + theme_minimal()
  p6 <- p6 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

  p7 <- ggplot(out[out$variable == "PiHelix",], aes(x = resid, y = value, fill = mut,))
  p7 <- p7 + geom_bar(stat = "identity", position = position_dodge(), width = 0.45)
  if(length(resids) >= 50) {
    p7 <- p7 + scale_x_continuous(breaks = c(resids[1], seq(10, resids[length(resids)], 10), resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  } else {
    p7 <- p7 + scale_x_continuous(breaks = c(resids[1]:resids[length(resids)]),
                                  expand = c(0,0), minor_breaks = resids)
  }
  p7 <- p7 + scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0.1), limits = c(0,100))
  p7 <- p7 + labs(x = "Residue Index", y = "Pi Helix (%)")
  p7 <- p7 + scale_fill_manual(values = colorPalet[c(color_number2, color_number1)])
  p7 <- p7 + theme_minimal()
  p7 <- p7 + theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 90, size = 18, face = "bold"),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 14),
                   axis.text.x = element_text(angle = 90, size = 10),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = "none",
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))

  p8 <- ggpubr::as_ggplot(ggpubr::get_legend(p1 + theme(legend.position = "top")))


  fig <- ggarrange(p8, p1, p2, p3, p4, p5, p6, p7, ncol = 1,
                   heights = c(1, rep(5, 7)))

  return(annotate_figure(fig, bottom = text_grob("Residue Index", face = "bold", size = 18)))
}

