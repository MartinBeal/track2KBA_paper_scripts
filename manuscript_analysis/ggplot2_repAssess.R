### ggplot of repAssess output ###
library(ggplot2)

# needs output of repAssess(bootTable=TRUE)
repr <- rep
# repr <- represent

reps <- repr[[2]]

Asymptote <- repr[[1]]$est_asym

summ <- reps %>%
  group_by(.data$SampleSize) %>%
  dplyr::summarise(
    meanPred = mean(na.omit(.data$pred)),
    sdInclude = sd(.data$InclusionRate))

yTemp <- c(
  summ$meanPred + Asymptote * summ$sdInclude, 
  rev(summ$meanPred - Asymptote * summ$sdInclude)
)
xTemp <- c(summ$SampleSize, rev(summ$SampleSize))

Temps <- data.frame(xTemp, yTemp)

repPlot <- ggplot() + 
  geom_polygon(aes(x=xTemp, y=yTemp), data=Temps, fill="gray85") +
  geom_point(aes(x=SampleSize, y=InclusionRate), data=reps, col="darkgray", size=0.75, alpha=0.15) +
  geom_line(aes(x=SampleSize, y=meanPred), data=summ, size=1.5) + 
  annotate("text", x=0, y=0.99, label=paste(round(repr[[1]]$out, 1), "%", sep=""), 
    size=6, col="grey30", adj=0) +
  theme(
    axis.text=element_text(size=14, color="black"),
    axis.title=element_text(size=16),
    panel.background=element_rect(fill="white", colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Inclusion") + xlab("Sample Size")
repPlot

## SAVE ## 
# turtles  
repPlot_gt <- repPlot

ggsave("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\figures\\sea_turtles\\rep_h2.18_n23_its800.png", plot = repPlot_gt, width=7, height=6.5)

# vultures 
repPlot_eg <- repPlot

ggsave("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\figures\\egyptian_vultures\\rep_h7.35_n15_its450.png", plot = repPlot_eg, width=7, height=6.5)

# COSH
ggsave("C:\\Users\\Martim Bill\\Documents\\annual_consistency\\figures\\cosh\\long_trips\\rep_h10_n64_its100.png", plot = repPlot, width=7, height=6.5)

# White storks
ggsave("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\figures\\white_storks\\all_represent_n76_h5_its200.png", plot = repPlot, width=7, height=6.5)
