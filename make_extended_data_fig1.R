suppressPackageStartupMessages({library(ggplot2)})
INK<-"#1A1A1A"; GREY1<-"#AAB2B8"; GREY2<-"#7E868D"
ex <- read.csv("exercise_contrast_data.csv")
ord <- c('Horvath1','Horvath2','Hannum','PCHorvath1','PCHorvath2','PCHannum','DNAmEMRAge','OMICmAge','PhenoAge','GrimAgeV1','GrimAgeV2','PCPhenoAge','PCGrimAge','SystemsAge','DunedinPoAm38','DunedinPACE')
ex$clock <- factor(ex$clock, levels=rev(ord))
all <- rbind(
  data.frame(clock=ex$clock, est=ex$sub1, lo=ex$sub1_lo, hi=ex$sub1_hi, p=NA, which="Low-carbohydrate arm: exercise minus none (n = 30 vs 30)"),
  data.frame(clock=ex$clock, est=ex$sub2, lo=ex$sub2_lo, hi=ex$sub2_hi, p=NA, which="Low-fat arm: exercise minus none (n = 30 vs 30)"),
  data.frame(clock=ex$clock, est=ex$diff, lo=ex$diff_lo, hi=ex$diff_hi, p=ex$diff_p, which="Pooled, inverse variance (95% CI)"))
all$which <- factor(all$which, levels=c("Low-carbohydrate arm: exercise minus none (n = 30 vs 30)",
  "Low-fat arm: exercise minus none (n = 30 vs 30)","Pooled, inverse variance (95% CI)"))
p <- ggplot(all, aes(est, clock, colour=which, shape=which, size=which, linewidth=which)) +
  geom_vline(xintercept=0, linewidth=.35, colour=INK) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0, position=position_dodge(width=.6)) +
  geom_point(fill="white", stroke=.7, position=position_dodge(width=.6)) +
  geom_text(data=all[all$which=="Pooled, inverse variance (95% CI)",], aes(x=.72, label=sprintf("P = %.2f", p)),
            size=2.2, hjust=0, family="Helvetica", colour=INK, show.legend=FALSE) +
  scale_colour_manual(values=c(GREY1, GREY2, INK), name=NULL) +
  scale_shape_manual(values=c(16, 16, 23), name=NULL) +
  scale_size_manual(values=c(1.5, 1.5, 2.3), name=NULL) +
  scale_linewidth_manual(values=c(.45, .45, .75), name=NULL) +
  scale_x_continuous(limits=c(-.75,.9), breaks=seq(-.6,.6,.2), labels=function(x) sprintf("%.1f",x)) +
  labs(x="Exercise minus no exercise, effect on epigenetic age (s.d. units), 95% CI", y=NULL) +
  theme_classic(base_size=8, base_family="Helvetica") +
  theme(axis.line=element_line(linewidth=.3, colour=INK), axis.ticks=element_line(linewidth=.3, colour=INK),
        axis.text=element_text(colour=INK), axis.text.y=element_text(size=7.4),
        axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major.y=element_line(colour="#EEEEEE", linewidth=.3),
        legend.position="top", legend.direction="vertical", legend.justification="left", legend.text=element_text(size=6.6),
        legend.key.size=unit(3,"mm"), legend.margin=margin(b=6))
ggsave("ExtendedDataFig1.pdf", p, width=5.2, height=5.6, device="pdf")
ggsave("ExtendedDataFig1.png", p, width=5.2, height=5.6, dpi=300)
cat("ED done\n")
