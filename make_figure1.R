# Figure 1, ggplot2 (R). Data from the article's own Supplementary Table 1.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

CDIS<-"#96393B"; CREF<-"#B08A2E"; CLON<-"#33536E"; CNS<-"#CBD1D6"; INK<-"#1A1A1A"; GREY<-"#8A929A"
fmt <- function(x) sprintf("%.1f", x)
fmt2 <- function(x) sprintf("%.2f", x)
base <- theme_classic(base_size=8, base_family="Helvetica") +
  theme(axis.line=element_line(linewidth=.3, colour=INK),
        axis.ticks=element_line(linewidth=.3, colour=INK),
        axis.text=element_text(colour=INK),
        plot.title=element_text(face="bold", size=8, hjust=0, margin=margin(b=2)),
        plot.title.position="plot",
        plot.margin=margin(2,4,2,2))

a <- read.csv("fig_data_a.csv")
XL <- c(-.37, .34)
block <- function(grp, title, show_x, annotate_ex=FALSE){
  d <- a[a$group==grp, ]
  d <- if(grp=="inc") d[order(d$mean, decreasing=TRUE), ] else d[order(d$mean), ]
  d$label <- factor(d$label, levels=rev(d$label))   # first row ends up at TOP
  d$cat <- factor(d$cat, levels=c("dis","ref","lon","ns"))
  face <- ifelse(d$cat[match(levels(d$label), d$label)]!="ns" |
                 d$notable[match(levels(d$label), d$label)]==1, "bold", "plain")
  col  <- with(d[match(levels(d$label), d$label), ],
               ifelse(cat=="dis", CDIS, ifelse(cat=="ref", CREF,
               ifelse(cat=="lon", INK, ifelse(notable==1, INK, GREY)))))
  ncol_ink <- with(d, ifelse(cat!="ns" | notable==1, INK, GREY))
  d$fillcat <- if(grp=="ns") factor("ns", levels=c("dis","ref","lon","ns")) else d$cat
  p <- ggplot(d, aes(mean, label, fill=fillcat)) +
    geom_col(width=.66) +
    geom_vline(xintercept=0, linewidth=.35, colour=INK) +
    geom_text(aes(x=.335, label=paste0("n = ", n)), hjust=1, size=1.95,
              colour=ncol_ink, family="Helvetica") +
    scale_fill_manual(values=c(dis=CDIS, ref=CREF, lon=INK, ns=CNS),
      breaks=c("dis","ref","lon","ns"),
      labels=c("Treatment of active disease",
               "Control or sham arm",
               "Other intervention","Not significant"), drop=FALSE,
      guide=guide_legend(ncol=2, byrow=TRUE)) +
    scale_x_continuous(limits=XL, breaks=seq(-.3,.3,.1), labels=fmt,
                       expand=expansion(mult=c(0,.01))) +
    labs(x=if(show_x) "Mean effect across 16 clocks (s.d. units)" else NULL,
         y=NULL, title=title) +
    base +
    theme(axis.text.y=element_text(size=6, face=face, colour=col),
          axis.line.y=element_blank(), axis.ticks.y=element_blank(),
          legend.position="none")
  if(!show_x) p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                             axis.line.x=element_blank())
  if(annotate_ex){
    ex <- d[d$label=="Exercise", ]
    p <- p + geom_text(data=ex, aes(x=-.365, y=label),
      label="0 of 108 biomarkers\nat nominal P < 0.05",
      size=2.0, hjust=0, lineheight=.92, fontface="italic", colour=INK, inherit.aes=FALSE)
  }
  p
}
pa1 <- block("dec", NULL, FALSE) +
  labs(subtitle="Decreased epigenetic age  (13 of 51)") +
  theme(plot.subtitle=element_text(face="bold", size=7.4, margin=margin(b=2)))
pa2 <- block("inc", NULL, FALSE) + labs(subtitle="Increased epigenetic age  (3 of 51)") +
  theme(plot.subtitle=element_text(face="bold", size=7.4, margin=margin(b=2)))
pa3 <- block("ns", NULL, TRUE, annotate_ex=TRUE) + labs(subtitle="Not significant after the correction  (35 of 51)") +
  theme(plot.subtitle=element_text(face="bold", size=7.4, margin=margin(b=2)))
legdat <- data.frame(
  x=c(.12,5.1,.12,5.1), y=c(1,1,0,0),
  cat=c("dis","ref","lon","ns"),
  lab=c("Treatment of active disease",
        "Control or sham arm",
        "Other intervention","Not significant"))
alegend <- ggplot(legdat) +
  geom_tile(aes(x, y, fill=cat), width=.18, height=.42, show.legend=FALSE) +
  geom_text(aes(x+.16, y, label=lab), hjust=0, size=2.3, family="Helvetica", colour=INK) +
  scale_fill_manual(values=c(dis=CDIS, ref=CREF, lon=INK, ns=CNS)) +
  scale_x_continuous(limits=c(-.05, 10.4), expand=c(0,0)) +
  scale_y_continuous(limits=c(-.45, 1.45), expand=c(0,0)) +
  labs(title="a") + theme_void(base_family="Helvetica") +
  theme(plot.title=element_text(face="bold", size=10, family="Helvetica", margin=margin(b=4)),
        plot.title.position="plot", plot.margin=margin(2,4,0,2))

# ---------------- Panel b, v2: three randomized trials, treated vs control vs difference ----------------
d2 <- read.csv("fig1b_v2_data.csv")
d2 <- d2[d2$clock %in% c("DunedinPACE","PCGrimAge"), ]
long <- do.call(rbind, lapply(c("treated","control","diff"), function(s){
  data.frame(trial=d2$trial, clock=d2$clock, series=s,
             est=d2[[s]], lo=d2[[paste0(s,"_lo")]], hi=d2[[paste0(s,"_hi")]],
             p=d2[[paste0(s,"_p")]], n_t=d2$n_treated, n_c=d2$n_control)
}))
long$series <- factor(long$series, levels=c("diff","control","treated"),
  labels=c("Treated minus control","Control arm, pre-post change","Treated arm, pre-post change"))
long$clock <- factor(long$clock, levels=c("PCGrimAge","DunedinPACE"))
long$plab <- ifelse(long$series=="Treated minus control",
                    ifelse(long$p<0.001, "P < 0.001", sprintf("P = %.2f", long$p)), "")
ttl <- c("DIRECT-PLUS"="DIRECT-PLUS: green Mediterranean vs guideline diet\n(n = 87 vs 88)",
         "CENTRAL"="CENTRAL: low-carbohydrate vs low-fat diet\n(n = 60 vs 60)",
         "Twin study"="Twin study: vegan vs omnivorous diet\n(n = 13 vs 12)")
XLB <- c(-.62, .60)
ARMS <- list("DIRECT-PLUS"=c("green Mediterranean","guideline diet"), "CENTRAL"=c("low-carbohydrate","low-fat"), "Twin study"=c("vegan","omnivorous"))
sub <- function(tr, show_x, title_letter=NULL){
  d <- long[long$trial==tr, ]
  d$armlab <- ifelse(d$series=="Treated arm, pre-post change", ARMS[[tr]][1],
              ifelse(d$series=="Control arm, pre-post change", ARMS[[tr]][2], ""))
  p <- ggplot(d, aes(est, clock, colour=series, shape=series)) +
    geom_label(aes(x=ifelse(hi > -.08 & hi < .03, .03, hi+.03), label=armlab), size=1.85, hjust=0, family="Helvetica", fill="white",
              label.size=0, label.padding=unit(.08,"lines"), label.r=unit(0,"lines"),
              position=position_dodge(width=.62), show.legend=FALSE) +
    geom_vline(xintercept=0, linewidth=.35, colour=INK) +
    geom_errorbarh(aes(xmin=lo, xmax=hi), height=0, linewidth=.55,
                   position=position_dodge(width=.62)) +
    geom_point(aes(size=series, stroke=series), fill="white", position=position_dodge(width=.62)) +
    scale_size_manual(values=c(2.4,1.9,1.9), guide="none") + scale_discrete_manual("stroke", values=c(1.0,.7,.7), guide="none") +
    geom_text(aes(x=.45, label=plab), size=1.9, hjust=0, family="Helvetica", fontface="bold",
              colour=INK, position=position_dodge(width=.62), show.legend=FALSE) +
    scale_colour_manual(values=c(INK, CREF, INK), drop=FALSE) +
    scale_shape_manual(values=c(23, 16, 16), drop=FALSE) +
    scale_x_continuous(limits=XLB, breaks=seq(-.6,.4,.2), labels=fmt) +
    labs(x=if(show_x) "Effect on epigenetic age (s.d. units), 95% CI" else NULL,
         y=NULL, subtitle=ttl[[tr]]) +
    base + theme(legend.position="none",
                 plot.subtitle=element_text(face="bold", size=7.2, lineheight=.95, margin=margin(b=2)),
                 axis.text.y=element_text(size=7.2), plot.margin=margin(2,4,0,2),
                 axis.line.y=element_blank(), axis.ticks.y=element_blank(),
                 panel.grid.major.y=element_blank())
  if(!is.null(title_letter)) p <- p + labs(title=title_letter) +
    theme(plot.title=element_text(face="bold", size=10, margin=margin(b=6)))
  if(!show_x) p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                             axis.line.x=element_blank())
  p
}
pb1 <- sub("DIRECT-PLUS", FALSE); pb2 <- sub("CENTRAL", FALSE); pb3 <- sub("Twin study", FALSE)

# ---------------- Panel c: CENTRAL exercise vs no exercise ----------------
ex <- read.csv("exercise_contrast_data.csv"); ex <- ex[ex$clock %in% c("DunedinPACE","PCGrimAge"), ]
exl <- do.call(rbind, lapply(c("treated","control","diff"), function(s){
  data.frame(clock=ex$clock, series=s, est=ex[[s]], lo=ex[[paste0(s,"_lo")]], hi=ex[[paste0(s,"_hi")]],
             p=if(s=="diff") ex$diff_p else NA)}))
exl$series <- factor(exl$series, levels=c("diff","control","treated"),
  labels=c("Treated minus control","Control arm, pre-post change","Treated arm, pre-post change"))
exl$clock <- factor(exl$clock, levels=c("PCGrimAge","DunedinPACE"))
exl$plab <- ifelse(exl$series=="Treated minus control", sprintf("P = %.2f", exl$p), "")
exl$armlab <- ifelse(exl$series=="Treated arm, pre-post change","exercise",
             ifelse(exl$series=="Control arm, pre-post change","no exercise",""))
pc <- ggplot(exl, aes(est, clock, colour=series, shape=series)) +
  geom_label(aes(x=ifelse(hi > -.08 & hi < .03, .03, hi+.03), label=armlab), size=1.85, hjust=0, family="Helvetica", fill="white",
            label.size=0, label.padding=unit(.08,"lines"), label.r=unit(0,"lines"),
            position=position_dodge(width=.62), show.legend=FALSE) +
  geom_vline(xintercept=0, linewidth=.35, colour=INK) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0, linewidth=.55, position=position_dodge(width=.62)) +
  geom_point(aes(size=series, stroke=series), fill="white", position=position_dodge(width=.62)) +
  scale_size_manual(values=c(2.4,1.9,1.9), guide="none") + scale_discrete_manual("stroke", values=c(1.0,.7,.7), guide="none") +
  geom_text(aes(x=.45, label=plab), size=1.9, hjust=0, family="Helvetica", fontface="bold", colour=INK,
            position=position_dodge(width=.62), show.legend=FALSE) +
  scale_colour_manual(values=c(INK, CREF, INK), drop=FALSE) +
  scale_shape_manual(values=c(23, 16, 16), drop=FALSE) +
  scale_x_continuous(limits=XLB, breaks=seq(-.6,.4,.2), labels=fmt) +
  labs(x="Effect on epigenetic age (s.d. units), 95% CI", y=NULL, title="c",
       subtitle="CENTRAL: exercise vs no exercise, pooled across diets\n(n = 60 vs 60)") +
  base + theme(legend.position="none",
               plot.title=element_text(face="bold", size=10, margin=margin(b=6)),
               plot.subtitle=element_text(face="bold", size=7.2, lineheight=.95, margin=margin(b=2)),
               axis.text.y=element_text(size=7.2), axis.line.y=element_blank(), axis.ticks.y=element_blank(),
               panel.grid.major.y=element_blank())

blegdat <- data.frame(x=c(.12,2.3,.12), y=c(1,1,0), s=c("t","c","d"),
  lab=c("Treated arm","Control arm","Treated minus control, 95% CI"))
blegend <- ggplot(blegdat) +
  geom_point(aes(x, y, colour=s, shape=s, size=s), fill="white", stroke=.8, show.legend=FALSE) +
  scale_size_manual(values=c(t=1.9,c=1.9,d=2.4)) +
  geom_text(aes(x+.25, y, label=lab), hjust=0, size=2.3, family="Helvetica", colour=INK) +
  scale_colour_manual(values=c(t=INK, c=CREF, d=INK)) + scale_shape_manual(values=c(t=16,c=16,d=23)) +
  scale_x_continuous(limits=c(-.05, 6), expand=c(0,0)) +
  scale_y_continuous(limits=c(-.6, 1.6), expand=c(0,0)) +
  labs(title="b") + theme_void(base_family="Helvetica") +
  theme(plot.title=element_text(face="bold", size=10, family="Helvetica", margin=margin(b=4)),
        plot.title.position="plot", plot.margin=margin(2,4,0,2))
left <- wrap_elements(full=ggplotGrob(alegend))/pa1/pa2/pa3 + plot_layout(heights=c(2.4,13,3,38))
right <- wrap_elements(full=ggplotGrob(blegend))/pb1/pb2/pb3/pc + plot_layout(heights=c(1.5,6,5,5,6.5))
fig <- (left | right) + plot_layout(widths=c(1.18, 1))
ggsave("Figure1.pdf", fig, width=7.2, height=8.6, device="pdf")
ggsave("Figure1.png", fig, width=7.2, height=8.6, dpi=300)
cat("done\n")
