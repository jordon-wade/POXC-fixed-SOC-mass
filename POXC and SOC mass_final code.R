library(raster)
library(emmeans)
library(rcompanion)
library(dplyr)
library(ggplot2)
library(lme4)
library(multcomp)
library(tidyr)
library(lmerTest)
library(ggtext)
library(ggpubr)

#### LOAD AND WRANGLE DATA ####
POXC.dat <- read.csv("Wade Margenot/POXC and SOC mass/Data/POXC and SOC data_for R.csv", header=TRUE)
POXC.dat$Unique.Soil.ID <- as.factor(POXC.dat$Unique.Soil.ID)
POXC.dat$SOC.mass.class <- as.factor(POXC.dat$SOC.mass.class)
levels(POXC.dat$SOC.mass.class)[levels(POXC.dat$SOC.mass.class)=="std"] <- "Standard (2.5g soil)"
levels(POXC.dat$SOC.mass.class)[levels(POXC.dat$SOC.mass.class)=="25mg"] <- "25 mg SOC"
levels(POXC.dat$SOC.mass.class)[levels(POXC.dat$SOC.mass.class)=="15mg"] <- "15 mg SOC"
POXC.dat$POXC.pct <- 100*(POXC.dat$POXC/(10000*POXC.dat$SOC.pct))
POXC.dat$mol.Mn.reduced <- 1000*POXC.dat$mmol.Mn.reduced
POXC.dat$mol.Mn.g.soil.sample <- (1000*POXC.dat$mmol.Mn.reduced)/(POXC.dat$Soil.mass)
POXC.dat$mol.Mn.mg.SOC.sample <- (1000*POXC.dat$mmol.Mn.reduced)/(POXC.dat$SOC.mass)
POXC.dat$mol.Mn.kg.soil <- (1000*POXC.dat$mmol.Mn.reduced)/(1000*POXC.dat$Soil.mass)

#### POXC: detectable by mass ####
# Note: using ALL data, not averaged by soil
### All SOC contents ###
Table.all <- as.matrix(table(POXC.dat$SOC.mass.class, POXC.dat$Detectable), header=TRUE, row.names=1)
Table.all

Table.pair.all <- pairwiseNominalIndependence(Table.all, compare = "row", fisher=TRUE, gtest=FALSE, chisq=FALSE, method="fdr", digits=3)
Table.pair.all
cldList(p.adj.Fisher ~ Comparison, data=Table.pair.all, threshold=0.1)

### "Low" SOC contents only (< 10%) ###
POXC.dat.lowSOC <- subset(POXC.dat, SOC.pct < 10)
Table.low <- as.matrix(table(POXC.dat.lowSOC$SOC.mass.class, POXC.dat.lowSOC$Detectable), header=TRUE, row.names=1)
Table.low

Table.pair.low <- pairwiseNominalIndependence(Table.low, compare = "row", fisher=TRUE, gtest=FALSE, chisq=FALSE, method="fdr", digits=3)
Table.pair.low
cldList(p.adj.Fisher ~ Comparison, data=Table.pair.low, threshold=0.1)


#### PRIMARY DATA ANALYSIS ####
## Eliminate values outside of theoretical range ##
POXC.data <- subset(POXC.dat, Detectable==1)
str(POXC.data)

## Average the (detectable) analytical replicates by soil*treatment combo ##
POXC.data.by.soil <- as.data.frame(na.omit(POXC.data %>% dplyr::group_by(Unique.Soil.ID, SOC.mass.class, SOC.pct, pH, Clay.pct, Soil.mass) %>% dplyr::summarise(POXC.avg=mean(POXC), CV.POXC=raster::cv(POXC), POXC.pct.avg=mean(POXC.pct), Mn.avg=mean(mol.Mn.reduced), CV.Mn=raster::cv(mol.Mn.reduced), avg.Mn.SOC=mean(mol.Mn.mg.SOC.sample), avg.Mn.mass=mean(mol.Mn.g.soil.sample), avg.mol.Mn.kg=mean(mol.Mn.kg.soil))))

### Table 1: Kolmogorov-Smirnov test results ###
df.15 <- POXC.data.by.soil[POXC.data.by.soil$SOC.mass.class=="15 mg SOC",]
df.25 <- POXC.data.by.soil[POXC.data.by.soil$SOC.mass.class=="25 mg SOC",]
df.std <- POXC.data.by.soil[POXC.data.by.soil$SOC.mass.class=="Standard (2.5g soil)",]

ks.test(df.15$avg.Mn.SOC, df.25$avg.Mn.SOC)
ks.test(df.15$avg.Mn.SOC, df.std$avg.Mn.SOC)
ks.test(df.25$avg.Mn.SOC, df.std$avg.Mn.SOC)

ks.test(df.15$avg.Mn.mass, df.25$avg.Mn.mass)
ks.test(df.15$avg.Mn.mass, df.std$avg.Mn.mass)
ks.test(df.25$avg.Mn.mass, df.std$avg.Mn.mass)

ks.test(df.15$avg.umol.Mn.kg, df.25$avg.umol.Mn.kg)
ks.test(df.15$avg.umol.Mn.kg, df.std$avg.umol.Mn.kg)
ks.test(df.25$avg.umol.Mn.kg, df.std$avg.umol.Mn.kg)

### Figure 1: absolute values of POXC and POXC/SOC ###
## First plot Mn/POXC absolute value data ##
Mn.data.wide <- POXC.data.by.soil[,c(1:3,14)] %>% dplyr::group_by(Unique.Soil.ID, SOC.pct) %>% pivot_wider(names_from=SOC.mass.class, values_from=avg.mol.Mn.kg)
Mn.data.wide$diff.15mg <- Mn.data.wide$`15 mg SOC` - Mn.data.wide$`Standard (2.5g soil)`
Mn.data.wide$diff.25mg <- Mn.data.wide$`25 mg SOC` - Mn.data.wide$`Standard (2.5g soil)`
Mn.data.diff <- Mn.data.wide[,c(1,2,6,7)] %>% pivot_longer(-c(Unique.Soil.ID, SOC.pct), names_to="Category", values_to="Difference")
Mn.data.diff$positive <- Mn.data.diff$Difference >= 0
Mn.data.diff$SOC.pct <- round(as.numeric(Mn.data.diff$SOC.pct), digits=2)
Mn.data.diff$Category <- as.factor(Mn.data.diff$Category)
levels(Mn.data.diff$Category)
levels(Mn.data.diff$Category) <- c("15mg SOC", "25mg SOC")

# How many are greater than the 2.5g (TRUE) and how many are less than (FALSE) for each SOC mass? #
with(subset(Mn.data.diff, SOC.pct<10), table(Category, positive))

Mn.diff.plot.main <- ggplot(data=Mn.data.diff, aes(group=Category)) + geom_point(aes(y=1000*Difference, x=SOC.pct, color=positive), alpha=0.8, size=2) + facet_grid(.~Category) + scale_color_manual(values=c("#DC3220", "#005AB5"), guide=FALSE) + geom_hline(yintercept=0, linetype="dashed") + scale_y_continuous(limits=c(-100, 400), sec.axis=sec_axis(~.*9, name=expression("Difference from Standard (mg POXC"~kg^{"-1"}~"soil)"))) + xlim(0,10) + xlab("SOC (%)") + ylab("Difference from Standard (\u03BCmol MnO<sub>4</sub><sup>-</sup> kg<sup>-1</sup> soil)") + theme_bw() + theme(axis.title.y=element_markdown())
Mn.diff.plot.main

Mn_get_inset <- function(df){
  p <- ggplot(data=df %>% drop_na(positive) %>% group_by(Unique.Soil.ID, Category) %>% dplyr::slice(1), aes(x=SOC.pct, color=positive)) + geom_point(aes(x=SOC.pct, y=1000*Difference), alpha=0.6, size=.75) + scale_x_continuous(breaks=c(0,15,30)) + scale_y_continuous(limits=c(-100 ,1500), breaks=c(0,750, 1500), sec.axis=sec_axis(~.*9)) + scale_color_manual(values=c("#DC3220", "#005AB5"), guide=FALSE) + theme_bw(base_size=9) + theme(panel.background = element_rect(fill="white"), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size=rel(1.1))) + geom_hline(yintercept=0, linetype="dashed")
  return(p)
}

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

Mn.insets <- Mn.data.diff %>% split(f = .$Category) %>% purrr::map(~annotation_custom2(grob=ggplotGrob(Mn_get_inset(.)), data=data.frame(Category=unique(.$Category)), xmin=0, xmax=9, ymin=150, ymax=400))

Mn.diff.plots <- Mn.diff.plot.main + Mn.insets
Mn.diff.plots

## POXC / SOC as % ##
POXC.pct.plot <- ggplot(data=POXC.data.by.soil, aes(x=SOC.mass.class, y=POXC.pct.avg, fill=SOC.mass.class)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.18, alpha=0.6) + scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) + theme_bw() + theme(axis.title.x = element_blank(), legend.position = "right") + ylab("POXC / SOC (%)") + guides(fill=FALSE) + scale_y_continuous(limits=c(-0.5,7.25)) + geom_text(x=1, y=-0.25, label="a", check_overlap = TRUE) + geom_text(x=2, y=-0.25, label="ab", check_overlap = TRUE) + geom_text(x=3, y=-0.25, label="b", check_overlap = TRUE)
POXC.pct.plot

## Combine them with ggpubr ##
Figure1 <- ggarrange(Mn.diff.plots, POXC.pct.plot, labels=c("a", "b"), heights=c(1, 0.8), widths = c(1,0.9) , ncol=1, nrow=2)
Figure1
ggsave("Figure 1.png", width=7, height=7) #export as PDF "Figure 2" at 7 x 7"

POXC.pct.lm <- lm(POXC.pct.avg ~ SOC.mass.class, data=POXC.data.by.soil)
plot(residuals(POXC.pct.lm))
shapiro.test(residuals(POXC.pct.lm))
summary(POXC.pct.lm)
anova(POXC.pct.lm)

em.pct.POXC <- emmeans(POXC.pct.lm, ~ SOC.mass.class) # can change to ".lm" too
em.pct.POXC
cld(em.pct.POXC, alpha = 0.05, reversed = TRUE, Letters = letters)

# Range by treatment #
POXC.data.by.soil %>% dplyr::group_by(SOC.mass.class) %>% dplyr::summarize(Min=min(POXC.pct.avg), Max=max(POXC.pct.avg), diff=(max(POXC.pct.avg)-min(POXC.pct.avg)))

### Figure 2: distribution of Mn reduction within sample ###
## 2a: Mn reduction across SOC content ##
Fig2.a <- ggplot(data=POXC.data.by.soil, aes(y=1000*Mn.avg, x=SOC.pct, color=SOC.mass.class)) + geom_point(alpha=0.6)  + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), guide=guide_legend(title=NULL)) + labs(y=expression("\u03BCmol"~MnO[4]^{"-"}~reduced), x="SOC (%)") + theme_bw() + geom_smooth(se=FALSE, method="loess", span=1.5) + theme(legend.position = "top", legend.title = element_blank(), legend.text=element_text(size=rel(1.1))) + guides(color = guide_legend(override.aes = list(size=rel(1.1))))
Fig2.a

## 2b and c: Mn reduction by SOC and soil mass of sample ##
B <- ddply(POXC.data.by.soil, "SOC.mass.class", summarise, grp.mean=mean(avg.Mn.SOC)) # group-level averages for 2b
C <- ddply(POXC.data.by.soil, "SOC.mass.class", summarise, grp.mean=mean(avg.Mn.mass)) # group-level averages for 2c
Fig2.b <- ggplot(POXC.data.by.soil, aes(x=avg.Mn.SOC, fill=SOC.mass.class)) + geom_density(color=NA, alpha=0.6) + theme_bw() + scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) + geom_segment(data=B, aes(x=grp.mean, y=0, xend=grp.mean, yend=Inf, color=SOC.mass.class), linetype="longdash", size=.6) + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) + labs(x=expression("\u03BCmol"~MnO[4]^{"-"}~reduced~mg^{-1}~SOC~"in"~sample), y="Density") + theme(legend.position="none")
Fig2.b
Fig2.c <- ggplot(POXC.data.by.soil, aes(x=avg.Mn.mass, fill=SOC.mass.class)) + geom_density(color=NA, alpha=0.6) + theme_bw() + scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) + geom_segment(data=C, aes(x=grp.mean, y=0, xend=grp.mean, yend=Inf, color=SOC.mass.class), linetype="longdash", size=.6) + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) + labs(x=expression("\u03BCmol"~MnO[4]^{"-"}~reduced~g^{-1}~soil~"in"~sample), y="Density") + theme(legend.position="none")
Fig2.c

Figure2.bc <- ggarrange(Fig2.b, Fig2.c, ncol=2, labels=c("b", "c"))
Figure2.abc <- ggarrange(Fig2.a, Figure2.bc, nrow=2, labels=c("a", ""))
Figure2.abc
ggsave("Figure 2.png", width=6.5, height=6) # export as PDF "Figure 2" at 8" x 6"

### Figure 3: Distribution of POXC values (as Mn) ###
## Distribution of umol Mn / kg soil ##
abs.Mn.plot <- ggplot(POXC.data.by.soil, aes(x=avg.umol.Mn.kg, fill=SOC.mass.class)) + geom_density(color=NA, alpha=0.6) + theme_bw() + scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), guide=guide_legend(title=NULL)) + labs(x=expression("POXC (\u03BCmol"~MnO[4]^{"-"}~reduced~kg^{-1}~"soil)"), y="Density") + theme(legend.key=element_blank())
abs.Mn.plot

## Kolmogorov-Smirnov plot of ECDF ##
ECDF.plot <- ggplot(data=POXC.data.by.soil, aes(x=avg.umol.Mn.kg, group=SOC.mass.class, color=SOC.mass.class)) + stat_ecdf(size=1) + theme_bw() +  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) + theme(legend.title=element_blank(), legend.position = "top") + labs(x=expression("POXC (\u03BCmol"~MnO[4]^{"-"}~reduced~kg^{-1}~"soil)"), y="Cumulative\nDistribution")
ECDF.plot

Figure3.ab <- ggpubr::ggarrange(abs.Mn.plot, ECDF.plot, nrow = 2, labels=c("a", "b"), heights=c(1.1, 0.8), widths=c(1,1), common.legend=TRUE, legend = "top", align="v")
Figure3.ab
ggsave("Figure 3.png", width=5.5, height=4.5) # export as PDF "Figure 3"

### VARIABILITY: by-soil basis ###
set.seed(12345)
By.soil.fit.lmer <- lmerTest::lmer(log(CV.POXC) ~ SOC.mass.class + (1|Unique.Soil.ID), data=POXC.data.by.soil, REML=TRUE)
plot(residuals(By.soil.fit.lmer))
shapiro.test(residuals(By.soil.fit.lmer))
summary(By.soil.fit.lmer)
anova(By.soil.fit.lmer, ddf="Satterthwaite")

Cov.lmer <- lmer(log(CV.POXC) ~ SOC.mass.class + SOC.pct + (1|Unique.Soil.ID), data=POXC.data.by.soil, REML=TRUE)
anova(Cov.lmer, By.soil.fit.lmer) # compare fits: no covariates help (p < 0.20)


emmeans(By.soil.fit.lmer,  ~  SOC.mass.class)
By.mass.soil <- emmeans(By.soil.fit.lmer,  ~  SOC.mass.class)
CLD(By.mass.soil, details = FALSE, sort = TRUE, reversed = TRUE, alpha = 0.01, Letters=letters)

## Mass effect on CV by treatment ##
cor.test(subset(POXC.data.by.soil, SOC.mass.class=="15 mg SOC")$Soil.mass, subset(POXC.data.by.soil, SOC.mass.class=="15 mg SOC")$CV.POXC, method="pearson")
cor.test(subset(POXC.data.by.soil, SOC.mass.class=="25 mg SOC")$Soil.mass, subset(POXC.data.by.soil, SOC.mass.class=="25 mg SOC")$CV.POXC, method="pearson")

## Mn reduced and soil pH ##
cor.test(subset(POXC.data.by.soil, SOC.mass.class=="15 mg SOC")$Mn.avg, subset(POXC.data.by.soil, SOC.mass.class=="15 mg SOC")$pH, method="pearson")
cor.test(subset(POXC.data.by.soil, SOC.mass.class=="25 mg SOC")$Mn.avg, subset(POXC.data.by.soil, SOC.mass.class=="25 mg SOC")$pH, method="pearson")

## Mn reduction range ##
POXC.data.by.soil %>% dplyr::group_by(SOC.mass.class) %>% dplyr::summarize(Mean=mean(Mn.avg), Min=min(Mn.avg), Max=max(Mn.avg), diff=(max(Mn.avg)-min(Mn.avg)))
