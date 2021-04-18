library("survival")
library("survminer")

setwd('E:/A_in_arev/Hannah_Melchiorre')
#check for difference in infection prevalence
#19/32 infected in p
#27/35 infected in ap
res <- prop.test(x = c(19, 27), n = c(32, 35))
res

#Import survival data
antib <- read.csv("survival.csv")

#fit model to generate survival curve
fit <- survfit(Surv(time, status) ~ treatment, data = antib)
print(fit)
summary(fit)

#ap, edae49
#p, d1495b
#a, 00BFC4
#c, 00B81F

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           palette = c('#00B81F','#00BFC4' ,'#edae49', '#d1495b'),
           risk.table = TRUE, # Add risk table
           break.x.by = 3,
           xlim = c(0, 9),
           legend.labs = c('Control', 'Antibiotic', 'Antibiotic + parasite', 'Parasite'),
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
            )

#Cox proportional hazard model
cph <- coxph(Surv(time, status) ~ treatment, data = antib)
cph
summary(cph)

#make pairwise comparisons
library(emmeans)
multzip <- emmeans(cph, list(pairwise ~ treatment), adjust = 'none')
multzip
