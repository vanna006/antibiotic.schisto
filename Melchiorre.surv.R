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

psurv_diff <- pairwise_survdiff(Surv(time, status) ~ treatment, data = antib, p.adjust.method = 'none')
psurv_diff


#Antibiotic vs Antibiotic + parasite
antib1 <- read.csv("survival1.csv")
fit <- survfit(Surv(time, status) ~ treatment, data = antib1)
print(fit)
summary(fit)

surv_diff <- survdiff(Surv(time, status) ~ treatment, data = antib1)
surv_diff

ggsurvplot(fit, data = antib1,
           xlab = 'Week', ylab = 'Proportion of hosts surviving',
           pval = F, conf.int = TRUE, font.x = 20, font.y = 20, font.legend = 20, font.tickslab = 20, 
           palette = c("#00BFC4", "#edae49"),
           legend.labs = c('Antibiotic', 'Antibiotic + parasite'),
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw() # Change ggplot2 theme
)


#Control vs parasite
antib2 <- read.csv("survival2.csv")
fit <- survfit(Surv(time, status) ~ treatment, data = antib2)
print(fit)
summary(fit)

surv_diff <- survdiff(Surv(time, status) ~ treatment, data = antib2)
surv_diff

ggsurvplot(fit, data = antib2,
           xlab = 'Week', ylab = 'Proportion of hosts surviving',
           pval = F, conf.int = TRUE, font.x = 20, font.y = 20, font.legend = 20, font.tickslab = 20, 
           palette = c("#00B81F", "#d1495b"),
           legend.labs = c('Control', 'Parasite'),
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw() # Change ggplot2 theme
)



#Time to event analysis for cercarial shedding between infection treatments.
shed <- read.csv('shed.csv')
fit <- survfit(Surv(time, status) ~ treatment, data = shed)
print(fit)
summary(fit)

surv_diff <- survdiff(Surv(time, status) ~ treatment, data = shed)
surv_diff

ggsurvplot(fit, data = shed,
           xlab = 'Week', ylab = 'Proportion of infected hosts producing parasites',
           pval = F, conf.int = TRUE, font.x = 20, font.y = 20, font.legend = 20, font.tickslab = 20, 
           palette = c("#edae49", "#d1495b"),
           legend.labs = c('Antibiotic + parasite', 'Parasite'),
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           fun = function(x){1-x} #make plot cumulative.
)
