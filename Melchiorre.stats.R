### stats for Melchiorre antibiotic/schisto study

setwd('E:/A_in_arev/Hannah_Melchiorre')

#load in cercarial data
byweek <- read.csv('byweek.csv')

#data is overdispersed
hist(byweek$cercs, breaks = 25)

#visualize cercarial shedding by week and treatment.
library("ggpubr")
ggboxplot(byweek, x = 'week', y = 'cercs', color = 'treat', palette = c("#d1495b", "#edae49"), xlab = 'Week', font.x = 40,
          font.y = 40, font.tickslab = 30, ylab = 'Average cercarial output', size = 1, width = 0.5, add = 'jitter', add.params = list(size = 2))

#for reproducibility
set.seed(123)

#import snail fecundity data
egg.week <- read.csv('egg.byweek.csv')

#data is overdispersed
hist(egg.week$eggs, breaks = 15)

#ap, edae49
#p, d1495b
#a, 00798c
#c, 66a182

####### for just number of cercs, not pres.abs
poshed <- subset(byweek, cercs > 0)
library(lme4)
cercsnb <- glmer.nb(cercs ~ treat*week + (week|ID), data=poshed)
summary(cercsnb)

####### pa
cpa <- read.csv('cerc.pa.csv')
gm3 <- glmer(pres ~ treat * week + (week | ID), data = cpa, family = binomial)
summary(gm3)

#solution using zero-inflated model of cercarial output
library(glmmTMB)
fit_zipoisson1 <- glmmTMB(cercs ~ treat * week + (1|ID), data = byweek, ziformula=~1, family=poisson)
summary(fit_zipoisson1)

#examine various distributions for the best fitting model
fit_hzip <- update(fit_zipoisson1,
                   ziformula = ~.,
                   data = byweek,
                   family = poisson(link="log"))
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom1,family=list(family="truncated_nbinom1",link="log"))

library(bbmle)
AICtab(fit_zipoisson1,fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)
#use best fitting model
summary(fit_zinbinom1)

library(effects)
plot(allEffects(fit_zinbinom1))

#create subsets for fecundity data comparisons
#create data for antibiotic + parasite vs parasite
p <- subset(egg.week, treat == "1p")
ap <- subset(egg.week, treat == "a1p")

p.ap4 <- subset(p.ap, week <= 4)
p.ap5 <- subset(p.ap, week >= 5)


dim(p)
dim(ap)
p.ap <- rbind(p, ap)
dim(p.ap)

p.ap3 <- subset(p.ap, week <= 3)
p.ap4 <- subset(p.ap, week >= 4)

####### for just number of eggs, not pres.abs
p.ap1 <- subset(p.ap, eggs > 0)
eggsnb <- glmer.nb(eggs ~ treat * week + (week|Group), data=p.ap1)
summary(eggsnb)

####### for just eggs pres.abs
epa <- read.csv('egg.pa.csv')
p <- subset(epa, treat == "1p")
ap <- subset(epa, treat == "a1p")
c <- subset(epa, treat == "0c")
a <- subset(epa, treat == "a")
p.ap <- rbind(p, ap)
c.a <- rbind(c, a)

gm1 <- glmer(pres ~ treat*week + (week | Group), data = p.ap, family = binomial)
summary(gm1)
gm2 <- glmer(pres ~ treat + week + (week | Group), data = c.a, family = binomial)
summary(gm2)
plot(allEffects(gm1))

#fecundity in parasite vs antibiotic + parasite for all weeks
fit_zipoisson3 <- glmmTMB(eggs ~ treat * week + (1|Group), data = p.ap, ziformula=~1, family=poisson)
summary(fit_zipoisson3)

#examine other distributions
fit_hzip <- update(fit_zipoisson3,
                   ziformula = ~.,
                   data = p.ap,
                   family = poisson(link="log"))
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom1, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_zipoisson3,fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)
#best fitting model
summary(fit_hnbinom2)

#examine trends in fecundity compensation window
egg.pap <- ggplot(p.ap, aes(x = as.factor(week), y = eggs, fill = treat)) + 
  geom_violin(trim = T) +
  labs(x="Week", y = "Egg masses", hjust=10)

egg.pap  +
  scale_fill_manual(values=c("#d1495b", "#edae49")) + theme_classic() +
  theme(text = element_text(size=40, color = 'black'), legend.title=element_blank())

plot(allEffects(fit_hnbinom2))

#ap, edae49
#p, d1495b
#a, 8d96a3
#c, 66a182

#create data for antibiotic vs control
c <- subset(egg.week, treat == "0c")
a <- subset(egg.week, treat == "a")

dim(c)
dim(a)
c.a <- rbind(c, a)
dim(c.a)

####### for just number of eggs, not pres.abs
c.a1 <- subset(c.a, eggs > 0)
eggsnb <- glmer.nb(eggs ~ treat + week + (week|Group), data=c.a1)
summary(eggsnb)

#antibiotic vs control for all weeks
fit_zipoisson4 <- glmmTMB(eggs ~ treat + week + (1|Group), data = c.a, ziformula=~1, family=poisson)
summary(fit_zipoisson4)

#examine other distirbutions
fit_hzip <- update(fit_zipoisson4,
                   ziformula = ~.,
                   data = c.a,
                   family = poisson(link="log"))
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_hzip, family=list(family="truncated_nbinom2",link="log"))
AICtab(fit_zipoisson4,fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)
#best fitting model with successful convergence
summary(fit_hnbinom2)

#visualize control versus antibiotic fecundity
plot(allEffects(fit_hnbinom2))

egg.ca <- ggplot(c.a, aes(x = treat, y = eggs, fill = treat)) + 
  geom_violin(trim = T) +
  labs(x="", y = "Egg masses (egg masses/snail/week)", hjust=10) +
  theme(text = element_text(size=50),
        axis.text.x = element_text(angle=90, hjust=1))
vioegg <- egg.ca +
  scale_fill_manual(values=c('#00B81F','#00BFC4')) + theme_classic() +
  theme(text = element_text(size=30), legend.title=element_blank(), legend.position = 'none')

egg.caweek <- ggplot(c.a, aes(x = as.numeric(week), y = eggs)) + 
  geom_smooth(method = lm) +
  labs(x="Week", y = "Egg masses (egg masses/snail/week)", hjust=10) +
  theme(text = element_text(size=50),
        axis.text.x = element_text(angle=90, hjust=1))
weekegg <- egg.caweek + theme_classic() +
  theme(text = element_text(size=30), legend.title=element_blank())

ggarrange(vioegg, weekegg, 
          heights = c(2, 2),
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          widths = c(3.5, 2.5))

####### interactions code from Spencer
library(interactions) #for interact_plot function

cercpro <- interact_plot( fit_zinbinom1, pred = week , modx = treat, #m7new is the model
                          main.title="",
                          line.thickness=2,
                          interval = T,
                          x.label = "Week", y.label = "Parasite production (cercariae/snail/hour)",
                          modx.labels=c("Parasite","Antibiotic + parasite")) + 
  scale_fill_manual(values=c("#d1495b", "#edae49")) +
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=14),
        axis.line=element_line(colour="black"),
        legend.position = 'top',
        legend.text = element_text(size = 18),
        legend.title=element_blank())

eggsint <- interact_plot( fit_hnbinom2, pred = week , modx = treat, #m7new is the model
                          main.title="",
                          line.thickness=2,
                          interval = T,
                          x.label = "Week", y.label = "Egg masses (egg masses/snail/week)",
                          modx.labels=c("Parasite","Antibiotic + parasite")) + 
  scale_fill_manual(values=c("#d1495b", "#edae49")) +
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=14),
        axis.line=element_line(colour="black"),
        legend.position = 'top',
        legend.text = element_text(size = 18),
        legend.title=element_blank())


