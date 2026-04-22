##This scrip uses the matched_did_data to conduct the dif-in-diff analysis

#figure 3 Difference in difference matched####
library(tidyverse)
require(lme4)
require(lmerTest)
require(sjPlot)
require(gridExtra)
require(ggpubr)
require(fixest)

library(broom.mixed)
library(readr)

library(insight)

set.seed(0)

#model test 
#sensitivity 

#read in matched data
matched.all.data <- read.csv('data_processed/matched_did_data.csv')

#relevel factirs
matched.all.data$group <- fct_relevel(matched.all.data$group, 'unburned', 'fire')
matched.all.data$fire.occ <- fct_relevel(matched.all.data$fire.occ,'before', 'after')

matched.all.data$group.fire.occ <- paste(matched.all.data$group, matched.all.data$fire.occ, sep = '-')

matched.all.data <- matched.all.data%>%
  filter(!is.na(sens))%>%
  mutate(group.num = ifelse(group == 'fire', 1, 0),
         fire.num = ifelse(fire.occ == 'after', 1, 0))

#DID tests sensitivity and efficeincy 
m.did.sens.match <-  lmer(data = matched.all.data, sens ~ fire.occ + group + fire.occ:group +   (1|pixel.id))
summary(m.did.sens.match)
confint(m.did.sens.match)


m.did.pue.match <-  lmer(data = matched.all.data, pue ~ fire.occ + group + fire.occ:group +   (1|pixel.id))
summary(m.did.pue.match)
confint(m.did.pue.match)



tab_model(m.did.sens.match, m.did.pue.match,
          pred.labels =  c('Intercept', 'Time (after)', 'Group (fire)', 'Time (after):Group (fire)'),
          dv.labels = c('Sensitivty', 'Efficiency' ))

### sensitivity did tests



m.did.sens.sensemakr <-  lm(data = matched.all.data, sens ~ group.num + fire.num + group.num:fire.num + map)

summary(m.did.sens.sensemakr)

emmeans::emmeans(m.did.sens.sensemakr, ~ group.num + fire.num + group.num:fire.num )


summary(m.did.sens.match)

m.did.pue.sensemakr <-  lm(data = matched.all.data, pue ~ group.num + fire.num + group.num:fire.num + map)

summary(m.did.pue.sensemakr)

emmeans::emmeans(m.did.sens.sensemakr, ~ group.num + fire.num + group.num:fire.num )


summary(m.did.sens.match)


emmeans::emmeans(m.did.sens.match, ~ fire.occ + group )
emmeans::emmeans(m.did.sens.match, pairwise~ fire.occ + group )
(223-214)/214
(208-161)/208

emmeans::emmeans(m.did.pue.match,  ~ fire.occ + group )

0.047/0.272
sensitivity <- sensemakr(model = m.did.sens.sensemakr, 
                         treatment = "group.num:fire.num",
                         benchmark_covariates = "map",
                         kd = 1:3)
summary(sensitivity)
plot(sensitivity, type = "contour")
m.did.sens.sensemakr$coefficients

#plot DID tests
#get model data 
plt.m.sens <- lmer(data = matched.all.data, sens ~ -1 + group.fire.occ +  (1|pixel.id))
sens.plt.m.output <- summary(plt.m.sens)$coefficients%>%
  cbind(as.data.frame(confint(plt.m.sens))[c(3:6),])%>%
  rename(low.ci = `2.5 %`,
         high.ci = `97.5 %`)%>%
  mutate(rwnms = str_sub(rownames(.), start = 15))%>%
  separate(rwnms, into = c('group', 'fire.occ'))

did.sens.plt <- ggplot(sens.plt.m.output, aes(x = fct_relevel(fire.occ,'before', 'after'), y = Estimate, color = group, shape = group))+
  # geom_point(aes(x = fire.occ, y = sens, color = group), alpha = 0.05)+
  #geom_jitter(alpha = 0.2)+
  geom_point(position = position_dodge(0.025), size = 2.5)+
  geom_line(aes(group = group), position = position_dodge(0.025))+
  geom_segment( aes(y = low.ci, yend = high.ci ),position = position_dodge(0.025))+
  scale_color_manual(values = c('tomato', "#0072B2"))+
  theme_bw(base_size = 13)+
  theme(legend.position = 'none')+
  xlab('')+
  ylab(expression("Sensitivity g m"^{-2}~"mm"^{-1}))

plt.m.pue <- lmer(data = matched.all.data, pue ~ -1 + group.fire.occ +  (1|pixel.id))
pue.plt.m.output <- summary(plt.m.pue)$coefficients%>%
  cbind(as.data.frame(confint(plt.m.pue))[c(3:6),])%>%
  rename(low.ci = `2.5 %`,
         high.ci = `97.5 %`)%>%
  mutate(rwnms = str_sub(rownames(.), start = 15))%>%
  separate(rwnms, into = c('group', 'fire.occ'))

did.pue.plt <- ggplot(pue.plt.m.output, aes(x = fct_relevel(fire.occ,'before', 'after'), y = Estimate, color = group, shape = group))+
  # geom_point(aes(x = fire.occ, y = sens, color = group), alpha = 0.05)+
  #geom_jitter(alpha = 0.2)+
  geom_line(aes(group = group), position = position_dodge(0.025))+
  geom_point(position = position_dodge(0.025), size = 2.5)+
  geom_segment( aes(y = low.ci, yend = high.ci ),position = position_dodge(0.025))+
  scale_color_manual( values = c('tomato', "#0072B2"))+
  theme_bw(base_size = 12)+
  theme(legend.position = 'none')+
  xlab('')+
  ylab(expression("Efficiency (PUE) g m"^{-2}~"mm"^{-1}))

did.lgd.plt <- ggplot(pue.plt.m.output, aes(x = fct_relevel(fire.occ,'before', 'after'), y = Estimate, color = group, shape = group))+
  # geom_point(aes(x = fire.occ, y = sens, color = group), alpha = 0.05)+
  #geom_jitter(alpha = 0.2)+
  geom_line(aes(group = group), position = position_dodge(0.025))+
  geom_point(position = position_dodge(0.025), size = 2.5)+
  geom_segment( aes(y = low.ci, yend = high.ci ),position = position_dodge(0.025))+
  scale_color_manual( values = c('tomato', "#0072B2"),  labels = c('Fire', 'Unburned'))+
  scale_shape( labels = c('Fire', 'Unburned'))+
  theme_bw(base_size = 13)+
  labs(color = 'Group', shape= 'Group')+
  theme(legend.position = 'top')+
  xlab('')+
  ylab('Sensitivity')


did.lgd <- get_legend(did.lgd.plt)



library(ggpubr)
library(patchwork)
grd.npp.match <- ggpubr::ggarrange(did.lgd,
                           ggarrange(did.sens.plt, did.pue.plt, ncol = 2, widths = c(3,3), labels = c('(a)','(b)'),
                                     font.label = list(face = "plain", size = 13, color = "black")),
                           nrow = 2, 
                           heights = c(0.4,3))

grd.npp.match
ggsave('didfig1.svg', grd.npp.match, width = 7, height = 4)

