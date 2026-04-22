#mediation

library(tidyverse)
library(lme4)
library(lmerTest)
library(mediation)
library(patchwork)

#1 estimate DiD for mediator 


matched.all.data <- read.csv('data/matched_pcoa_data.csv')%>%
  mutate(s.sens  = scale(sens))%>%
 # filter(s.sens < 3)%>%
  unite(group.fire.occ, group, fire.occ, remove = F)%>%
  mutate(fire.occ = fct_relevel(fire.occ, c( 'before', 'after')),
         group.fire.occ = fct_relevel(group.fire.occ, c('unburned_before', 'unburned_after', 'fire_before', 'fire_after')),
         fire.num = ifelse(fire.occ == 'before', 0,1),
         group.num = ifelse(group == 'unburned', 0, 1),
         firexgroup = fire.num * group.num,
         s.annual = as.numeric(scale(m.annual)))%>%
  dplyr::select(sens, pue, fire.occ, fire.num, group.num, firexgroup, m.annual, s.annual, pixel.id, map, V2, V1, group)%>%
  na.omit()


m.sens.out.v2 <- glmer(data = matched.all.data, sens ~ fire.num + group.num + firexgroup + V2 + (1|pixel.id), family = gaussian() )
summary(m.sens.out.v2)

m.sens.med.v2 <- glmer(data = matched.all.data, V2 ~ fire.num + group.num + firexgroup  + (1|pixel.id), family  = gaussian() )
summary(m.sens.med.v2)


med.out <- mediate( m.sens.med.v2, m.sens.out.v2,
                    treat = "firexgroup",
                    mediator = "V2",
                    boot = F, sims = 100)
summary(med.out)


m.sens.out.ann <- glmer(data = matched.all.data, sens ~ fire.num + group.num + firexgroup + s.annual +(1|pixel.id), family = gaussian() )
summary(m.sens.out.ann)

m.sens.med.ann <- glmer(data = matched.all.data, s.annual ~ fire.num + group.num + firexgroup  +  (1|pixel.id), family  = gaussian() )
summary(m.sens.med.ann)
#generate instrumental variable


med.out.annual <- mediate( m.sens.med.ann, m.sens.out.ann,
                           treat = "firexgroup",
                           mediator = "s.annual",
                           covariates = 'map',
                          # robustSE = TRUE,
                           sims = 100)
summary(med.out.annual)

outcome_model <- glmer(data = matched.all.data, sens ~ fire.num + group.num + firexgroup * V2 + firexgroup*map +(1|pixel.id), family = gaussian() )
summary(m.sens.out.ann)
p

mediation_model <- glmer(data = matched.all.data, s.annual ~ fire.num + group.num + firexgroup  +  firexgroup*map + (1|pixel.id), family  = gaussian() )
summary(m.sens.med.ann)
summary(med.out.annual)


m.sens.out.v2 <- glmer(data = matched.all.data, sens ~ fire.num + group.num + firexgroup * V2 + firexgroup*map + (1|pixel.id), family = gaussian() )
summary(m.sens.out)



m.sens.med.v2 <- glmer(data = matched.all.data, V2 ~ fire.num + group.num + firexgroup +  firexgroup*map +(1|pixel.id), family  = gaussian() )
summary(m.sens.med.v2)


med.out.v2 <- mediate( m.sens.med.v2, m.sens.out.v2,
                       treat = "firexgroup",
                       mediator = "V2",
                       covariates = 'map',
                       #robustSE = TRUE,
                       sims = 100)
#med.out.v2.map <- med.out.v2
summary(med.out.v2)



##try with CMAverse
library(CMAverse)

###

m.sens.out.v1v2 <- glm(data = matched.all.data, sens ~ fire.num + group.num + firexgroup + V2  + V1  ,family = gaussian() )
summary(m.sens.out.v1v2)

m.sens.med.v2 <- glm(data = matched.all.data, V2 ~ fire.num + group.num + firexgroup , family  = gaussian() )
summary(m.sens.med.v2 )

m.sens.med.v1 <- glm(data = matched.all.data, V1 ~ fire.num + group.num + firexgroup  , family  = gaussian() )
summary(m.sens.med.v1 )

cma.sens <- cmest(data = matched.all.data, 
                  model = 'gformula',
                  estimation = 'imputation',
                  inference =  'bootstrap',
                  outcome = 'sens',
                  exposure = 'firexgroup',
                  mediator = c('V1', 'V2'),
                  EMint = T,
                  mval = list(mean(matched.all.data$V1), mean(matched.all.data$V1)),
                  mreg = list(m.sens.med.v1, m.sens.med.v2),
                  yreg = m.sens.out.v1v2,
                  basec = c('fire.num', 'group.num'),
                  nboot = 1000, 
                  astar = 0,
                  a = 1,
                  full = T)

summary(cma.sens)
cma.sens.sens <- cmsens(cma.sens)

#################


####################
library(CMAverse)

outcome_model_annual_lm <- lm(data = matched.all.data, sens ~ fire.num + group.num + firexgroup * s.annual + firexgroup* map   )
mediation_model_annual_lm <- lm(data = matched.all.data, s.annual ~ fire.num + group.num + firexgroup * map )

cma.sens.ann <- cmest(data = matched.all.data, 
                       model = 'gformula',
                       estimation = 'imputation',
                       inference =  'bootstrap',
                       outcome = 'sens',
                       exposure = 'firexgroup',
                       mediator = c('s.annual'),
                       EMint = T,
                       mval = list(mean(matched.all.data$s.annual)),
                       mreg = list(mediation_model_annual_lm),
                       yreg = outcome_model_annual_lm,
                       basec = c('fire.num', 'group.num', 'map'),
                       nboot = 1000,
                       full = T)
summary(cma.sens.ann)


outcome_model_V2_lm <- lm(data = matched.all.data, sens ~ fire.num + group.num + firexgroup * V2  +firexgroup *map)
mediation_model_V2_lm <- lm(data = matched.all.data, V2 ~ fire.num + group.num + firexgroup *map )

#manual_cmest_decomp_single_boot(outcome_model_V2_lm, mediation_model_V2_lm,
     #                           data = matched.all.data, exposure = "firexgroup", mediator = "V2")

cma.sens.V2 <- cmest(data = matched.all.data, 
                           model = 'gformula',
                           estimation = 'imputation',
                           inference =  'bootstrap',
                           outcome = 'sens',
                           exposure = 'firexgroup',
                           mediator = c('V2'),
                          EMint = T,
                          mval = list(mean(matched.all.data$V2)),
                          mreg = list(mediation_model_V2_lm),
                          yreg = outcome_model_V2_lm,
                          basec = c('fire.num', 'group.num', 'map'),
                          nboot = 1000,
                          full = T)

summary( cma.sens.V2 )

###v1
outcome_model_V1_lm <- lm(data = matched.all.data, sens ~ fire.num + group.num + firexgroup * V1  +  firexgroup * map)
mediation_model_V1_lm <- lm(data = matched.all.data, V1 ~ fire.num + group.num + firexgroup * map  )

cma.sens.V1 <- cmest(data = matched.all.data, 
                          model = 'gformula',
                          estimation = 'imputation',
                          inference =  'bootstrap',
                          outcome = 'sens',
                          exposure = 'firexgroup',
                          mediator = c('V1'),
                          EMint = T,
                          mval = list(mean(matched.all.data$V1)),
                          mreg = list(mediation_model_V1_lm),
                          yreg = outcome_model_V1_lm,
                          basec = c('fire.num', 'group.num', 'map'),
                          nboot = 1000,
                          full = T)

summary( cma.sens.V1 )

##compare to cmest
#first run an lm 
outcome_model_lm <- lm(data = matched.all.data, sens ~ fire.num + group.num + firexgroup * V2  + firexgroup * V1 + firexgroup * map  )
#plot(outcome_model_lm)

summary(outcome_model_lm)
require(sjPlot)
outcome_model_lm_pue
tab_model(outcome_model_lm, outcome_model_lm_pue,
          pred.labels =  c('Intercept', 'Time', 'Group (fire)','TimexGroup', 'PCOA Axis 2', 'PCOA Axis 1', 'MAP', 'Axis 2 *Time* Group (fire)',  'Axis 1 *Time* Group (fire)', 'MAP *Time* Group (fire)'),
          dv.labels = c('Sensitivty', 'Efficiency'))


mediation_model_v2_lm <- lm(data = matched.all.data, V2 ~ fire.num + group.num + firexgroup + firexgroup *map)

#check fit
#plot(mediation_model_v2_lm)

mediation_model_v1_lm <- lm(data = matched.all.data, V1 ~ fire.num + group.num + firexgroup + firexgroup *map )



cma.sens.v1v2<- cmest(data = matched.all.data, 
                  model = 'gformula',
                  estimation = 'imputation',
                  inference =  'bootstrap',
                  outcome = 'sens',
                  exposure = 'firexgroup',
                  mediator = c('V1', 'V2'),
                  EMint = T,
                  mval = list(mean(matched.all.data$V1), mean(matched.all.data$V2)),
                  mreg = list(mediation_model_v1_lm, mediation_model_v2_lm),
                  yreg = outcome_model_lm,
                  basec = c('fire.num', 'group.num', 'map'),
                  nboot = 100,
                  full = T)

summary(cma.sens.v1v2)

#########################
sens.v1v2.df <- data.frame(pe = cma.sens.v1v2$effect.pe, ci.l = cma.sens.v1v2$effect.ci.low, ci.high = cma.sens.v1v2$effect.ci.high)%>%
  mutate(model = 'PCoA1+PCoA2')%>%
  rownames_to_column(var = "Estimate")



sensv2.df <- data.frame(pe = cma.sens.V2$effect.pe, ci.l = cma.sens.V2$effect.ci.low, ci.high = cma.sens.V2$effect.ci.high)%>%
  mutate(model = 'PCoA2')%>%
  rownames_to_column(var = "Estimate")

sensv1.df <- data.frame(pe = cma.sens.V1$effect.pe, ci.l = cma.sens.V1$effect.ci.low, ci.high = cma.sens.V1$effect.ci.high)%>%
  mutate(model = 'PCoA1')%>%
  rownames_to_column(var = "Estimate")


sensann.df <- data.frame(pe = cma.sens.ann$effect.pe, ci.l = cma.sens.ann$effect.ci.low, ci.high = cma.sens.ann$effect.ci.high)%>%
  mutate(model = '%Annual')%>%
  rownames_to_column(var = "Estimate")

sens.meds.df <- bind_rows(sens.v1v2.df, sensv1.df)%>%
  bind_rows(sensv2.df)%>%
  #bind_rows(sensann.df)%>%
  mutate(Estimate = toupper(fct_relevel(Estimate, 'te', 'tnie', 'pnde')))

write.csv(sens.meds.df , 'sens_med.csv')
###################PUE################

outcome_model_annual_lm_pue <- lm(data = matched.all.data, pue ~ fire.num + group.num + firexgroup * s.annual + firexgroup* map   )


mediation_model_annual_lm_pue <- lm(data = matched.all.data, s.annual ~ fire.num + group.num + firexgroup * map )



#manual_cmest_decomp_single_boot(outcome_model_annual_lm, mediation_model_annual_lm,
#                                data = matched.all.data, exposure = "firexgroup", mediator = "s.annual")

cma.pue.ann <- cmest(data = matched.all.data, 
                      model = 'gformula',
                      estimation = 'imputation',
                      inference =  'bootstrap',
                      outcome = 'pue',
                      exposure = 'firexgroup',
                      mediator = c('s.annual'),
                      EMint = T,
                      mval = list(mean(matched.all.data$s.annual)),
                      mreg = list(mediation_model_annual_lm_pue),
                      yreg = outcome_model_annual_lm_pue,
                      basec = c('fire.num', 'group.num', 'map'),
                      nboot = 1000,
                      full = T)
summary(cma.pue.ann)


outcome_model_V2_lm_pue <- lm(data = matched.all.data, pue ~ fire.num + group.num + firexgroup * V2  +firexgroup *map)
mediation_model_V2_lm_pue <- lm(data = matched.all.data, V2 ~ fire.num + group.num + firexgroup *map )

#

cma.pue.V2 <- cmest(data = matched.all.data, 
                     model = 'gformula',
                     estimation = 'imputation',
                     inference =  'bootstrap',
                     outcome = 'pue',
                     exposure = 'firexgroup',
                     mediator = c('V2'),
                     EMint = T,
                     mval = list(mean(matched.all.data$V2)),
                     mreg = list(mediation_model_V2_lm_pue),
                     yreg = outcome_model_V2_lm_pue,
                     basec = c('fire.num', 'group.num', 'map'),
                     nboot = 1000,
                     full = T)

summary( cma.pue.V2 )

###v1
outcome_model_V1_lm_pue <- lm(data = matched.all.data, pue ~ fire.num + group.num + firexgroup * V1  +  firexgroup * map)
mediation_model_V1_lm_pue <- lm(data = matched.all.data, V1 ~ fire.num + group.num + firexgroup * map  )


cma.pue.V1 <- cmest(data = matched.all.data, 
                     model = 'gformula',
                     estimation = 'imputation',
                     inference =  'bootstrap',
                     outcome = 'pue',
                     exposure = 'firexgroup',
                     mediator = c('V1'),
                     EMint = T,
                     mval = list(mean(matched.all.data$V1)),
                     mreg = list(mediation_model_V1_lm_pue),
                     yreg = outcome_model_V1_lm_pue,
                     basec = c('fire.num', 'group.num', 'map'),
                     nboot = 1000,
                     full = T)

summary( cma.pue.V1 )

##compare to cmest
#first run an lm 
outcome_model_lm_pue <- lm(data = matched.all.data, pue ~ fire.num + group.num + firexgroup * V2  + firexgroup * V1 + firexgroup *map  )

mediation_model_v2_lm_pue <- lm(data = matched.all.data, V2 ~ fire.num + group.num + firexgroup + firexgroup *map)

mediation_model_v1_lm_pue <- lm(data = matched.all.data, V1 ~ fire.num + group.num + firexgroup + firexgroup *map )



cma.pue.v1v2<- cmest(data = matched.all.data, 
                      model = 'gformula',
                      estimation = 'imputation',
                      inference =  'bootstrap',
                      outcome = 'pue',
                      exposure = 'firexgroup',
                      mediator = c('V1', 'V2'),
                      EMint = T,
                      mval = list(mean(matched.all.data$V1), mean(matched.all.data$V2)),
                      mreg = list(mediation_model_v1_lm_pue, mediation_model_v2_lm_pue),
                      yreg = outcome_model_lm_pue,
                      basec = c('fire.num', 'group.num', 'map'),
                      nboot = 100,
                      full = T)

summary(cma.pue.v1v2)

#########################
pue.v1v2.df <- data.frame(pe = cma.pue.v1v2$effect.pe, ci.l = cma.pue.v1v2$effect.ci.low, ci.high = cma.pue.v1v2$effect.ci.high)%>%
  mutate(model = 'PCoA1+PCoA2')%>%
  rownames_to_column(var = "Estimate")


puev2.df <- data.frame(pe = cma.pue.V2$effect.pe, ci.l = cma.pue.V2$effect.ci.low, ci.high = cma.pue.V2$effect.ci.high)%>%
  mutate(model = 'PCoA2')%>%
  rownames_to_column(var = "Estimate")

puev1.df <- data.frame(pe = cma.pue.V1$effect.pe, ci.l = cma.pue.V1$effect.ci.low, ci.high = cma.pue.V1$effect.ci.high)%>%
  mutate(model = 'PCoA1')%>%
  rownames_to_column(var = "Estimate")

cma.pue.ann
pueann.df <- data.frame(pe = cma.pue.ann$effect.pe, ci.l = cma.pue.ann$effect.ci.low, ci.high = cma.pue.ann$effect.ci.high)%>%
  mutate(model = '%Annual')%>%
  rownames_to_column(var = "Estimate")

pue.meds.df <- bind_rows(pue.v1v2.df, puev1.df)%>%
  bind_rows(puev2.df)%>%
  #bind_rows(pueann.df)%>%
  mutate(Estimate = toupper(Estimate))
write.csv(pue.meds.df,'pue_meds.csv')
paths.pue <- pue.meds.df%>%
  filter(Estimate %in% c('TE', 'PNDE', 'TNIE'))%>%
  mutate(Estimate = fct_relevel(Estimate, 'TE', 'TNIE', 'PNDE'))
pm.pue <- pue.meds.df%>%
  filter(Estimate %in% c('PM'))

decomp.pue <- ggplot(paths.pue, aes(x = Estimate, y = pe, color = fct_relevel(model, c('PCoA1', 'PCoA2', 'PCoA1+PCoA2'))))+
  geom_point(position = position_dodge(0.5), size = 3)+
  geom_segment(aes(y = ci.l, yend = ci.high),position = position_dodge(0.5), linewidth = 1)+
  labs(color = 'Model', y = expression("Effect estimate g m"^{-2}~"mm"^{-1}), x = 'Decomposition parameter')+
  guides(color = 'none')+
  scale_color_manual(values = c("#E6AB02","#4682B4", "#40E0D0", "#608090"), labels = c('PCoA1 (Annual - Shrub)', 'PCoA2 (Bare  - Perennial) (Annual - Shrub)', 'PCoA1+PCoA2'))+ ggtitle("PUE") +
  theme_bw(base_size = 13)+
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      vjust = -6
    )
  )
pm.pue.plt <- ggplot(pm.pue, aes(x = Estimate, y = pe, fill = fct_relevel(model, c('PCoA1', 'PCoA2', 'PCoA1+PCoA2'))))+
  geom_bar(stat = 'identity', position = position_dodge(1))+
  geom_segment(aes(y = ci.l, yend = ci.high),position = position_dodge(1), linewidth = 1)+
  labs(fill = 'Mediators',  y = '% Mediated ',, x = 'Proportion mediated')+
  scale_fill_manual(values = c("#E6AB02","#4682B4", "#40E0D0", "#608090"), labels = c('PCoA1 (Annual - Shrub)', 'PCoA2 (Bare  - Perennial)', 'PCoA1+PCoA2'))+
  theme_bw(base_size = 13)



paths.sens <- sens.meds.df%>%
  filter(Estimate %in% c('TE', 'PNDE', 'TNIE'))%>%
  mutate(Estimate = fct_relevel(Estimate, 'TE', 'TNIE', 'PNDE'))
pm.sens <- sens.meds.df%>%
  filter(Estimate %in% c('PM'))

decomp.sens <- ggplot(paths.sens, aes(x = Estimate, y = pe, color = fct_relevel(model, c('PCoA1', 'PCoA2', 'PCoA1+PCoA2'))))+
  geom_point(position = position_dodge(0.5), size = 3)+
  geom_segment(aes(y = ci.l, yend = ci.high),position = position_dodge(0.5), linewidth = 1)+
  labs(color = 'Model', y =  expression("Effect estimate g m"^{-2}~"mm"^{-1}), x = 'Decomposition parameter')+
  guides(color = 'none')+
  scale_color_manual(values = c("#E6AB02","#4682B4", "#40E0D0", "#608090"), labels = c('PCoA1 (Annual - Shrub)', 'PCoA2 (Bare  - Perennial)', 'PCoA1+PCoA2'))+ 
  ggtitle("Sensitivity") +
  theme_bw(base_size = 13)+
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      vjust = -6
    )
  )


pm.sens.plt <- ggplot(pm.sens, aes(x = Estimate, y = pe, fill = fct_relevel(model, c('PCoA1', 'PCoA2', 'PCoA1+PCoA2'))))+
  geom_bar(stat = 'identity', position = position_dodge(1))+
  geom_segment(aes(y = ci.l, yend = ci.high),position = position_dodge(1), linewidth = 1)+
  labs(fill = 'Mediators',  y = '% Mediated', x = 'Proportion mediated')+
  scale_fill_manual(values = c("#E6AB02","#4682B4", "#40E0D0", "#608090"), labels = c('PCoA1 (Annual - Shrub)', 'PCoA2 (Bare  - Perennial)', 'PCoA1+PCoA2'))+
  theme_bw(base_size = 13)


layout = '
aaab
cccd'
library(patchwork)
decomp.sens + pm.sens.plt + decomp.pue + pm.pue.plt + plot_layout( design = layout, guides = 'collect', axes = 'collect') + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

##simplified mediation figure

decomp.sens <- ggplot(paths.sens, aes(x = Estimate, y = pe, color = fct_relevel(model, c('PCoA1+PCoA2','PCoA1', 'PCoA2', '%Annual'))))+
  geom_point(position = position_dodge(0.5), size = 3)+
  geom_segment(aes(y = ci.l, yend = ci.high),position = position_dodge(0.5), linewidth = 1)+
  labs(color = 'Model', y =  expression("Effect estimate g m"^{-2}~"mm"^{-1}), x = 'Decomposition parameter')+
  guides(color = 'none')+
  scale_color_manual(values = c("#E6AB02","#4682B4", "#40E0D0", "#608090"))+ 
  ggtitle("Sensitivity") +
  theme_bw(base_size = 13)+
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      vjust = -6
    )
  )


pm.sens.plt <- ggplot(pm.sens, aes(x = Estimate, y = pe, fill = fct_relevel(model, c('PCoA1+PCoA2','PCoA1', 'PCoA2'))))+
  geom_bar(stat = 'identity', position = position_dodge(1))+
  geom_segment(aes(y = ci.l, yend = ci.high),position = position_dodge(1), linewidth = 1)+
  labs(fill = 'Mediators',  y = '% Mediated', x = 'Proportion mediated')+
  scale_fill_manual(values = c("#E6AB02","#4682B4", "#40E0D0", "#608090"))+
  theme_bw(base_size = 13)


layout = '
aaab
cccd'
library(patchwork)
decomp.sens + pm.sens.plt + decomp.pue + pm.pue.plt + plot_layout( design = layout, guides = 'collect', axes = 'collect') + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

##################################
###########sensitivity analysis##########
summary(cma.sens.v1v2)
cmsens(cma.sens.v1v2)$evalues%>%
  as_tibble()
cmsens(cma.pue.v1v2)
sens.sens <- cmsens(cma.sens.v1v2)$evalues%>%
  as.data.frame()%>%
  rownames_to_column(var = 'Param')

pue.sens <- cmsens(cma.pue.v1v2)$evalues%>%
  as.data.frame()%>%
  rownames_to_column(var = 'Param')

senses <- bind_rows( sens.sens, pue.sens)

write.csv(senses, 'senses.csv')
0.027615 - 0.045651


##########################
#predict model 
outcome_model_lm_pue <- lm(data = matched.all.data, pue ~ fire.num + group.num + firexgroup * V2  + firexgroup * V1 + firexgroup *map  )
pue.outcome.v2.plt <- effects::predictorEffect('V2', outcome_model_lm_pue, xlevels = list(firexgroup = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = V2, y = fit, color = as.factor(firexgroup)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(firexgroup) ), alpha = 0.3)+ 
  labs(y = expression("PUE g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'PCoA2')+
  scale_color_manual(values = c("#0072B2", "#D55E00"))+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  theme_bw()


pue.outcome.v1.plt <- effects::predictorEffect('V1', outcome_model_lm_pue, xlevels = list(firexgroup = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = V1, y = fit, color = as.factor(firexgroup)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(firexgroup) ), alpha = 0.3)+ 
  labs(y = expression("PUE g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'PCoA1')+
  scale_color_manual(values = c("#0072B2", "#D55E00"))+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  theme_bw()

sens.outcome.v1.plt <- effects::predictorEffect('V1', outcome_model_lm, xlevels = list(firexgroup = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = V1, y = fit, color = as.factor(firexgroup)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(firexgroup)), alpha = 0.3)+
  labs(y = expression("Sens g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'PCoA1')+
  scale_color_manual(values = c("#0072B2", "#D55E00"))+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  theme_bw()

sens.outcome.v2.plt <- effects::predictorEffect('V2', outcome_model_lm, xlevels = list(firexgroup = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = V2, y = fit, color = as.factor(firexgroup)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(firexgroup)), alpha = 0.3)+
  labs(y = expression("Sens g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'PCoA2')+
  scale_color_manual(values = c("#0072B2", "#D55E00"))+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  theme_bw()

sens.outcome.map <- effects::predictorEffect('map', outcome_model_lm, xlevels = list(firexgroup = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = map, y = fit, color = as.factor(firexgroup)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(firexgroup)), alpha = 0.3)+
  labs(y = expression("Sens g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'MAP')+
  scale_color_manual(values = c("#0072B2", "#D55E00"))+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  theme_bw()


pue.outcome.map <- effects::predictorEffect('map', outcome_model_lm_pue, xlevels = list(firexgroup = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = map, y = fit, color = as.factor(firexgroup)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(firexgroup)), alpha = 0.3)+
  scale_color_manual(values = c("#0072B2", "#D55E00"))+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  labs(y = expression("PUE g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'MAP')+
  theme_bw()
sens.outcome.v1.plt + sens.outcome.v2.plt + sens.outcome.map +pue.outcome.v1.plt + pue.outcome.v2.plt + pue.outcome.map + plot_layout(guides = 'collect', axis_titles = 'collect') + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')
###############
summary(cma.sens.comp.rb)






##
m.sens.out.v2 <- glm(data = matched.all.data, sens ~ fire.num + group.num + firexgroup + V2 , family = gaussian() )
summary(m.sens.out.v2)



m.sens.med.v2 <- glm(data = matched.all.data, V2 ~ fire.num + group.num + firexgroup   , family  = gaussian() )
summary(m.sens.med.v2)

med.out.v2 <- mediate( m.sens.med.v2, m.sens.out.v2,
                       treat = "firexgroup",
                       mediator = "V2",
                       boot = F, sims = 100)
summary(med.out.v2)
cma.sens.v2 <- cmest(data = matched.all.data, 
                     model = 'gformula',
                     estimation = 'imputation',
                     inference =  'bootstrap',
                     outcome = 'sens',
                     exposure = 'firexgroup',
                     mediator = c('V2'),
                     EMint = T,
                     mval = list(mean(matched.all.data$V2)),
                     mreg = list(m.sens.med.v2),
                     yreg = m.sens.out.v2,
                     basec = c('fire.num', 'group.num'),
                     nboot = 1000, 
                     astar = 0,
                     a = 1,
                     full = T)
summary(cma.sens.v2)


m.sens.out.v1 <- glm(data = matched.all.data, sens ~ fire.num + group.num + firexgroup + V1 , family = gaussian() )
summary(m.sens.out)



m.sens.med.v1 <- glm(data = matched.all.data, V1 ~ fire.num + group.num + firexgroup  , family  = gaussian() )
summary(m.sens.med.ann)
med.out.v1 <- mediate( m.sens.med.v1, m.sens.out.v1,
                       treat = "firexgroup",
                       mediator = "V1",
                       boot = F, sims = 100)
summary(med.out.v1)

cma.sens.v1 <- cmest(data = matched.all.data, 
                     model = 'gformula',
                     estimation = 'imputation',
                     inference =  'bootstrap',
                     outcome = 'sens',
                     exposure = 'firexgroup',
                     mediator = c('V1'),
                     EMint = F,
                     mval = list(mean(matched.all.data$V1)),
                     mreg = list(m.sens.med.v1),
                     yreg = m.sens.out.v1,
                     basec = c('fire.num', 'group.num'),
                     nboot = 1000, 
                     astar = 0,
                     a = 1,
                     full = T)
summary(cma.sens.v1)

##annual mediation model 
m.sens.out.ann <- glm(data = matched.all.data, sens ~ fire.num + group.num + firexgroup + m.annual + firexgroup*map, family = gaussian() )
summary(m.sens.out)

m.sens.med.ann <- glm(data = matched.all.data, m.annual ~ fire.num + group.num + firexgroup*map  , family  = gaussian() )
summary(m.sens.med.ann)
med.out.ann <- mediate( m.sens.med.ann, m.sens.out.ann,
                        treat = "firexgroup",
                        mediator = "m.annual",
                        boot = F, sims = 100)
summary(med.out.ann)
cma.sens.ann <- cmest(data = matched.all.data, 
                      model = 'gformula',
                      estimation = 'imputation',
                      inference =  'bootstrap',
                      outcome = 'sens',
                      exposure = 'firexgroup',
                      mediator = c('m.annual'),
                      EMint = T,
                      mval = list(mean(matched.all.data$m.annual)),
                      mreg = list(m.sens.med.ann),
                      yreg = m.sens.out.ann,
                      basec = c('fire.num', 'group.num', 'map'),
                      nboot = 1000, 
                      astar = 0,
                      a = 1,
                      full = T)
#cmest()
summary(cma.sens.ann)


####### Elevation and MAP moderator analysis ##########

matched.unrestored <- read.csv('matched_unrestored.csv')
unique(matched.unrestored$group)
matched.all.data.pred <- read.csv('matched_pcoa_data.csv')%>%
  mutate(s.sens  = scale(sens))%>%
  # filter(s.sens < 3)%>%
  unite(group.fire.occ, group, fire.occ, remove = F)%>%
 # merge(subclass.id, by = c('group', 'pixel.id', 'pix'))%>%
  filter(pixel.id %in% matched.unrestored$pixel.id)%>%
  mutate(fire.occ = fct_relevel(fire.occ, c( 'before', 'after')),
         group.fire.occ = fct_relevel(group.fire.occ, c('unburned_before', 'unburned_after', 'fire_before', 'fire_after')),
         fire.num = ifelse(fire.occ == 'before', 0,1),
         group.num = ifelse(group == 'unburned', 0, 1),
         firexgroup = fire.num * group.num,
         s.annual = as.numeric(scale(m.annual)))%>%
  dplyr::select(sens, pue, fire.occ, fire.num, group.num, firexgroup, m.annual, s.annual, pixel.id, map, V2, V1, elev)%>%
  na.omit()%>%
  group_by(pixel.id)%>%
  mutate(initial.v1 = V1[fire.occ == 'before'],
         initial.v2 = V2[fire.occ == 'before'])%>%
  rename(fire.eff = firexgroup)

matched.all.data.pred%>%
  filter(is.na(V1))
  
V1.10 <- quantile(matched.all.data.pred$initial.v1, .05)
V1.90 <- quantile(matched.all.data.pred$initial.v1, .95)
V2.10 <- quantile(matched.all.data.pred$initial.v2, .05)
V2.90 <- quantile(matched.all.data.pred$initial.v2, .95)
elev.10 <- quantile(matched.all.data.pred$elev, .05)
elev.90 <- quantile(matched.all.data.pred$elev, .95)

matched.unrest.pred.reduced <- matched.all.data.pred%>%
  filter(V1 > V1.10 & V1 < V1.90)%>%
  filter(V2 > V2.10 & V2 < V2.90)%>%
  filter(elev > elev.10 & elev < elev.90)

#matched.all.data.pred<- matched.unrest.pred.reduced



#######now compare elev and map using initial comp##########

library(sjPlot)



sens_lm_elev <- lmer(data = matched.unrest.pred.reduced, sens ~ fire.num + group.num + fire.eff * initial.v1  + fire.eff * initial.v2 +  fire.eff*elev + (1|pixel.id) )
summary(sens_lm_elev)

sens_elev_nocomp <- lmer(data = matched.unrest.pred.reduced, sens ~ fire.num + group.num + elev*fire.eff + (1|pixel.id) )
summary(sens_elev_nocomp)
AIC(sens_elev_nocomp)


sens_elev_trend_nocomp <- emmeans::emtrends(sens_elev_nocomp, var = 'elev', pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'no comp')

nocomp.elev.eff <- effects::predictorEffect('elev',sens_elev_nocomp , xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  mutate(model = 'No comp.')


sens.elev.pvalue <- data.frame(model = c('Initial comp.', 'No comp.'),
                              p.val = c(summary(sens_lm_elev)$coefficients[10,5], summary(sens_elev_nocomp)$coefficients[6,5]))%>%
  mutate(sig = ifelse(p.val < 0.001, '***', 
                      ifelse(p.val < 0.01, '**', 
                             ifelse(p.val < 0.05, "*", ''))))


Sens.elev <- effects::predictorEffect('elev', sens_lm_elev, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  mutate(model = 'Initial comp.')%>%
  bind_rows(nocomp.elev.eff)%>%
  ggplot( aes(x = elev, y = fit, color = as.factor(fire.eff)))+
  geom_ribbon( aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff)), alpha = 0.2)+
  geom_line(linewidth = 1)+
  geom_text(data = sens.elev.pvalue, aes(x = 1650, y = 0.3, label = sig), inherit.aes = F, size = 7.5)+
   # geom_point(data =  matched.unrest.pred.reduced , aes(y = sens), shape =21)+
  #ylim(c(-0.1, 0.34))+
  labs(y = '', fill = 'Fire effect', color = 'Fire effect', x = 'Elevation')+
  coord_cartesian(ylim = c(0.0, 0.35))+
  scale_color_manual(values = c("#0072B2", 'tomato'), guide = 'none')+
  scale_fill_manual(values = c("#0072B2", 'tomato'), guide = 'none')+
  theme_bw()+
  facet_wrap(~model)


Sens.elev 

##sens map

sens_lm_map <-  lmer(data = matched.unrest.pred.reduced, sens ~ fire.num + group.num + fire.eff * initial.v1  + fire.eff * initial.v2 +  fire.eff*map + (1|pixel.id) )
summary(sens_lm_map)

sens_map_nocomp <-  lmer(data = matched.unrest.pred.reduced, sens ~ fire.num + group.num + fire.eff*map + (1|pixel.id) )
summary(sens_map_nocomp)

nocomp.map.eff <- effects::predictorEffect('map', sens_map_nocomp, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  mutate(model = 'No comp.')

sens.map.pvalue <- data.frame(model = c('Initial comp.', 'No comp.'),
                               p.val = c(summary(sens_lm_map)$coefficients[10,5], summary(sens_map_nocomp)$coefficients[6,5]))%>%
  mutate(sig = ifelse(p.val < 0.001, '***', 
                      ifelse(p.val < 0.01, '**', 
                             ifelse(p.val < 0.05, "*", ''))))

Sens.map <- effects::predictorEffect('map', sens_lm_map, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  mutate(model = 'Initial comp.')%>%
  bind_rows(nocomp.map.eff)%>%
  ggplot( aes(x = map, y = fit, color = as.factor(fire.eff)))+
  geom_ribbon( aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff)), alpha = 0.2)+
  geom_line( linewidth = 1)+
  labs(y = expression("Sens g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'Mean annual precip (mm)')+
  geom_text(data = sens.map.pvalue, aes(x = 400, y = 0.3, label = sig), inherit.aes = F, size = 7.5)+
  ylim(c(-0.21, 0.40))+
  scale_color_manual(values = c("#0072B2", 'tomato'))+
  scale_fill_manual(values = c("#0072B2", 'tomato'))+
  theme_bw()+
  theme(
    legend.position = c(0.25, 0.4),  # Correct argument
    legend.justification = c("right", "top"),  # Anchors legend to that corner
    legend.background = element_rect(fill = "white", color = "black")  # Optional styling
  )+
  facet_wrap(~model)


Sens.map 

#pue map
pue_lm_map <- lmer(data = matched.unrest.pred.reduced, pue ~ fire.num + group.num + fire.eff * initial.v1  + fire.eff * initial.v2 +  fire.eff*map + (1|pixel.id) )
plot(pue_lm_map)
summary(pue_lm_map)

pue_map_nocomp <- lmer(data = matched.unrest.pred.reduced, pue ~ fire.num + group.num + fire.eff*map + (1|pixel.id) )

summary(pue_map_nocomp  )

pue.nocomp.map.eff <- effects::predictorEffect('map', pue_map_nocomp, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  mutate(model = 'No comp.')


pue.map.pvalue <- data.frame(model = c('Initial comp.', 'No comp.'),
                             p.val = c(summary(pue_lm_map)$coefficients[10,5], summary(pue_map_nocomp)$coefficients[6,5]))%>%
  mutate(sig = ifelse(p.val < 0.001, '***', 
                  ifelse(p.val < 0.01, '**', 
                         ifelse(p.val < 0.05, "*", ''))))


PUE.map <- effects::predictorEffect('map', pue_lm_map, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  mutate(model = 'Initial comp.')%>%
  bind_rows(pue.nocomp.map.eff)%>%
  ggplot( aes(x = map, y = fit, color = as.factor(fire.eff)))+
    geom_line( linewidth = 1)+#inital comp 
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff)), alpha = 0.2)+
  geom_text(data = pue.map.pvalue, aes(x = 400, y = 0.35, label = sig), inherit.aes = F, size = 7.5)+
  scale_color_manual(values =  c( "#0072B2", 'tomato'), guide = 'none')+
  scale_fill_manual(values = c("#0072B2", 'tomato'), guide = 'none')+
  labs(y = expression("PUE g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'Mean annual precip (mm)')+
  #geom_line( data = data.frame(map = NA, fit = NA, comp.group = factor(c("No comp.", "Initial comp."))), inherit.aes = F, aes(x = map, y = fit, linetype = comp.group),
  # show.legend = TRUE) +
  #scale_linetype_manual( name = "Model",values = c("Initial comp." = "solid",   "No comp." = "dashed"))+
  ylim(c(0.05,0.38))+
  theme_bw()+
  #theme(
  #  legend.position = c(0.41, 0.47),  # Correct argument
  #  legend.justification = c("right", "top"),  # Anchors legend to that corner
  #  legend.background = element_rect(fill = "white", color = "black")  # Optional styling
  #)+
  facet_wrap(~model)
PUE.map

#pue elev

pue_lm_elev <- lmer(data = matched.unrest.pred.reduced, pue ~ fire.num + group.num + fire.eff * initial.v1  + fire.eff * initial.v2 +  fire.eff*elev + (1|pixel.id) )
plot(pue_lm_map)
AIC(pue_lm_map)
summary(pue_lm_elev)

pue_elev_nocomp <- lmer(data = matched.unrest.pred.reduced, pue ~ fire.num + group.num + fire.eff*elev + (1|pixel.id) )

summary(pue_elev_nocomp)

pue.nocomp.elev.eff <- effects::predictorEffect('elev', pue_elev_nocomp, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  mutate(model = 'No comp.')

pue.elev.pvalue <- data.frame(model = c('Initial comp.', 'No comp.'),
                             p.val = c(summary(pue_lm_elev)$coefficients[10,5], summary(pue_elev_nocomp)$coefficients[6,5]))%>%
  mutate(sig = ifelse(p.val < 0.001, '***', 
                      ifelse(p.val < 0.01, '**', 
                             ifelse(p.val < 0.05, "*", ''))))


PUE.elev <- effects::predictorEffect('elev', pue_lm_elev, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  mutate(model = 'Initial comp.')%>%
  bind_rows(pue.nocomp.elev.eff)%>%
  ggplot( aes(x = elev, y = fit, color = as.factor(fire.eff)))+

  geom_line( linewidth = 1)+#no comp effect
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff)), alpha = 0.2)+#no comp effect
  
 # geom_ribbon(data = pue.nocomp.elev.eff, aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff)), alpha = 0.1, linewidth  = 0.15)+
  #geom_line(data = pue.nocomp.elev.eff,  linetype = 'dashed', linewidth = 0.75)+
  geom_text(data = pue.elev.pvalue, aes(x = 1650, y = 0.35, label = sig), inherit.aes = F, size = 7.5)+
  scale_color_manual(values = c("#0072B2", 'tomato'), guide = 'none')+
  scale_fill_manual(values = c("#0072B2", 'tomato'), guide = 'none')+
  #scale_linetype_manual(values = c('dashed', 'Soild'), name = 'Model', guide = guide_legend(override.aes=aes(linetype=NA)))+
  ylim(c(0.19,0.41))+
  labs(y = '', fill = 'Fire effect', color = 'Fire effect', x = 'Elevation')+
  theme_bw()+
  facet_wrap(~model)


map_elev_trends <- Sens.map + Sens.elev + PUE.map + PUE.elev + plot_layout( axes = 'collect') + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

map_elev_trends
ggsave( 'map_elev_fig_rot.svg', map_elev_trends, width = 9.68, height = 6, units = 'in')

###########

sens_elev_trend <- emmeans::emtrends(sens_lm_elev, var = 'elev',pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'comp')

emmeans::emtrends(sens_lm_elev, var = 'elev',~fire.eff)%>%
  emmeans::contrast(method = 'trt.vs.ctrl', adjust = 'none')%>%
  as.data.frame()%>%
  mutate(covariates = 'comp')

sens_elev_trend_nocomp <- emmeans::emtrends(sens_elev_nocomp, var = 'elev', pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'no comp')

sens_map_trend <- emmeans::emtrends(sens_lm_map, var = 'map', pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'comp')


sens_map_trend_nocomp <- emmeans::emtrends(sens_map_nocomp, var = 'map', pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'no comp')

pue_elev_trend <- emmeans::emtrends(pue_lm_elev, var = 'elev', pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'comp')

pue_elev_cont <- emmeans::emtrends(pue_lm_elev, var = 'elev', pairwise~fire.eff)%>%
  emmeans::contrast(method = 'trt.vs.ctrl')%>%
  as.data.frame()%>%
  mutate(covariates = 'comp')

pue_elev_trend_nocomp <- emmeans::emtrends(pue_elev_nocomp, var = 'elev', pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'no comp')

pue_elev_cont_nocomp <- emmeans::emtrends(pue_elev_nocomp, var = 'elev', pairwise~fire.eff)%>%
  emmeans::contrast(method = 'trt.vs.ctrl')%>%
  as.data.frame()%>%
  mutate(covariates = 'no comp')


emmeans::emtrends(pue_elev_nocomp, var = 'elev', pairwise~fire.eff)%>%
  emmeans::contrast(method = 'trt.vs.ctrl')
pue_map_trend <- emmeans::emtrends(pue_lm_map, var = 'map', pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'comp')

pue_map_trend_nocomp <- emmeans::emtrends(pue_map_nocomp, var = 'map', pairwise~fire.eff)%>%
  as.data.frame()%>%
  mutate(covariates = 'no comp')

emmeans::emtrends(pue_map_nocomp, var = 'map', ~fire.eff)%>%
  emmeans::contrast(method = 'trt.vs.ctrl')

map_elev_trends <- sens_elev_trend %>%
  bind_rows(sens_elev_trend_nocomp)%>%
  bind_rows(sens_map_trend)%>%
  bind_rows(sens_map_trend_nocomp)%>%
  bind_rows(pue_elev_trend) %>%
  bind_rows(pue_elev_trend_nocomp)%>%
  bind_rows(pue_map_trend)%>%
  bind_rows(pue_map_trend_nocomp)%>%
  as_tibble()

write.csv(map_elev_trends, 'map_elev_trends.csv')

######plot relationships with initial cover########

pue_initial_map <- lmer(data = matched.unrest.pred.reduced, pue ~ fire.num + group.num + fire.eff * initial.v1  + fire.eff * initial.v2 +  fire.eff*map + (1|pixel.id) )
sens_initial_map <- lmer(data = matched.unrest.pred.reduced, sens ~ fire.num + group.num + fire.eff * initial.v1  + fire.eff * initial.v2 +  fire.eff*map + (1|pixel.id) )

#PUE initial v2
pue.initial.v2.plt <- effects::predictorEffect('initial.v2', pue_initial_map, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = initial.v2, y = fit, color = as.factor(fire.eff)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff) ), alpha = 0.3)+ 
  labs(y = expression("PUE g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'Initial PCoA2')+
  ylim(c(0.13,0.48))+
  scale_color_manual(values = c('steelblue', 'tomato'))+
  scale_fill_manual(values = c('steelblue', 'tomato'))+
  theme_bw()+
  theme(axis.title.y = element_blank())

#PUE initial v1
pue.initial.v1.plt <- effects::predictorEffect('initial.v1', pue_initial_map, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = initial.v1, y = fit, color = as.factor(fire.eff)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff) ), alpha = 0.3)+ 
  labs(y = expression("PUE g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'Initial PCoA1')+
  ylim(c(0.13,0.48))+
  scale_color_manual(values = c('steelblue', 'tomato'))+
  scale_fill_manual(values = c('steelblue', 'tomato'))+
  theme_bw()

#Sens initial v2
sens.initial.v2.plt <- effects::predictorEffect('initial.v2', sens_initial_map, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = initial.v2, y = fit, color = as.factor(fire.eff)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff) ), alpha = 0.3)+ 
  labs(y = expression("Sens g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'Initial PCoA2')+
  ylim(c(0.07, 0.37))+
  scale_color_manual(values = c('steelblue', 'tomato'))+
  scale_fill_manual(values = c('steelblue', 'tomato'))+
  theme_bw()+
  theme(axis.title.y = element_blank())

#Sens initial v1
sens.initial.v1.plt <- effects::predictorEffect('initial.v1', sens_initial_map, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = initial.v1, y = fit, color = as.factor(fire.eff)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff) ), alpha = 0.3)+ 
  labs(y = expression("Sens g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'Initial PCoA1')+
  ylim(c(0.07, 0.37))+
  scale_color_manual(values = c('steelblue', 'tomato'))+
  scale_fill_manual(values = c('steelblue', 'tomato'))+
  theme_bw()

###intial cover plot
sens.initial.v1.plt + sens.initial.v2.plt +pue.initial.v1.plt + pue.initial.v2.plt +  plot_layout(guides = 'collect', axes = 'collect')+ plot_annotation(tag_levels = 'a')&
  theme(plot.margin = margin(1, 1, 1, 1),
        plot.tag.position = c(0.05, 0.95),
        plot.tag = element_text(size = 14, face = "bold")
  )


########
?plot_layout
sens.outcome.map <- effects::predictorEffect('map', outcome_model_lm, xlevels = list(firexgroup = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = map, y = fit, color = as.factor(firexgroup)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(firexgroup)), alpha = 0.3)+
  labs(y = expression("Sens g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'MAP')+
  scale_color_manual(values = c('steelblue', 'tomato'))+
  scale_fill_manual(values = c('steelblue', 'tomato'))+
  theme_bw()


pue.outcome.map <- effects::predictorEffect('map', outcome_model_lm_pue, xlevels = list(firexgroup = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = map, y = fit, color = as.factor(firexgroup)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(firexgroup)), alpha = 0.3)+
  scale_color_manual(values = c('steelblue', 'tomato'))+
  scale_fill_manual(values = c('steelblue', 'tomato'))+
  labs(y = expression("PUE g m"^{-2}~"mm"^{-1}), fill = 'Fire effect', color = 'Fire effect', x = 'MAP')+
  theme_bw()


#extra
outcome_model_lmer_elev <- lmer(data = matched.unrest.pred.reduced, sens ~ fire.num + group.num + fire.eff * initial.v1  + fire.eff * initial.v2 +  fire.eff*elev + (1|pixel.id) )
summary(outcome_model_lmer_elev)
effects::predictorEffect('elev', outcome_model_lmer_elev, xlevels = list(fire.eff= 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = elev, y = fit, color = as.factor(fire.eff)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff) ), alpha = 0.3)


outcome_model_lmer_elev_pue <- lmer(data = matched.all.data.pred, pue ~ fire.num + group.num + fire.eff * initial.v1  + fire.eff * initial.v2 +  fire.eff*elev + (1|pixel.id) )
summary(outcome_model_lmer_elev_pue)

effects::predictorEffect('elev', outcome_model_lmer_elev_pue, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = elev, y = fit, color = as.factor(fire.eff)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff) ), alpha = 0.3)

effects::predictorEffect('initial.v1', outcome_model_lmer_elev_pue, xlevels = list(fire.eff = 2))%>%
  as.data.frame()%>%
  ggplot( aes(x = initial.v1, y = fit, color = as.factor(fire.eff)))+
  geom_line()+
  geom_ribbon(aes(ymin =lower, ymax = upper, fill = as.factor(fire.eff) ), alpha = 0.3)

