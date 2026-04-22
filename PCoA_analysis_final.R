##this script requires the processed matched dataset
#This co

library(tidyverse)
library(lme4)
library(lmerTest)
library(gridExtra)
library(ggpubr)
library(vegan)

set.seed(0)

#model test 
#sensitivity 


matched.all.data <- read.csv('matched_did_data.csv')
subset <- matched.all.data%>%
  ungroup()%>%
  dplyr::select(subclass)%>%
  distinct()%>%
  sample_n(912)

subset.data <- matched.all.data%>%
  filter(subclass %in% subset$subclass)

mycolors =  c("#56B4E9", "#003366", "tomato","#7F0000" )

##figure 3 c-e PCoA 
set.seed(123)
match.mds.df <- subset.data%>%
  ungroup()%>%
  #mutate(m.bare = m.bare + m.litter)%>%
  dplyr::select('m.annual', 'm.perennial', 'm.bare', 'm.shrub')

### PERMANOVA to test burn vs unburned pre fire cover

matched.all.data$fire.occ <- as.factor(matched.all.data$fire.occ)
matched.all.data$group <- as.factor(matched.all.data$group)

#permanova test takes a long time
match.permanova.2 <- adonis2(match.mds.df ~ fire.occ + group + fire.occ:group, data = subset.data, method = 'euclid', by = 'terms', permutations = 999)
match.permanova.2

match.permanova.2%>%
  as.data.frame()%>%
  write.csv('cover_change_permanova.csv')
#match.permanova

#run pcoa
dist.m <- vegdist(match.mds.df, method = 'euclid')
match.pcoa<- cmdscale(dist.m, eig = T)


match.pcoa.df <- as.data.frame(match.pcoa$points)%>%
  bind_cols(matched.all.data)%>%
  unite( 'groups.all', fire.occ:group, remove = F)%>%
  mutate(fire.occ )



annual.fit <- envfit(match.pcoa~ match.mds.df$m.annual)

pcoa.points <- match.pcoa$points

rotate_ordination <- function(scores, fit, axes = c(1, 2)) {
  # scores: matrix or data frame of ordination scores (rows = samples)
  # fit:    envfit result (vegan::envfit)
  # axes:   which two axes to rotate (default: 1 and 2)
  
  # Extract the vector to align with axis 1
  v <- fit$vectors$arrows[axes]
  
  if (any(is.na(v))) {
    stop("The envfit object does not contain a vector for the specified axes.")
  }
  
  # Angle of the vector relative to x-axis
  theta <- atan2(v[2], v[1])
  
  # To align the vector with the x-axis, rotate by -theta
  R <- matrix(c(cos(-theta), -sin(-theta),
                sin(-theta),  cos(-theta)), nrow = 2)
  
  # Apply rotation to the selected axes
  rotated <- scores
  rotated[, axes] <- as.matrix(scores[, axes]) %*% R
  
  list(
    scores = rotated,
    rotation_matrix = R,
    angle = -theta,
    original_vector = v
  )
}

rot.points <- rotate_ordination(pcoa.points, annual.fit)


match.pcoa.df <- rot.points$scores%>%
  as.data.frame()%>%
  bind_cols(matched.all.data)%>%
  unite( 'groups.all', fire.occ:group, remove = F)%>%
  mutate(fire.occ )


write.csv(match.pcoa.df,'matched_pcoa_data.csv')




#extract cover loading vecotrs
coverscores <- as.data.frame(wascores(rot.points$scores, match.mds.df, expand = TRUE))
rownames(coverscores) <- c('Annual', 'Perennial', 'Bare', 'Shrub')


#plot pcoa in elipse  geom_point(aes(color = groups.all), alpha  = 0.3)+
scale_color_manual(values = mycolors, name = 'Group', breaks = c(  'before_unburned','after_unburned', 'before_fire','after_fire'), labels = c('Unburned: Before','Unburned: After',  'Burned: Before', 'Burned: After'))+
  stat_ellipse(aes(color = groups.all), level = 0.9, size = 1)+
  xlab('Axis 1')+
  ylab('Axis 2')+
  theme_bw(base_size = 12)+
  theme(legend.position = 'none')

match.ellipse.pcoa.unburned <- match.pcoa.df%>%
  filter(group == 'unburned')%>%
  ggplot(aes(V1,V2))+
  geom_point(aes(color = fire.occ), shape = 17,  alpha  = 0.3)+
  stat_ellipse(aes(color = fire.occ, linetype = fire.occ),  level = 0.95, linewidth = 2)+
  scale_color_manual(values = mycolors[c(1,2)], name = 'Group', breaks = c(  'before','after'), labels = c('Unburned: Before','Unburned: After'))+
  #scale_fill_manual(values = mycolors[c(1,2)], name = 'Group', breaks = c(  'before','after'), labels = c('Unburned: Before','Unburned: After'))+
  #scale_shape_manual(values = c(21,22))+
  #scale_linetype_manual(values = c('dashed', 'solid'))+
  xlab('Axis 1')+
  ylab('Axis 2')+
  theme_bw(base_size = 12)+
  coord_cartesian(xlim = c(-30,40), ylim = c( -50,30))+
  theme(legend.position = 'none')

match.ellipse.pcoa.burned <- match.pcoa.df%>%
  filter(group == 'fire')%>%
  ggplot(aes(V1,V2))+
  geom_point(aes(color = fire.occ), shape = 16, alpha  = 0.3)+
  stat_ellipse(aes(color = fire.occ, linetype = fire.occ),  level = 0.95, linewidth = 2)+
  scale_color_manual(values = mycolors[c(3,4)], name = 'Group', breaks = c(  'before','after'), labels = c('Burned: Before','Burned: After'))+
 # scale_linetype_manual(values = c('dashed', 'solid'))+
  xlab('Axis 1')+
  ylab('Axis 2')+
  coord_cartesian(xlim = c(-30,40), ylim = c( -50,30))+
  theme_bw(base_size = 12)+
  theme(legend.position = 'none')



## add sens and pue vectors 
response.df <- matched.all.data%>%
  ungroup()%>%
  dplyr::select(sens, pue, elev)
respfit <- vegan::scores(envfit(rot.points, response.df, na.rm = T), 'vectors')
respfit <- as.data.frame(respfit)

vegfit <- vegan::scores(envfit(rot.points, match.mds.df, na.rm = T), 'vectors')
vegfit <- as.data.frame(vegfit)
colnames(vegfit) <- c('V1', 'V2')


match.envplot <- ggplot(match.pcoa.df, aes(V1,V2))+
  geom_point(aes(color = groups.all), alpha  = 0.15, color = 'grey60')+
  scale_color_manual(values = mycolors, name = 'Group', breaks = c(  'before_unburned','after_unburned', 'before_fire','after_fire'), labels = c('Unburned: Before','Unburned: After',  'Burned: Before', 'Burned: After'))+
  #geom_segment(data = coverscores, aes(x =0, y = 0, xend = V1/1.1, yend =V2), size = 0.75)+
  #geom_text(data = coverscores, aes(x =V1, y = V2, label = c('Annual', 'Perennial', 'Bare', 'Shrub')), color = 'red', size = 4)+
  geom_segment(data = vegfit, aes(x =0, y = 0, xend = V1*20, yend =V2*20 ), size = 0.75, arrow = arrow(length = unit(0.03, "npc")))+
  geom_text(data = vegfit, aes(x =V1*23, y = V2*23, label = c('Annual', 'Perennial', 'Bare', 'Shrub')), color = 'black', size = 4)+
  labs(x = 'Axis 1', y = 'Axis 2')+
  theme_bw(base_size = 12)+
  coord_cartesian(xlim = c(-30,40), ylim = c( -50,30))+
  theme(legend.position = 'none')

c(
   "#0072B2",
  "#56B4E9",
  "#D55E00",
  "#E69F00",
  "#BB0000",
  '#E41A1C',
  
  
  "#CC3311"
)


lgd.plt <- ggplot(match.pcoa.df, aes(V1,V2))+
  geom_point(aes(color = groups.all, shape =groups.all), alpha  = 0.2)+
  stat_ellipse(aes(color = groups.all, linetype = groups.all), level = 0.9, size = 1)+
  scale_color_manual(values = mycolors, name = 'Group', breaks = c(  'before_unburned','after_unburned', 'before_fire','after_fire'), labels = c('Unburned: Before','Unburned: After',  'Fire: Before', 'Fire: After'))+
  
  scale_linetype_manual(values = c( 'solid','dashed', 'solid', 'dashed'), name = 'Group', breaks = c(  'before_unburned','after_unburned', 'before_fire','after_fire'), labels = c('Unburned: Before','Unburned: After',  'Fire: Before', 'Fire: After'))+
  scale_shape_manual(values = c(16,16,17,17), name = 'Group', breaks = c(  'before_unburned','after_unburned', 'before_fire','after_fire'), labels = c('Unburned: Before','Unburned: After',  'Fire: Before', 'Fire: After'))+
  guides(color=guide_legend(ncol=2))+
  geom_text(data = coverscores, aes(x =V1, y = V2, label = rownames(coverscores)))+
  xlab('Axis 1')+
  ylab('Axis 2')+
  theme_bw(base_size = 12)+
  theme(legend.position = 'top')

lgd <- get_legend(lgd.plt)

#cover dids
#bare
plt.m.bare <- lmer(data = matched.all.data, m.bare ~ -1 + group.fire.occ +  (1|pixel.id))
bare.plt.m.output <- summary(plt.m.bare)$coefficients%>%
  cbind(as.data.frame(confint(plt.m.bare))[c(3:6),])%>%
  rename(low.ci = `2.5 %`,
         high.ci = `97.5 %`)%>%
  mutate(rwnms = str_sub(rownames(.), start = 15))%>%
  separate(rwnms, into = c('group', 'fire.occ'))

bare.did.plt <- ggplot(bare.plt.m.output, aes(x = fct_relevel(fire.occ,'before', 'after'), y = Estimate, color = group))+
  # geom_point(aes(x = fire.occ, y = sens, color = group), alpha = 0.05)+
  #geom_jitter(alpha = 0.2)+
  geom_point(aes(shape = group), position = position_dodge(0.025), size = 2.5)+
  geom_line(aes(group = group), position = position_dodge(0.025))+
  geom_segment( aes(y = low.ci, yend = high.ci ),position = position_dodge(0.025))+
  scale_color_manual(values = c('tomato', "#0072B2"))+
  theme_bw(base_size = 13)+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('Bare ground %')

# bare.did.plt <- ggplot(matched.all.data, aes(x = fct_relevel(fire.occ,'before', 'after'), y = m.bare, group = group))+
#   stat_summary(fun.data = mean_cl_normal, size = 0.25, aes(color = group),position = position_dodge(0.025))+
#   stat_summary(fun =mean, geom="line", aes(x = fire.occ, y = m.bare, color = group),position = position_dodge(0.025))+
#   scale_color_manual(values = c( 'tomato',"#0072B2"))+
#   theme_bw(base_size = 12)+
#   theme(legend.position = 'none')+
#   xlab('')+
#   ylab('Bare ground')


plt.m.annual <- lmer(data = matched.all.data, m.annual ~ -1 + group.fire.occ +  (1|pixel.id))
annual.plt.m.output <- summary(plt.m.annual)$coefficients%>%
  cbind(as.data.frame(confint(plt.m.annual))[c(3:6),])%>%
  rename(low.ci = `2.5 %`,
         high.ci = `97.5 %`)%>%
  mutate(rwnms = str_sub(rownames(.), start = 15))%>%
  separate(rwnms, into = c('group', 'fire.occ'))
annual.did.plt <- ggplot(annual.plt.m.output, aes(x = fct_relevel(fire.occ,'before', 'after'), y = Estimate, color = group))+
  # geom_point(aes(x = fire.occ, y = sens, color = group), alpha = 0.05)+
  #geom_jitter(alpha = 0.2)+
  geom_point(aes(shape = group), position = position_dodge(0.025), size = 2.5)+
  geom_line(aes(group = group), position = position_dodge(0.025))+
  geom_segment( aes(y = low.ci, yend = high.ci ),position = position_dodge(0.025))+
  scale_color_manual(values = c('tomato', "#0072B2"))+
  theme_bw(base_size = 13)+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('Annual herb cover %')
 # ggplot(matched.all.data, aes(x = fct_relevel(fire.occ,'before', 'after'), y = m.annual , group = group))+
 #  stat_summary(fun.data = mean_cl_normal, size = 0.25, aes(color = group, shape = group),position = position_dodge(0.025))+
 #  stat_summary(fun =mean, geom="line", aes(x = fire.occ, y = m.annual, color = group),position = position_dodge(0.025))+
 #  scale_color_manual(values = c('tomato', "#0072B2"))+
 #  theme_bw(base_size = 12)+
 #  theme(legend.position = 'none')+
 #  xlab('')+
 #  ylab('Annual herb cover')


plt.m.perennial <- lmer(data = matched.all.data, m.perennial ~ -1 + group.fire.occ +  (1|pixel.id))
perennial.plt.m.output <- summary(plt.m.perennial)$coefficients%>%
  cbind(as.data.frame(confint(plt.m.perennial))[c(3:6),])%>%
  rename(low.ci = `2.5 %`,
         high.ci = `97.5 %`)%>%
  mutate(rwnms = str_sub(rownames(.), start = 15))%>%
  separate(rwnms, into = c('group', 'fire.occ'))
perennial.did.plt <- ggplot(perennial.plt.m.output, aes(x = fct_relevel(fire.occ,'before', 'after'), y = Estimate, color = group))+
  # geom_point(aes(x = fire.occ, y = sens, color = group), alpha = 0.05)+
  #geom_jitter(alpha = 0.2)+
  geom_point(aes(shape = group), position = position_dodge(0.025), size = 2.5)+
  geom_line(aes(group = group), position = position_dodge(0.025))+
  geom_segment( aes(y = low.ci, yend = high.ci ),position = position_dodge(0.025))+
  scale_color_manual(values = c('tomato', "#0072B2"))+
  theme_bw(base_size = 13)+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('Perennial herb cover %')
# perennial.did.plt <- ggplot(matched.all.data, aes(x = fct_relevel(fire.occ,'before', 'after'), y = m.perennial , group = group))+
#   stat_summary(fun.data = mean_cl_normal, size = 0.25, aes(color = group),position = position_dodge(0.025))+
#   stat_summary(fun =mean, geom="line", aes(x = fire.occ, y = m.perennial, color = group),position = position_dodge(0.025))+
#   scale_color_manual(values =  c('tomato', "#0072B2"))+
#   theme_bw(base_size = 12)+
#   theme(legend.position = 'none')+
#   xlab('')+
#   ylab('Perennial herb cover')

plt.m.shrub <- lmer(data = matched.all.data, m.shrub ~ -1 + group.fire.occ +  (1|pixel.id))
shrub.plt.m.output <- summary(plt.m.shrub)$coefficients%>%
  cbind(as.data.frame(confint(plt.m.shrub))[c(3:6),])%>%
  rename(low.ci = `2.5 %`,
         high.ci = `97.5 %`)%>%
  mutate(rwnms = str_sub(rownames(.), start = 15))%>%
  separate(rwnms, into = c('group', 'fire.occ'))
shrub.did.plt <- ggplot(shrub.plt.m.output, aes(x = fct_relevel(fire.occ,'before', 'after'), y = Estimate, color = group))+
  # geom_point(aes(x = fire.occ, y = sens, color = group), alpha = 0.05)+
  #geom_jitter(alpha = 0.2)+
  geom_point(aes(shape = group), position = position_dodge(0.025), size = 2.5)+
  geom_line(aes(group = group), position = position_dodge(0.025))+
  geom_segment( aes(y = low.ci, yend = high.ci ),position = position_dodge(0.025))+
  scale_color_manual(values = c('tomato', "#0072B2"))+
  theme_bw(base_size = 13)+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('Shrub cover %')
# shrub.did.plt <- ggplot(matched.all.data, aes(x = fct_relevel(fire.occ,'before', 'after'), y = m.shrub , group = group))+
#   stat_summary(fun.data = mean_cl_normal, size = 0.25, aes(color = group),position = position_dodge(0.025))+
#   stat_summary(fun = mean, geom="line", aes(x = fire.occ, y = m.shrub, color = group),position = position_dodge(0.025))+
#   scale_color_manual(values =  c('tomato', "#0072B2"))+
#   theme_bw(base_size = 12)+
#   theme(legend.position = 'none')+
#   xlab('')+
#   ylab('Shrub cover')

lgd.plt.cov <- ggplot(annual.plt.m.output, aes(x = fct_relevel(fire.occ,'before', 'after'), y = Estimate, color = group))+
  # geom_point(aes(x = fire.occ, y = sens, color = group), alpha = 0.05)+
  #geom_jitter(alpha = 0.2)+
  geom_point(aes(shape = group), position = position_dodge(0.025), size = 2.5)+
  geom_line(aes(group = group), position = position_dodge(0.025))+
  geom_segment( aes(y = low.ci, yend = high.ci ),position = position_dodge(0.025))+
  scale_color_manual(values = c('tomato', "#0072B2"), name = 'Group', labels = c('Fire', 'Unburned'))+
  scale_shape( name = 'Group', labels = c('Fire', 'Unburned'))+
  theme_bw(base_size = 13)+
  theme(legend.position = 'top')+
  xlab('')+
  ylab('Annual herb cover %')
cov.lgd <- get_legend(lgd.plt.cov)
matchedpcoas <-ggarrange( lgd,
                          ggarrange(match.envplot, match.ellipse.pcoa.unburned, match.ellipse.pcoa.burned, nrow = 1, widths = c(3, 3, 3), labels = c('(a)', '(b)', '(c)'),
                                    font.label = list(face = "plain", size = 13, color = "black")),
                          cov.lgd,
                          ggarrange(bare.did.plt, annual.did.plt, perennial.did.plt, shrub.did.plt, nrow = 1, widths = c(2.25, 2.25, 2.25, 2.25), labels = c('(d)', '(e)', '(f)', '(g)'),
                                    font.label = list(face = "plain", size = 13, color = "black"),
                                    
                                    label.x = 0.25,
                                    label.y = 0.95),
                          nrow = 4,
                          heights = c(0.6, 3,0.4, 3) )
matchedpcoas
ggsave('matchpcoas_rot.svg', matchedpcoas, width = 9, height = 7)


#get mean changes
firedid.v1 <- lmer( V1 ~ fire.occ + group + fire.occ:group + (1|pixel.id), data = match.pcoa.df)
firedid.v2 <- lmer( V2 ~ fire.occ + group + fire.occ:group + (1|pixel.id), data = match.pcoa.df)

summary(firedid.v1)
summary(firedid.v2)

dif.pcoa.df <- match.pcoa.df%>%
  group_by(pixel.id, group)%>%
  summarize(dif.V1 = V1[fire.occ == 'after'] - V1[fire.occ == 'before'],
            dif.V2 = V2[fire.occ == 'after'] - V1[fire.occ == 'before'])%>%
  ungroup()

dif.V1.lm <- lm(data = dif.pcoa.df, dif.V1 ~ group)
summary(dif.V1.lm)
confint(dif.V1.lm)

dif.V2.lm <- lm(data = dif.pcoa.df, dif.V2 ~ group)
summary(dif.V2.lm)
confint(dif.V2.lm)
