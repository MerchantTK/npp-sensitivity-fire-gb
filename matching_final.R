#This script takes the cleanned data and applies a matching algorithm to match pixels

#read in masked full data 
did.data.all <- read.csv('data/all_data.csv')



did.data.all%>%
  filter(is.na(pue))
##desctriptive statistics of data 

means <- did.data.all%>%
  group_by(pixel.id)%>%
  summarize(m.mm = mean(mm))
mean(means$m.mm)
max(means$m.mm)
min(means$m.mm)

means <- did.data.all%>%
  group_by(year)%>%
  summarize(m.mm = mean(mm))
mean(means$m.mm)
max(means$m.mm)
min(means$m.mm)

histogram(did.data.all$npp)

##remove anomalous PUE values 

max.pue 
did.data.all <- did.data.all%>% 
  filter(!(pue > (mean(pue) + sd(pue)*3) | pue < (mean(pue) - sd(pue)*3)))%>%
  group_by(pixel.id, fire.occ)%>%
  mutate(n.yrs = n())%>%
  filter(n.yrs > 5)

histogram(sqrt(did.data.all$pue))

#summarize did data - Create dataset with two values for each pixel 
did.sum.data <- did.data.all %>% 
  group_by(pixelxocc, pixel.id, ID, fire.occ, group, elev, fire.yr)%>%
  mutate(npp = npp/10)%>%
  summarise(
    sens = slopefunct(mm, npp),
    int = intfunct(mm, npp),
    pue = mean(npp)/mean(mm),
    m.npp = mean(npp),
    map = mean(mm),
    m.bare = mean(bare),
    m.shrub = mean(shrub),
    m.annual = mean(annual),
    m.perennial = mean(perennial),
    m.litter = mean(litter)
  ) %>%
  #filter(pue < 15)%>% 
  separate(pixel.id, sep = '-',into= c('fire', 'pix'), remove = F)
did.sum.data$pue

#check
missing.data <- did.sum.data%>%
  group_by(pixel.id)%>%
  summarize(nobs = n())%>%
  filter(nobs < 2)

#### matching #####
#install.packages('MatchIt')
library(MatchIt)
#set up dataframe
pre.match <- did.sum.data%>%
  left_join(all.coords, by = c('pixel.id'))%>%
  #left_join(vpd.data, by = c('pixel.id', 'fire.occ'))%>%
  filter(!(pixel.id %in% missing.data$pixel.id))%>%
  filter(fire.occ == 'before')%>%
  mutate(group.num = as.numeric(factor(group, levels = c( 'unburned', 'fire'), ordered = T)) - 1)%>%
  ungroup()



pre.match
# group_by(group.num)%>%
#sample_n(956)
##pre fire t test
summary(lm(data = pre.match, sens ~ -1+ group))
confint(lm(data = pre.match, sens ~ -1 + group))
summary(lm(data = pre.match, pue ~ -1+ group))
confint(lm(data = pre.match, pue ~ -1  + group))


#perform matching on pre-fire cover and envrionmental values
set.seed(123)
matched.out <- matchit(group.num ~ elev + m.npp + map + m.bare + m.shrub + m.annual + m.perennial , data = pre.match, method = 'nearest', distance = 'glm', caliper = 0.05)

summary(matched.out)
#plot(matched.out,type = "jitter",)

matched.out$match.matrix
#identify matched pixels
matched.id <- match.data(matched.out)%>%
  dplyr::select(pixel.id, subclass)
subclass.id <-  match.data(matched.out)%>%
  dplyr::select(pixel.id, subclass, pix, group)
matched.data <- match.data(matched.out)

fire.pair <- match.data(matched.out)%>%
  dplyr::select(subclass, fire.yr, group)%>%
  unique()%>%
  filter(group == 'fire')%>%
  rename(fire.yr.pair = fire.yr)%>%
  dplyr::select(subclass, fire.yr.pair)



matched.id <- match.data(matched.out)%>%
  dplyr::select(pixel.id, subclass)%>%
  unique()%>%
  merge(fire.pair, by = c('subclass'))


#select before and after points for each pixel
matched.all.data <- did.sum.data%>%
  left_join(matched.id, by = 'pixel.id')%>%
  filter( pixel.id %in% matched.data$pixel.id)


matched.all.yrs <- did.data.all%>%
  filter( pixel.id %in% matched.data$pixel.id)%>%
  merge(matched.id, by = 'pixel.id')%>%
  filter(!is.na(subclass))%>%
  mutate(fire.occ = ifelse( year < fire.yr.pair, 'before', ifelse(year > fire.yr.pair + 4 , 'after', 'recovery' )))%>%
  filter(fire.occ != 'recovery')

  
matched.all.data <- matched.all.yrs%>% 
  group_by(pixelxocc, pixel.id, ID, fire.occ, group, elev, fire.yr)%>%
  mutate(npp = npp/10)%>%
  summarise(
    sens = slopefunct(mm, npp),
    int = intfunct(mm, npp),
    pue = mean(npp)/mean(mm),
    m.npp = mean(npp),
    map = mean(mm),
    m.bare = mean(bare),
    m.shrub = mean(shrub),
    m.annual = mean(annual),
    m.perennial = mean(perennial),
    m.litter = mean(litter)
  ) %>%
  #filter(pue < 15)%>% 
  separate(pixel.id, sep = '-',into= c('fire', 'pix'), remove = F)


matched.all.data$group <- fct_relevel(matched.all.data$group, 'unburned', 'fire')
matched.all.data$fire.occ <- fct_relevel(matched.all.data$fire.occ,'before', 'after')

matched.all.data$group.fire.occ <- paste(matched.all.data$group, matched.all.data$fire.occ, sep = '-')

matched.all.data <- matched.all.data%>%
  filter(!is.na(sens))%>%
  mutate(group.num = ifelse(group == 'fire', 1, 0),
         fire.num = ifelse(fire.occ == 'after', 1, 0))

#save matched data

write.csv(matched.all.data, 'data/matched_did_data.csv')

#create unmatched data of the same size
unmatch <- did.sum.data%>%
  filter(fire.occ == 'before')%>%
  mutate(group.num = as.numeric(factor(group, levels = c('unburned', 'fire'), ordered = T)) - 1)%>%
  ungroup()%>%
  group_by(group.num)%>%
  sample_n(917)
###



