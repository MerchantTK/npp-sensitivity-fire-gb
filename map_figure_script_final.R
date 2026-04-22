#map figure with historgrams 
#update 1/29/26

##figure 2 : distrubutions pre and post match#####

require(tidyverse)
require(rasterVis)
require(terra)
require(tidyterra)
require(maps)
require(USAboundaries)
library(stars)
library(sf)

#
matched.all.data <- read.csv('data/matched_pcoa_data.csv')
did.sum.data <- read.csv( 'data/all_did_sum_data.csv')

#create unmatched data of the same size
unmatch <- did.sum.data%>%
  filter(fire.occ == 'before')%>%
  mutate(group.num = as.numeric(factor(group, levels = c('unburned', 'fire'), ordered = T)) - 1)%>%
  ungroup()%>%
  group_by(group.num)%>%
  sample_n(956)


# Read the Shapefile
study_area <- st_read("HC_NVpoly/HC_NVpoly.shp")
gs_class_sf <- st_read("focal_class_sf.shp")



boundary <-st_as_sfc(st_bbox(study_area))
study_area_gs <- st_intersection( gs_class_sf, boundary)
sim_study_gs <- st_transform(study_area_gs, "ESRI:54009")%>%
  st_simplify(dTolerance = 300)%>%
  st_transform('EPSG:4326')%>%
  st_make_valid()

focal_fires_sf <- st_read('mtbs_perimeter_data/mtbs_perims_DD.shp')%>%
  st_transform(4326)%>%
  filter(Year %in% c(1995:2002))%>%
  st_intersection(boundary)

sim_study_gs_full <- sim_study_gs%>%
  st_intersection(study_area)
ggplot()+
  geom_sf(data = sim_study_gs_full,fill = 'grey70', color = 'white')+
  theme_minimal()

rest.treats.hc <- st_intersection( st_make_valid(rest.treats), sim_study_gs)
focal_fires_gs_sf <- st_intersection( focal_fires_sf, sim_study_gs)
focal.fire.rest <- st_intersection( rest.treats.hc, sim_study_gs)

rest.treats.full <- st_intersection( st_make_valid(rest.treats), sim_study_gs_full)
fire.rest.full <-  st_intersection(rest.treats.full, focal_fires_sf)

focal_fires_gs_full <- st_intersection( focal_fires_sf, sim_study_gs_full)
focal.fire.rest.full <- st_intersection( rest.treats.hc, focal_fires_gs_full)

#supplemental figure 
rest.fig <- ggplot()+
  geom_sf(data = sim_study_gs_full,fill = 'grey70', color = 'white')+
  geom_sf(data = focal_fires_gs_full, fill ='grey30', color = 'transparent')+
  geom_sf(data = focal.fire.rest.full, fill ='red', color = 'transparent', alpha = 9)+
  theme_minimal()





study_area_fires <- ggplot(study_area_gs)+
  geom_sf(fill = 'white')+
  geom_sf(data = focal_fires_sf, fill = 'red', alpha = 0.5)+
  #geom_sf(data = indian_lands_study, fill = 'grey30')+
  ggspatial::annotation_scale(location = 'tr', text_cex = 1.2, pad_y = unit(-0.5, "cm"), pad_x = unit(0.5, "cm"))+
  ggspatial::annotation_north_arrow(height = unit(0.7, 'cm'), width = unit(0.6, 'cm'), pad_y = unit(0.8, "cm"), pad_x = unit(0.7, "cm"))+
  #labs(caption = 'Fire boundaries \n 1995-2002')+
  labs(title = 'Fire boundaries')+
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust=0.5 ))
  #theme(plot.caption = element_text(hjust=0.5, size=rel(1.4) ))

#ggsave('study_area_fires.png', study_area_fires, width = 2, height = 4)



# create raster stacks for each cover type then extract mean 
#annual stack
ldf<-list.files(path = 'cover/', pattern = '\\.tif$')

all86 <- terra::rast(paste0('cover/',ldf[1]))
all86mask<- mask(all86, sim_study_gs_full)


ann86plt <- ggplot()+
  geom_spatraster(data = all86mask, aes(fill = cover1986_1))+
  #scale_fill_viridis_c(option = 'turbo')+
  #scale_fill_distiller(palette = 'YlGnBu', trans = 'reverse', )+
  #scale_fill_gradient2(low = 'blue', mid = 'yellow', high = 'green', midpoint = 30)+
  scale_fill_gradient2(low = '#D3E6FF', mid = "#328CEC", high = '#0F63BD', midpoint = 40, na.value = 'transparent')+
  labs(title = 'Annual', fill = '%Cover')+
  guides(fill = guide_colorbar(reverse = TRUE))+
  theme_void(base_size = 11)+
  coord_sf(clip = "off")+
  theme(legend.position = "bottom",
        plot.margin = margin(5.5, 5.5, 25, 5.5))
  
  
#coord_sf(clip = "off") & theme(plot.margin = margin(5.5, 5.5, 25, 5.5))

per86plt <- ggplot()+
  geom_spatraster(data = all86mask, aes(fill = cover1986_4))+
  #scale_fill_viridis_c(option = 'turbo')+
  #scale_fill_distiller(palette = 'YlGnBu', trans = 'reverse', )+
  #scale_fill_gradient2(low = 'steelblue', mid = 'yellow', high = 'green', midpoint = 30)+
  scale_fill_gradient2(low = '#D3E6FF', mid = "#328CEC", high = '#0F63BD', midpoint =  40, na.value = 'transparent')+
  labs(title = 'Perennial', fill = '%Cover')+
  guides(fill = guide_colorbar(reverse = TRUE))+
  coord_sf(clip = "off")+
  theme_void(base_size = 11)+
  theme(legend.position = "bottom",
        plot.margin = margin(5.5, 5.5, 25, 5.5))


bare86plt <- ggplot()+
  geom_spatraster(data = all86mask, aes(fill = cover1986_2))+
  #scale_fill_viridis_c(option = 'turbo')+
  # scale_fill_distiller(palette = 'YlGnBu', trans = 'reverse', )+
  scale_fill_gradient2(low = '#D3E6FF', mid = "#328CEC", high = '#0F63BD', midpoint =  50, na.value = 'transparent')+
  coord_sf(clip = "off")+
  labs(title = 'Bare', fill = '%Cover')+
  guides(fill = guide_colorbar(reverse = TRUE))+
  theme_void(base_size = 11)+
  theme(legend.position = "bottom",
        plot.margin = margin(5.5, 5.5, 25, 5.5))


shrub86plt <- ggplot()+
  geom_spatraster(data = all86mask, aes(fill = cover1986_5))+
  
  #scale_fill_viridis_c(option = 'turbo')+
  #scale_fill_distiller(palette = 'YlGnBu', trans = 'reverse', )+
  scale_fill_gradient2(low = '#D3E6FF', mid = "#328CEC", high = '#0F63BD', midpoint =  30, na.value = 'transparent')+
  labs(title = 'Shrub', fill = '%Cover')+
  coord_sf(clip = "off")+
  guides(fill = guide_colorbar(reverse = TRUE))+
  theme_void(base_size = 11)+
  theme(legend.position = "bottom",
        plot.margin = margin(5.5, 5.5, 25, 5.5))



temp.elev.1 <- raster::raster('nv_dem/USGS_1_n41w118_20130911.tif')
temp.elev.2 <- raster::raster('nv_dem/USGS_1_n41w119_20130911.tif')
temp.elev.3 <- raster::raster('nv_dem/USGS_1_n42w118_20130911.tif')
temp.elev.4 <- raster::raster('nv_dem/USGS_1_n42w119_20130911.tif')

elev <- merge(temp.elev.1, temp.elev.2, temp.elev.3, temp.elev.4)

elevrast <- terra::rast(elev)
crs(elevrast) <- 'EPSG:4326'
elevmask <- mask(elevrast, sim_study_gs_full)

elevplt <- ggplot()+
  geom_spatraster(data = elevrast, aes(fill = layer))+
  scale_fill_gradient2(low = '#D9EAFC', mid = "#328CEC", high = '#0F63BD', midpoint =  2100, na.value = 'transparent')+
  #scale_fill_distiller(palette = 'YlGnBu', trans = 'reverse', )+
  labs(title = 'Elevation', fill = 'Elev. (m)')+
  guides(fill = guide_colorbar(reverse = TRUE))+
  theme_void(base_size = 11)+
  theme(legend.position = "bottom")

ann86plt
per86plt
shrub86plt
bare86plt
elevplt

###### unmatched and matched distributions 
un_ann <- unmatch%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = m.annual, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Pre-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('prematch_ann.png', width = 2, height = 1.5)


mat_ann <- matched.all.data%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = m.annual, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Post-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())


#ggsave('postmatch_ann.png',  width = 2, height = 1.5)

#bare distribution
un_bare <- unmatch%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = m.bare, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Pre-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('prematch_bare.png', width = 2, height = 1.5)



mat_bare <- matched.all.data%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = m.bare, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Post-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('postmatch_bare.png', width = 2, height = 1.5)

#perennial

un_per <- unmatch%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = m.perennial, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Pre-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 12)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('prematch_perennial.png',  width = 2, height = 1.5)



mat_per <- matched.all.data%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = m.perennial, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Post-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('postmatch_perennial.png',  width = 2, height = 1.5)

#shrub
un_shrub <- unmatch%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = m.shrub, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Pre-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('prematch_shrub.png',  width = 2, height = 1.5)



mat_shrub <- matched.all.data%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = m.shrub, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Post-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('postmatch_shrub.png', width = 2, height = 1.5)


#elevation

un_elev <- unmatch%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = elev, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Pre-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('prematch_elev.png',  width = 2, height = 1.5)




mat_elev <- matched.all.data%>%
  filter(fire.occ == 'before')%>%
  ggplot()+
  geom_density(aes(x = elev, fill = group, color = group), alpha = 0.7, size = 1)+
  scale_color_manual(values = c(  'tomato1', 'steelblue'))+
  scale_fill_manual(values = c(  'tomato1', 'steelblue'))+
  labs(x = 'Post-match')+
  guides(color = 'none', fill = 'none')+
  theme_classic(base_size = 11)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
#ggsave('postmatch_elev.png', width = 2, height = 1.5)



#figure 1 arrange
require(gridExtra)
require(ggpubr)
require(patchwork)
ggarrange(un_ann, un_per, un_bare, un_shrub, un_elev, study_area_fires, ann86plt, per86plt, bare86plt, shrub86plt, elevplt,  mat_ann, mat_per, mat_bare, mat_shrub, mat_elev, ncol = 5, nrow = 3, heights = c(1, 2, 1), font.label=list(color="black",size=17))

layout1 <- c('###aaa##ddd#######
              ###aaabbdddee#####
              pppaaaccdddffmmm##
              ppp##########mmmnn
              pppggg##jjj##mmmoo
              ###ggghhjjjkk#####
              ###gggiijjjll#####')
matching.maps.ann <- (ann86plt) + (un_ann / mat_ann) + plot_layout(widths = c(3, 2))+  plot_annotation(
  title = 'Annual cover',
  theme = theme(plot.title = element_text(size = 15, hjust = 0.5))
)
matching.maps.per <- (per86plt) + (un_per / mat_per)+ plot_layout(widths = c(3, 2))+  plot_annotation(
  title = 'Perennial cover',
  theme = theme(plot.title = element_text(size = 15, hjust = 0.5))
)

matching.maps.shrub <- (shrub86plt) + (un_shrub / mat_shrub)+ plot_layout(widths = c(3, 2))+  plot_annotation(
  title = 'Shrub cover',
  theme = theme(plot.title = element_text(size = 15, hjust = 0.5))
)


matching.maps.shrub <- (bare86plt) + (un_bare / mat_bare)+ plot_layout(widths = c(3, 2))+  plot_annotation(
  title = 'Bare cover',
  theme = theme(plot.title = element_text(size = 15, hjust = 0.5))
)

matching.maps.elev <- (elevplt) + (un_elev / mat_elev)+ plot_layout(widths = c(3, 2))+  plot_annotation(
  title = 'Elevation',
  theme = theme(plot.title = element_text(size = 15, hjust = 0.5))
)

matching.maps.ann + matching.maps.per

 free(ann86plt) + un_ann + mat_ann + free(per86plt) + un_per + mat_per + free(shrub86plt) + un_shrub  + mat_shrub + free(bare86plt) +  un_bare + mat_bare  + free(elevplt) +  un_elev + mat_elev  +  study_area_fires + plot_layout(design = layout1)

  un_elev  +    elevplt +   mat_elev + 

matching.maps <- un_ann + un_per + un_bare + un_shrub + un_elev  +  ann86plt + per86plt + bare86plt +  shrub86plt + elevplt +mat_ann + mat_per +mat_bare + mat_shrub + mat_elev + plot_layout(nrow = 3, ncol = 5, heights = c(1,3,1))


study_area_fires_elev <- ggplot(data = elevmask)+
  geom_spatraster(data = elevmask, aes(fill = layer))+
  geom_sf(data = focal_fires_sf, fill = 'red', color = 'red', linewidth = 0.5, alpha = 0.05)+
  #geom_sf(data = indian_lands_study, fill = 'grey30')+
  ggspatial::annotation_scale(location = 'tr', text_cex = 1.2, pad_y = unit(-0.5, "cm"), pad_x = unit(0.5, "cm"))+
  ggspatial::annotation_north_arrow(height = unit(0.5, 'cm'), width = unit(0.4, 'cm'), pad_y = unit(0.5, "cm"), pad_x = unit(0.4, "cm"))+
  scale_fill_gradient2(low = '#D9EAFC', mid = "#328CEC", high = '#0F63BD', midpoint =  2000, na.value = 'transparent')+
  #labs(caption = 'Fire boundaries \n 1995-2002')+
  labs(title = 'Elevation', fill = 'Elev. (m)')+
  theme_minimal(base_size = 11)+
  #theme_void(base_size = 11)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

layout2 <- c('######aaa##ddd##
              mmmm##aaabbdddee
              mmmmnnaaaccdddff
              mmmmooggg##jjj##
              mmmm##ggghhjjjkk
              ######gggiijjjll')


map.fig <-  free(ann86plt) + un_ann + mat_ann + free(per86plt) + un_per + mat_per + free(shrub86plt) + un_shrub  + mat_shrub + free(bare86plt) +  un_bare + mat_bare  + free(study_area_fires_elev) +  un_elev + mat_elev  +   plot_layout(design = layout2)


#ggsave('figures/cover_map_mask.svg', map.fig, width = 9, height = 6)
