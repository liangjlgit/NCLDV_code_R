########working dir###############
setwd("C:\\Users\\fengsw\\Desktop\\")
library(tidyverse)
library(reshape2)
library(ggpubr)
library(vegan)
library(vroom)
library(stringr)
library(ggplot2)
library(basicTrendline)
library(patchwork)
library(igraph)
library(ggraph)
library(readxl)
library(albersusa)
library(biscale)
library(sf)
library(hrbrthemes)
library(ggtext)
library(scatterpie)
library(ggsci)
library(cowplot)
library(ggspatial)
library(ggprism)
############################### R codes for Figure 1a #######################################################
china_data <- read.csv('all_pie2.csv' ,header = 1) #pie file
china_data_pro <- china_data %>% sf::st_as_sf(coords = c("lon", "lat"),crs = 4326)%>%
  st_transform(crs =  "+proj=laea +lat_0=40 +lon_0=104")%>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1], lat = sf::st_coordinates(.)[,2])%>% 
  select(farmland ,forest ,grass ,gobi ,tailings ,regin ,radius  ,lon,lat) %>% 
  data.frame(stringsAsFactors = F) 
china <- sf::read_sf(china_shp)
nine_line <- sf::read_sf(nine)
p <- ggplot()+
  geom_sf(data = china,fill="#99FFFF",size=.125,color="black") + 
  geom_sf(data = nine_line,size=.125) + 
  coord_sf(crs =  "+proj=laea +lat_0=40 +lon_0=104",  xlim = c(-2500000,2877844),ylim = c(-2387082,1654989)) +
  scatterpie::geom_scatterpie(data=china_data_pro, aes(x = lon,y = lat,r=radius*30000),
                              cols=c('farmland','forest','grass','gobi','tailings'),
                              alpha=0.8) +
  annotation_north_arrow(location = "tl", which_north = "false",
                         style = north_arrow_fancy_orienteering) +
  theme(plot.subtitle = element_text(size = 27, 
                                     colour = "red3", hjust = 0.5), axis.line = element_line(colour = "orangered",  size = 1, linetype = "solid"), axis.ticks = element_line(colour = "orangered"), 
        axis.text = element_text(colour = "orangered" ,size = 20), 
        axis.text.x = element_text(colour = "orangered" ,size = 20), 
        axis.text.y = element_text(colour = "orangered" ,size = 20), 
        panel.background = element_rect(fill = "azure2" )) +
  #scale_fill_ucscgb(name="Habitat")+#ggsci
  scale_fill_prism(name="Habitat")+
  labs(x = "Longitude", y = "Latitude", col = "orangered", size = 30,
       fill = "Habitat", subtitle = "Sampling points distribution of  NCLDV in China")+
  theme(axis.ticks = element_line(size = 0.6), 
        axis.title = element_text(size = 22, 
                                  colour = "orangered")) +theme(plot.subtitle = element_text(size = 32)) + 
								  theme(legend.text = element_text(size = 13),  legend.title = element_text(size = 15))
p2 <- p+annotation_scale(location = "bl") +
  geom_scatterpie_legend(china_data_pro$radius*30000,x = -2200000,y=-1800000,n = 3 ,
                         labeller=function(x) x%/%3000)
p2
china <- st_read(china_shp)
pp  <- ggplot() + theme(axis.text.x = element_text(size = 8 ,angle = 60), 
                        axis.text.y = element_text(size = 8))+
  geom_sf(data = china, colour = "black", fill="#99FFFF")+
  coord_sf(xlim = c(117131.4,2115095), ylim = c(-4028017,-1877844),
           crs = "+proj=laea +lat_0=40 +lon_0=104")+
  theme(aspect.ratio = 1.25,
        panel.border = element_rect(fill = NA, colour = "#525252"),
        plot.margin = unit(c(0,0,0,0),"mm")) + theme(axis.text.x = element_text(colour = "black"), 
                                                     axis.text.y = element_text(colour = "black")) +
  theme(axis.text.x = element_text(vjust = 0.5))
aa <- ggdraw() + 
  draw_plot(p2) +
  draw_plot(pp, x = 0.675, y = 0.08, width = 0.1, height = 0.3)
aa
#ggsave(width  =13 ,height   = 7,"123555.pdf")
##################################R codes for Figure 2a,2b, 2c,2d,2e ##############################################################################
#Figure 2a
dt = data.frame(A = c(2,3,12,14), B = c('Asfa','Coco','pan','un'))
dt = dt[order(dt$A, decreasing = TRUE),] 
myLabel = as.vector(dt$B)   
myLabel = paste(myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)        ", sep = "")   
pp <- ggplot(dt, aes(x = "", y = A, fill = B)) + 
  geom_bar(stat = "identity", width = 1,color = 'white' ,alpha = 1 ) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "left") +
  scale_fill_manual(values =c("#dba5b4","#81dd81",'#a3c2ef' ,"#d5a27c"))+
  theme_void()+
  # scale_fill_discrete(breaks = dt$B, labels = myLabel) + 
  theme(axis.text.x = element_blank()) +theme(panel.grid=element_blank()) +    ## 去掉白色圆框和中间的坐标线
  theme(panel.border=element_blank())  +labs(title = NULL, x = NULL, y = NULL)
ppp <- pp + theme(legend.position = "left")
ppp
#Figure 2b 
library(UpSetR) 
a <- read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\ven.csv' ,header = 1 ,row.names = 1) #share data
P<-  upset(data.frame(a), nsets = 5, nintersects = 50, mb.ratio = c(0.7, 0.33),main.bar.color = '#7CC17D',
           matrix.color = '#7CC17D', #
           shade.color = "#f7f7fa",
           matrix.dot.alpha = 0.3,
           queries = list(
             list(query=intersects, params=list("farmland", "forest"), color="#9ECE9B", active=T),
             list(query=intersects, params=list("farmland",  "tailings"), color="#9ECE9B", active=T),
             list(query=intersects, params=list("farmland", "grass" ), color="#9ECE9B", active=T),
             list(query=intersects, params=list("farmland",  "gobi"), color="#9ECE9B", active=T),
             list(query=intersects, params=list("forest" , "tailings"), color="#9ECE9B", active=T),
             list(query=intersects, params=list("forest" , "grass"), color="#9ECE9B", active=T),
             list(query=intersects, params=list("forest" , "gobi"), color="#9ECE9B", active=T),
             list(query=intersects, params=list("tailings", "grass"), color="#9ECE9B", active=T),
             list(query=intersects, params=list("tailings", "gobi"), color="#9ECE9B", active=T),
             list(query=intersects, params=list( "grass", "gobi"), color="#9ECE9B", active=T),
             list(query=intersects, params=list( "farmland" ), color="#C1DDBB", active=T),
             list(query=intersects, params=list( "gobi"), color="#C1DDBB", active=T),
             list(query=intersects, params=list( "grass"), color="#C1DDBB", active=T),
             list(query=intersects, params=list(  "tailings"), color="#C1DDBB", active=T),
             list(query=intersects, params=list( "forest" ), color="#C1DDBB", active=T),
             list(query=intersects, params=list( "forest" , "tailings", "grass", "gobi"), color="#279149", active=T),
             list(query=intersects, params=list("farmland", "forest" ,  "grass", "gobi"), color="#279149", active=T),
             list(query=intersects, params=list("farmland", "forest" , "tailings",  "gobi"), color="#279149", active=T),
             list(query=intersects, params=list("farmland", "forest" , "tailings", "grass"), color="#279149", active=T),
             list(query=intersects, params=list("farmland", "forest" , "tailings", "grass", "gobi"), color="#104928", active=T)),
           mainbar.y.label = "Count of Intersection",  # y 
           number.angles = 0,                     #
           sets.x.label = "The number of NCLDV \nPolBs in habitats",  # x 
           text.scale = c(4, 3.5, 3, 3.5, 3.5,3.5),sets.bar.color=c('#B0BEC5','#EF9A9A','#E1BEE7','#7986CB','#42A5F5'),shade.alpha=0.5,
		   order.by = c("freq"),  point.size =6,
           line.size = 1,                       #
           decreasing = c(T ,F))
#Figure 2c
#env  range
da <- read.csv('D:\\typora dir\\3.environment_range\\physical_all_321_2.csv',header = 1 ,row.names = 1)
max(da)
scale(da)
library('vegan')
library('dplyr')
decostand(da,'max') %>% summary() #
daa <- decostand(da,'range')
DB <- data.frame(rowMeans(daa)) #
DAA <- decostand(DB,'range')

vroom("new_nclav_coverage/nor_coverage.csv", col_names = T) %>%  data.frame() -> da_ncldv
rownames(da_ncldv) <- da_ncldv$...1
da_ncldv %>% dplyr::select(-1) %>% t() %>% data.frame()  %>% mutate(roo = rowSums(.)) -> da_ncldv_sum
da_ncldv_sum %>% filter(roo > 0) %>% select(-roo) -> da_ncldv_t
#write.csv(da_ncldv_t,"env_range_ncldv.csv") #
dd <-  read.csv('env_range_ncldv.csv',header = 1 ,row.names = 1) #
#
fa <- NULL
for (i in c(2:1876)){
  daa <- dd[which(dd[i] > 0), ]
  #d.geo <- distm(daa[1:2] ,fun = distHaversine)/1000
  #dis_geo <- as.data.frame(as.vector(as.dist(d.geo)))
  ranges <- max(daa[1])-min(daa[1])
  mean_env <- mean(as.matrix(daa[1])) #
  mean_cov <-mean(as.matrix(daa[i])) #mean coverage
  sum_cov <-sum(as.matrix(daa[i])) #sum covergae
  occ <- dim(daa[1])[1] #
  all <- cbind(colnames(dd)[i] ,occ ,ranges,mean_env,mean_cov,sum_cov)#
  fa <- rbind(fa ,all )#
}
fa #
#write.csv(fa,"scale_mean_all2.csv")
ggplot(add_habs_da, aes(NU_habs , log(da$mean_cov), color = factor(da$ranges)))+
  geom_jitter(size = 3 ,shape = 19 ,alpha = 0.75,width = 0.2)+
  scale_colour_manual(values = colorRampPalette(c("#c8e6c9", "#2e7d32"))(1020))+
  stat_smooth(method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 28)+scale_y_continuous( limits=c(-2.5, 5))+ theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif",  colour = "black"), axis.text.x = element_text(family = "serif",  colour = "black"), axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = NA)) +labs(x = "Number of habitats (r = -0.374, P < 0.001)", y = "Mean NCLDV abundance (ln)" ) -> p_fig2c
p_fig2c
#Figure 2d
ggplot(data = da ,aes(log(da$occ), log(da$mean_cov), color = factor(da$ranges)))+
  geom_point(size = 4 ,shape = 19 ,alpha = 0.75)+
  scale_colour_manual(values = colorRampPalette(c("#b2ebf2", "#00838f"))(1020))+
  stat_smooth(method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 28)+scale_y_continuous( limits=c(-2.5, 5))+ theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif",  colour = "black"), axis.text.x = element_text(family = "serif",  colour = "black"), axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = NA)) +labs(x = "Number of sampling sites (r = -0.5, P < 0.001)", y = "Mean NCLDV abundance (ln)" ) -> p_fig2d
p_fig2d
#Figure 2e
ggplot(data = da ,aes(decostand(log(da$ranges+1),'max'), log(da$mean_cov), color = factor(da$ranges)))+
  geom_point(size = 4 ,shape = 19 ,alpha = 0.75)+
  scale_colour_manual(values = colorRampPalette(c("#BA68C8", "#7B1FA2"))(1140))+
  stat_smooth(method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 28)+scale_y_continuous( limits=c(-2.5, 5))+ theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"),
                                            plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif",  colour = "black"), axis.text.x = element_text(family = "serif",  colour = "black"), axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = NA)) +labs(x = "Environmental range (r = -0.412, P < 0.001)", y = "Mean NCLDV abundance (ln)" ) -> p_fig2e
trendline(decostand(log(da$ranges+1),'max'), log(da$mean_cov), model="line2P", ePos.x = "topright", 
          summary=TRUE, eDigit=2)
##################################R codes for Figure 3 ##############################################################################
#Figure 3a
ncldv_world <- read.csv("new_nclav_coverage/bases_mapping_fix_global.csv", header = T)
ncldv_world %>% select(1,2,3,4) %>% mutate(hab=str_replace(ID, '.*farmland.*', "Farmland")) %>% mutate(hab=str_replace(hab, '.*forest.*', "Forest")) %>% 
  mutate(hab=str_replace(hab, '.*Gr.*', "Grassland"))  %>% mutate(hab=str_replace(hab, '.*Go.*', "Gobi desert"))  %>%  mutate(hab=str_replace(hab, '.*gr.*', "Grassland")) %>% 
  mutate(hab=str_replace(hab, '.*ta.*', "Mine wasteland"))  %>% mutate(hab=str_replace(hab, '.*D_.*', "Mine wasteland")) %>% group_by(family, hab) %>% count() %>% 
  select(2,1,3) %>% arrange(hab) -> national_count  #map to world NO. of hab and tax
ncldv_world <- read.csv("new_nclav_coverage/bases_mapping_fix_global.csv", header = T)
ncldv_world[,5:292]
ncldv_world[,5:292][ncldv_world[,5:292] < 150] <- 0 
ncldv_world[,5:292][ncldv_world[,5:292] > 0] <- 1 
rownames(ncldv_world) <- ncldv_world$ID
ncldv_world[,5:292] %>% t() %>% data.frame() -> ncldv_world_t
rownames(ncldv_world_t) 
df <- read.csv('E:\\1.aliyun\\ncldv\\0_topsoil\\clean_data_2land2world\\1.csv',header = F)
merge(df,ncldv_world_t,by.x = 'V1', by.y = 0)  %>% select(2:576) %>% group_by(V2) %>%  summarise(across(where(is.numeric), sum)) %>% 
  melt() %>% select(2,1,3) %>% filter(value>0) %>% merge(ncldv_world[,1:4], by.x = "variable", by.y = "ID") %>% select(2,6)  %>% group_by(family,V2) %>%summarise(Nu=n()) -> world_count
colnames(world_count) <- c('hab','family','n')
library(sankeywheel)
rbind(national_count, world_count)  -> na_wo_map
colnames(na_wo_map) <- c('from','to','weight') 
sankeywheel(from = na_wo_map$from,theme = "grid",
            to = na_wo_map$to,
            weight = na_wo_map$weight,
            type = "sankey",
            width = "100%",height = "30%")  
national_count %>% select(hab,n)%>% group_by(hab)%>% summarise(sum(n))
world_count %>% select(family,n)%>% group_by(family)%>% summarise(sum(n))
national_count %>% select(family,n)%>% group_by(family)%>% summarise(sum(n))
#Figure3b
cov_na_world <- read.csv("new_nclav_coverage/coverage_map_world_out.csv",header = T,row.names  = 1)
cov_na_world[cov_na_world < 0.07] <- 0 #
cov_na_world[cov_na_world > 0.06] <- 1 
cov_na_world %>% mutate(summ=rowSums(.) ) %>% filter(summ>0) %>% select(-summ) %>% t() %>% data.frame() %>% mutate(summm = rowSums(.)) %>% filter(summm>0) %>%
  select(-summm) %>%  mutate(site = rownames(.)) %>%  rowwise() %>%  mutate(NU_habs = sum(across(where(is.numeric)))) %>% select(site, NU_habs) %>% 
  merge(read.csv("E:\\1.aliyun\\ncldv\\0_topsoil\\clean_data_2land2world\\dat_world2.csv",header = 1), by.x= 'site', by.y= 'variable') %>% 
  select(1,3,4,2,5) -> No_plob_map_world
world <- map_data("world")
library(ggnewscale)
library(RColorBrewer)
cols<-brewer.pal(8, "YlOrRd")
ggplot() +
  geom_map(data = world, map = world, aes(x=long, y=lat, map_id = region),color = "#d3d3d3", fill = "#d3d3d3", size = 0.2) +theme_void()+
  geom_hline(yintercept = 0, col = "gray90")+
  geom_jitter( data = No_plob_map_world,aes(x=long, y=lat,color =Habitats, fill =Habitats),  size =log(No_plob_map_world$NU_habs+1)*4,alpha = 1,shape= 21 , width = 0) +
  scale_color_manual(values = c('#006064','#dd2c00','#006400','#3d5afe')) + scale_fill_manual(values = c('#006064','#dd2c00','#006400','#3d5afe'))+
  new_scale_color()+  new_scale_fill() +
  #geom_jitter( data = da,aes(x=long, y=lat, color =Habitats,fill =Habitats),  size = log(da$nu)*1,alpha = 1 ,shape=24, width = 0) + 
  scale_color_manual(values = c('#1e90ff','#dd2c00','#ffab00','#006400','#6d4c41'))+
  scale_fill_manual(values =  c('#1e90ff','#dd2c00','#ffab00','#006400','#6d4c41'))
No_plob_map_world$NU_habs %>% max()
No_plob_map_world$NU_habs %>% min()
##################################R codes for Figure 4 ##############################################################################
NCLDV_tax <- read.csv("new_nclav_coverage/ncldv_taxo_rmbasal.csv", header = 1,row.names = 1,na.strings = "")
NCLDV_tax[is.na(NCLDV_tax)] <- "Unassigned"
source("D:/typora dir/tran_tb.R")
#farmland
NCLDV_tax %>% select(3) %>% merge(ncldv_fa, by = 0) %>%  filter(X.1!="Unassigned") %>% group_by(X.1) %>%  summarise(across(where(is.numeric), sum))  %>% data.frame() %>% transpose_df() -> abun_fa
merge(da_env_fa[,1:2], abun_fa, by=0) %>% merge(rich_fa_18s_ncldv,by.x = 'Row.names', by.y = 'site' ) %>% select(3:14, 16) -> abun_fa2
abun_fa2[,2:13] <- as.numeric(unlist(abun_fa2[,2:13]))
aggregate(abun_fa2[,2:13],by=list(abun_fa2$LAT),FUN=mean) -> mean_fa
mean_fa[2:12] <- mean_fa[2:12]/rowSums(mean_fa[2:12])*100
rownames(mean_fa) <- mean_fa$Group.1
richness2 <- data.frame(rich = rep(mean_fa$rich_ncldv,each = 11))
pal2_fa <- c("#5DB1DDFF" ,"#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF" ,"#BA6338FF" ,"#004ae9" ,"#802268FF", "#6BD76BFF", "#D595A7FF", "#924822FF")#, "#837B8DFF"
cut1 = 1.5
t(mean_fa[,2:12]) %>% melt() %>% cbind(richness2)%>% set_names('family','hab','value','richness') -> fa_abun
fa_abun %>% ggplot() +
  geom_bar(stat='identity', aes(factor(hab),value , fill = family), width = 0.9 ,alpha =0.88) +
  labs(x = ' ', y = ' ' ) +
  theme(axis.text = element_text(size = 16 ,family="serif"), axis.title = element_text(size = 20,family="serif")) +
  scale_fill_manual(values = pal2_fa)+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.background = element_blank()) +
  theme(axis.text.x=  element_text(colour = 'black' ,angle = 0 ,size =10 ,vjust = 0.5), axis.text.y=  element_text(colour = 'black' ,angle = 0 ,size = 12)) +
  theme(axis.line = element_line(colour = "black"),  legend.spacing.x = unit(1.0, 'cm') , panel.border = element_blank())+
  theme(axis.text.x = element_text(size = 15),  axis.text.y = element_text(size = 20))+
  theme(axis.ticks = element_line(linetype = "blank"), 
        axis.text = element_text(colour = "gray3"))+ theme(axis.text = element_text(hjust = 0.75,vjust = 0.25), axis.text.x = element_text(vjust = 0.75,angle = 90)) +
  theme(axis.title = element_text(family = "sans", size = 15, face = "bold"), plot.title = element_text(size = 22,  face = "bold", colour = "blue", hjust = 0.5)) +labs(title = "Farmland")+
  geom_point(aes(x = factor(hab),y = richness/1.5) ,shape=16,size = 5 ,colour = 'white' ,alpha = 0.5)+
  scale_y_continuous(expand = c(0, 0),sec.axis = sec_axis(~.*cut1))+
  theme(legend.position = "none")+
  theme(plot.subtitle = element_text(size = 14, colour = "gray14", hjust = 1), plot.caption = element_text(face = "bold", hjust = 1.25), axis.ticks = element_line(colour = "gray14", linetype = "solid")) +
  labs(y = " Relative abundance (%)")+theme(text=element_text(family = 'serif')) + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) -> p1
  #forest
NCLDV_tax %>% select(3) %>% merge(ncldv_fo, by = 0) %>%  filter(X.1!="Unassigned") %>% group_by(X.1) %>%  summarise(across(where(is.numeric), sum))  %>% data.frame() %>% transpose_df() -> abun_fo
merge(da_env_fo[,1:2], abun_fo, by=0) %>% merge(rich_fo_18s_ncldv,by.x = 'Row.names', by.y = 'site' ) %>% select(3:13, 15) -> abun_fo2
abun_fo2[,2:12] <- as.numeric(unlist(abun_fo2[,2:12]))
aggregate(abun_fo2[,2:12],by=list(abun_fo2$LAT),FUN=mean) -> mean_fo
mean_fo[2:11] <- mean_fo[2:11]/rowSums(mean_fo[2:11])*100
rownames(mean_fo) <- mean_fo$Group.1
richness2 <- data.frame(rich = rep(mean_fo$rich_ncldv,each = 10))
library(ggsci)
library(scales)
show_col(pal2)
pal2_fo <- c("#5DB1DDFF" ,"#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF" ,"#BA6338FF" ,"#004ae9" ,"#802268FF", "#6BD76BFF", "#D595A7FF")#, "#924822FF", ,"#837B8DFF"
cut1 = 1.5
t(mean_fo[,2:11]) %>% melt() %>% cbind(richness2)%>% set_names('family','hab','value','richness') ->fo_abun
fo_abun %>% ggplot() +
  geom_bar(stat='identity', aes(factor(hab),value , fill = family), width = 0.9 ,alpha =0.88) +
  labs(x = ' ', y = ' ' ) +
  theme(axis.text = element_text(size = 16 ,family="serif"), axis.title = element_text(size = 20,family="serif")) +
  scale_fill_manual(values = pal2_fo)+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.background = element_blank()) +
  theme(axis.text.x=  element_text(colour = 'black' ,angle = 0 ,size =10 ,vjust = 0.5), axis.text.y=  element_text(colour = 'black' ,angle = 0 ,size = 12)) +
  theme(axis.line = element_line(colour = "black"),  legend.spacing.x = unit(1.0, 'cm') , panel.border = element_blank())+
  theme(axis.text.x = element_text(size = 15),  axis.text.y = element_text(size = 20))+
  theme(axis.ticks = element_line(linetype = "blank"), 
        axis.text = element_text(colour = "gray3"))+ theme(axis.text = element_text(hjust = 0.75,vjust = 0.25), axis.text.x = element_text(vjust = 0.75,angle = 90)) +
  theme(axis.title = element_text(family = "sans", size = 15, face = "bold"), plot.title = element_text(size = 22,  face = "bold", colour = "blue", hjust = 0.5)) +labs(title = "Forest")+
  geom_point(aes(x = factor(hab),y = richness/1.5) ,shape=16,size = 5 ,colour = 'white' ,alpha = 0.5)+
  scale_y_continuous(expand = c(0, 0),sec.axis = sec_axis(~.*cut1))+
  theme(legend.position = "none")+
  theme(plot.subtitle = element_text(size = 14, colour = "gray14", hjust = 1), plot.caption = element_text(face = "bold", hjust = 1.25), axis.ticks = element_line(colour = "gray14", linetype = "solid")) +
  labs(y = " Relative abundance (%)")+theme(text=element_text(family = 'serif')) + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) -> p2
#grassland
NCLDV_tax %>% select(3) %>% merge(ncldv_gr, by = 0) %>%  filter(X.1!="Unassigned") %>% group_by(X.1) %>%  summarise(across(where(is.numeric), sum))  %>% data.frame() %>% transpose_df() -> abun_gr
merge(da_env_gr[,1:2], abun_gr, by=0) %>% merge(rich_gr_18s_ncldv,by.x = 'Row.names', by.y = 'site' ) %>% select(3: (dim(abun_gr)[2]+3),  (dim(abun_gr)[2]+5)) -> abun_gr2
abun_gr2[,2: (dim(abun_gr)[2]+2)] <- as.numeric(unlist(abun_gr2[,2:(dim(abun_gr)[2]+2)]))
aggregate(abun_gr2[,2:(dim(abun_gr)[2]+2)],by=list(abun_gr2$LAT),FUN=mean) -> mean_gr
mean_gr[2:(dim(abun_gr)[2]+1)] <- mean_gr[2:(dim(abun_gr)[2]+1)]/rowSums(mean_gr[2:(dim(abun_gr)[2]+1)])*100
rownames(mean_gr) <- mean_gr$Group.1
richness2 <- data.frame(rich = rep(mean_gr$rich_ncldv,each = (dim(abun_gr)[2])))
pal2_gr <- c("#5DB1DDFF" ,"#CE3D32FF", "#749B58FF",  "#466983FF" ,"#004ae9" ,"#802268FF", "#6BD76BFF", "#D595A7FF")#, "#837B8DFF"
cut1 = 1.5
t(mean_gr[,2:(dim(abun_gr)[2]+1)]) %>% melt() %>% cbind(richness2)%>% set_names('family','hab','value','richness') -> gr_abun
gr_abun %>% ggplot() +
  geom_bar(stat='identity', aes(factor(hab),value , fill = family), width = 0.9 ,alpha =0.88) +
  labs(x = ' ', y = ' ' ) +
  theme(axis.text = element_text(size = 16 ,family="serif"), axis.title = element_text(size = 20,family="serif")) +
  scale_fill_manual(values = pal2_gr)+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.background = element_blank()) +
  theme(axis.text.x=  element_text(colour = 'black' ,angle = 0 ,size =10 ,vjust = 0.5), axis.text.y=  element_text(colour = 'black' ,angle = 0 ,size = 12)) +
  theme(axis.line = element_line(colour = "black"),  legend.spacing.x = unit(1.0, 'cm') , panel.border = element_blank())+
  theme(axis.text.x = element_text(size = 15),  axis.text.y = element_text(size = 20))+
  theme(axis.ticks = element_line(linetype = "blank"), 
        axis.text = element_text(colour = "gray3"))+ theme(axis.text = element_text(hjust = 0.75,vjust = 0.25), axis.text.x = element_text(vjust = 0.75,angle = 90)) +
  theme(axis.title = element_text(family = "sans", size = 15, face = "bold"), plot.title = element_text(size = 22,  face = "bold", colour = "blue", hjust = 0.5)) +labs(title = "Grassland")+
  geom_point(aes(x = factor(hab),y = richness/1.5) ,shape=16,size = 5 ,colour = 'white' ,alpha = 0.5)+
  scale_y_continuous(expand = c(0, 0),sec.axis = sec_axis(~.*cut1))+
  theme(legend.position = "none")+
  theme(plot.subtitle = element_text(size = 14, colour = "gray14", hjust = 1), plot.caption = element_text(face = "bold", hjust = 1.25), axis.ticks = element_line(colour = "gray14", linetype = "solid")) +
  labs(y = " Relative abundance (%)")+theme(text=element_text(family = 'serif')) + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) -> p3
#gobi desert
NCLDV_tax %>% select(3) %>% merge(ncldv_go, by = 0) %>%  filter(X.1!="Unassigned") %>% group_by(X.1) %>%  summarise(across(where(is.numeric), sum))  %>% data.frame() %>% transpose_df() -> abun_go
merge(da_env_go[,1:2], abun_go, by=0) %>% merge(rich_go_18s_ncldv,by.x = 'Row.names', by.y = 'site' ) %>% select(3: (dim(abun_go)[2]+3),  (dim(abun_go)[2]+5)) -> abun_go2
abun_go2[,2: (dim(abun_go)[2]+2)] <- as.numeric(unlist(abun_go2[,2:(dim(abun_go)[2]+2)]))
aggregate(abun_go2[,2:(dim(abun_go)[2]+2)],by=list(abun_go2$LAT),FUN=mean) -> mean_go
mean_go[2:(dim(abun_go)[2]+1)] <- mean_go[2:(dim(abun_go)[2]+1)]/rowSums(mean_go[2:(dim(abun_go)[2]+1)])*100
rownames(mean_go) <- mean_go$Group.1
richness2 <- data.frame(rich = rep(mean_go$rich_ncldv,each = (dim(abun_go)[2])))
pal2_go <- c("#5DB1DDFF" ,"#CE3D32FF", "#749B58FF",  "#466983FF" ,"#004ae9" , "#6BD76BFF")#, "#837B8DFF"
cut1 = 1.5
t(mean_go[,2:(dim(abun_go)[2]+1)]) %>% melt() %>% cbind(richness2)%>% set_names('family','hab','value','richness') -> go_abun
go_abun %>% ggplot() +
  geom_bar(stat='identity', aes(factor(hab),value , fill = family), width = 0.9 ,alpha =0.88) +
  labs(x = ' ', y = ' ' ) +
  theme(axis.text = element_text(size = 16 ,family="serif"), axis.title = element_text(size = 20,family="serif")) +
  scale_fill_manual(values = pal2_go)+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.background = element_blank()) +
  theme(axis.text.x=  element_text(colour = 'black' ,angle = 0 ,size =10 ,vjust = 0.5), axis.text.y=  element_text(colour = 'black' ,angle = 0 ,size = 12)) +
  theme(axis.line = element_line(colour = "black"),  legend.spacing.x = unit(1.0, 'cm') , panel.border = element_blank())+
  theme(axis.text.x = element_text(size = 15),  axis.text.y = element_text(size = 20))+
  theme(axis.ticks = element_line(linetype = "blank"), 
        axis.text = element_text(colour = "gray3"))+ theme(axis.text = element_text(hjust = 0.75,vjust = 0.25), axis.text.x = element_text(vjust = 0.75,angle = 90)) +
  theme(axis.title = element_text(family = "sans", size = 15, face = "bold"), plot.title = element_text(size = 22,  face = "bold", colour = "blue", hjust = 0.5)) +labs(title = "Gobi desert")+
  geom_point(aes(x = factor(hab),y = richness/1.5) ,shape=16,size = 5 ,colour = 'white' ,alpha = 0.5)+
  scale_y_continuous(expand = c(0, 0),sec.axis = sec_axis(~.*cut1))+
  theme(legend.position = "none")+
  theme(plot.subtitle = element_text(size = 14, colour = "gray14", hjust = 1), plot.caption = element_text(face = "bold", hjust = 1.25), axis.ticks = element_line(colour = "gray14", linetype = "solid")) +
  labs(y = " Relative abundance (%)")+theme(text=element_text(family = 'serif')) + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) -> p4
#mine wasteland
NCLDV_tax %>% select(3) %>% merge(ncldv_mine, by = 0) %>%  filter(X.1!="Unassigned") %>% group_by(X.1) %>%  summarise(across(where(is.numeric), sum))  %>% data.frame() %>% transpose_df() -> abun_mine
merge(da_env_mine[,1:2], abun_mine, by=0) %>% merge(rich_mine_18s_ncldv,by.x = 'Row.names', by.y = 'site' ) %>% select(3: (dim(abun_mine)[2]+3),  (dim(abun_mine)[2]+5)) -> abun_mine2
abun_mine2[,2: (dim(abun_mine)[2]+2)] <- as.numeric(unlist(abun_mine2[,2:(dim(abun_mine)[2]+2)]))
aggregate(abun_mine2[,2:(dim(abun_mine)[2]+2)],by=list(abun_mine2$LAT),FUN=mean) -> mean_mine
mean_mine[2:(dim(abun_mine)[2]+1)] <- mean_mine[2:(dim(abun_mine)[2]+1)]/rowSums(mean_mine[2:(dim(abun_mine)[2]+1)])*100
rownames(mean_mine) <- mean_mine$Group.1
richness2 <- data.frame(rich = rep(mean_mine$rich_ncldv,each = (dim(abun_mine)[2])))
pal2_mine<- c("#5DB1DDFF" ,"#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF" ,"#004ae9" , "#6BD76BFF", "#D595A7FF")#, "#837B8DFF"
cut1 = 1.5
t(mean_mine[,2:(dim(abun_mine)[2]+1)]) %>% melt() %>% cbind(richness2)%>% set_names('family','hab','value','richness') -> mine_abun
mine_abun %>% ggplot() +
  geom_bar(stat='identity', aes(factor(hab),value , fill = family), width = 0.9 ,alpha =0.88) +
  labs(x = ' ', y = ' ' ) +
  theme(axis.text = element_text(size = 16 ,family="serif"), axis.title = element_text(size = 20,family="serif")) +
  scale_fill_manual(values = pal2_mine)+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.background = element_blank()) +
  theme(axis.text.x=  element_text(colour = 'black' ,angle = 0 ,size =10 ,vjust = 0.5), axis.text.y=  element_text(colour = 'black' ,angle = 0 ,size = 12)) +
  theme(axis.line = element_line(colour = "black"),  legend.spacing.x = unit(1.0, 'cm') , panel.border = element_blank())+
  theme(axis.text.x = element_text(size = 15),  axis.text.y = element_text(size = 20))+
  theme(axis.ticks = element_line(linetype = "blank"), 
        axis.text = element_text(colour = "gray3"))+ theme(axis.text = element_text(hjust = 0.75,vjust = 0.25), axis.text.x = element_text(vjust = 0.75,angle = 90)) +
  theme(axis.title = element_text(family = "sans", size = 15, face = "bold"), plot.title = element_text(size = 22,  face = "bold", colour = "blue", hjust = 0.5)) +labs(title = "Mine wasteland")+
  geom_point(aes(x = factor(hab),y = richness/1.5) ,shape=16,size = 5 ,colour = 'white' ,alpha = 0.5)+
  scale_y_continuous(expand = c(0, 0),sec.axis = sec_axis(~.*cut1))+
  theme(legend.position = "none")+
  theme(plot.subtitle = element_text(size = 14, colour = "gray14", hjust = 1), plot.caption = element_text(face = "bold", hjust = 1.25), axis.ticks = element_line(colour = "gray14", linetype = "solid")) +
  labs(y = " Relative abundance (%)")+theme(text=element_text(family = 'serif')) + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) -> p5
#combined
library(patchwork)
layout <- 'ABBB 
CDDD 
EFFF'
p4+p1+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  p3+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  p2+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  p4+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  p5+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  plot_layout(design = layout)
##################################R codes for Figure 5 ##############################################################################
#richness
rich_fa_18s_ncldv %>% mutate(hab= rep(each=dim(rich_fa_18s_ncldv)[1], "Farmland")) %>% rbind( mutate(rich_fo_18s_ncldv,hab= rep(each=dim(rich_fo_18s_ncldv)[1], "Forest"))) %>% 
  rbind( mutate(rich_gr_18s_ncldv,hab= rep(each=dim(rich_gr_18s_ncldv)[1], "Grassland"))) %>% rbind( mutate(rich_go_18s_ncldv,hab= rep(each=dim(rich_go_18s_ncldv)[1], "Gobi desert")))%>% 
  rbind( mutate(rich_mine_18s_ncldv,hab= rep(each=dim(rich_mine_18s_ncldv)[1], "Mine wasteland"))) -> rich_18s_ncldv_combine
rich_18s_ncldv_combine$hab  <- factor(rich_18s_ncldv_combine$hab,levels = c('Farmland','Forest','Grassland','Gobi desert','Mine wasteland'))

rich_18s_ncldv_combine %>%  ggboxplot( x="hab", y="rich_ncldv", color = "hab", bxp.errorbar = T, outlier.shape = NA,
                                                palette = c("#3B4992FF" ,"#CC3333", "#008B45FF" ,"#631879FF" ,"#008280FF") ,alpha = 0.7,
                                                add = "jitter", shape=1 ,size = 1.5)+
  theme_pubr(base_size = 20)+theme(text = element_text(family = 'serif'))+theme_test(base_size = 20 )+theme(text = element_text(family = 'serif'))+ theme(legend.position = "none")+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = "white")) +labs(x = "", y = "The number of phylotypes")+
  theme(axis.title = element_text(size = 30, vjust = 1)) + 
  theme(axis.title = element_text(vjust = 1), axis.text = element_text(hjust = 1), axis.text.x = element_text(vjust = 1,hjust = 1 ,angle = 30), axis.text.y = element_text(vjust = 0))+theme(axis.text = element_text(face = "bold")) +labs(x = NULL)+
  theme(text = element_text(family = 'serif'))+  ylim(c(-10,160))  -> p_rich_ncldv
#randomForest, Figure 5b-5f
library(randomForest)
#install.packages(c('rfPermute' ,'A3'))
library(rfPermute)
library(A3)
library(ggplot2)

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\fa_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
#set.seed(2755)
randomForest(rich_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
set.seed(2755)
rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> fa_RF_result
fa.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
fa.pval
fa_RF_result2 <- cbind(fa_RF_result,row.names(fa_RF_result))
fa_RF_result2$`row.names(fa_RF_result)` <- factor(fa_RF_result2$`row.names(fa_RF_result)` ,level = fa_RF_result2$`row.names(fa_RF_result)`)
ggplot(fa_RF_result2) +
  aes(x = `row.names(fa_RF_result)`, weight = X.IncMSE) +
  geom_bar(fill = "#3B4992FF" ,color = 'black',alpha = 0.85 ,width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=c(0,18)) +
  theme_test(base_size = 20 ) + 
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), panel.background = element_rect(fill = "white")) +
  labs(x = "", y = "Increased in MSE (%)") +
  theme(axis.title = element_text(size = 30, vjust = 1)) + 
  theme(axis.title = element_text(vjust = 1), axis.text = element_text(hjust = 1), axis.text.x = element_text(vjust = 1,hjust = 1 ,angle = 30), axis.text.y = element_text(vjust = 0)) +
  theme(axis.text = element_text(face = "bold")) +labs(x = NULL) +
  theme(text = element_text(family = 'serif')) -> p_fa_RF 

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\fo_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(12213)
randomForest(rich_ncldv~., data = tmp_RF[, c(1:7,9:10)], importance = TRUE) %>%importance( type = 1) 
set.seed(12213)
rfPermute(rich_ncldv~., data = tmp_RF[, c(1:7,9:10)], na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> fo_RF_result
fo.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
fo.pval
fo_RF_result2 <- cbind(fo_RF_result,row.names(fo_RF_result))
fo_RF_result2$`row.names(fo_RF_result)` <- factor(fo_RF_result2$`row.names(fo_RF_result)` ,level = fo_RF_result2$`row.names(fo_RF_result)`)
ggplot(fo_RF_result2) +
  aes(x = `row.names(fo_RF_result)`, weight = X.IncMSE) +
  geom_bar(fill = "#CC3333" ,color = 'black',alpha = 0.8 ,width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=c(0,max(fo_RF_result2$X.IncMSE+2))) +
  theme_test(base_size = 20 ) + 
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), panel.background = element_rect(fill = "white")) +
  labs(x = "", y = "Increased in MSE (%)") +
  theme(axis.title = element_text(size = 30, vjust = 1)) + 
  theme(axis.title = element_text(vjust = 1), axis.text = element_text(hjust = 1), axis.text.x = element_text(vjust = 1,hjust = 1 ,angle = 30), axis.text.y = element_text(vjust = 0)) +
  theme(axis.text = element_text(face = "bold")) +labs(x = NULL) +
  theme(text = element_text(family = 'serif')) -> p_fo_RF 
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\gr_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
#set.seed(64499)
randomForest(rich_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
set.seed(64499)
rfPermute(rich_ncldv~., data = tmp_RF,na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> gr_RF_result
gr.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
gr.pval
gr_RF_result2 <- cbind(gr_RF_result,row.names(gr_RF_result))
gr_RF_result2$`row.names(gr_RF_result)` <- factor(gr_RF_result2$`row.names(gr_RF_result)` ,level = gr_RF_result2$`row.names(gr_RF_result)`)
ggplot(gr_RF_result2) +
  aes(x = `row.names(gr_RF_result)`, weight = X.IncMSE) +
  geom_bar(fill = "#008B45FF" ,color = 'black',alpha = 0.75 ,width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=c(0,max(gr_RF_result2$X.IncMSE+2))) +
  theme_test(base_size = 20 ) + 
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), panel.background = element_rect(fill = "white")) +
  labs(x = "", y = "Increased in MSE (%)") +
  theme(axis.title = element_text(size = 30, vjust = 1)) + 
  theme(axis.title = element_text(vjust = 1), axis.text = element_text(hjust = 1), axis.text.x = element_text(vjust = 1,hjust = 1 ,angle = 30), axis.text.y = element_text(vjust = 0)) +
  theme(axis.text = element_text(face = "bold")) +labs(x = NULL) +
  theme(text = element_text(family = 'serif')) -> p_gr_RF 
#devtools::install_github('ericarcher/swfscMisc')
library('swfscMisc')
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\go_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(10079)
randomForest(rich_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
set.seed(10079)
rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>%  data.frame() -> go_RF_result
go.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
go.pval
set.seed(123)
otu_forest.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
otu_forest.pval
go_RF_result2 <- cbind(go_RF_result,row.names(go_RF_result))
go_RF_result2$`row.names(go_RF_result)` <- factor(go_RF_result2$`row.names(go_RF_result)` ,level = go_RF_result2$`row.names(go_RF_result)`)
ggplot(go_RF_result2) +
  aes(x = `row.names(go_RF_result)`, weight = X.IncMSE) +
  geom_bar(fill = "#631879FF" ,color = 'black',alpha = 0.75 ,width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=c(0,max(go_RF_result2$X.IncMSE+2))) +
  theme_test(base_size = 20 ) + 
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), panel.background = element_rect(fill = "white")) +
  labs(x = "", y = "Increased in MSE (%)") +
  theme(axis.title = element_text(size = 30, vjust = 1)) + 
  theme(axis.title = element_text(vjust = 1), axis.text = element_text(hjust = 1), axis.text.x = element_text(vjust = 1,hjust = 1 ,angle = 30), axis.text.y = element_text(vjust = 0)) +
  theme(axis.text = element_text(face = "bold")) +labs(x = NULL) +
  theme(text = element_text(family = 'serif')) -> p_go_RF 
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\ta_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(10079)
randomForest(rich_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
set.seed(10079)
rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>%  data.frame() -> mine_RF_result
set.seed(123)
mine.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
mine.pval
mine_RF_result2 <- cbind(mine_RF_result,row.names(mine_RF_result))
mine_RF_result2$`row.names(mine_RF_result)` <- factor(mine_RF_result2$`row.names(mine_RF_result)` ,level = mine_RF_result2$`row.names(mine_RF_result)`)
ggplot(mine_RF_result2) +
  aes(x = `row.names(mine_RF_result)`, weight = X.IncMSE) +
  geom_bar(fill = "#008280FF" ,color = 'black',alpha = 0.75 ,width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=c(0,max(mine_RF_result2$X.IncMSE+2))) +
  theme_test(base_size = 20 ) + 
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), panel.background = element_rect(fill = "white")) +
  labs(x = "", y = "Increased in MSE (%)") +
  theme(axis.title = element_text(size = 30, vjust = 1)) + 
  theme(axis.title = element_text(vjust = 1), axis.text = element_text(hjust = 1), axis.text.x = element_text(vjust = 1,hjust = 1 ,angle = 30), axis.text.y = element_text(vjust = 0)) +
  theme(axis.text = element_text(face = "bold")) +labs(x = NULL) +
  theme(text = element_text(family = 'serif')) -> p_mine_RF 
library(patchwork)
p_rich_ncldv +p_fa_RF+p_fo_RF+p_gr_RF+p_go_RF+p_mine_RF
##################################R codes for Figure 6 ##############################################################################
#Figure 6a-e
library(vegan)
library(geosphere)
library(ggpubr)
#fa
dim(merge_18s_env_ncldv) 
fa_da_18S_t %>% merge(merge_fa_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
#18S
merge_18s_env_ncldv[,2:(dim(fa_da_18S_t)[2]+1)] -> vpa_fa_18s
vegdist(vpa_fa_18s,method = "bray") %>%  cmdscale( k = (nrow(vpa_fa_18s) - 1), eig = TRUE)  -> pcoafung
pcoafung$points[,1:5] #
pcoa_exp <- pcoafung$eig/sum(pcoafung$eig) #
sum(round(100*pcoa_exp[1:50], 2)) #
vpa_fa_18s_pcoa <- pcoafung$points[,1:50]  #
#geo
merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+2):(dim(fa_da_18S_t)[2]+4)] -> vpa_fa_geo
#clima
merge_18s_env_ncldv[(dim(fa_da_18S_t)[2]+5)] -> vpa_fa_clima
#phy_che
merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+6):(dim(fa_da_18S_t)[2]+15)] -> vpa_fa_phyche
#ncldv
merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+16):dim(merge_18s_env_ncldv)[2]] -> vpa_fa_nclav
#vpa_fa
mod_fa <- varpart(vpa_fa_nclav ,vpa_fa_18s_pcoa , vpa_fa_geo, vpa_fa_clima, vpa_fa_phyche ,transfo="hel")
plot(mod_fa ,digits = 1 ,bg = c("#d50000", "#70f3ff", "#ffff00","#6200ea") ,cex = 1.5 ,alpha=70,
     adj = 0.6 ,lwd = 1 ,lty =1,fg = 'black',pch = 9,family = "serif", Xnames = c( NA))
#fo
dim(fo_da_18S_t) 
fo_da_18S_t %>% merge(merge_fo_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
#18S
merge_18s_env_ncldv[,2:(dim(fo_da_18S_t)[2]+1)] -> vpa_fo_18s
vegdist(vpa_fo_18s,method = "bray") %>%  cmdscale( k = (nrow(vpa_fo_18s) - 1), eig = TRUE)  -> pcoafung
pcoafung$points[,1:5] #
pcoa_exp <- pcoafung$eig/sum(pcoafung$eig) #
sum(round(100*pcoa_exp[1:50], 2)) #
vpa_fo_18s_pcoa <- pcoafung$points[,1:50]  #
#geo
merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+2):(dim(fo_da_18S_t)[2]+4)] -> vpa_fo_geo
#clima
merge_18s_env_ncldv[(dim(fo_da_18S_t)[2]+5)] -> vpa_fo_clima
#phy_che
merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+6):(dim(fo_da_18S_t)[2]+15)] -> vpa_fo_phyche
#ncldv
merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+16):dim(merge_18s_env_ncldv)[2]] -> vpa_fo_nclav
#vpa_fo
mod_fo <- varpart(vpa_fo_nclav ,vpa_fo_18s_pcoa , vpa_fo_geo, vpa_fo_clima, vpa_fo_phyche ,transfo="hel")
plot(mod_fo ,digits = 1 ,bg = c("#d50000", "#70f3ff", "#ffff00","#6200ea") ,cex = 1.5 ,alpha=70,
     adj = 0.6 ,lwd = 1 ,lty =1,fg = 'black',pch = 9,family = "serif", Xnames = c( NA))
#gr
dim(gr_da_18S_t) 
gr_da_18S_t %>% merge(merge_gr_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
#18S
merge_18s_env_ncldv[,2:(dim(gr_da_18S_t)[2]+1)] -> vpa_gr_18s
vegdist(vpa_gr_18s,method = "bray") %>%  cmdscale( k = (nrow(vpa_gr_18s) - 1), eig = TRUE)  -> pcoafung
pcoafung$points[,1:5] #
pcoa_exp <- pcoafung$eig/sum(pcoafung$eig) #
sum(round(100*pcoa_exp[1:11], 2)) #
vpa_gr_18s_pcoa <- pcoafung$points[,1:11]  #
#geo
merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+2):(dim(gr_da_18S_t)[2]+4)] -> vpa_gr_geo
#clima
merge_18s_env_ncldv[(dim(gr_da_18S_t)[2]+5)] -> vpa_gr_clima
#phy_che
merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+6):(dim(gr_da_18S_t)[2]+15)] -> vpa_gr_phyche
#ncldv
merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+16):dim(merge_18s_env_ncldv)[2]] -> vpa_gr_nclav
#vpa_fo
mod_gr <- varpart(vpa_gr_nclav ,vpa_gr_18s_pcoa , vpa_gr_geo, vpa_gr_clima, vpa_gr_phyche ,transfo="hel")
plot(mod_gr ,digits = 1 ,bg = c("#d50000", "#70f3ff", "#ffff00","#6200ea") ,cex = 1.5 ,alpha=70,
     adj = 0.6 ,lwd = 1 ,lty =1,fg = 'black',pch = 9,family = "serif", Xnames = c( NA))
#go
dim(go_da_18S_t) 
go_da_18S_t %>% merge(merge_go_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
#18S
merge_18s_env_ncldv[,2:(dim(go_da_18S_t)[2]+1)] -> vpa_go_18s
vegdist(vpa_go_18s,method = "bray") %>%  cmdscale( k = (nrow(vpa_go_18s) - 1), eig = TRUE)  -> pcoafung
pcoafung$points[,1:5] #
pcoa_exp <- pcoafung$eig/sum(pcoafung$eig) #
sum(round(100*pcoa_exp[1:6], 2)) #
vpa_go_18s_pcoa <- pcoafung$points[,1:4]  #
#geo
merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+2):(dim(go_da_18S_t)[2]+4)] -> vpa_go_geo
#clima
merge_18s_env_ncldv[(dim(go_da_18S_t)[2]+5)] -> vpa_go_clima
#phy_che
merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+6):(dim(go_da_18S_t)[2]+15)] -> vpa_go_phyche
da_geo_edu <- vegdist(vpa_go_phyche,method = "euclidean")
pcoan<- cmdscale(da_geo_edu, k = (10 - 1), eig = TRUE)
pcoan$points[,1:3] #
pcoa_expn <- pcoan$eig/sum(pcoan$eig) #
sum(round(100*pcoa_expn[1:1], 2)) #
n_fa_pcoa <- pcoan$points[,1:1]
#ncldv
merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+16):dim(merge_18s_env_ncldv)[2]] -> vpa_go_nclav
#vpa_fo
mod_go <- varpart(vpa_go_nclav ,vpa_go_18s_pcoa , vpa_go_geo, vpa_go_clima, n_fa_pcoa ,transfo="hel")
plot(mod_go ,digits = 1 ,bg = c("#d50000", "#70f3ff", "#ffff00","#6200ea") ,cex = 1.5 ,alpha=70,
     adj = 0.6 ,lwd = 1 ,lty =1,fg = 'black',pch = 9,family = "serif", Xnames = c( NA))
#mine
dim(mine_da_18S_t) 
mine_da_18S_t %>% merge(merge_mine_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
#18S
merge_18s_env_ncldv[,2:(dim(mine_da_18S_t)[2]+1)] -> vpa_mine_18s
vegdist(vpa_mine_18s,method = "bray") %>%  cmdscale( k = (nrow(vpa_mine_18s) - 1), eig = TRUE)  -> pcoafung
pcoafung$points[,1:5] #
pcoa_exp <- pcoafung$eig/sum(pcoafung$eig) #
sum(round(100*pcoa_exp[1:50], 2)) #
vpa_mine_18s_pcoa <- pcoafung$points[,1:50]  #
#geo
merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+2):(dim(mine_da_18S_t)[2]+4)] -> vpa_mine_geo
#clima
merge_18s_env_ncldv[(dim(mine_da_18S_t)[2]+5)] -> vpa_mine_clima
#phy_che
merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+6):(dim(mine_da_18S_t)[2]+15)] -> vpa_mine_phyche
#ncldv
merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+16):dim(merge_18s_env_ncldv)[2]] -> vpa_mine_nclav
#vpa_fo
mod_mine <- varpart(vpa_mine_nclav ,vpa_mine_18s_pcoa , vpa_mine_geo, vpa_mine_clima, vpa_mine_phyche ,transfo="hel")
plot(mod_mine ,digits = 1 ,bg = c("#d50000", "#70f3ff", "#ffff00","#6200ea") ,cex = 1.5 ,alpha=70,
     adj = 0.6 ,lwd = 1 ,lty =1,fg = 'black',pch = 9,family = "serif", Xnames = c( NA))
#7:7

#Figure 6f-g
#control dis
rbind(partial_mantel_fa_result, partial_mantel_fo_result, partial_mantel_gr_result, partial_mantel_go_result, partial_mantel_mine_result) %>% data.frame()  -> mantel_re_habs
mantel_re_habs[,3:4] <- as.data.frame(lapply(mantel_re_habs[,3:4],as.numeric))
mantel_re_habs %>% mutate(sig = case_when(X4 <= 0.001 ~'***' , (X4 > 0.001 & X4 <= 0.01) ~'**', (X4 > 0.01 & X4 <= 0.05) ~'*')) %>% mutate(X2= str_replace(X2, "WC", "Moisture")) %>% mutate(X2= str_replace(X2, ".*EXCa.*", "EX-Ca")) %>% mutate(X2= str_replace(X2, ".*Clays.*", "Clay")) -> mantel_re_habs_sig
mantel_re_habs_sig$X2 <- factor(mantel_re_habs_sig$X2, levels = c('ALT','CEC','Clay','EC','EX-Ca','MAP','Moisture','pH','TC','TK','TN','TP','Eukaryota'))
mantel_re_habs_sig$group <- factor(mantel_re_habs_sig$group, levels = rev(c('Farmland','Forest','Grassland','Gobi desert','Mine wasteland')))
mantel_re_habs_sig -> mantel_re_habs_sig_con_dis
ggplot(mantel_re_habs_sig_con_dis) +
  geom_point(aes(x=group ,y = X2 ,size=X3),color = '#f47f7e'  ,alpha =0.9)+
  scale_shape_manual(values = 19) +
  scale_size_continuous(range=c(5,18))+ #
  theme(legend.position = "none")+
  theme_bw(base_family = "serif" ,base_size = 20)+
  theme(axis.text.x = element_text(colour = "black", angle = 30 ,hjust = 1), axis.text.y = element_text(colour = "black", angle = 0 ,hjust = 1, vjust = 0))+
  geom_text(aes(x=group, y = X2,label = sig), family = 'serif' ,color = 'black' ,size =8)+coord_flip()+ theme(legend.position = "none") -> p_partial_mantel_habs
#control Eu
rbind(partial2_mantel_fa_result, partial2_mantel_fo_result, partial2_mantel_gr_result, partial2_mantel_go_result, partial2_mantel_mine_result) %>% data.frame()  -> mantel_re_habs
mantel_re_habs[,3:4] <- as.data.frame(lapply(mantel_re_habs[,3:4],as.numeric))
mantel_re_habs %>% mutate(sig = case_when(X4 <= 0.001 ~'***' , (X4 > 0.001 & X4 <= 0.01) ~'**', (X4 > 0.01 & X4 <= 0.05) ~'*')) %>% mutate(X2= str_replace(X2, "WC", "Moisture")) %>% mutate(X2= str_replace(X2, ".*EXCa.*", "EX-Ca")) %>% mutate(X2= str_replace(X2, ".*Clays.*", "Clay")) -> mantel_re_habs_sig
mantel_re_habs_sig$X2 <- factor(mantel_re_habs_sig$X2, levels = c('ALT','CEC','Clay','EC','EX-Ca','MAP','Moisture','pH','TC','TK','TN','TP','Distance'))
mantel_re_habs_sig$group <- factor(mantel_re_habs_sig$group, levels = rev(c('Farmland','Forest','Grassland','Gobi desert','Mine wasteland')))
ggplot(mantel_re_habs_sig) +
  geom_point(aes(x=group ,y = X2 ,size=X3),color = '#41c9df'  ,alpha =0.9)+
  scale_shape_manual(values = 19) +
  scale_size_continuous(range=c(4,12))+ #
  theme(legend.position = "none")+
  theme_bw(base_family = "serif" ,base_size = 20)+
  theme(axis.text.x = element_text(colour = "black", angle = 30 ,hjust = 1), axis.text.y = element_text(colour = "black", angle = 0 ,hjust = 1, vjust = 0))+
  geom_text(aes(x=group, y = X2,label = sig), family = 'serif' ,color = 'black' ,size =8)+coord_flip()+ theme(legend.position = "none") -> p_partial2_mantel_habs
p_partial_mantel_habs+theme(axis.text.y = element_blank() )+labs(x = NULL, y = NULL)+p_partial2_mantel_habs+labs(x = NULL, y = NULL)
##################################R codes for Figure 7 ##############################################################################
#farmland
merge_fa[,2:(dim(ncldv_fa_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fa[,2:(dim(ncldv_fa_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%
  filter(col>8)  %>% select(-col) -> merge_fa_naldv_select 
merge_fa[,(dim(ncldv_fa_t)[2]+2):dim(merge_fa)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fa[,(dim(ncldv_fa_t)[2]+2):dim(merge_fa)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>8) %>% select(-col)  -> merge_fa_18S_select 
library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(merge_fa_naldv_select), y=t(merge_fa_18S_select),type = "spearman") 
net <- CoDF(cor_ls$r, cor_ls$P) 
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_fa
dasign_fa <- read.csv('18S/fa_net_padj10%_desi.csv' ,header = 1 )#
aaa <- net_fa[order(net_fa$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))
nodeattrib_root_combine$indicgroup <- 0
all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)
V(all_root_net)$degree <- degree(as.matrix(all_root_net))
factor(dasign_fa$treat)
color_pal1 <- c('#ffff00','#01579b','#90caf9','#673ab7','#bdbdbd','#388e3c','#f50057','#26a69a',   
               '#795548','#ffd180','#ff8125')
all_root_net -> faaaa
p <- ggraph(faaaa, layout = 'kk', maxiter =20000) +#
  geom_node_voronoi(fill = 'black')+
  geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + 
   geom_edge_density(fill = '#FEF2F7' ,n=600, ) +
  geom_node_point(aes(size =V(faaaa)$degree, shape =dasign_fa$treat2,color=dasign_fa$treat), alpha =1 ) +  
   scale_size_continuous(range = c(4.5, 12) )+ 
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
   scale_color_manual(values = color_pal1)+
  geom_node_text(aes(label=dasign_fa$group),color="#F2FBFC",check_overlap = F,size = 5)+
  theme_void() +  
  expand_limits(x = c(-0.25, 0.25), y = c(-0.5, 0.5))  
p
p + theme(legend.position = "none") -> p_fa
p_fa
# forest
dim(merge_fo)
merge_fo[,2:(dim(ncldv_fo_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fo[,2:(dim(ncldv_fo_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>7)  %>% select(-col) -> merge_fo_naldv_select 
merge_fo[,(dim(ncldv_fo_t)[2]+2):dim(merge_fo)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fo[,(dim(ncldv_fo_t)[2]+2):dim(merge_fo)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>7) %>% select(-col)  -> merge_fo_18S_select 
cor_ls <- rcorr(x=t(merge_fo_naldv_select), y=t(merge_fo_18S_select),type = "spearman") 
net <- CoDF(cor_ls$r, cor_ls$P) 
net$padj <- p.adjust(net$p, method="BH")
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_fo
dasign_fo <- read.csv('18S/fo_net_padj10%_desi.csv' ,header = 1 )#
library(ggraph)
library(igraph)
aaa <- net_fo[order(net_fo$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))
nodeattrib_root_combine$indicgroup <- 0
all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)
V(all_root_net)$degree <- degree(as.matrix(all_root_net))
factor(dasign_fo$treat)
color_pal2 <- c('#ffff00',     '#01579b',     '#90caf9',   '#673ab7',          '#bdbdbd',        '#388e3c',  '#f50057',      '#26a69a',   '#795548',   '#ffd180',     "#595800",     '#ff8125', "#F5515D")
all_root_net -> foooo
p <- ggraph(foooo, layout = 'kk', maxiter =20000) +#
  geom_node_voronoi(fill = 'black')+
    geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + #
    geom_edge_density(fill = '#FEF2F7' ,n=600 ) +
  geom_node_point(aes(size =V(foooo)$degree, shape =dasign_fo$treat2,color=dasign_fo$treat), alpha =1 ) +  #
   scale_size_continuous(range = c(4.5, 12) )+ #
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
   scale_color_manual(values = color_pal2)+
  geom_node_text(aes(label=dasign_fo$group),color="#F2FBFC",size = 5)+
  theme_void() +  #
  expand_limits(x = c(-0.25, 0.25), y = c(-0.5, 0.5))  
p
p + theme(legend.position = "none") -> p_fo
p_fo
#grass
merge_gr[,2:(dim(ncldv_gr_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_gr[,2:(dim(ncldv_gr_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%
  filter(col>7)  %>% select(-col) -> merge_gr_naldv_select #

merge_gr[,(dim(ncldv_gr_t)[2]+2):dim(merge_gr)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_gr[,(dim(ncldv_gr_t)[2]+2):dim(merge_gr)[2]], 2, FUN =function(x) sum(x>0)))  %>%
  filter(col>8) %>% select(-col)  -> merge_gr_18S_select #

cor_ls <- rcorr(x=t(merge_gr_naldv_select), y=t(merge_gr_18S_select),type = "spearman") 
net <- CoDF(cor_ls$r, cor_ls$P) 
net$padj <- p.adjust(net$p, method="BH")
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_gr
dasign_gr <- read.csv('18S/gr_net_padj10%_desi.csv' ,header = 1 )#
aaa <- net_gr[order(net_gr$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))
nodeattrib_root_combine$indicgroup <- 0
all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)
V(all_root_net)$degree <- degree(as.matrix(all_root_net))
factor(dasign_gr$treat)
color_pal3 <- c('#ffff00',     '#01579b',                    '#673ab7',          '#bdbdbd',        '#388e3c',  '#f50057',      '#26a69a',                '#ffd180',     "#595800",     '#ff8125'    )
all_root_net ->grrrr
p <- ggraph(grrrr, layout = 'kk', maxiter =10000) +#
  geom_node_voronoi(fill = 'black')+
  geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + #
  geom_edge_density(fill = '#FEF2F7' ,n=600 ) +
  geom_node_point(aes(size =V(grrrr)$degree, shape =dasign_gr$treat2,color=dasign_gr$treat), alpha =1 ) +  #
   scale_size_continuous(range = c(4.5, 12) )+ #
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
  scale_color_manual(values = color_pal3)+
  geom_node_text(aes(label=dasign_gr$group),color="#F2FBFC", size = 5)+
  theme_void() +  #
  expand_limits(x = c(-0.5, 0.5), y = c(-0.5, 0.5)) 
p
p + theme(legend.position = "none") -> p_gr
p_gr
#gobi
dim(merge_go)
merge_go[,2:(dim(ncldv_go_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_go[,2:(dim(ncldv_go_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%
  filter(col>3)  %>% select(-col) -> merge_go_naldv_select #

merge_go[,(dim(ncldv_go_t)[2]+2):dim(merge_go)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_go[,(dim(ncldv_go_t)[2]+2):dim(merge_go)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>6) %>% select(-col)  -> merge_go_18S_select 
cor_ls <- rcorr(x=t(merge_go_naldv_select), y=t(merge_go_18S_select),type = "spearman") #
net <- CoDF(cor_ls$r, cor_ls$P) #
net$padj <- p.adjust(net$p, method="BH")
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 #
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_go
dasign_go <- read.csv('18S/go_net_padj10%_desi.csv' ,header = 1 )#
aaa <- net_go[order(net_go$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))
nodeattrib_root_combine$indicgroup <- 0
all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)
V(all_root_net)$degree <- degree(as.matrix(all_root_net))
color_pal4 <- c(   '#388e3c',  '#f50057',    "#595800",     '#ff8125', "#F5515D"               )
all_root_net -> goooo
p <- ggraph(goooo, layout = 'kk', maxiter =10000) +#
  geom_node_voronoi(fill = 'black')+
  geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + 
  geom_edge_density(fill = '#FEF2F7' ,n=600 ) +
  geom_node_point(aes(size =V(goooo)$degree, shape =dasign_go$treat2,color=dasign_go$treat), alpha =1 ) +  #，
  scale_size_continuous(range = c(4.5, 6) )+ #
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
  scale_color_manual(values = color_pal4)+
  geom_node_text(aes(label=dasign_go$group),color="#F2FBFC",size = 5)+
  theme_void() +  #
  expand_limits(x = c(-0.5, 0.5), y = c(-0.5, 0.5))  
p
p + theme(legend.position = "none") -> p_go
p_go
#mine
merge_mine[,2:(dim(ncldv_mine_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_mine[,2:(dim(ncldv_mine_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>7)  %>% select(-col) -> merge_mine_naldv_select # 
merge_mine[,(dim(ncldv_mine_t)[2]+2):dim(merge_mine)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_mine[,(dim(ncldv_mine_t)[2]+2):dim(merge_mine)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>10) %>% select(-col)  -> merge_mine_18S_select #
cor_ls <- rcorr(x=t(merge_mine_naldv_select), y=t(merge_mine_18S_select),type = "spearman") 
net <- CoDF(cor_ls$r, cor_ls$P) 
net$padj <- p.adjust(net$p, method="BH")
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_mine
dasign_mine <- read.csv('18S/mine_net_padj10%_desi.csv' ,header = 1 )#
aaa <- net_mine[order(net_mine$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))
nodeattrib_root_combine$indicgroup <- 0
all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)
V(all_root_net)$degree <- degree(as.matrix(all_root_net))
color_pal5 <- c('#ffff00',     '#01579b',     '#673ab7','#bdbdbd',        '#388e3c',  '#f50057',      '#26a69a',     '#795548',   "#595800",     '#ff8125', "#F5515D")
all_root_net -> mineeee
p <- ggraph(mineeee, layout = 'kk', maxiter =20000) +#
  geom_node_voronoi(fill = 'black')+
   geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + #
   geom_edge_density(fill = '#FEF2F7' ,n=600 ) +
  geom_node_point(aes(size =V(mineeee)$degree, shape =dasign_mine$treat2,color=dasign_mine$treat), alpha =1 ) +  #
  scale_size_continuous(range = c(4.5, 16) )+ #
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
   scale_color_manual(values = color_pal5)+ 
  geom_node_text(aes(label=dasign_mine$group),color="#F2FBFC", size = 5)+
  theme_void() +  #
  expand_limits(x = c(-0.5, 0.5), y = c(-0.5, 0.5))  
p
p + theme(legend.position = "none") -> p_mine
p_mine

library(patchwork)
p_fa+p_fo+p_gr+p_go+p_mine #16：24
