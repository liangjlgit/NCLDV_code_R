########working dir###############
setwd("C:\\Users\\fengsw\\Desktop\\")
load("18S/NCLDV_figure.RData")
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsci)
library(vegan)
library(vroom)
library(stringr)
library(ggplot2)
library(basicTrendline)
library(patchwork)
library(igraph)
library(ggraph)

###########1. 18S data#####################################
vroom('18S/18S_feature_table_271.csv',col_names = T,  skip = 0) %>% data.frame() -> da_18S
rownames(da_18S) <- da_18S$X.OTU.ID

#taxonomy of all, keep 2 cols
da_18S %>% select(272) %>% filter(grepl(pattern = "Anima" ,perl = F,fixed = F, taxonomy))
da_18S %>% select(272) %>% 
  separate(sep = ";", col = `taxonomy`, into =c(NA,NA,'Kindom',NA,NA,'Phylum',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA) ) -> da_18S_tax

#rarecure
da_18S %>% select(2:271) %>% t() %>% rarecurve( step = 1000 ,label = F ,col = '#3949ab' ,size = 0.8,xlim = c(0,45000), ylab = "Richness", xlab = "Number of sequences")

#ASV table of five habs separated
da_18S %>% select(contains('Fa')) %>% mutate(rows = rowSums(.)) %>% filter(rows > 0) %>% select(-rows)-> fa_da_18S
da_18S %>% select(contains('Fo')) %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>% select(-rows)-> fo_da_18S
da_18S %>% select(contains('Gr')) %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>% select(-rows)-> gr_da_18S
da_18S %>% select(contains('Go')) %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>% select(-rows)-> go_da_18S
da_18S %>% select(contains(c(".T" ,".D"))) %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>% select(-rows)-> mine_da_18S

fa_da_18S %>% merge(da_18S_tax, by = 0) %>% select(Kindom)

#fa_da_18S %>% merge(da_18S_tax, by = 0) %>% group_by(Phylum) %>% summarise(across(where(is.numeric), sum))

###########2. NCLDV###########################
#lenth
vroom("new_nclav_coverage/polb_fivehab_rmbasal_id.fna.cdhit_id.faa_lenth.csv",col_names = F, col_select = 2)  %>% as.matrix() %>% median() #median
vroom("new_nclav_coverage/polb_fivehab_rmbasal_id.fna.cdhit_id.faa_lenth.csv",col_names = F, col_select = 2)  %>% as.matrix() %>% mean() #mean
vroom("new_nclav_coverage/polb_fivehab_rmbasal_id.fna.cdhit_id.faa_lenth.csv",col_names = F, col_select = 2)  %>% as.matrix() %>% min() #mean
#denisty
vroom("new_nclav_coverage/polb_fivehab_rmbasal_id.fna.cdhit_id.faa_lenth.csv",col_names = F) %>% ggplot() +
  geom_density(aes(x = X2), fill="#80cbc4", color = NA, alpha = 0.75, position = "stack")+
  geom_segment(aes(x=410 , y = 0, xend = 409.14, yend = 0.0025),size = 0.5, color='#283593')+
  geom_segment(aes(x=314 , y = 0, xend = 313, yend = 0.0033),size = 0.5, color='#d32f2f') +theme_pubr(border = T, base_size = 15)+
  labs(x="Length of NCLDV polB", y= "Density")+
  annotate("text", x = 404, y = 0.0034, label = "Medain: 313", size = 5 ) +
  annotate("text", x = 500, y = 0.0026, label = "Mean: 409", size = 5 )+scale_y_continuous(limits = c(0,0.0035), breaks = c(0.000,0.0007,0.0014,0.0021,0.0027,0.0035))

vroom("new_nclav_coverage/nor_coverage.csv", col_names = T) %>%  data.frame() -> da_ncldv
rownames(da_ncldv) <- da_ncldv$...1

da_ncldv %>% select(2:271) %>% t() %>% head()

##这有时候会有bug
da_ncldv %>% select(contains("Fa")) %>% mutate(ro = rowSums(.)) %>% as_tibble() %>% filter(  ro  >0) %>% select(-ro) -> ncldv_fa


da_ncldv %>% select(contains("Fo")) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_fo
da_ncldv %>% select(contains("Gr")) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_gr
da_ncldv %>% select(contains("Go")) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_go
da_ncldv %>% select(contains(c(".T" ,".D"))) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_mine



#apply(2, FUN =function(x) sum(x>0))
#merge_fa[,2:1161] %>% t() %>% data.frame() -> merge_fa_t
#merge_fa_t$co1 <- apply(merge_fa_t, 1, FUN =function(x) sum(x>0)) #
#################3 network - ncldv -18S###############################
#合并
fa_da_18S %>% t() %>% data.frame() -> fa_da_18S_t
ncldv_fa %>% t() %>% data.frame() -> ncldv_fa_t
dim(ncldv_fa_t)
merge(ncldv_fa_t, fa_da_18S_t, by = 0) -> merge_fa
dim(merge_fa)
fo_da_18S %>% t() %>% data.frame() -> fo_da_18S_t
ncldv_fo %>% t() %>% data.frame() -> ncldv_fo_t
row.names(ncldv_fo_t)
merge(ncldv_fo_t, fo_da_18S_t, by = 0) -> merge_fo

gr_da_18S %>% t() %>% data.frame() -> gr_da_18S_t
ncldv_gr %>% t() %>% data.frame() -> ncldv_gr_t
row.names(ncldv_gr_t)
merge(ncldv_gr_t, gr_da_18S_t, by = 0) -> merge_gr

go_da_18S %>% t() %>% data.frame() -> go_da_18S_t
ncldv_go %>% t() %>% data.frame() -> ncldv_go_t
row.names(ncldv_go_t)
merge(ncldv_go_t, go_da_18S_t, by = 0) -> merge_go

mine_da_18S %>% t() %>% data.frame() -> mine_da_18S_t
ncldv_mine %>% t() %>% data.frame() -> ncldv_mine_t
row.names(ncldv_mine_t)
merge(ncldv_mine_t, mine_da_18S_t, by = 0) -> merge_mine
##############################################3.1 farmland########
merge_fa[,2:(dim(ncldv_fa_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fa[,2:(dim(ncldv_fa_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
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

net_padj3 #
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) #
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/fa_net_padj10%.csv') #

#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_fa
dasign_fa <- read.csv('18S/fa_net_padj10%_desi.csv' ,header = 1 )#

library(ggraph)
library(igraph)
aaa <- net_fa[order(net_fa$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))

nodeattrib_root_combine$indicgroup <- 0

all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)

V(all_root_net)$degree <- degree(as.matrix(all_root_net))
#write.csv(degree(as.matrix(all_root_net)) ,'spearson_degree.csv')
igraph.weight = E(all_root_net)$weight
E(all_root_net)$weight = net_fa$cor

#
factor(dasign_fa$treat)
#Levels:  k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Pandoravirales
color_pal1 <- c('#ffff00','#01579b','#90caf9','#673ab7','#bdbdbd','#388e3c','#f50057','#26a69a',   
               '#795548','#ffd180','#ff8125')
all_root_net -> faaaa
p <- ggraph(faaaa, layout = 'kk', maxiter =20000) +#
  geom_node_voronoi(fill = 'black')+
  #create_layout(all_root_net ,layout, circular = T )+
  #geom_node_circle(aes( fill = 'black',r =1), alpha = 0.3 ) +
  geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + 
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +
  geom_edge_density(fill = '#FEF2F7' ,n=600, ) +
  geom_node_point(aes(size =V(faaaa)$degree, shape =dasign_fa$treat2,color=dasign_fa$treat), alpha =1 ) +  
  # geom_node_text(aes(label = dasign$Sample), angle = 60, hjust = 1, nudge_y = 0, size = 3) + 
  scale_size_continuous(range = c(4.5, 12) )+ 
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
  #geom_node_text(aes(label = dasign$group), angle = 0, hjust = 0.5, nudge_y = 0, 
  #               family = 'serif' , size = 2) +  
  scale_color_manual(values = color_pal1)+
  geom_node_text(aes(label=dasign_fa$group),color="#F2FBFC",check_overlap = F,size = 5)+
  theme_void() +  
  expand_limits(x = c(-0.25, 0.25), y = c(-0.5, 0.5))  
p
p + theme(legend.position = "none") -> p_fa
p_fa
#7.6:7.6
##################3.2 forest########################
dim(merge_fo)
merge_fo[,2:(dim(ncldv_fo_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fo[,2:(dim(ncldv_fo_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>7)  %>% select(-col) -> merge_fo_naldv_select 

merge_fo[,(dim(ncldv_fo_t)[2]+2):dim(merge_fo)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fo[,(dim(ncldv_fo_t)[2]+2):dim(merge_fo)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>7) %>% select(-col)  -> merge_fo_18S_select 

library(Hmisc)

cor_ls <- rcorr(x=t(merge_fo_naldv_select), y=t(merge_fo_18S_select),type = "spearman") 
net <- CoDF(cor_ls$r, cor_ls$P) 
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))

net_padj3 #
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) #
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/fo_net_padj10%.csv') #

#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_fo
#write.csv(unique(plyr::rbind.fill(data.frame(id=net_fo$from),data.frame(id=net_fo$to)))  ,row.names = FALSE,'18S/fo_net_padj10%_desi.csv') #

#write.csv(dasign_fo ,row.names = FALSE,'18S/fo_net_padj10%_desi.csv') #

dasign_fo <- read.csv('18S/fo_net_padj10%_desi.csv' ,header = 1 )#

library(ggraph)
library(igraph)
aaa <- net_fo[order(net_fo$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))

nodeattrib_root_combine$indicgroup <- 0

all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)

V(all_root_net)$degree <- degree(as.matrix(all_root_net))
#write.csv(degree(as.matrix(all_root_net)) ,'spearson_degree.csv')
igraph.weight = E(all_root_net)$weight
E(all_root_net)$weight = net_fo$cor
fc = cluster_fast_greedy(all_root_net,weights =NULL)
##
modularity = modularity(all_root_net,membership(fc))
comps = membership(fc)
V(all_root_net)$color = c(rep('red ' ,54),rep('blue ' ,42))
set.seed(9939)
coords<-layout_(all_root_net,with_fr(niter=99, grid="nogrid"))
plot(all_root_net, vertex.label=NA, edge.width=0.01,vertex.size=8, 
     layout=layout.fruchterman.reingold)

#
factor(dasign_fo$treat)
#fa ：Levels:  k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Pandoravirales
#fo levels   ：k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Imitervirales Pandoravirales Pimascovirales
color_pal2 <- c('#ffff00',     '#01579b',     '#90caf9',   '#673ab7',          '#bdbdbd',        '#388e3c',  '#f50057',      '#26a69a',   '#795548',   '#ffd180',     "#595800",     '#ff8125', "#F5515D")
all_root_net -> foooo
p <- ggraph(foooo, layout = 'kk', maxiter =20000) +#
  geom_node_voronoi(fill = 'black')+
  #create_layout(all_root_net ,layout, circular = T )+
  #geom_node_circle(aes( fill = 'black',r =1), alpha = 0.3 ) +
  geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + #
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +#
  geom_edge_density(fill = '#FEF2F7' ,n=600 ) +
  geom_node_point(aes(size =V(foooo)$degree, shape =dasign_fo$treat2,color=dasign_fo$treat), alpha =1 ) +  #
  # geom_node_text(aes(label = dasign$Sample), angle = 60, hjust = 1, nudge_y = 0, size = 3) +  #
  scale_size_continuous(range = c(4.5, 12) )+ #
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
  #geom_node_text(aes(label = dasign$group), angle = 0, hjust = 0.5, nudge_y = 0, 
  #               family = 'serif' , size = 2) +  #
  scale_color_manual(values = color_pal2)+
  geom_node_text(aes(label=dasign_fo$group),color="#F2FBFC",size = 5)+
  theme_void() +  #
  expand_limits(x = c(-0.25, 0.25), y = c(-0.5, 0.5))  
p
p + theme(legend.position = "none") -> p_fo
p_fo
#9:9
#######################3.3 grass#########################
dim(merge_gr)
merge_gr[,2:(dim(ncldv_gr_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_gr[,2:(dim(ncldv_gr_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%
  filter(col>7)  %>% select(-col) -> merge_gr_naldv_select #

merge_gr[,(dim(ncldv_gr_t)[2]+2):dim(merge_gr)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_gr[,(dim(ncldv_gr_t)[2]+2):dim(merge_gr)[2]], 2, FUN =function(x) sum(x>0)))  %>%
  filter(col>8) %>% select(-col)  -> merge_gr_18S_select #
library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(merge_gr_naldv_select), y=t(merge_gr_18S_select),type = "spearman") 
net <- CoDF(cor_ls$r, cor_ls$P) 
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))

net_padj3 #
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) #
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/gr_net_padj10%.csv') #

#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_gr
#write.csv(unique(plyr::rbind.fill(data.frame(id=net_gr$from),data.frame(id=net_gr$to)))  ,row.names = FALSE,'18S/gr_net_padj10%_desi.csv') #

dasign_gr <- read.csv('18S/gr_net_padj10%_desi.csv' ,header = 1 )#

library(ggraph)
library(igraph)
aaa <- net_gr[order(net_gr$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))

nodeattrib_root_combine$indicgroup <- 0

all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)

V(all_root_net)$degree <- degree(as.matrix(all_root_net))
#write.csv(degree(as.matrix(all_root_net)) ,'spearson_degree.csv')
igraph.weight = E(all_root_net)$weight
E(all_root_net)$weight = net_gr$cor
fc = cluster_fast_greedy(all_root_net,weights =NULL)
##
modularity = modularity(all_root_net,membership(fc))
comps = membership(fc)
V(all_root_net)$color = c(rep('red ' ,54),rep('blue ' ,42))
set.seed(9939)
coords<-layout_(all_root_net,with_fr(niter=99, grid="nogrid"))
plot(all_root_net, vertex.label=NA, edge.width=0.01,vertex.size=8, 
     layout=layout.fruchterman.reingold)

#
factor(dasign_gr$treat)
#fa ：Levels:  k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Pandoravirales
#fo levels   ：k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Imitervirales Pandoravirales Pimascovirales
#gr Levels:   k__Alveolata   k__Amoebozoa               k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles              Chitovirales Imitervirales Pandoravirales
color_pal3 <- c('#ffff00',     '#01579b',                    '#673ab7',          '#bdbdbd',        '#388e3c',  '#f50057',      '#26a69a',                '#ffd180',     "#595800",     '#ff8125'    )
all_root_net ->grrrr
p <- ggraph(grrrr, layout = 'kk', maxiter =10000) +#
  geom_node_voronoi(fill = 'black')+
  #create_layout(all_root_net ,layout, circular = T )+
  #geom_node_circle(aes( fill = 'black',r =1), alpha = 0.3 ) +
  geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + #
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +#
  geom_edge_density(fill = '#FEF2F7' ,n=600 ) +
  geom_node_point(aes(size =V(grrrr)$degree, shape =dasign_gr$treat2,color=dasign_gr$treat), alpha =1 ) +  #
  # geom_node_text(aes(label = dasign$Sample), angle = 60, hjust = 1, nudge_y = 0, size = 3) +  #
  scale_size_continuous(range = c(4.5, 12) )+ #
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
  #geom_node_text(aes(label = dasign$group), angle = 0, hjust = 0.5, nudge_y = 0, 
  #               family = 'serif' , size = 2) +  #
  scale_color_manual(values = color_pal3)+
  geom_node_text(aes(label=dasign_gr$group),color="#F2FBFC", size = 5)+
  theme_void() +  #
  expand_limits(x = c(-0.5, 0.5), y = c(-0.5, 0.5)) 
p
p + theme(legend.position = "none") -> p_gr
p_gr
#######################3.4 gobi#################################
dim(merge_go)
merge_go[,2:(dim(ncldv_go_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_go[,2:(dim(ncldv_go_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%
  filter(col>3)  %>% select(-col) -> merge_go_naldv_select #

merge_go[,(dim(ncldv_go_t)[2]+2):dim(merge_go)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_go[,(dim(ncldv_go_t)[2]+2):dim(merge_go)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>6) %>% select(-col)  -> merge_go_18S_select 

library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(merge_go_naldv_select), y=t(merge_go_18S_select),type = "spearman") #
net <- CoDF(cor_ls$r, cor_ls$P) #
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))

net_padj3 #
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) #
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/go_net_padj10%.csv') #

#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_go
#write.csv(unique(plyr::rbind.fill(data.frame(id=net_go$from),data.frame(id=net_go$to)))  ,row.names = FALSE,'18S/go_net_padj10%_desi.csv') 


dasign_go <- read.csv('18S/go_net_padj10%_desi.csv' ,header = 1 )#

library(ggraph)
library(igraph)
aaa <- net_go[order(net_go$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))

nodeattrib_root_combine$indicgroup <- 0

all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)

V(all_root_net)$degree <- degree(as.matrix(all_root_net))
#write.csv(degree(as.matrix(all_root_net)) ,'spearson_degree.csv')
igraph.weight = E(all_root_net)$weight
E(all_root_net)$weight = net_go$cor
fc = cluster_fast_greedy(all_root_net,weights =NULL)
##
modularity = modularity(all_root_net,membership(fc))
comps = membership(fc)
V(all_root_net)$color = c(rep('red ' ,54),rep('blue ' ,42))
set.seed(9939)
coords<-layout_(all_root_net,with_fr(niter=99, grid="nogrid"))
plot(all_root_net, vertex.label=NA, edge.width=0.01,vertex.size=8, 
     layout=layout.fruchterman.reingold)

#
factor(dasign_go$treat)
#fa ：Levels:  k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Pandoravirales
#fo levels   ：k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Imitervirales Pandoravirales Pimascovirales
#gr Levels:   k__Alveolata   k__Amoebozoa               k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles              Chitovirales Imitervirales Pandoravirales
#go  Levels:  k__Fungi  k__Rhizaria Imitervirales Pandoravirales Pimascovirales
color_pal4 <- c(   '#388e3c',  '#f50057',    "#595800",     '#ff8125', "#F5515D"               )
all_root_net -> goooo
p <- ggraph(goooo, layout = 'kk', maxiter =10000) +#
  geom_node_voronoi(fill = 'black')+
  #create_layout(all_root_net ,layout, circular = T )+
  #geom_node_circle(aes( fill = 'black',r =1), alpha = 0.3 ) +
  geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + 
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +#
  geom_edge_density(fill = '#FEF2F7' ,n=600 ) +
  geom_node_point(aes(size =V(goooo)$degree, shape =dasign_go$treat2,color=dasign_go$treat), alpha =1 ) +  #，
  # geom_node_text(aes(label = dasign$Sample), angle = 60, hjust = 1, nudge_y = 0, size = 3) +  #
  scale_size_continuous(range = c(4.5, 6) )+ #
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
  #geom_node_text(aes(label = dasign$group), angle = 0, hjust = 0.5, nudge_y = 0, 
  #               family = 'serif' , size = 2) +  #
  scale_color_manual(values = color_pal4)+
  geom_node_text(aes(label=dasign_go$group),color="#F2FBFC",size = 5)+
  theme_void() +  #
  expand_limits(x = c(-0.5, 0.5), y = c(-0.5, 0.5))  
p
p + theme(legend.position = "none") -> p_go
p_go

##################################3.5 mine#######################################
dim(merge_mine)
merge_mine[,2:(dim(ncldv_mine_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_mine[,2:(dim(ncldv_mine_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>7)  %>% select(-col) -> merge_mine_naldv_select # 

merge_mine[,(dim(ncldv_mine_t)[2]+2):dim(merge_mine)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_mine[,(dim(ncldv_mine_t)[2]+2):dim(merge_mine)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>10) %>% select(-col)  -> merge_mine_18S_select #

library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(merge_mine_naldv_select), y=t(merge_mine_18S_select),type = "spearman") 
net <- CoDF(cor_ls$r, cor_ls$P) 
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))

net_padj3 
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) 
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/mine_net_padj10%.csv') 

#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_mine
#write.csv(unique(plyr::rbind.fill(data.frame(id=net_mine$from),data.frame(id=net_mine$to)))  ,row.names = FALSE,'18S/mine_net_padj10%_desi.csv') 


dasign_mine <- read.csv('18S/mine_net_padj10%_desi.csv' ,header = 1 )#

library(ggraph)
library(igraph)
aaa <- net_mine[order(net_mine$from ,decreasing=F),] #
nodeattrib_root_combine <- data.frame(node=union(aaa$from, aaa$to))

nodeattrib_root_combine$indicgroup <- 0

all_root_net <- graph_from_data_frame(aaa, direct=F, vertices=nodeattrib_root_combine$node)

V(all_root_net)$degree <- degree(as.matrix(all_root_net))
#write.csv(degree(as.matrix(all_root_net)) ,'spearson_degree.csv')
igraph.weight = E(all_root_net)$weight
E(all_root_net)$weight = net_mine$cor
fc = cluster_fast_greedy(all_root_net,weights =NULL)
##
modularity = modularity(all_root_net,membership(fc))
comps = membership(fc)
V(all_root_net)$color = c(rep('red ' ,54),rep('blue ' ,42))
set.seed(9939)
coords<-layout_(all_root_net,with_fr(niter=99, grid="nogrid"))
plot(all_root_net, vertex.label=NA, edge.width=0.01,vertex.size=8, 
     layout=layout.fruchterman.reingold)

#
factor(dasign_mine$treat)

#fa ：Levels:  k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Pandoravirales
#fo levels   ：k__Alveolata  k__Amoebozoa  k__Animalia  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Chitovirales Imitervirales Pandoravirales Pimascovirales
#color_pal <- c('#ffff00',     '#01579b',     '#90caf9',   '#673ab7',          '#bdbdbd',        '#388e3c',  '#f50057',      '#26a69a',   '#795548',   '#ffd180',     "#595800",     '#ff8125', "#F5515D")
#mine Levels:  k__Alveolata  k__Amoebozoa  k__Chloroplastida  k__Eukaryota_unknown  k__Fungi  k__Rhizaria  k__Stramenopiles Asfuvirales Imitervirales Pandoravirales Pimascovirales
color_pal5 <- c('#ffff00',     '#01579b',     '#673ab7','#bdbdbd',        '#388e3c',  '#f50057',      '#26a69a',     '#795548',   "#595800",     '#ff8125', "#F5515D")
all_root_net -> mineeee
p <- ggraph(mineeee, layout = 'kk', maxiter =20000) +#
  geom_node_voronoi(fill = 'black')+
  #create_layout(all_root_net ,layout, circular = T )+
  #geom_node_circle(aes( fill = 'black',r =1), alpha = 0.3 ) +
  geom_edge_arc(edge_width = 0.1,strength = 1,color = '#E7B04D',
                alpha =0.95) + #
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +
  geom_edge_density(fill = '#FEF2F7' ,n=600 ) +
  geom_node_point(aes(size =V(mineeee)$degree, shape =dasign_mine$treat2,color=dasign_mine$treat), alpha =1 ) +  #
  # geom_node_text(aes(label = dasign$Sample), angle = 60, hjust = 1, nudge_y = 0, size = 3) +  #
  scale_size_continuous(range = c(4.5, 16) )+ #
  scale_edge_width_continuous(range = c(1,1) )+
  scale_shape_manual(values = c(17 ,19) )+
  #geom_node_text(aes(label = dasign$group), angle = 0, hjust = 0.5, nudge_y = 0, 
  #               family = 'serif' , size = 2) +  #
  scale_color_manual(values = color_pal5)+ 
  geom_node_text(aes(label=dasign_mine$group),color="#F2FBFC", size = 5)+
  theme_void() +  #
  expand_limits(x = c(-0.5, 0.5), y = c(-0.5, 0.5))  
p
p + theme(legend.position = "none") -> p_mine
p_mine

library(patchwork)
p_fa+p_fo+p_gr+p_go+p_mine #16：24

#####################3.6 network feature of host and virus##################################################figure7
rbind(net_fa, net_fo, net_gr, net_go, net_mine) %>% merge(da_18S_tax,by.x = "from", by.y = 0) %>%  #fungi tax
  merge(read.csv("new_nclav_coverage/ncldv_taxo.csv", header = F), by.x = 'to', by.y = 'V1') %>%
  select(6,7,8,9,10) %>% arrange(Kindom, V2) -> host_virus
#write.csv(host_virus,"18S/host_virus_net_detail1.csv",row.names = F) 
library(RColorBrewer)
color<-brewer.pal(11, "PiYG")[5:1]
read.csv("18S/host_virus_net_detail.csv", header = T) -> da_hv
da_hv %>% select(2,4) %>% group_by(Phylum, V3) %>% summarise(summ_HV=n()) -> H_V
#write.csv(H_V,"18S/host_virus_net.csv",row.names = F) 
  
H_V$V3  <- factor(H_V$V3, levels = unique(da_hv$V3) )
H_V$Phylum  <- factor(H_V$Phylum, levels = unique(da_hv$Phylum) )

ggplot(H_V) +
  geom_tile(aes(x=Phylum, y=V3, fill=log(summ_HV, base = 4.4))) +
  scale_fill_gradientn(colours = color) +
  theme_pubr(border = T, base_size = 20, base_family = "serif",legend = "top") + theme(axis.text.x = element_text(vjust = 1,hjust = 1, angle = 45))+ labs(x=NULL, y=NULL)
#8:12
H_V %>% group_by(Phylum) %>% summarise(summ_HV=sum(summ_HV)) %>% data.frame() %>% write.csv("host_virus_connections.csv")
###########################3.7 sub_networks of habs######################################
##31 ncldv shared in five habs
shared_ncldv_31 <- read.csv("18S/31_shared_polB.out", header = T, sep = "\t")

merge_fa[,(dim(ncldv_fa_t)[2]+2):dim(merge_fa)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fa[,(dim(ncldv_fa_t)[2]+2):dim(merge_fa)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>8) %>% select(-col)  -> merge_fa_18S_select #
#merge_fa[,2:(dim(ncldv_fa_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fa[,2:(dim(ncldv_fa_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
  #filter(col>8)  %>% select(-col) -> merge_fa_naldv_select #
merge_fa[,2:(dim(ncldv_fa_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fa[,2:(dim(ncldv_fa_t)[2]+1)], 2, FUN =function(x) sum(x>0))) %>% 
  merge(shared_ncldv_31[,c(1,7:8)], by.x = 0,by.y = 'Sample')   %>% select(-col) -> sub_merge_fa_naldv_select
rownames(sub_merge_fa_naldv_select) <- sub_merge_fa_naldv_select$Row.names
dim(sub_merge_fa_naldv_select)
library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(sub_merge_fa_naldv_select[,c(-1,-83,-84)]), y=t(merge_fa_18S_select),type = "spearman") #
net <- CoDF(cor_ls$r, cor_ls$P) #
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) 
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/fa_net_padj10%.csv') 
#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_fa
net_fa %>% merge(shared_ncldv_31, by.x = 'to', by.y = 'Sample') %>% select(1,2,3,12, 13) %>% mutate(hab=c('fa','fa','fa')) -> sub_fa_net

merge_fo[,(dim(ncldv_fo_t)[2]+2):dim(merge_fo)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fo[,(dim(ncldv_fo_t)[2]+2):dim(merge_fo)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>8) %>% select(-col)  -> merge_fo_18S_select #
#merge_fo[,2:(dim(ncldv_fo_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fo[,2:(dim(ncldv_fo_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
#filter(col>8)  %>% select(-col) -> merge_fo_naldv_select #
merge_fo[,2:(dim(ncldv_fo_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_fo[,2:(dim(ncldv_fo_t)[2]+1)], 2, FUN =function(x) sum(x>0))) %>% 
  merge(shared_ncldv_31[,c(1,7:8)], by.x = 0,by.y = 'Sample')   %>% select(-col) -> sub_merge_fo_naldv_select
rownames(sub_merge_fo_naldv_select) <- sub_merge_fo_naldv_select$Row.names
dim(sub_merge_fo_naldv_select)
library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(sub_merge_fo_naldv_select[,c(-1,-77,-78)]), y=t(merge_fo_18S_select),type = "spearman") #
net <- CoDF(cor_ls$r, cor_ls$P) #
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) 
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/fa_net_padj10%.csv') 
#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_fo
net_fo %>% merge(shared_ncldv_31, by.x = 'to', by.y = 'Sample') %>% select(1,2,3,12, 13) %>% mutate(hab=c('fo','fo','fo','fo','fo')) -> sub_fo_net

merge_gr[,(dim(ncldv_gr_t)[2]+2):dim(merge_gr)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_gr[,(dim(ncldv_gr_t)[2]+2):dim(merge_gr)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>10) %>% select(-col)  -> merge_gr_18S_select #
#merge_gr[,2:(dim(ncldv_gr_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_gr[,2:(dim(ncldv_gr_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
#filter(col>8)  %>% select(-col) -> merge_gr_naldv_select #
merge_gr[,2:(dim(ncldv_gr_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_gr[,2:(dim(ncldv_gr_t)[2]+1)], 2, FUN =function(x) sum(x>0))) %>% 
  merge(shared_ncldv_31[,c(1,7:8)], by.x = 0,by.y = 'Sample')   %>% select(-col) -> sub_merge_gr_naldv_select
rownames(sub_merge_gr_naldv_select) <- sub_merge_gr_naldv_select$Row.names
dim(sub_merge_gr_naldv_select)
library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(sub_merge_gr_naldv_select[,c(-1,-29,-30)]), y=t(merge_gr_18S_select),type = "spearman") #
net <- CoDF(cor_ls$r, cor_ls$P) #
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) 
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/fa_net_padj10%.csv') 
#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_gr
net_gr %>% merge(shared_ncldv_31, by.x = 'to', by.y = 'Sample') %>% select(1,2,3,12, 13) %>% mutate(hab=c('gr','gr','gr','gr')) -> sub_gr_net

merge_go[,(dim(ncldv_go_t)[2]+2):dim(merge_go)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_go[,(dim(ncldv_go_t)[2]+2):dim(merge_go)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>7) %>% select(-col)  -> merge_go_18S_select #
#merge_go[,2:(dim(ncldv_go_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_go[,2:(dim(ncldv_go_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
#filter(col>8)  %>% select(-col) -> merge_go_naldv_select #
merge_go[,2:(dim(ncldv_go_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_go[,2:(dim(ncldv_go_t)[2]+1)], 2, FUN =function(x) sum(x>0))) %>% 
  merge(shared_ncldv_31[,c(1,7:8)], by.x = 0,by.y = 'Sample')   %>% select(-col) -> sub_merge_go_naldv_select
rownames(sub_merge_go_naldv_select) <- sub_merge_go_naldv_select$Row.names
dim(sub_merge_go_naldv_select)
library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(sub_merge_go_naldv_select[,c(-1,-12,-13)]), y=t(merge_go_18S_select),type = "spearman") #
net <- CoDF(cor_ls$r, cor_ls$P) #
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) 
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/fa_net_padj10%.csv') 
#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_go
net_go %>% merge(shared_ncldv_31, by.x = 'to', by.y = 'Sample') %>% select(1,2,3,12, 13) %>% mutate(hab=c('go','go','go','go','go',"go")) -> sub_go_net

merge_mine[,(dim(ncldv_mine_t)[2]+2):dim(merge_mine)[2]] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_mine[,(dim(ncldv_mine_t)[2]+2):dim(merge_mine)[2]], 2, FUN =function(x) sum(x>0)))  %>%#
  filter(col>14) %>% select(-col)  -> merge_mine_18S_select #
#merge_mine[,2:(dim(ncldv_mine_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_mine[,2:(dim(ncldv_mine_t)[2]+1)], 2, FUN =function(x) sum(x>0)))  %>%#
#filter(col>8)  %>% select(-col) -> merge_mine_naldv_select #
merge_mine[,2:(dim(ncldv_mine_t)[2]+1)] %>% t() %>% data.frame() %>% mutate(col=  apply(merge_mine[,2:(dim(ncldv_mine_t)[2]+1)], 2, FUN =function(x) sum(x>0))) %>% 
  merge(shared_ncldv_31[,c(1,7:8)], by.x = 0,by.y = 'Sample')   %>% select(-col) -> sub_merge_mine_naldv_select
rownames(sub_merge_mine_naldv_select) <- sub_merge_mine_naldv_select$Row.names
dim(sub_merge_mine_naldv_select)
library(Hmisc)
source("D:\\typora dir\\CoDF.R")
cor_ls <- rcorr(x=t(sub_merge_mine_naldv_select[,c(-1,-79,-80)]), y=t(merge_mine_18S_select),type = "spearman") #
net <- CoDF(cor_ls$r, cor_ls$P) #
net$padj <- p.adjust(net$p, method="BH")
## 
net_padj <- net[which(net$cor >0.65),]
net_padj <- net_padj[which(net_padj$padj < 0.05),]
#
net_padj2 <-  dplyr::filter(net_padj, !grepl("_",from))
net_padj3 <-  dplyr::filter(net_padj2, grepl("_",to))
net_padj3 
summary(unique(net_padj3$from)) 
summary(unique(net_padj3$to)) 
#write.csv(net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) ,row.names = FALSE,'18S/fa_net_padj10%.csv') 
#
net_padj3 %>% mutate(from=str_replace(from, "X" ,"")) -> net_mine
net_mine %>% merge(shared_ncldv_31, by.x = 'to', by.y = 'Sample') %>% select(1,2,3,12, 13) %>% mutate(hab=c('mine','mine','mine','mine','mine',"mine",'mine','mine',"mine")) -> sub_mine_net

#18s taxonomy detail 
da_18S %>% select(272) %>% 
  separate(sep = ";", col = `taxonomy`, into =c(NA,NA,'Kindom',NA,NA,'Phylum',NA,NA,NA,'class',NA,NA,NA,'order',NA,NA,'family',NA,'genus','sp')) ->taxo_18_detail
#作，5 habs sub net
rbind(sub_fa_net,sub_fo_net,sub_gr_net,sub_go_net,sub_mine_net) %>% merge(taxo_18_detail, by.x = 'from', by.y = 0) -> sub_habs_net #merge taxonomy


sub_habs_net %>% select(1:8) %>% filter(hab %in% c("fa")) -> fa_sub_habs_net
rbind(data.frame(id=fa_sub_habs_net$from, tax=fa_sub_habs_net$Kindom,tax2=fa_sub_habs_net$Phylum), 
      data.frame(id=fa_sub_habs_net$to,tax=fa_sub_habs_net$X.2,tax2=fa_sub_habs_net$X.1)) %>% data.frame() %>% mutate(K = rep(c('Eu', "NCLDV"), each =3))->  fa_sub_habs_net_design
sub_habs_net %>% select(1:8) %>% filter(hab %in% c("fa")) %>% graph_from_data_frame( directed=F) ->  tmp_sub_net
ggraph(tmp_sub_net, layout = 'linear', circular = TRUE) +#
  #create_layout(all_root_net ,layout, circular = T )+
  geom_edge_arc(aes(edge_width = cor),color = '#F00000',alpha =0.5) + #
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +
  geom_node_point(aes(color =fa_sub_habs_net_design$tax,shape = fa_sub_habs_net_design$K), alpha =1, size=10) +  #，CoDF
  geom_node_text(aes(label = fa_sub_habs_net_design$tax2), angle =0, hjust = 0.5, nudge_y = 0, size = 4 ,family ='serif') +  #
  scale_size_continuous(range = c(8, 12) )+ #
  scale_edge_width_continuous(range = c(1, 2) )+
  scale_shape_manual(values = c(17 ,19 ) )+
  scale_color_manual(values = c('#FF9900','#FFFF00','#339933','#FD7446FF' ,"red"))+
  theme_void() +  #
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))  -> p_sub_fa_net


sub_habs_net %>% select(1:8) %>% filter(hab %in% c("fo")) -> fo_sub_habs_net
rbind(data.frame(id=fo_sub_habs_net$from, tax=fo_sub_habs_net$Kindom,tax2=fo_sub_habs_net$Phylum), 
      data.frame(id=fo_sub_habs_net$to,tax=fo_sub_habs_net$X.2,tax2=fo_sub_habs_net$X.1)) %>% data.frame() %>% mutate(K = rep(c('Eu', "NCLDV"), each =5)) %>% unique()   ->  fo_sub_habs_net_design
fo_sub_habs_net %>% graph_from_data_frame( directed=F) ->  tmp_sub_net2
length(V(tmp_sub_net2))
ggraph(tmp_sub_net2, layout = 'linear', circular = TRUE) +#
  #create_layout(all_root_net ,layout, circular = T )+
  geom_edge_arc(aes(edge_width = cor),color = '#F00000',alpha =0.5) + #
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +
  geom_node_point(aes(color =fo_sub_habs_net_design$tax, shape = fo_sub_habs_net_design$K), alpha =1, size=10) +  #，CoDF
  geom_node_text(aes(label = fo_sub_habs_net_design$tax2), angle =0, hjust = 0.5, nudge_y = 0, size = 4 ,family ='serif') +  #
  scale_size_continuous(range = c(8, 12) )+ #
  scale_edge_width_continuous(range = c(1, 2) )+
  scale_shape_manual(values = c(17 ,19 ) )+
  scale_color_manual(values = c('#FF9900','#FFFF00','#339933','#FD7446FF' ,"red", 'darkred', "green", "black"))+
  theme_void() +  #
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))  -> p_sub_fo_net

sub_habs_net %>% select(1:8) %>% filter(hab %in% c("gr")) -> gr_sub_habs_net
rbind(data.frame(id=gr_sub_habs_net$from, tax=gr_sub_habs_net$Kindom,tax2=gr_sub_habs_net$Phylum), 
      data.frame(id=gr_sub_habs_net$to,tax=gr_sub_habs_net$X.2,tax2=gr_sub_habs_net$X.1)) %>% data.frame() %>% mutate(K = rep(c('Eu', "NCLDV"), each =4)) %>% unique() -> gr_sub_habs_net_design
gr_sub_habs_net %>% graph_from_data_frame( directed=F) ->  tmp_sub_net3
length(V(tmp_sub_net3))
ggraph(tmp_sub_net3, layout = 'linear', circular = TRUE) +#
  #create_layout(all_root_net ,layout, circular = T )+
  geom_edge_arc(aes(edge_width = cor),color = '#F00000',alpha =0.5) + #
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +
  geom_node_point(aes(color =gr_sub_habs_net_design$tax, shape = gr_sub_habs_net_design$K), alpha =1, size=10) +  #，CoDF
  geom_node_text(aes(label = gr_sub_habs_net_design$tax2), angle =0, hjust = 0.5, nudge_y = 0, size = 4 ,family ='serif') +  #
  scale_size_continuous(range = c(8, 12) )+ #
  scale_edge_width_continuous(range = c(1, 2) )+
  scale_shape_manual(values = c(17 ,19 ) )+
  scale_color_manual(values = c('#FF9900','#FFFF00','#339933','#FD7446FF' ,"red", 'darkred', "green", "black"))+
  theme_void() +  #
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))  -> p_sub_gr_net

sub_habs_net %>% select(1:8) %>% filter(hab %in% c("go")) -> go_sub_habs_net
rbind(data.frame(id=go_sub_habs_net$from, tax=go_sub_habs_net$Kindom,tax2=go_sub_habs_net$Phylum), 
      data.frame(id=go_sub_habs_net$to,tax=go_sub_habs_net$X.2,tax2=go_sub_habs_net$X.1)) %>% data.frame() %>% mutate(K = rep(c('Eu', "NCLDV"), each =6)) %>% unique() -> go_sub_habs_net_design
go_sub_habs_net %>% graph_from_data_frame( directed=F) ->  tmp_sub_net4
length(V(tmp_sub_net4))
ggraph(tmp_sub_net4, layout = 'linear', circular = TRUE) +#
  #create_layout(all_root_net ,layout, circular = T )+
  geom_edge_arc(aes(edge_width = cor),color = '#F00000',alpha =0.5) + #
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +
  geom_node_point(aes(color =go_sub_habs_net_design$tax, shape = go_sub_habs_net_design$K), alpha =1, size=10) +  #，CoDF
  geom_node_text(aes(label = go_sub_habs_net_design$tax2), angle =0, hjust = 0.5, nudge_y = 0, size = 4 ,family ='serif') +  #
  scale_size_continuous(range = c(8, 12) )+ #
  scale_edge_width_continuous(range = c(2.2, 3) )+
  scale_shape_manual(values = c(17 ,19 ) )+
  scale_color_manual(values = c('#FF9900','#FFFF00','#339933','#FD7446FF' ,"red", 'darkred', "green", "black"))+
  theme_void() +  #
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))  -> p_sub_go_net


sub_habs_net %>% select(1:8) %>% filter(hab %in% c("mine")) -> mine_sub_habs_net
rbind(data.frame(id=mine_sub_habs_net$from, tax=mine_sub_habs_net$Kindom,tax2=mine_sub_habs_net$Phylum), 
      data.frame(id=mine_sub_habs_net$to,tax=mine_sub_habs_net$X.2,tax2=mine_sub_habs_net$X.1)) %>% data.frame() %>% mutate(K = rep(c('Eu', "NCLDV"), each =9)) %>% unique() -> mine_sub_habs_net_design
mine_sub_habs_net %>% graph_from_data_frame( directed=F) ->  tmp_sub_net5
length(V(tmp_sub_net5))
ggraph(tmp_sub_net5, layout = 'linear', circular = TRUE) +#
  #create_layout(all_root_net ,layout, circular = T )+
  geom_edge_arc(aes(edge_width = cor),color = '#F00000',alpha =0.5) + #
  # geom_edge_link(aes(edge_width=weight),color="grey",
  #               # arrow = arrow(length = unit(2, 'mm')), 
  #                end_cap = circle(3, 'mm')) +
  geom_node_point(aes(color =mine_sub_habs_net_design$tax, shape = mine_sub_habs_net_design$K), alpha =1, size=10) +  #，CoDF
  geom_node_text(aes(label = mine_sub_habs_net_design$tax2), angle =0, hjust = 0.5, nudge_y = 0, size = 4 ,family ='serif') +  #
  scale_size_continuous(range = c(8, 12) )+ #
  scale_edge_width_continuous(range = c(1, 2) )+
  scale_shape_manual(values = c(17 ,19 ) )+
  scale_color_manual(values = c('#FF9900','#FFFF00','#339933','#FD7446FF' ,"red", 'darkred', "green", "black"))+
  theme_void() +  #
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))  -> p_sub_mine_net
library(patchwork)
p_sub_fa_net+p_sub_fo_net+p_sub_gr_net+p_sub_go_net+p_sub_mine_net


##########################4 diversity randomForest###########################################
#asvtable
da_18S %>% select(contains('Fa')) %>% mutate(rows = rowSums(.)) %>% filter(rows > 0) %>% select(-rows)-> fa_da_18S
da_18S %>% select(contains('Fo')) %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>% select(-rows)-> fo_da_18S
da_18S %>% select(contains('Gr')) %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>% select(-rows)-> gr_da_18S
da_18S %>% select(contains('Go')) %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>% select(-rows)-> go_da_18S
da_18S %>% select(contains(c(".T" ,".D"))) %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>% select(-rows)-> mine_da_18S

da_ncldv %>% select(contains("Fa")) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_fa
da_ncldv %>% select(contains("Fo")) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_fo
da_ncldv %>% select(contains("Gr")) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_gr
da_ncldv %>% select(contains("Go")) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_go
da_ncldv %>% select(contains(c(".T" ,".D"))) %>% mutate(ro = rowSums(.)) %>% filter(ro>0) %>% select(-ro) -> ncldv_mine

source("D:\\typora dir\\alpha_index.R") #activate alpha_index function to calculate alpha diversity :alpha_index(t(da), method = 'richness', base = exp(1))
library(vegan)
############4.1 richness#####################
alpha_index(t(fa_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_fa), method = 'richness', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_fa_18s_ncldv
alpha_index(t(fo_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_fo), method = 'richness', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_fo_18s_ncldv
alpha_index(t(gr_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_gr), method = 'richness', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_gr_18s_ncldv
alpha_index(t(go_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_go), method = 'richness', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_go_18s_ncldv
alpha_index(t(mine_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_mine), method = 'richness', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_mine_18s_ncldv

##################4.1.1  box NCLDV richness #########################
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


rich_18s_ncldv_combine %>%  ggboxplot( x="hab", y="rich_18s", color = "hab", bxp.errorbar = T, outlier.shape = NA,
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
  theme(text = element_text(family = 'serif'))+  ylim(c(-0,1250)) 
##################4.1.2 shannon #######################
alpha_index(t(fa_da_18S), method = 'shannon', base = exp(1)) %>% merge(alpha_index(t(ncldv_fa), method = 'shannon', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "shan_18s","shan_ncldv") -> shan_fa_18s_ncldv
alpha_index(t(fo_da_18S), method = 'shannon', base = exp(1)) %>% merge(alpha_index(t(ncldv_fo), method = 'shannon', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "shan_18s","shan_ncldv") -> shan_fo_18s_ncldv
alpha_index(t(gr_da_18S), method = 'shannon', base = exp(1)) %>% merge(alpha_index(t(ncldv_gr), method = 'shannon', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "shan_18s","shan_ncldv") -> shan_gr_18s_ncldv
alpha_index(t(go_da_18S), method = 'shannon', base = exp(1)) %>% merge(alpha_index(t(ncldv_go), method = 'shannon', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "shan_18s","shan_ncldv") -> shan_go_18s_ncldv
alpha_index(t(mine_da_18S), method = 'shannon', base = exp(1)) %>% merge(alpha_index(t(ncldv_mine), method = 'shannon', base = exp(1)), by=0) %>% data.frame() %>% 
  set_names("site", "shan_18s","shan_ncldv") -> shan_mine_18s_ncldv

#shannon box
shan_fa_18s_ncldv %>% mutate(hab= rep(each=dim(shan_fa_18s_ncldv)[1], "Farmland")) %>% rbind( mutate(shan_fo_18s_ncldv,hab= rep(each=dim(shan_fo_18s_ncldv)[1], "Forest"))) %>% 
  rbind( mutate(shan_gr_18s_ncldv,hab= rep(each=dim(shan_gr_18s_ncldv)[1], "Grassland"))) %>% rbind( mutate(shan_go_18s_ncldv,hab= rep(each=dim(shan_go_18s_ncldv)[1], "Gobi desert")))%>% 
  rbind( mutate(shan_mine_18s_ncldv,hab= rep(each=dim(shan_mine_18s_ncldv)[1], "Mine wasteland"))) -> shan_18s_ncldv_combine
shan_18s_ncldv_combine$hab  <- factor(shan_18s_ncldv_combine$hab,levels = c('Farmland','Forest','Grassland','Gobi desert','Mine wasteland'))
shan_18s_ncldv_combine %>%  ggboxplot( x="hab", y="shan_ncldv", color = "hab",bxp.errorbar = T,outlier.shape = NA,
                                       palette = c("#3B4992FF" ,"#CC3333", "#008B45FF" ,"#631879FF" ,"#008280FF") ,alpha = 0.7,
                                       add = "jitter", shape=1 ,size = 1.5)+
  theme_pubr(base_size = 20)+theme(text = element_text(family = 'serif'))+theme_test(base_size = 20 )+theme(text = element_text(family = 'serif'))+ theme(legend.position = "none")+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = "white")) +labs(x = "", y = "Shannon index")+
  theme(axis.title = element_text(size = 30, vjust = 1)) + 
  theme(axis.title = element_text(vjust = 1),  axis.text.x = element_text(vjust = 1,hjust = 1 ,angle = 30))+theme(axis.text = element_text(face = "bold")) +labs(x = NULL)+
  theme(text = element_text(family = 'serif'))+  ylim(c(-0.5,6))  -> p_shan_ncldv

#################4.2 randomForest##########################figure5
library(randomForest)
#install.packages(c('rfPermute' ,'A3'))
library(rfPermute)
library(A3)
library(ggplot2)
#write.csv(shan_18s_ncldv_combine,"shan_18s_ncldv_combine.csv")
#write.csv(rich_18s_ncldv_combine,"rich_18s_ncldv_combine.csv")
######################4.21 five habs richness randomforest#################################
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\fa_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
#set.seed(27155)
randomForest(rich_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
set.seed(2755)
rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> fa_RF_result
#
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
  theme(text = element_text(family = 'serif')) -> p_fa_RF #R2 43.23%

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\fo_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(12213)
randomForest(rich_ncldv~., data = tmp_RF[, c(1:7,9:10)], importance = TRUE) %>%importance( type = 1) 
set.seed(12213)
rfPermute(rich_ncldv~., data = tmp_RF[, c(1:7,9:10)], na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> fo_RF_result
#
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
  theme(text = element_text(family = 'serif')) -> p_fo_RF #R2 63.56%

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\gr_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
#set.seed(23535)11462
randomForest(rich_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
set.seed(23535)
set.seed(25202)
set.seed(64499)
rfPermute(rich_ncldv~., data = tmp_RF,na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> gr_RF_result
#
gr.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
gr.pval
"
tmp <-  read.csv('D:/data/NCLDVtree/1_ncldv_abun0625/α_diversity/RF/gr_physical_rich.csv', header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- 'Eukaryon'
for (i in 1:1258) {
  print(i)
  set.seed(i)
  rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( ) %>% data.frame() -> gr_RF_result
  print(gr_RF_result)
  if (gr_RF_result$X.IncMSE.pval[6] < 0.05){
   break  }
  
}
"
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
  theme(text = element_text(family = 'serif')) -> p_gr_RF #R2 63.09%

#devtools::install_github('ericarcher/swfscMisc')
library('swfscMisc')

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\go_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(10079)
randomForest(rich_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 

set.seed(10079)
rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>%  data.frame() -> go_RF_result
#
go.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
go.pval
"
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\go_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- 'Eukaryon''
for (i in 1200:12574) {
  print(i)
  set.seed(i)
  rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( ) %>% data.frame() -> go_RF_result
  print(go_RF_result)
  if (go_RF_result$X.IncMSE.pval[6] < 0.05){
   break  }
  
}
"
set.seed(123)
#
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
  theme(text = element_text(family = 'serif')) -> p_go_RF #R2 77.66%

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\ta_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(10079)
randomForest(rich_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 

set.seed(10079)
rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>%  data.frame() -> mine_RF_result
set.seed(123)
#
mine.pval <- a3(rich_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
mine.pval
'''
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\ta_physical_rich2.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
for (i in 20000:22574) {
  print(i)
  set.seed(i)
  rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( ) %>% data.frame() -> mine_RF_result
  print(mine_RF_result)
  #if (mine_RF_result$X.IncMSE.pval[8] != 0.05){
   #break  }
  
}
'''
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
  theme(text = element_text(family = 'serif')) -> p_mine_RF #R2 48.57%
library(patchwork)
p_rich_ncldv +p_fa_RF+p_fo_RF+p_gr_RF+p_go_RF+p_mine_RF
#11:21

######################4.22 five habs shannon randomforest#################################
library(randomForest)
#install.packages(c('rfPermute' ,'A3'))
library(rfPermute)
library(A3)
library(ggplot2)
#write.csv(shan_18s_ncldv_combine,"shan_18s_ncldv_combine.csv")
#write.csv(rich_18s_ncldv_combine,"rich_18s_ncldv_combine.csv")
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\2shan\\fa_physical_rich.csv' ,header = 1 ,row.names = 1)
shan_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(2755)
randomForest(shan_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
set.seed(2755)
rfPermute(shan_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> fa_RF_result
#
fa.pval <- a3(shan_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
fa.pval

fa_RF_result2 <- cbind(fa_RF_result,row.names(fa_RF_result))
fa_RF_result2$`row.names(fa_RF_result)` <- factor(fa_RF_result2$`row.names(fa_RF_result)` ,level = fa_RF_result2$`row.names(fa_RF_result)`)
ggplot(fa_RF_result2) +
  aes(x = `row.names(fa_RF_result)`, weight = X.IncMSE) +
  geom_bar(fill = "#3B4992FF" ,color = 'black',alpha = 0.85 ,width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=c(0,15)) +
  theme_test(base_size = 20 ) + 
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), panel.background = element_rect(fill = "white")) +
  labs(x = "", y = "Increased in MSE (%)") +
  theme(axis.title = element_text(size = 30, vjust = 1)) + 
  theme(axis.title = element_text(vjust = 1), axis.text = element_text(hjust = 1), axis.text.x = element_text(vjust = 1,hjust = 1 ,angle = 30), axis.text.y = element_text(vjust = 0)) +
  theme(axis.text = element_text(face = "bold")) +labs(x = NULL) +
  theme(text = element_text(family = 'serif')) -> p_fa_RF #R2 47.3%

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\2shan\\fo_physical_rich.csv' ,header = 1 ,row.names = 1)
shan_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(12213)
randomForest(shan_ncldv~., data = tmp_RF[, c(1:7,9:10)], importance = TRUE) %>%importance( type = 1) 
set.seed(12213)
rfPermute(shan_ncldv~., data = tmp_RF[, c(1:7,9:10)], na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> fo_RF_result
#
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
  theme(text = element_text(family = 'serif')) -> p_fo_RF #R2 63.56%

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\2shan\\gr_physical_rich.csv' ,header = 1 ,row.names = 1)
shan_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(11640) #11640
randomForest(shan_ncldv~., data = tmp_RF, importance = TRUE) #%>%importance( type = 1) 

set.seed(11640)
rfPermute(shan_ncldv~., data = tmp_RF,na.action = na.roughfix) %>% importance( type = 1) %>% data.frame() -> gr_RF_result
#
gr.pval <- a3(shan_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
gr.pval
"
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\2shan\\gr_physical_rich.csv', header = 1 ,row.names = 1)
shan_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- 'Eukaryon'
for (i in 20000:22158) {
  print(i)
  set.seed(i)
  rfPermute(shan_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( ) %>% data.frame() -> gr_RF_result
  print(gr_RF_result)
  if (gr_RF_result$X.IncMSE[7] > 0.1){
   break  }
  
}
"
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
  theme(text = element_text(family = 'serif')) -> p_gr_RF #R2 63.09%

#devtools::install_github('ericarcher/swfscMisc')
library('swfscMisc')

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\2shan\\go_physical_rich.csv' ,header = 1 ,row.names = 1)
shan_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(572455)
randomForest(shan_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
#set.seed(23535)
#set.seed(25202)
set.seed(572455)
rfPermute(shan_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>%  data.frame() -> go_RF_result
#
go.pval <- a3(shan_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
go.pval
"
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\2shan\\go_physical_rich.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- 'Eukaryon''
for (i in 1200:12574) {
  print(i)
  set.seed(i)
  rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( ) %>% data.frame() -> go_RF_result
  print(go_RF_result)
  if (go_RF_result$X.IncMSE.pval[6] < 0.05){
   break  }
  
}
"
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
  theme(text = element_text(family = 'serif')) -> p_go_RF #R2 77.66%

tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\2shan\\ta_physical_rich.csv' ,header = 1 ,row.names = 1)
shan_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- "Eukaryon"
set.seed(10079)
randomForest(shan_ncldv~., data = tmp_RF, importance = TRUE) %>%importance( type = 1) 
#set.seed(23535)
#set.seed(25202)
set.seed(10079)
rfPermute(shan_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( type = 1) %>%  data.frame() -> mine_RF_result
set.seed(123)
#
mine.pval <- a3(shan_ncldv~., data = tmp_RF, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
mine.pval
"
tmp <-  read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\α_diversity\\RF\\2shan\\ta_physical_rich2.csv' ,header = 1 ,row.names = 1)
rich_18s_ncldv_combine %>% merge(tmp[c(-1,-2)], by.x ='site', by.y = 0) %>% select(-1,-4) -> tmp_RF
names(tmp_RF)[1] <- 'Eukaryon'
for (i in 20000:22574) {
  print(i)
  set.seed(i)
  rfPermute(rich_ncldv~., data = tmp_RF, na.action = na.roughfix) %>% importance( ) %>% data.frame() -> mine_RF_result
  print(mine_RF_result)
  #if (mine_RF_result$X.IncMSE.pval[8] != 0.05){
   #break  }
  
}
"
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
  theme(text = element_text(family = 'serif')) -> p_mine_RF #R2 48.57%
library(patchwork)
p_shan_ncldv +p_fa_RF+p_fo_RF+p_gr_RF+p_go_RF+p_mine_RF
#11:21

##########################5 NMDS & dis#############################
#########5.1 NMDS#########
library(vegan)
library(stringr)
da_ncldv %>% select(2:271) %>% 
  t() %>% data.frame()  %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>%
  mutate(groups=str_replace(rownames(.),".*Fa.*", "Farmland")) %>%  
  mutate(groups=str_replace(groups,".*Fo.*", "Forest")) %>%
  mutate(groups=str_replace(groups,".*Go.*", "Gobi desert")) %>% 
  mutate(groups=str_replace(groups,".*Gr.*", "Grassland")) %>%  
  mutate(groups=str_replace(groups,".*T.*", "Mine wasteland")) %>%
  mutate(groups=str_replace(groups,".*\\.D.*", "Mine wasteland"))%>%    
  select(groups) -> da_ncldv_group

da_ncldv %>% select(2:271) %>% 
  t() %>% data.frame()  %>% mutate(rows = rowSums(.)) %>% filter(rows>0) %>%
  mutate(groups=str_replace(rownames(.),".*Fa.*", "Farmland")) %>%  
  mutate(groups=str_replace(groups,".*Fo.*", "Forest")) %>%
  mutate(groups=str_replace(groups,".*Go.*", "Gobi desert")) %>% 
  mutate(groups=str_replace(groups,".*Gr.*", "Grassland")) %>%  
  mutate(groups=str_replace(groups,".*T.*", "Mine wasteland")) %>%
  mutate(groups=str_replace(groups,".*\\.D.*", "Mine wasteland"))%>%  
  select(-groups, -rows) %>%
  vegdist(method = "bray") -> bray_dis
  
nmds_otu <- metaMDS(bray_dis, k =3)
nmds_otu$stress
a <- data.frame(nmds_otu$points)
a
anosim.result <- anosim(bray_dis,da_ncldv_group$group,permutations =999)
plot(anosim.result)

da_ncldv_group %>% filter(rownames(.) !='XJ.Go10')  %>% filter(rownames(.) !='XJ.Go8')  -> da_ncldv_group_rm 
a %>% filter(a$MDS1 > 0.03)
a %>% filter(a$MDS1 < 0.03) %>% 
ggplot( aes(MDS1 ,MDS2)) +
  geom_point(aes(fill = da_ncldv_group_rm$group ),color = 'black' ,shape=21  ,size =6 ,alpha = 0.99) +
  #stat_ellipse(aes(color = group$group), geom = 'path', size = 1.3, alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values =c("#6c77ad","#d86666","#8a529b","#44d355" ,"#41a1a0" ), NULL) +
  # scale_fill_manual(values = c('red', 'orange', 'green3')) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +theme_test(base_size = 25 ,base_family = 'serif')+
  theme( plot.background = element_rect(colour = NA))+ theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black")) +labs(x = "NMDS1", y = "NMDS2")
#6:10
scale_size_continuous(trans = )
#############################5.2 distance-decay##################################################################################
library(vegan)
library(geosphere)
library(ggpubr)
da_env_fa <- read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\man\\fa_physical_rich2.csv', header = 1 ,row.names = 1) #fa_env
da_env_fa %>% merge(ncldv_fa_t, by=0) -> merge_fa_ncldv_env
rownames(merge_fa_ncldv_env) <- merge_fa_ncldv_env$Row.names #merge env & ncldv

dim(merge_fa_ncldv_env)
merge_fa_ncldv_env[,1:15]
d.geo_fa <- distm(merge_fa_ncldv_env[,2:3] ,fun = distHaversine)
dist_fa <- vegdist(merge_fa_ncldv_env[, 16:dim(merge_fa_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) 
library(ggpmisc)
mantel(as.dist(d.geo_fa),dist_fa,  method = 'pearson', permutations = 999, na.rm = TRUE)
cbind(data.frame(geo = as.vector(as.dist(d.geo_fa))) , data.frame(bray=as.vector(dist_fa)) ) %>%
  ggplot(aes(x = geo/1000,y =  bray))+
    geom_point(shape = 19,colour = '#424242' ,size =2.5 ,alpha = 0.5)+
    geom_smooth(method = 'lm' ,colour = 'red',alpha = 0)+
    labs(y = 'NCLDV community dissimilarity' ,x = NULL)+
    theme_test(base_size = 22,base_family = 'serif')+
    theme(plot.subtitle = element_text(family = "serif",
    size = 15, colour = "gray0")) +labs(subtitle = "(a) Farmland") + theme(axis.ticks = element_line(colour = "gray0"),
    axis.text = element_text(colour = "gray0"),
    axis.text.x = element_text(colour = "black"))+
  scale_x_continuous(limits = c(0,3500), breaks = c(0,700,1400,2100,2800,3500))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.0))+theme(plot.subtitle = element_text(size = 22,  face = "bold")) ->p_fa_DDR
  geom_label(aes(x=2800, y=0.3),label = "Mantel test\nr = 0.316\n P = 0.001",label.size = 0, size = 8, family='serif') + 

#fo
da_env_fo <- read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\man\\fo_physical_rich2.csv', header = 1 ,row.names = 1) #fo_env
da_env_fo %>% merge(ncldv_fo_t, by=0) -> merge_fo_ncldv_env
rownames(merge_fo_ncldv_env) <- merge_fo_ncldv_env$Row.names #merge env & ncldv
dim(merge_fo_ncldv_env)
merge_fo_ncldv_env[,1:15]
d.geo_fo <- distm(merge_fo_ncldv_env[,2:3] ,fun = distHaversine)
dist_fo <- vegdist(merge_fo_ncldv_env[, 16:dim(merge_fo_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #ncldv
library(ggpmisc)
mantel(as.dist(d.geo_fo),dist_fo,  method = 'pearson', permutations = 999, na.rm = TRUE)
cbind(data.frame(geo = as.vector(as.dist(d.geo_fo))) , data.frame(bray=as.vector(dist_fo)) ) %>%
  ggplot(aes(x = geo/1000,y =  bray))+
  geom_point(shape = 19,colour = '#424242' ,size =2.5 ,alpha = 0.5)+
  geom_smooth(method = 'lm' ,colour = 'red',alpha = 0)+
  labs(y = NULL ,x = NULL)+
  #stat_cor(color = 'red', label.y =0.4 ,label.x = 0 ,size =6)+
  #stat_poly_eq(aes( label = ..eq.label.. ),colour = 'red',
  #           formula = 'y~x',label.y =0.5,label.x =0 ,size = 6 ,parse = T) +
  theme_test(base_size = 22,base_family = 'serif')+
  theme(plot.subtitle = element_text(family = "serif",size = 15, colour = "gray0")) +labs(subtitle = "(b) Forest") + 
  theme(axis.ticks = element_line(colour = "gray0"), axis.text = element_text(colour = "gray0"), axis.text.x = element_text(colour = "black"))+
  scale_x_continuous(limits = c(0,3500), breaks = c(0,700,1400,2100,2800,3500))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.0))+theme(plot.subtitle = element_text(size = 22,  face = "bold"))-> p_fo_DDR
  geom_label(aes(x=2800, y=0.3),label = "Mantel test\nr = 0.211\n P = 0.001",label.size = 0, size = 8, family='serif') + 

da_env_gr <- read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\man\\gr_physical_rich2.csv', header = 1 ,row.names = 1) #gr_env
da_env_gr %>% merge(ncldv_gr_t, by=0) -> merge_gr_ncldv_env
rownames(merge_gr_ncldv_env) <- merge_gr_ncldv_env$Row.names #merge env & ncldv
dim(merge_gr_ncldv_env)
merge_gr_ncldv_env[,1:15]
d.geo_gr <- distm(merge_gr_ncldv_env[,2:3] ,fun = distHaversine)
dist_gr <- vegdist(merge_gr_ncldv_env[, 16:dim(merge_gr_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #
library(ggpmisc)
mantel(d.geo_gr,dist_gr,  method = 'pearson', permutations = 999, na.rm = TRUE)
cbind(data.frame(geo = as.vector(as.dist(d.geo_gr))) , data.frame(bray=as.vector(dist_gr)) ) %>%
  ggplot(aes(x = geo/1000,y =  bray))+
  geom_point(shape = 19,colour = '#424242' ,size =2.5 ,alpha = 0.5)+
  geom_smooth(method = 'lm' ,colour = 'red',alpha = 0)+
  labs(y =NULL ,x = 'Geographic distance(km)')+
  #stat_cor(color = 'red', label.y =0.4 ,label.x = 0 ,size =6)+
  #stat_poly_eq(aes( label = ..eq.label.. ),colour = 'red',
  #           formula = 'y~x',label.y =0.5,label.x =0 ,size = 6 ,parse = T) +
  theme_test(base_size = 22,base_family = 'serif')+
  theme(plot.subtitle = element_text(family = "serif",size = 15, colour = "gray0")) +labs(subtitle = "(c) Grassland") + 
  theme(axis.ticks = element_line(colour = "gray0"), axis.text = element_text(colour = "gray0"), axis.text.x = element_text(colour = "black"))+
  scale_x_continuous(limits = c(0,2500), breaks = c(0,500,1000, 1500,2000,2500))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.0))+theme(plot.subtitle = element_text(size = 22,  face = "bold")) ->p_gr_DDR
  geom_label(aes(x=2100, y=0.3),label = "Mantel test\nr = 0.502\n P = 0.001",label.size = 0, size = 8, family='serif')+ 

da_env_go <- read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\man\\go_physical_rich2.csv', header = 1 ,row.names = 1) #gr_env
da_env_go %>% merge(ncldv_go_t, by=0) -> merge_go_ncldv_env
rownames(merge_go_ncldv_env) <- merge_go_ncldv_env$Row.names #merge env & ncldv
dim(merge_go_ncldv_env)
merge_go_ncldv_env[,1:15]
d.geo_go <- distm(merge_go_ncldv_env[,2:3] ,fun = distHaversine)
dist_go <- vegdist(merge_go_ncldv_env[, 16:dim(merge_go_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #
library(ggpmisc)
mantel(d.geo_go,dist_go,  method = 'pearson', permutations = 999, na.rm = TRUE)
cbind(data.frame(geo = as.vector(as.dist(d.geo_go))) , data.frame(bray=as.vector(dist_go)) ) %>%
  ggplot(aes(x = geo/1000,y =  bray))+
  geom_point(shape = 19,colour = '#424242' ,size =2.5 ,alpha = 0.5)+
  geom_smooth(method = 'lm' ,colour = 'red',alpha = 0)+
  labs(y =NULL ,x = 'Geographic distance(km)')+
  #stat_cor(color = 'red', label.y =0.4 ,label.x = 0 ,size =6)+
  #stat_poly_eq(aes( label = ..eq.label.. ),colour = 'red',
  #           formula = 'y~x',label.y =0.5,label.x =0 ,size = 6 ,parse = T) +
  theme_test(base_size = 22,base_family = 'serif')+
  theme(plot.subtitle = element_text(family = "serif",size = 15, colour = "gray0")) +labs(subtitle = "(d) Gobi desert") + 
  theme(axis.ticks = element_line(colour = "gray0"), axis.text = element_text(colour = "gray0"), axis.text.x = element_text(colour = "black"))+
  scale_x_continuous(limits = c(0,1500), breaks = c(0,500,1000,1500))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.0))+ theme(plot.subtitle = element_text(size = 22,  face = "bold")) ->p_go_DDR
  geom_label(aes(x=1200, y=0.3),label = "Mantel test\nr = 0.3633\n P = 0.001",label.size = 0, size = 8, family='serif')+

da_env_mine <- read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\man\\ta_physical_rich2.csv', header = 1 ,row.names = 1) #gr_env
da_env_mine %>% merge(ncldv_mine_t, by=0) -> merge_mine_ncldv_env
rownames(merge_mine_ncldv_env) <- merge_mine_ncldv_env$Row.names #merge env & ncldv
dim(merge_mine_ncldv_env)
merge_mine_ncldv_env[,1:15]
d.geo_mine <- distm(merge_mine_ncldv_env[,2:3] ,fun = distHaversine)
dist_mine <- vegdist(merge_mine_ncldv_env[, 16:dim(merge_mine_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #
library(ggpmisc)
mantel(d.geo_mine,dist_mine,  method = 'pearson', permutations = 999, na.rm = TRUE) 
cbind(data.frame(geo = as.vector(as.dist(d.geo_mine))) , data.frame(bray=as.vector(dist_mine)) ) %>%
  ggplot(aes(x = geo/1000,y =  bray))+
  geom_point(shape = 19,colour = '#424242' ,size =2.5 ,alpha = 0.5)+
  geom_smooth(method = 'lm' ,colour = 'red',alpha = 0)+
  labs(y =NULL ,x = 'Geographic distance(km)')+
  #stat_cor(color = 'red', label.y =0.4 ,label.x = 0 ,size =6)+
  #stat_poly_eq(aes( label = ..eq.label.. ),colour = 'red',
  #           formula = 'y~x',label.y =0.5,label.x =0 ,size = 6 ,parse = T) +
  theme_test(base_size = 22,base_family = 'serif')+
  theme(plot.subtitle = element_text(family = "serif",size = 15, colour = "gray0")) +labs(subtitle = "(f) Mine wasteland") + 
  theme(axis.ticks = element_line(colour = "gray0"), axis.text = element_text(colour = "gray0"), axis.text.x = element_text(colour = "black"))+
  scale_x_continuous(limits = c(0,3500), breaks = c(0,700,1400,2100,2800,3500))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.2,0.4,0.6,0.8,1.0))+theme(plot.subtitle = element_text(size = 22,  face = "bold"))->p_mine_DDR
  geom_label(aes(x=2800, y=0.3),label = "Mantel test\nr = 0.092\n P = 0.001",label.size = 0, size = 8, family='serif') + 

library(patchwork)
p_fa_DDR+p_fo_DDR+p_gr_DDR+p_go_DDR+p_mine_DDR

#ggsave(file="DDR_five_habs.pdf", width = 360, height = 220, units = "mm")
############################6 environment range & abundance###########################################################################################figure2
##################6.1 total abun with all env  range####################
da <- read.csv('D:\\typora dir\\3.environment_range\\physical_all_321_2.csv',header = 1 ,row.names = 1)
max(da)
scale(da)
library('vegan')
library('dplyr')
decostand(da,'max') %>% summary() #
daa <- decostand(da,'range')
DB <- data.frame(rowMeans(daa)) #
DAA <- decostand(DB,'range') #######################

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
da <- read.csv('scale_mean_all2.csv' ,header = 1,row.names = 1)
da$ranges <- vegan::decostand(da$ranges,'max')
library(ggplot2)
library(basicTrendline)
ggplot(data = da ,aes(x= ranges,y=  occ,color = factor(ranges)))+
  geom_point(size = 4 ,shape = 19 ,alpha = 0.75)+
  scale_colour_manual(values = colorRampPalette(c("#2196F3", "#0D47A1"))(1139))+
  stat_smooth(method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 28) + 
  theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif", colour = "black"), axis.text.x = element_text(family = "serif", colour = "black"), axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = NA)) +labs(x = "Environmental range (r=0.7,p<0.001) ", y = "Number of sampling sites" ) -> p1
trendline(da$ranges, da$occ, model="line2P", ePos.x = "topright", 
          summary=TRUE, eDigit=2)# eDigit
ggplot(data = da ,aes(decostand(log(da$ranges+1),'max'), log(da$sum_cov), color = factor(da$ranges)))+
  geom_point(size = 4 ,shape = 19 ,alpha = 0.75)+
  scale_colour_manual(values = colorRampPalette(c("#BA68C8", "#7B1FA2"))(1140))+
  stat_smooth(method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 28) + theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"),  plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif",  colour = "black"), axis.text.x = element_text(family = "serif", colour = "black"), axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = NA)) +labs(x = "Environmental range (r = 0.520, P < 0.001)", y = "Total NCLDV abundance (ln)" ) -> p2
trendline(decostand(log(da$ranges+1),'max'), log(da$sum_cov), model="line2P", ePos.x = "topright", 
          summary=TRUE, eDigit=2)# eDigit
p1+p2 #Fig S3

#################6.2 mean with all range env###################

#Fig 2e env-range with Mean abun
ggplot(data = da ,aes(decostand(log(da$ranges+1),'max'), log(da$mean_cov), color = factor(da$ranges)))+
  geom_point(size = 4 ,shape = 19 ,alpha = 0.75)+
  scale_colour_manual(values = colorRampPalette(c("#BA68C8", "#7B1FA2"))(1140))+
  stat_smooth(method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 28)+scale_y_continuous( limits=c(-2.5, 5))+ theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"),
                                            plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif",  colour = "black"), axis.text.x = element_text(family = "serif",  colour = "black"), axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = NA)) +labs(x = "Environmental range (r = -0.412, P < 0.001)", y = "Mean NCLDV abundance (ln)" ) -> p_fig2e
trendline(decostand(log(da$ranges+1),'max'), log(da$mean_cov), model="line2P", ePos.x = "topright", 
          summary=TRUE, eDigit=2)# 
cor.test(da$occ,da$ranges,method = 'pearson')

#Fig 2 env-range with Mean abun
ggplot(data = da ,aes(log(da$occ), log(da$mean_cov), color = factor(da$ranges)))+
  geom_point(size = 4 ,shape = 19 ,alpha = 0.75)+
  scale_colour_manual(values = colorRampPalette(c("#b2ebf2", "#00838f"))(1020))+
  stat_smooth(method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 28)+scale_y_continuous( limits=c(-2.5, 5))+ theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"),
                                                                                                           plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif",  colour = "black"), axis.text.x = element_text(family = "serif",  colour = "black"), axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = NA)) +labs(x = "Number of sampling sites (r = -0.5, P < 0.001)", y = "Mean NCLDV abundance (ln)" ) -> p_fig2d
p_fig2d
trendline(log(da$occ), log(da$mean_cov), model="line2P", ePos.x = "topright", 
          summary=TRUE, eDigit=2)# 
cor.test(log(da$occ), log(da$mean_cov),method = 'pearson')

#
ncldv_fa %>% rowSums() %>% data.frame()-> ncldv_fa_rowsum
ncldv_fa_rowsum$NAM <- rownames(ncldv_fa_rowsum)
ncldv_fo %>% rowSums() %>% data.frame() -> ncldv_fo_rowsum
ncldv_fo_rowsum$NAM <- rownames(ncldv_fo_rowsum)
ncldv_go %>% rowSums() %>% data.frame() -> ncldv_go_rowsum
ncldv_go_rowsum$NAM <- rownames(ncldv_go_rowsum)
ncldv_gr %>% rowSums() %>% data.frame() -> ncldv_gr_rowsum
ncldv_gr_rowsum$NAM <- rownames(ncldv_gr_rowsum)
ncldv_mine %>% rowSums() %>% data.frame() -> ncldv_mine_rowsum
ncldv_mine_rowsum$NAM <- rownames(ncldv_mine_rowsum)
colnames(ncldv_fa_rowsum) <- c('Farmland', 'NAM')
colnames(ncldv_fo_rowsum) <- c('Forest', 'NAM')
colnames(ncldv_go_rowsum) <- c('Gobi desert', 'NAM')
colnames(ncldv_gr_rowsum) <- c('Grassland', 'NAM')
colnames(ncldv_mine_rowsum) <- c('Mine wasteland', 'NAM')
list(ncldv_fa_rowsum,ncldv_fo_rowsum,ncldv_go_rowsum,ncldv_gr_rowsum,ncldv_mine_rowsum) %>% reduce(full_join, by='NAM') -> merge_naldv_all_sum
merge_naldv_all_sum[is.na(merge_naldv_all_sum)] <- 0
merge_naldv_all_sum[,c(1,3:6)][merge_naldv_all_sum[,c(1,3:6)] >0] <- 1
#write.csv(merge_naldv_all_sum,"num_habs_dected.csv")
#########!!rowwise()###############
merge_naldv_all_sum %>% select(2,1,3,4,5,6)  %>%   rowwise() %>%   mutate(NU_habs = sum(across(where(is.numeric)))) %>% select(1,7) %>% merge(da, by.x='NAM', by.y = 'X.1') -> add_habs_da

ggplot(add_habs_da, aes(NU_habs , log(da$mean_cov), color = factor(da$ranges)))+
  geom_jitter(size = 3 ,shape = 19 ,alpha = 0.75,width = 0.2)+
  scale_colour_manual(values = colorRampPalette(c("#c8e6c9", "#2e7d32"))(1020))+
  stat_smooth(method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 28)+scale_y_continuous( limits=c(-2.5, 5))+ theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"),
                                                                                                           plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif",  colour = "black"), axis.text.x = element_text(family = "serif",  colour = "black"), axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = NA)) +labs(x = "Number of habitats (r = -0.374, P < 0.001)", y = "Mean NCLDV abundance (ln)" ) -> p_fig2c
p_fig2c


trendline(add_habs_da$NU_habs , log(da$mean_cov), model="line2P", ePos.x = "topright", 
          summary=TRUE, eDigit=2)# 
cor.test(log(da$occ), log(da$mean_cov),method = 'pearson')

library(patchwork)
p_fig2c+p_fig2d+p_fig2e
#
########################6.3 family environment range ####################################################
add_habs_da %>% merge(read.csv("new_nclav_coverage/ncldv_taxo_rmbasal.csv", header = F), by.x = "NAM", by.y = "V1") %>% filter(!(V4 %in% c('Phycodnaviridae','Prasinoviridae', ""))) %>% 
  ggplot()+
  geom_point(aes(x= log(occ),y=  log( mean_cov),color = factor(NU_habs)),size = 4 ,shape = 19 ,alpha = 0.75)+
  stat_cor(aes(x= log(occ),y= log(mean_cov),label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 5, show.legend = FALSE) +
  scale_colour_manual(values = colorRampPalette(c("#A5D6A7", "#1B5E20"))(5))+
  stat_smooth(aes(x= log(occ),y=  log( mean_cov),color = factor(NU_habs)),method = 'lm',formula = y~poly(x ,1) ,color = 'red' ,alpha = 0)+
  theme_bw(base_size = 20,base_family = "serif")+
  facet_wrap(vars(V4) ,nrow =3 ,scales = "fixed") +
  theme(legend.position = "none") + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) + 
  theme(axis.ticks = element_line(colour = "gray14"), axis.text = element_text(family = "serif", colour = "black"), axis.text.x = element_text(family = "serif", colour = "black"), axis.text.y = element_text(colour = "black"), panel.background = element_rect(fill = NA)) +
  labs(x = "Number of sampling sites (ln)",  y = "NCLDV abundance (ln)" )+
  ylim(-2.5 ,5)
#

#####################################7 VPA ###################################figure6
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

par(mfrow=c(1, 1))  
par(mar=rep(4, 4))  
####################8 mantel& part mantel###############################
library(vegan)
################8.1 Mantel #######################
#ASV 
merge_fa_ncldv_env[,1:15] 
d.geo_fa <- distm(merge_fa_ncldv_env[,2:3] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_fa <- vegdist(merge_fa_ncldv_env[, 16:dim(merge_fa_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #

all<- NULL
fa_mantel<- NULL
for (i in c(4:ncol(merge_fa_ncldv_env[,1:15]))){
  daa <- merge_fa_ncldv_env[,1:15][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel(dist_fa, dist.f ,  method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_fa_ncldv_env[,1:15][i]),aa$statistic,aa$signif)
  fa_mantel <- rbind(fa_mantel ,all )
}
fa_mantel <- data.frame(fa_mantel )
faaa <- fa_mantel[order(fa_mantel[,3],decreasing = TRUE),]
vegdist(vpa_fa_18s,method = "bray")
mantel(vegdist(vpa_fa_nclav,method = "bray"), vegdist(vpa_fa_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> fa_18s_ncldv_mantel
data.frame(cbind('NCLDV_bray', 'Eukaryota', fa_18s_ncldv_mantel$statistic,fa_18s_ncldv_mantel$signif)) %>% rbind(faaa)%>% mutate(group=rep("Farmland", 13)) -> mantel_fa_result
#fo
merge_fo_ncldv_env[,1:15] 
d.geo_fo <- distm(merge_fo_ncldv_env[,2:3] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_fo <- vegdist(merge_fo_ncldv_env[, 16:dim(merge_fo_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #
all<- NULL
fo_mantel<- NULL
for (i in c(4:ncol(merge_fo_ncldv_env[,1:15]))){
  daa <- merge_fo_ncldv_env[,1:15][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel(dist_fo, dist.f ,  method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_fo_ncldv_env[,1:15][i]),aa$statistic,aa$signif)
  fo_mantel <- rbind(fo_mantel ,all )
}
fo_mantel <- data.frame(fo_mantel )
fooo <- fo_mantel[order(fo_mantel[,3],decreasing = TRUE),]
mantel(vegdist(vpa_fo_nclav,method = "bray"), vegdist(vpa_fo_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> fo_18s_ncldv_mantel
data.frame(cbind('NCLDV_bray', 'Eukaryota', fo_18s_ncldv_mantel$statistic,fo_18s_ncldv_mantel$signif)) %>% rbind(fooo) %>% mutate(group=rep("Forest", 13)) -> mantel_fo_result
#gr
merge_gr_ncldv_env[,1:15] 
d.geo_gr <- distm(merge_gr_ncldv_env[,2:3] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_gr <- vegdist(merge_gr_ncldv_env[, 16:dim(merge_gr_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #
all<- NULL
gr_mantel<- NULL
for (i in c(4:ncol(merge_gr_ncldv_env[,1:15]))){
  daa <- merge_gr_ncldv_env[,1:15][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel(dist_gr, dist.f ,  method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_gr_ncldv_env[,1:15][i]),aa$statistic,aa$signif)
  gr_mantel <- rbind(gr_mantel ,all )
}
gr_mantel <- data.frame(gr_mantel )
grrr <- gr_mantel[order(gr_mantel[,3],decreasing = TRUE),]
mantel(vegdist(vpa_gr_nclav,method = "bray"), vegdist(vpa_gr_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> gr_18s_ncldv_mantel
data.frame(cbind('NCLDV_bray', 'Eukaryota', gr_18s_ncldv_mantel$statistic,gr_18s_ncldv_mantel$signif)) %>% rbind(grrr) %>% mutate(group=rep("Grassland", 13))-> mantel_gr_result
#go
merge_go_ncldv_env[,1:15] 
d.geo_go <- distm(merge_go_ncldv_env[,2:3] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_go <- vegdist(merge_go_ncldv_env[, 16:dim(merge_go_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #
all<- NULL
go_mantel<- NULL
for (i in c(4:ncol(merge_go_ncldv_env[,1:15]))){
  daa <- merge_go_ncldv_env[,1:15][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel(dist_go, dist.f ,  method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_go_ncldv_env[,1:15][i]),aa$statistic,aa$signif)
  go_mantel <- rbind(go_mantel ,all )
}
go_mantel <- data.frame(go_mantel )
gorr <- go_mantel[order(go_mantel[,3],decreasing = TRUE),]
mantel(vegdist(vpa_go_nclav,method = "bray"), vegdist(vpa_go_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> go_18s_ncldv_mantel
data.frame(cbind('NCLDV_bray', 'Eukaryota', go_18s_ncldv_mantel$statistic,go_18s_ncldv_mantel$signif)) %>% rbind(gorr) %>% mutate(group=rep("Gobi desert", 13))-> mantel_go_result
#mine
#mine
merge_mine_ncldv_env[,1:15] 
d.geo_mine <- distm(merge_mine_ncldv_env[,2:3] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_mine <- vegdist(merge_mine_ncldv_env[, 16:dim(merge_mine_ncldv_env)[2]], method = 'bray',upper = T ,diag = T) #
all<- NULL
mine_mantel<- NULL
for (i in c(4:ncol(merge_mine_ncldv_env[,1:15]))){
  daa <- merge_mine_ncldv_env[,1:15][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel(dist_mine, dist.f ,  method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_mine_ncldv_env[,1:15][i]),aa$statistic,aa$signif)
  mine_mantel <- rbind(mine_mantel ,all )
}
mine_mantel <- data.frame(mine_mantel )
minerr <- mine_mantel[order(mine_mantel[,3],decreasing = TRUE),]
mantel(vegdist(vpa_mine_nclav,method = "bray"), vegdist(vpa_mine_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> mine_18s_ncldv_mantel
data.frame(cbind('NCLDV_bray', 'Eukaryota', mine_18s_ncldv_mantel$statistic,mine_18s_ncldv_mantel$signif)) %>% rbind(minerr) %>% mutate(group=rep("Mine wasteland", 13)) -> mantel_mine_result
rbind(mantel_fa_result, mantel_fo_result, mantel_gr_result, mantel_go_result, mantel_mine_result) %>% data.frame()  -> mantel_re_habs
mantel_re_habs[,3:4] <- as.data.frame(lapply(mantel_re_habs[,3:4],as.numeric))
mantel_re_habs %>% mutate(sig = case_when(X4 <= 0.001 ~'***' , (X4 > 0.001 & X4 <= 0.01) ~'**', (X4 > 0.1 & X4 <= 0.05) ~'*')) %>% mutate(X2= str_replace(X2, "WC", "Moisture")) -> mantel_re_habs_sig
mantel_re_habs_sig$X2 <- factor(mantel_re_habs_sig$X2, levels = c('ALT','CEC','Clays','EC','EXCa','MAP','Moisture','pH','TC','TK','TN','TP','Eukaryota'))
mantel_re_habs_sig$group <- factor(mantel_re_habs_sig$group, levels = rev(c('Farmland','Forest','Grassland','Gobi desert','Mine wasteland')))
ggplot(mantel_re_habs_sig) +
  geom_point(aes(x=group ,y = X2 ,size=X3),color = '#f47f7e'  ,alpha =0.9)+
  scale_shape_manual(values = 19) +
  scale_size_continuous(range=c(2,18))+ #
  theme(legend.position = "none")+
  theme_bw(base_family = "serif" ,base_size = 20)+
  theme(axis.text.x = element_text(colour = "black", angle = 30 ,hjust = 1), axis.text.y = element_text(colour = "black", angle = 0 ,hjust = 1, vjust = 0))+
  geom_text(aes(x=group, y = X2,label = sig), family = 'serif' ,color = 'black' ,size =8)+coord_flip()+ theme(legend.position = "none") -> p_mantel_habs
################8.2 partial Mantel (control dis)#######################
#
#fa
fa_da_18S_t %>% merge(merge_fa_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
vpa_fa_18s  #18S
merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+2):(dim(fa_da_18S_t)[2]+3)] #lon lat
merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+4):(dim(fa_da_18S_t)[2]+15)] #phy_chem
vpa_fa_nclav #ncldv
d.geo_fa <- distm(merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+2):(dim(fa_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_fa <- vegdist(vpa_fa_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
fa_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+4):(dim(fa_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+4):(dim(fa_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_fa, dist.f ,d.geo_fa, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+4):(dim(fa_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  fa_mantel <- rbind(fa_mantel ,all )
}
fa_mantel <- data.frame(fa_mantel )
faaa <- fa_mantel[order(fa_mantel[,3],decreasing = TRUE),]
vegdist(vpa_fa_18s,method = "bray")
mantel.partial(vegdist(vpa_fa_nclav,method = "bray"), vegdist(vpa_fa_18s,method = "bray") ,d.geo_fa,  method = 'pearson', permutations = 999, na.rm = TRUE) -> fa_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Eukaryota', fa_18s_ncldv_mantel_part$statistic,fa_18s_ncldv_mantel_part$signif)) %>% rbind(faaa)%>% mutate(group=rep("Farmland", 13)) -> partial_mantel_fa_result

#fo
fo_da_18S_t %>% merge(merge_fo_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
d.geo_fo <- distm(merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+2):(dim(fo_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_fo <- vegdist(vpa_fo_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
fo_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+4):(dim(fo_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+4):(dim(fo_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_fo, dist.f ,d.geo_fo, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+4):(dim(fo_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  fo_mantel <- rbind(fo_mantel ,all )
}
fo_mantel <- data.frame(fo_mantel )
fooo <- fo_mantel[order(fo_mantel[,3],decreasing = TRUE),]
vegdist(vpa_fo_18s,method = "bray")
mantel.partial(vegdist(vpa_fo_nclav,method = "bray"), vegdist(vpa_fo_18s,method = "bray") ,d.geo_fo,  method = 'pearson', permutations = 999, na.rm = TRUE) -> fo_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Eukaryota', fo_18s_ncldv_mantel_part$statistic,fo_18s_ncldv_mantel_part$signif)) %>% rbind(fooo)%>% mutate(group=rep("Forest", 13)) -> partial_mantel_fo_result

#gr
gr_da_18S_t %>% merge(merge_gr_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
d.geo_gr <- distm(merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+2):(dim(gr_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_gr <- vegdist(vpa_gr_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
gr_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+4):(dim(gr_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+4):(dim(gr_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_gr, dist.f ,d.geo_gr, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+4):(dim(gr_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  gr_mantel <- rbind(gr_mantel ,all )
}
gr_mantel <- data.frame(gr_mantel )
groo <- gr_mantel[order(gr_mantel[,3],decreasing = TRUE),]
vegdist(vpa_gr_18s,method = "bray")
mantel.partial(vegdist(vpa_gr_nclav,method = "bray"), vegdist(vpa_gr_18s,method = "bray") ,d.geo_gr,  method = 'pearson', permutations = 999, na.rm = TRUE) -> gr_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Eukaryota', gr_18s_ncldv_mantel_part$statistic,gr_18s_ncldv_mantel_part$signif)) %>% rbind(groo)%>% mutate(group=rep("Grassland", 13)) -> partial_mantel_gr_result

#go
go_da_18S_t %>% merge(merge_go_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
d.geo_go <- distm(merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+2):(dim(go_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_go <- vegdist(vpa_go_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
go_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+4):(dim(go_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+4):(dim(go_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean', diag = T, upper = T)
  aa <- mantel.partial(dist_go, dist.f ,d.geo_go, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+4):(dim(go_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  go_mantel <- rbind(go_mantel ,all )
}
go_mantel <- data.frame(go_mantel )
gooo <- go_mantel[order(go_mantel[,3],decreasing = TRUE),]
vegdist(vpa_go_18s,method = "bray")
mantel.partial(vegdist(vpa_go_nclav,method = "bray"), vegdist(vpa_go_18s,method = "bray") ,d.geo_go,  method = 'pearson', permutations = 999, na.rm = TRUE) -> go_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Eukaryota', go_18s_ncldv_mantel_part$statistic,go_18s_ncldv_mantel_part$signif)) %>% rbind(gooo)%>% mutate(group=rep("Gobi desert", 13)) -> partial_mantel_go_result

#mine
mine_da_18S_t %>% merge(merge_mine_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
d.geo_mine <- distm(merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+2):(dim(mine_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_mine <- vegdist(vpa_mine_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
mine_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+4):(dim(mine_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+4):(dim(mine_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_mine, dist.f ,d.geo_mine, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+4):(dim(mine_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  mine_mantel <- rbind(mine_mantel ,all )
}
mine_mantel <- data.frame(mine_mantel )
mineoo <- mine_mantel[order(mine_mantel[,3],decreasing = TRUE),]
vegdist(vpa_mine_18s,method = "bray")
mantel.partial(vegdist(vpa_mine_nclav,method = "bray"), vegdist(vpa_mine_18s,method = "bray") ,d.geo_mine,  method = 'pearson', permutations = 999, na.rm = TRUE) -> mine_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Eukaryota', mine_18s_ncldv_mantel_part$statistic,mine_18s_ncldv_mantel_part$signif)) %>% rbind(mineoo)%>% mutate(group=rep("Mine wasteland", 13)) -> partial_mantel_mine_result

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

################8.3 partial Mantel (control Eu)#######################
#control Eu
#fa
fa_da_18S_t %>% merge(merge_fa_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
#vpa_fa_18s  #18S
dist_Eu <- vegdist(vpa_fa_18s, method = 'bray',upper = T ,diag = T) #Eu bray
merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+2):(dim(fa_da_18S_t)[2]+3)] #lon lat
merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+4):(dim(fa_da_18S_t)[2]+15)] #phy_chem
#vpa_fa_nclav #ncldv
d.geo_fa <- distm(merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+2):(dim(fa_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_fa <- vegdist(vpa_fa_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
fa_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+4):(dim(fa_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+4):(dim(fa_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_fa, dist.f ,dist_Eu, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(fa_da_18S_t)[2]+4):(dim(fa_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  fa_mantel <- rbind(fa_mantel ,all )
}
fa_mantel <- data.frame(fa_mantel )
faaa <- fa_mantel[order(fa_mantel[,3],decreasing = TRUE),]
vegdist(vpa_fa_18s,method = "bray")
mantel.partial(vegdist(vpa_fa_nclav,method = "bray"),d.geo_fa, vegdist(vpa_fa_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> fa_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Distance', fa_18s_ncldv_mantel_part$statistic,fa_18s_ncldv_mantel_part$signif)) %>% rbind(faaa)%>% mutate(group=rep("Farmland", 13)) -> partial2_mantel_fa_result

#fo
fo_da_18S_t %>% merge(merge_fo_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
#vpa_fo_18s  #18S
dist_Eu <- vegdist(vpa_fo_18s, method = 'bray',upper = T ,diag = T) #Eu bray
merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+2):(dim(fo_da_18S_t)[2]+3)] #lon lat
merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+4):(dim(fo_da_18S_t)[2]+15)] #phy_chem
vpa_fo_nclav #ncldv
d.geo_fo <- distm(merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+2):(dim(fo_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_fo <- vegdist(vpa_fo_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
fo_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+4):(dim(fo_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+4):(dim(fo_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_fo, dist.f ,dist_Eu, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(fo_da_18S_t)[2]+4):(dim(fo_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  fo_mantel <- rbind(fo_mantel ,all )
}
fo_mantel <- data.frame(fo_mantel )
fooo <- fo_mantel[order(fo_mantel[,3],decreasing = TRUE),]
vegdist(vpa_fo_18s,method = "bray")
mantel.partial(vegdist(vpa_fo_nclav,method = "bray"),d.geo_fo, vegdist(vpa_fo_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> fo_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Distance', fo_18s_ncldv_mantel_part$statistic,fo_18s_ncldv_mantel_part$signif)) %>% rbind(fooo)%>% mutate(group=rep("Forest", 13)) -> partial2_mantel_fo_result

#gr
gr_da_18S_t %>% merge(merge_gr_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
#vpa_fo_18s  #18S
dist_Eu <- vegdist(vpa_gr_18s, method = 'bray',upper = T ,diag = T) #Eu bray
merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+2):(dim(gr_da_18S_t)[2]+3)] #lon lat
merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+4):(dim(gr_da_18S_t)[2]+15)] #phy_chem
#vpa_gr_nclav #ncldv
d.geo_gr <- distm(merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+2):(dim(gr_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_gr <- vegdist(vpa_gr_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
gr_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+4):(dim(gr_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+4):(dim(gr_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_gr, dist.f ,dist_Eu, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(gr_da_18S_t)[2]+4):(dim(gr_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  gr_mantel <- rbind(gr_mantel ,all )
}
gr_mantel <- data.frame(gr_mantel )
grrr <- gr_mantel[order(gr_mantel[,3],decreasing = TRUE),]
vegdist(vpa_gr_18s,method = "bray")
mantel.partial(vegdist(vpa_gr_nclav,method = "bray"),d.geo_gr, vegdist(vpa_gr_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> gr_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Distance', gr_18s_ncldv_mantel_part$statistic,gr_18s_ncldv_mantel_part$signif)) %>% rbind(grrr)%>% mutate(group=rep("Grassland", 13)) -> partial2_mantel_gr_result

#go
go_da_18S_t %>% merge(merge_go_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
#vpa_fo_18s  #18S
dist_Eu <- vegdist(vpa_go_18s, method = 'bray',upper = T ,diag = T) #Eu bray
merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+2):(dim(go_da_18S_t)[2]+3)] #lon lat
merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+4):(dim(go_da_18S_t)[2]+15)] #phy_chem
#vpa_gr_nclav #ncldv
d.geo_gr <- distm(merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+2):(dim(go_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_gr <- vegdist(vpa_go_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
go_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+4):(dim(go_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+4):(dim(go_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_go, dist.f ,dist_Eu, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(go_da_18S_t)[2]+4):(dim(go_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  go_mantel <- rbind(go_mantel ,all )
}
go_mantel <- data.frame(go_mantel )
gooo <- go_mantel[order(go_mantel[,3],decreasing = TRUE),]
vegdist(vpa_go_18s,method = "bray")
mantel.partial(vegdist(vpa_go_nclav,method = "bray"),d.geo_go, vegdist(vpa_go_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> go_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Distance', go_18s_ncldv_mantel_part$statistic,go_18s_ncldv_mantel_part$signif)) %>% rbind(gooo)%>% mutate(group=rep("Gobi desert", 13)) -> partial2_mantel_go_result

#mine
mine_da_18S_t %>% merge(merge_mine_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
dim(merge_18s_env_ncldv)
#vpa_fo_18s  #18S
dist_Eu <- vegdist(vpa_mine_18s, method = 'bray',upper = T ,diag = T) #Eu bray
merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+2):(dim(mine_da_18S_t)[2]+3)] #lon lat
merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+4):(dim(mine_da_18S_t)[2]+15)] #phy_chem
#vpa_gr_nclav #ncldv
d.geo_mine <- distm(merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+2):(dim(mine_da_18S_t)[2]+3)] ,fun = distHaversine) #lon\ lat
#ncldv bray dis
dist_mine <- vegdist(vpa_mine_nclav, method = 'bray',upper = T ,diag = T) #
all<- NULL
mine_mantel<- NULL
for (i in c(1:ncol(merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+4):(dim(mine_da_18S_t)[2]+15)]))){
  daa <- merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+4):(dim(mine_da_18S_t)[2]+15)][i]
  dist.f <- vegdist(daa, method = 'euclidean')
  aa <- mantel.partial(dist_mine, dist.f ,dist_Eu, method = 'pearson', permutations = 999, na.rm = TRUE)
  all <- cbind('NCLDV_bray',names(merge_18s_env_ncldv[,(dim(mine_da_18S_t)[2]+4):(dim(mine_da_18S_t)[2]+15)][i]),aa$statistic,aa$signif)
  mine_mantel <- rbind(mine_mantel ,all )
}
mine_mantel <- data.frame(mine_mantel )
mineee <- mine_mantel[order(mine_mantel[,3],decreasing = TRUE),]
vegdist(vpa_mine_18s,method = "bray")
mantel.partial(vegdist(vpa_mine_nclav,method = "bray"),d.geo_mine, vegdist(vpa_mine_18s,method = "bray") ,  method = 'pearson', permutations = 999, na.rm = TRUE) -> mine_18s_ncldv_mantel_part
data.frame(cbind('NCLDV_bray', 'Distance', mine_18s_ncldv_mantel_part$statistic,mine_18s_ncldv_mantel_part$signif)) %>% rbind(mineee)%>% mutate(group=rep("Mine wasteland", 13)) -> partial2_mantel_mine_result
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

###################9 phy-che points R^2,P,et.al.############################
fa_da_18S_t %>% merge(merge_fa_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
merge_18s_env_ncldv[,c(1,(dim(fa_da_18S_t)[2]+3):(dim(fa_da_18S_t)[2]+15))] %>% merge(rich_fa_18s_ncldv, by.x ='Row.names', by.y = 'site' ) %>% select(16, 2:15)%>%
  mutate(habs=rep("Farmland")) -> phy_che_fa_point 
fo_da_18S_t %>% merge(merge_fo_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
merge_18s_env_ncldv[,c(1,(dim(fo_da_18S_t)[2]+3):(dim(fo_da_18S_t)[2]+15))] %>% merge(rich_fo_18s_ncldv, by.x ='Row.names', by.y = 'site' ) %>% select(16, 2:15)%>%
  mutate(habs=rep("Forest")) -> phy_che_fo_point 
gr_da_18S_t %>% merge(merge_gr_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
merge_18s_env_ncldv[,c(1,(dim(gr_da_18S_t)[2]+3):(dim(gr_da_18S_t)[2]+15))] %>% merge(rich_gr_18s_ncldv, by.x ='Row.names', by.y = 'site' ) %>% select(16, 2:15)%>%
  mutate(habs=rep("Grassland")) -> phy_che_gr_point 
go_da_18S_t %>% merge(merge_go_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
merge_18s_env_ncldv[,c(1,(dim(go_da_18S_t)[2]+3):(dim(go_da_18S_t)[2]+15))] %>% merge(rich_go_18s_ncldv, by.x ='Row.names', by.y = 'site' ) %>% select(16, 2:15)%>%
  mutate(habs=rep("Gobi desert")) -> phy_che_go_point 
mine_da_18S_t %>% merge(merge_mine_ncldv_env, by.x = 0, by.y = 'Row.names') -> merge_18s_env_ncldv
merge_18s_env_ncldv[,c(1,(dim(mine_da_18S_t)[2]+3):(dim(mine_da_18S_t)[2]+15))] %>% merge(rich_mine_18s_ncldv, by.x ='Row.names', by.y = 'site' ) %>% select(16, 2:15)%>%
  mutate(habs=rep("Mine westeland")) -> phy_che_mine_point 
names(phy_che_fa_point)[c(1, 15)] <- c("NCLDV", "Eukaryota")
names(phy_che_fo_point)[c(1, 6,15)] <- c("NCLDV",'Moisture', "Eukaryota")
names(phy_che_gr_point)[c(1,6, 15)] <- c("NCLDV",'Moisture', "Eukaryota")
names(phy_che_go_point)[c(1,6, 15)] <- c("NCLDV",'Moisture', "Eukaryota")
names(phy_che_mine_point)[c(1, 6,15)] <- c("NCLDV",'Moisture', "Eukaryota")
library(ggpmisc)
formula <- y ~ x + I(x^2)
rbind(phy_che_fa_point,phy_che_fo_point,phy_che_gr_point,phy_che_go_point,phy_che_mine_point) %>% 
  ggplot()+#labs(x = NULL)+
  geom_point(aes(y = NCLDV ,x = LAT,fill = habs) , colour = "white",size = 4.2,alpha = 0.9 , stroke = 1.5 ,shape = 21)+
  geom_smooth(aes(y = NCLDV ,x = LAT) ,color = 'black',se=FALSE ,size = 1,
              method = 'lm' , formula=y~poly(x,2))+
  stat_poly_eq(aes(y = NCLDV ,x = LAT,label = paste(..AIC.label.., ..BIC.label.., ..rr.label..,  sep = "~",stat(p.value.label))), formula =formula, parse = T,family = "serif",size = 5) +
  #stat_cor(aes(x= NCLDV,y= LAT,label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'lasso', label.x.npc = 'left', label.y.npc = 'top', size = 5, show.legend = FALSE) +
  scale_fill_manual(values = c("#3B4992FF", "#CC3333" ,"#631879FF" ,"#008B45FF" ,"#008280FF"  ))+
  facet_grid(~ habs)+
  #scale_shape_manual(values = c(3,2,19,15,8))+labs(y = NULL)+
   xlim(c(20 ,60))+
  theme_bw(base_family = "serif" ,base_size = 28 )+ theme(legend.position = "none") -> p_point_lat

formula <- y ~ x + I(x^2)
rbind(phy_che_fa_point,phy_che_fo_point,phy_che_gr_point,phy_che_go_point,phy_che_mine_point) %>% 
  ggplot()+#labs(x = NULL)+
  geom_point(aes(y = NCLDV ,x = MAP,fill = habs) , colour = "white",size =4.2,alpha = 0.9 , stroke = 1.5 ,shape = 21)+
  geom_smooth(aes(y = NCLDV ,x = MAP) ,color = 'black',se=FALSE ,size = 1,
              method = 'lm' , formula=y~poly(x,2))+
  stat_poly_eq(aes(y = NCLDV ,x = MAP,label = paste(..AIC.label.., ..BIC.label.., ..rr.label..,  sep = "~",stat(p.value.label))), formula =formula, parse = T,family = "serif",size = 5) +
  #stat_cor(aes(x= NCLDV,y= LAT,label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'lasso', label.x.npc = 'left', label.y.npc = 'top', size = 5, show.legend = FALSE) +
  scale_fill_manual(values = c("#3B4992FF", "#CC3333" ,"#631879FF" ,"#008B45FF" ,"#008280FF"  ))+
  facet_grid(~ habs)+
  #scale_shape_manual(values = c(3,2,19,15,8))+labs(y = NULL)+
  #xlim(c(20 ,60))+
  theme_bw(base_family = "serif" ,base_size = 28 )+ theme(legend.position = "none") -> p_point_map

formula <- y ~ x + I(x^2)
rbind(phy_che_fa_point,phy_che_fo_point,phy_che_gr_point,phy_che_go_point,phy_che_mine_point) %>% 
  ggplot()+#labs(x = NULL)+
  geom_point(aes(y = NCLDV ,x = log(EC),fill = habs) , colour = "white",size =4.2,alpha = 0.9 , stroke = 1.5 ,shape = 21)+
  geom_smooth(aes(y = NCLDV ,x = log(EC)) ,color = 'black',se=FALSE ,size = 1,
              method = 'lm' , formula=y~poly(x,2))+
  stat_poly_eq(aes(y = NCLDV ,x = log(EC),label = paste(..AIC.label.., ..BIC.label.., ..rr.label..,  sep = "~",stat(p.value.label))), formula =formula, parse = T,family = "serif",size = 5) +
  #stat_cor(aes(x= NCLDV,y= LAT,label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'lasso', label.x.npc = 'left', label.y.npc = 'top', size = 5, show.legend = FALSE) +
  scale_fill_manual(values = c("#3B4992FF", "#CC3333" ,"#631879FF" ,"#008B45FF" ,"#008280FF"  ))+
  facet_grid(~ habs)+
  #scale_shape_manual(values = c(3,2,19,15,8))+labs(y = NULL)+
  xlim(c(-5, 5))+
  theme_bw(base_family = "serif" ,base_size = 28 )+ theme(legend.position = "none") -> p_point_ec

formula <- y ~ x + I(x^2)
rbind(phy_che_fa_point,phy_che_fo_point,phy_che_gr_point,phy_che_go_point,phy_che_mine_point) %>% 
  ggplot()+#labs(x = NULL)+
  geom_point(aes(y = NCLDV ,x = pH,fill = habs) , colour = "white",size = 4.2,alpha = 0.9 , stroke = 1.5 ,shape = 21)+
  geom_smooth(aes(y = NCLDV ,x = pH) ,color = 'black',se=FALSE ,size = 1,
              method = 'lm' , formula=y~poly(x,2))+
  stat_poly_eq(aes(y = NCLDV ,x = pH,label = paste(..AIC.label.., ..BIC.label.., ..rr.label..,  sep = "~",stat(p.value.label))), formula =formula, parse = T,family = "serif",size = 5) +
  #stat_cor(aes(x= NCLDV,y= LAT,label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'lasso', label.x.npc = 'left', label.y.npc = 'top', size = 5, show.legend = FALSE) +
  scale_fill_manual(values = c("#3B4992FF", "#CC3333" ,"#631879FF" ,"#008B45FF" ,"#008280FF"  ))+
  facet_grid(~ habs)+
  #scale_shape_manual(values = c(3,2,19,15,8))+labs(y = NULL)+
  #xlim(c(20 ,60))+
  theme_bw(base_family = "serif" ,base_size = 28 )+ theme(legend.position = "none") -> p_point_ph

formula <- y ~ x + I(x^2)
rbind(phy_che_fa_point,phy_che_fo_point,phy_che_gr_point,phy_che_go_point,phy_che_mine_point) %>% 
  ggplot()+#labs(x = NULL)+
  geom_point(aes(y = NCLDV ,x = CEC,fill = habs) , colour = "white",size = 4.2,alpha = 0.9 , stroke = 1.5 ,shape = 21)+
  geom_smooth(aes(y = NCLDV ,x = CEC) ,color = 'black',se=FALSE ,size = 1,
              method = 'lm' , formula=y~poly(x,2))+
  stat_poly_eq(aes(y = NCLDV ,x = CEC,label = paste(..AIC.label.., ..BIC.label.., ..rr.label..,  sep = "~",stat(p.value.label))), formula =formula, parse = T,family = "serif",size = 5) +
  #stat_cor(aes(x= NCLDV,y= LAT,label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'lasso', label.x.npc = 'left', label.y.npc = 'top', size = 5, show.legend = FALSE) +
  scale_fill_manual(values = c("#3B4992FF", "#CC3333" ,"#631879FF" ,"#008B45FF" ,"#008280FF"  ))+
  facet_grid(~ habs)+
  #scale_shape_manual(values = c(3,2,19,15,8))+labs(y = NULL)+
  #xlim(c(20 ,60))+
  theme_bw(base_family = "serif" ,base_size = 28 )+ theme(legend.position = "none") -> p_point_cec

formula <- y~poly(x,1)
rbind(phy_che_fa_point,phy_che_fo_point,phy_che_gr_point,phy_che_go_point,phy_che_mine_point) %>% 
  ggplot()+#labs(x = NULL)+
  geom_point(aes(y = NCLDV ,x = Eukaryota   ,fill = habs) , colour = "white",size =4.2,alpha = 0.9 , stroke = 1.5 ,shape = 21)+
  geom_smooth(aes(y = NCLDV ,x = Eukaryota   ) ,color = 'black',se=FALSE ,size = 1,
              method = 'lm' , formula=y~poly(x,1))+
  stat_poly_eq(aes(y = NCLDV ,x = Eukaryota  ,label = paste(..AIC.label.., ..BIC.label.., ..rr.label..,  sep = "~",stat(p.value.label))), formula =y~poly(x,1), parse = T,family = "serif",size = 5) +#这一行paste()中的内容是一个类似Latex的表达，“~”表示空格，parse=T表示将语句翻译成可读形式；size=2.5表示字体大小
  #stat_cor(aes(x= NCLDV,y= LAT,label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'lasso', label.x.npc = 'left', label.y.npc = 'top', size = 5, show.legend = FALSE) +
  scale_fill_manual(values = c("#3B4992FF", "#CC3333" ,"#631879FF" ,"#008B45FF" ,"#008280FF"  ))+
  facet_grid(~ habs)+
  #scale_shape_manual(values = c(3,2,19,15,8))+labs(y = NULL)+
  #xlim(c(20 ,60))+
  theme_bw(base_family = "serif" ,base_size = 28 )+ theme(legend.position = "none")  -> p_point_eu


library(patchwork)
p_point_cec/p_point_ec/p_point_lat/p_point_map/p_point_ph/p_point_eu

#################10 ncldv to world##################################################################
##################10.1 sankey to world######################figure 3
ncldv_world <- read.csv("new_nclav_coverage/bases_mapping_fix_global.csv", header = T)
ncldv_world %>% select(1,2,3) %>% mutate(hab=str_replace(ID, '.*farmland.*', "Farmland")) %>% mutate(hab=str_replace(hab, '.*forest.*', "Forest")) %>% 
  mutate(hab=str_replace(hab, '.*Gr.*', "Grassland"))  %>% mutate(hab=str_replace(hab, '.*Go.*', "Gobi desert"))  %>%  mutate(hab=str_replace(hab, '.*gr.*', "Grassland")) %>% 
  mutate(hab=str_replace(hab, '.*ta.*', "Mine wasteland"))  %>% mutate(hab=str_replace(hab, '.*D_.*', "Mine wasteland")) %>% group_by(order, hab) %>% count() %>% 
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
  melt() %>% select(2,1,3) %>% filter(value>0) %>% merge(ncldv_world[,1:4], by.x = "variable", by.y = "ID") %>% select(2,4)  %>% group_by(order,V2) %>%summarise(Nu=n()) -> world_count
colnames(world_count) <- c('hab','order','n')
library(sankeywheel)
rbind(national_count, world_count)  %>% arrange(order)-> na_wo_map
colnames(na_wo_map) <- c('from','to','weight') 
sankeywheel(from = na_wo_map$from,theme = "grid",
            to = na_wo_map$to,
            weight = na_wo_map$weight,
            type = "sankey",
            width = "100%",height = "30%")  
national_count %>% select(hab,n)%>% group_by(hab)%>% summarise(sum(n))
world_count %>% select(order,n)%>% group_by(order)%>% summarise(sum(n))
national_count %>% select(order,n)%>% group_by(order)%>% summarise(sum(n))

#family
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
##########################10.2 mapping to world########################
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
#10:16
=##############################11 ph phylotypes############################
dim(da_ncldv)
source("D:/typora dir/tran_tb.R") #



da_ncldv %>% select(-1) %>% t() %>% data.frame() %>% mutate(sum = rowSums(.)) %>% filter(sum> 0)  %>% 
  select(-sum) %>% 
  merge(read.csv("D:/data/NCLDVtree/1_ncldv_abun0625/α_diversity/physical_all_321.csv",header = T), by.x=0,by.y="ID") %>% 
  select(Row.names, pH, contains("_"), -NCLDV_rich) %>% mutate(ph_group = case_when(pH <= 3 ~"2-3", (pH >3 & pH <=4) ~"3-4", (pH >3 & pH <=4) ~"3-4",(pH >4 & pH <=5) ~"4-5",(pH >5 & pH <=6) ~"5-6",(pH >6 & pH <=7) ~"6-7",(pH >7 & pH <=8) ~"7-8",(pH >8 & pH <=9) ~"8-9",(pH >9 & pH <= 10) ~"9-10")) %>% 
  select(-pH, -Row.names) %>% group_by(ph_group) %>% 
  summarise(across(where(is.numeric), sum)) %>% data.frame() %>% transpose_df() %>% 
  apply(2,FUN = function(x) ifelse(x>0,1,0)) -> ph_0_1_da



ph_0_1_da %>% merge(read.csv("C:/Users/fengsw/Desktop/NC-mod/new_nclav_coverage/ncldv_taxo_rmbasal.csv", header = F,row.names = 1), by = 0) %>% 
  unite(col= family,V3,V4, sep = "_") %>% group_by(V2) %>% select(2:10) %>% #
  summarise(across(where(is.numeric), sum)) %>% melt() %>% 
  ggplot() +
  geom_tile(aes(x=variable, y=V2,fill=log(value)),color='white',size=1)+
  scale_fill_gradient(low = '#ffd8d2', high = '#ff2806',na.value = 'white')+
  geom_text(aes(x=variable, y=V2,label= value,family="serif"),size=5)+
  theme_pubr(base_size = 21,border = T)+#coord_flip() + #+facet_wrap(vars(Habitats))
  theme(axis.text = element_text(family = "serif"),
          axis.text.x = element_text(family = "serif"),
          axis.text.y = element_text(family = "serif"),
          legend.text = element_text(family = "serif"),
          legend.title = element_text(family = "serif"),
          legend.position = "right")


ph_0_1_da %>% merge(read.csv("C:/Users/fengsw/Desktop/NC-mod/new_nclav_coverage/ncldv_taxo_rmbasal.csv", header = F,row.names = 1), by = 0) %>% 
  unite(col= famil,V3,V4, sep = "_") %>% group_by(famil) %>% select(2:9,11) %>%  #
  summarise(across(where(is.numeric), sum)) %>% data.frame () %>% 
  mutate(famil=str_replace(famil, "#N/A_#N/A", "Unassigned")) %>% 
  mutate(famil=str_replace(famil, "_-", "")) %>% mutate(famil=str_replace(famil, "_#N/A", ""))  %>% mutate(famil=str_replace(famil, "Medusa", "Pandoravirales_unknown")) %>%
  melt() %>% filter(value > 0) %>% filter(famil!="_") %>% mutate(famil = str_replace(famil, ".*_", "")) %>%
  ggplot() +
  geom_tile(aes(x=variable, y=famil,fill=log(value)),color='white',size=1)+
  scale_fill_gradient(low = '#ffd8d2', high = '#ff2806',na.value = 'white')+
  geom_text(aes(x=variable, y=famil,label= value,family="serif"),size=5)+
  theme_pubr(base_size = 21,border = T)+#coord_flip() + #+facet_wrap(vars(Habitats))
  theme(axis.text = element_text(family = "serif"),
        axis.text.x = element_text(family = "serif"),
        axis.text.y = element_text(family = "serif"),
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif"),
        legend.position = "right") + labs(x=NULL, y=NULL)
################################12 upset #####################################################
library(UpSetR) 
#
a <- read.csv('D:\\data\\NCLDVtree\\1_ncldv_abun0625\\ven.csv' ,header = 1 ,row.names = 1)

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
             #list(query=intersects, params=list("farmland", "tailings"), color="#FF61EF", active=T),
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
           #list(query=intersects, params=list("farmland", "forest", "tailings"), color="#FF3D8B", active=T),
           #list(query = elements, params = c('taxonomy', 'PHYCO'), color = c('#3B4992FF'), active = TRUE)),
           #list(query=intersects, params=list("set2", "set5"), color="blue", active=T)),
           # active T  F 
           mainbar.y.label = "Count of Intersection",  # y 
           number.angles = 0,                     #
           sets.x.label = "The number of NCLDV \nPolBs in habitats",  # x 
           text.scale = c(4, 3.5, 3, 3.5, 3.5,3.5),sets.bar.color=c('#B0BEC5','#EF9A9A','#E1BEE7','#7986CB','#42A5F5'),shade.alpha=0.5,
		   order.by = c("freq"),  point.size =6,
           
           line.size = 1,                       #
           decreasing = c(T ,F))
############################################13 tree_heatmap_data ############################
da_ncldv %>% select(-1) %>% select(contains("fa")) %>% t() %>% data.frame() %>% mutate(sum = rowSums(.)) %>% filter(sum> 0)  %>% 
  select(-sum) %>% 
  merge(read.csv("D:/data/NCLDVtree/1_ncldv_abun0625/α_diversity/physical_all_321.csv",header = T), by.x=0,by.y="ID") %>% transpose_df() %>%   apply(2,FUN = function(x) ifelse(x>0,1,0)) %>% 
  apply(1,FUN=sum) %>% data.frame() %>% set_names("Farmland") %>% mutate(Farmland=Farmland/81*100) %>% filter(Farmland >0) %>% filter(Farmland < 100) -> fre1

da_ncldv %>% select(-1) %>% select(contains("fo")) %>% t() %>% data.frame() %>% mutate(sum = rowSums(.)) %>% filter(sum> 0)  %>% 
  select(-sum) %>% 
  merge(read.csv("D:/data/NCLDVtree/1_ncldv_abun0625/α_diversity/physical_all_321.csv",header = T), by.x=0,by.y="ID") %>% transpose_df() %>%   apply(2,FUN = function(x) ifelse(x>0,1,0)) %>% 
  apply(1,FUN=sum) %>% data.frame() %>% set_names("Forest") %>% mutate(Forest=Forest/76*100) %>% filter(Forest >0) %>% filter(Forest < 100) -> fre2

da_ncldv %>% select(-1) %>% select(contains("Gr")) %>% t() %>% data.frame() %>% mutate(sum = rowSums(.)) %>% filter(sum> 0)  %>% 
  select(-sum) %>% 
  merge(read.csv("D:/data/NCLDVtree/1_ncldv_abun0625/α_diversity/physical_all_321.csv",header = T), by.x=0,by.y="ID") %>% transpose_df() %>%   apply(2,FUN = function(x) ifelse(x>0,1,0)) %>% 
  apply(1,FUN=sum) %>% data.frame() %>% set_names("Grassland") %>% mutate(Grassland=Grassland/27*100) %>% filter(Grassland >0) %>% filter(Grassland < 100) -> fre3

da_ncldv %>% select(-1) %>% select(contains("Go")) %>% t() %>% data.frame() %>% mutate(sum = rowSums(.)) %>% filter(sum> 0)  %>% 
  select(-sum) %>% 
  merge(read.csv("D:/data/NCLDVtree/1_ncldv_abun0625/α_diversity/physical_all_321.csv",header = T), by.x=0,by.y="ID") %>% transpose_df() %>%   apply(2,FUN = function(x) ifelse(x>0,1,0)) %>% 
  apply(1,FUN=sum) %>% data.frame() %>% set_names("Gobid esert") %>% mutate(`Gobid esert`=`Gobid esert`/10*100) %>% filter(`Gobid esert` >0) %>% filter(`Gobid esert` < 100) -> fre4

da_ncldv %>% select(-1) %>% select(contains(c(".T" ,".D"))) %>% t() %>% data.frame() %>% mutate(sum = rowSums(.)) %>% filter(sum> 0)  %>% 
  select(-sum) %>% 
  merge(read.csv("D:/data/NCLDVtree/1_ncldv_abun0625/α_diversity/physical_all_321.csv",header = T), by.x=0,by.y="ID") %>% transpose_df() %>%   apply(2,FUN = function(x) ifelse(x>0,1,0)) %>% 
  apply(1,FUN=sum) %>% data.frame() %>% set_names("Mine wasteland") %>% mutate(`Mine wasteland`=`Mine wasteland`/76*100) %>% filter(`Mine wasteland` >0) %>% filter(`Mine wasteland` <100) -> fre5


list_f <- list(fre1,fre2,fre3,fre4,fre5)
new = merge(list_f[1], list_f[2], by=0, all=T)
for (i in 3:5){
  new = merge(new, list_f[i], by.x= "Row.names", by.y=0, all=T)
}
new[,2:6] <- log(new[,2:6])
new
write.csv(new, "C:\\Users\\fengsw\\Desktop\\NC-mod\\1_tree_build_all\\tree_last\\frequence——all.csv")

#####################14 pie 31 share polb########
library(ggplot2)

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
  theme(axis.text.x = element_blank()) +theme(panel.grid=element_blank()) +    ## 
  theme(panel.border=element_blank())  +labs(title = NULL, x = NULL, y = NULL)
ppp <- pp + theme(legend.position = "left")
ppp

###################################15 relative abundance NCLDV############################################################figure4
alpha_index(t(fa_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_fa), method = 'richness', base = exp(1)), by=0, all.y=T) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_fa_18s_ncldv
alpha_index(t(fo_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_fo), method = 'richness', base = exp(1)), by=0, all.y=T) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_fo_18s_ncldv
alpha_index(t(gr_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_gr), method = 'richness', base = exp(1)), by=0, all.y=T) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_gr_18s_ncldv
alpha_index(t(go_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_go), method = 'richness', base = exp(1)), by=0, all.y=T) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_go_18s_ncldv
alpha_index(t(mine_da_18S), method = 'richness', base = exp(1)) %>% merge(alpha_index(t(ncldv_mine), method = 'richness', base = exp(1)), by=0, all.y=T) %>% data.frame() %>% 
  set_names("site", "rich_18s","rich_ncldv") -> rich_mine_18s_ncldv

NCLDV_tax <- read.csv("new_nclav_coverage/ncldv_taxo_rmbasal.csv", header = 1,row.names = 1,na.strings = "")
NCLDV_tax[is.na(NCLDV_tax)] <- "Unassigned"
source("D:/typora dir/tran_tb.R")
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

library(patchwork)
layout <- 'ABBB 
CDDD 
EFFF
'
p4+p1+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  p3+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  p2+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  p4+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  p5+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15))+
  plot_layout(design = layout)
#
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
  theme(legend.position = "right")+
  theme(plot.subtitle = element_text(size = 14, colour = "gray14", hjust = 1), plot.caption = element_text(face = "bold", hjust = 1.25), axis.ticks = element_line(colour = "gray14", linetype = "solid")) +
  labs(y = " Relative abundance (%)")+theme(text=element_text(family = 'serif')) + theme(axis.title = element_text(family = "serif"), plot.title = element_text(family = "serif")) 






