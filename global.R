library(ggplot2)
library(plotly)
library(data.table)
library(shiny)
library(plotly)

base <- 'Data/'

express_vals <- readRDS(paste0(base,'express_vals.rds'))
umap <- readRDS(paste0(base,'umap.rds'))
diffmap <- readRDS(paste0(base,'diffmap.rds'))
gene.names<-colnames(express_vals)

is_df<-readRDS(paste0(base,'is_df.rds'))
tr_tablef<-readRDS(paste0(base,'tr_tablef.rds'))
### LOAD data ###

plot.std.col <- function(df, xname, yname, title){
  
  ggplot(df)+
    geom_point(size=2,aes(x=x, y=y, color=z))+
    
    xlab(xname)+
    ylab(yname)+
    ggtitle(title)+
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          plot.margin=unit(c(1,1,1.5,1.2),"cm"),
          legend.text=element_text(size=12),#size of legend
          legend.title=element_text(size=7),
        
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(lineheight=.8, face="bold",size=10)) +
    scale_fill_gradient(low = "#032927",high ="#f7fdb1",space='Lab', guide = "colourbar", aesthetics = "colour") +
    labs(colour = "natural.log")
  
}
plot_2gene<-function(gene.1,gene.2){
  data.1<-unlist(express_vals[[gene.1]])
  data.2<-unlist(express_vals[[gene.2]])
  data.1.norm<-(data.1-min(data.1))/(max(data.1)-min(data.1))
  data.2.norm<-(data.2-min(data.2))/(max(data.2)-min(data.2))
  col.g1<-c(103,39,112)/255#c(30/255,144/255,255/255)
  col.g2<-c(col2rgb("darkorange"))/255   #c(105,154,51)/255#c(1,0,0)
  n<-c(1,1,1)
  f<-c(0,0.961,0.961) #cyan
  A=col.g1[1]-n[1]
  B=col.g2[1]-n[1]
  C=n[1]
  D=col.g1[2]-n[2]
  E=col.g2[2]-n[2]
  Ff=n[2]
  G=col.g1[3]-n[3]
  H=col.g2[3]-n[3]
  I=n[3]
  alpha<-f[1]-A-B-C
  beta<-f[2]-D-E-Ff
  gamma<-f[3]-G-H-I
  data.norm<-cbind(data.1.norm, data.2.norm)
  z1<-apply(data.norm, 1, function(x) max(c(0,A*x[1]+B*x[2]+C+alpha*x[1]*x[2])))
  z2<-apply(data.norm, 1, function(x) max(c(0,D*x[1]+E*x[2]+Ff+beta*x[1]*x[2])))
  z3<-apply(data.norm, 1, function(x) max(c(0,G*x[1]+H*x[2]+I+gamma*x[1]*x[2])))
  df.3<-data.frame(x=umap$X0, y=umap$X1, 
                   z1=z1,#0.5*(1+data.1.norm),
                   z2=z2,#0.5-max(c(0.5*(1+data.1.norm), 0.5*(1+data.2.norm))),#rep(0,length(data.1.norm))
                   z3=z3)
  
  p3<-ggplot()+
    xlab("")+
    ylab("")+
    theme(plot.title = element_text(size=2, face="bold"))+
    geom_point(data=df.3, aes(x=x,y=y),size=5.1,pch=21 ,fill="black")+
    geom_point(data=df.3, aes(x=x, y=y,color=rgb(z1,z2, z3))   ,size=5)+
    #  scale_color_gradient(low="gray50", high="red")+
    #theme_bw() + 
    theme(panel.background = element_rect(fill = 'grey80', colour = 'white'),
          panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_color_identity()+
    theme(
      plot.margin=unit(c(1,1,1.5,1.2),"cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.text=element_text(size=25),#size of legend
      legend.title=element_text(size=25))
  p3
  
}
plot_2genecolor<-function(gene.1,gene.2){
  data.1<-unlist(express_vals[[gene.1]])
  data.2<-unlist(express_vals[[gene.2]])
  data.1.norm<-(data.1-min(data.1))/(max(data.1)-min(data.1))
  data.2.norm<-(data.2-min(data.2))/(max(data.2)-min(data.2))
  col.g1<-c(103,39,112)/255#c(30/255,144/255,255/255)
  col.g2<-c(col2rgb("darkorange"))/255   #c(105,154,51)/255#c(1,0,0)
  n<-c(1,1,1)
  f<-c(0,0.961,0.961) #cyan
  A=col.g1[1]-n[1]
  B=col.g2[1]-n[1]
  C=n[1]
  D=col.g1[2]-n[2]
  E=col.g2[2]-n[2]
  Ff=n[2]
  G=col.g1[3]-n[3]
  H=col.g2[3]-n[3]
  I=n[3]
  alpha<-f[1]-A-B-C
  beta<-f[2]-D-E-Ff
  gamma<-f[3]-G-H-I
  data.norm<-cbind(data.1.norm, data.2.norm)
  z1<-apply(data.norm, 1, function(x) max(c(0,A*x[1]+B*x[2]+C+alpha*x[1]*x[2])))
  z2<-apply(data.norm, 1, function(x) max(c(0,D*x[1]+E*x[2]+Ff+beta*x[1]*x[2])))
  z3<-apply(data.norm, 1, function(x) max(c(0,G*x[1]+H*x[2]+I+gamma*x[1]*x[2])))
  df.3<-data.frame(x=umap$X0, y=umap$X1, 
                   z1=z1,#0.5*(1+data.1.norm),
                   z2=z2,#0.5-max(c(0.5*(1+data.1.norm), 0.5*(1+data.2.norm))),#rep(0,length(data.1.norm))
                   z3=z3)
  x<-seq(0,1,0.1)
  y<-seq(0,1,0.1)
  m<-expand.grid(x,y)
  res<-t(apply(m, 1, function(x) c(max(c(0,A*x[1]+B*x[2]+C+alpha*x[1]*x[2])), 
                                   max(c(0,D*x[1]+E*x[2]+Ff+beta*x[1]*x[2])),
                                   max(c(0,G*x[1]+H*x[2]+I+gamma*x[1]*x[2]))))
  )
  id<-which(res>1,arr.ind=T)
  res[id]<-1
  df<-data.frame(m=m,col=res)
  
  p4<-ggplot(df, aes(x=m.Var1,y=m.Var2,fill=rgb(col.1,col.2,col.3)))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
    xlab(gene.1)+
    ylab(gene.2)+
    geom_raster()+
    scale_fill_identity()
  p4
  
  }
 
plot_cluster<-function(){
  df_scatter<-data.frame(x=umap$X0,
                         y=umap$X1,
                         cluster=umap$cluster_id)
  ggplot(df_scatter, aes(x=x, y=y,  color=cluster)) +
    geom_point() +
    theme_classic()+scale_colour_manual(values=c('#4E9C6D',"#84584E",'#58BDCC','#B6BD6E','#EF8635','#3C76AF',
                                                 '#C53933', '#9D48F3','#D57DBF','#B3C6E5','#f4aeae'))+
    xlab("")+
    ylab("")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          legend.text=element_text(size=12),
          axis.ticks.y=element_blank())+
    ggtitle("Cluster annotation")+
    guides(color = guide_legend(override.aes = list(size=5)))
    
  
}

boxplot_cluster<-function(df_boxp,gene){
  
  ggplot(df_boxp, aes(x=cluster, y=z,color=cluster)) + 
    geom_boxplot()+scale_colour_manual(values=c('#4E9C6D',"#84584E",'#58BDCC','#B6BD6E','#EF8635','#3C76AF',
                                                 '#C53933', '#9D48F3','#D57DBF','#B3C6E5','#f4aeae'))+ 
    xlab("")+
    ylab("Nat. log. expr")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(legend.text=element_text(size=12),axis.text.x  = element_text( size=12))
  
}


plot_3d<-function(){
  df_diffmap<-data.frame(x=diffmap$X0,
                         y=diffmap$X1, z = diffmap$X2,
                         cluster=umap$cluster_id)
  fig <- plot_ly(df_diffmap, x = ~x, y = ~y, z = ~z, color = ~cluster, colors = c('#4E9C6D',"#84584E",'#58BDCC','#B6BD6E','#EF8635','#3C76AF',
                                                                                  '#C53933', '#9D48F3','#D57DBF','#B3C6E5','#f4aeae'))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'DC1'),
                                     yaxis = list(title = 'DC2'),
                                     zaxis = list(title = 'DC3')))
  
  fig
  
}

plot_iso<-function(g) {
  selectedd<-as.vector(tr_tablef[which(tr_tablef$Gene.name==g),]$Transcript.stable.ID)
  select_gn<-as.vector(tr_tablef[which(tr_tablef$Gene.name==g),]$Transcript.name)
  #is_dfx<-is_df[order(is_df$cluster),]
  iso_sel<-is_df[selectedd]
  colnames(iso_sel)<-select_gn
  iso_sel$cluster<-is_df$cluster
  iso_sel<-aggregate(.~cluster,data=iso_sel, FUN=mean)
  rownames(iso_sel)<-iso_sel$cluster
  iso_sel$cluster<-NULL
  p<-plot_ly(x=colnames(iso_sel),y=rownames(iso_sel),z = as.matrix(iso_sel),type='heatmap',zauto=FALSE,zmin=0,zsmooth='none',colorscale='Portland') %>%
    layout(yaxis = list(autorange = "reversed"))%>% colorbar(title = "Mean log TPM")
  p
}
