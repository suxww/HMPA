rm(list = ls())
library(impute)
library(limma)
library(readxl)

setwd("/Users/suxinwan/Desktop/cptac_reanalysis/UCEC/results/02p_exp/")

file="sample_raw_case.csv" 

#read file
rt=read.table(file,sep=",",header=T,check.names=F,row.names=1)
rt=as.matrix(rt)


#data<-rt[which(rowMeans(!is.na(rt)) > 0.85), which(colMeans(!is.na(rt)) > 0.85)]

# Calculate the proportion of missing values per row
missing_percentages_row <- rowMeans(is.na(rt))

# Found row > 80% missing values
rows_to_delete <- which(missing_percentages_row > 0.85)

# Deletes the specified row
rt_select <- rt[-rows_to_delete, ]

#mat=impute.knn(rt)

write.table(rt_select,file="sample_select_raw.csv",sep=",",quote=F,col.names=T,row.names = T)

mat = impute.knn(rt_select ,k = 10, rowmax = 0.85, colmax = 0.85, maxp = 1500, rng.seed=362436069)

rt=mat$data
write.table(rt,file="mat_sample_raw.csv",sep=",",quote=F,col.names=T,row.names = T)


library(pheatmap)
fig1<-pheatmap(rt,
                     main=(""),
                     cluster_rows = T,cluster_cols= F,
                     clustering_distance_rows="correlation",clustering_distance_cols="correlation",
                     #color=colorRampPalette(c("blue", "white","red"))(100),
                     fontsize_row=6,scale="row",border_color=NA,
                     treeheight_row = 60,
                     treeheight_col = 100)


# 2.热图
## 2.2 进行scale，但是不设定最大最小值
library(pheatmap)
library(extrafont)
file="mat_sample_raw_group.csv"
#读取文件
rt=read.table(file,sep=",",header=T,check.names=F,row.names=1)

t <- t(scale(t(rt)))
p1 <- pheatmap(t,
               #annotation_col = col,
               annotation_legend = T,
               border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
               cluster_rows = T, #对行聚类
               cluster_cols = T, #队列聚类
               show_colnames = F, #是否显示列名
               show_rownames = F #是否显示行名
)

## 2.3 进行scale后设定最大最小值的情况 
table(abs(t)>2)
t[t>=2]=2
t[t<=-2]=-2
table(is.na(t))

annotation_col = data.frame(
  group = c(rep("NAT",59),rep("Tumor",65))
)
row.names(annotation_col) <- colnames(t)

groupcolor <- c("#9ccfe6","#e0897e") 
names(groupcolor) <- c("NAT","Tumor") 
ann_colors <- list(group=groupcolor)

p2 <- pheatmap(t,
               annotation_col = annotation_col,
               annotation_legend = T,
               color=colorRampPalette(c("#005f81","white","#b03d26"))(100),
               border_color = F,
               cluster_rows = T, 
               cluster_cols = T, 
               show_colnames = F, 
               show_rownames = F, 
               annotation_colors = ann_colors,
               fontfamily= "serif"
)

##heatmap group
#heatmap
#rt=rt[apply(rt,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
colorList[[gene]]=c("Low"="blue", "High"="red")
j=0
for(cli in colnames(rt2[,3:ncol(rt2)])){
  cliLength=length(levels(factor(rt2[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(rt2[,cli]))
  cliCol["unknow"]="grey75"
  colorList[[cli]]=cliCol
}

#plot
ha=HeatmapAnnotation(df=rt2[,c(1,3:6)], col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt2))
Hm=Heatmap(zero_row_mat, top_annotation=ha)

#save
pdf(file=paste0(gene,"_heatmap.pdf"), width=7, height=5)
draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
