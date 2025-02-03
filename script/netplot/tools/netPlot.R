library(igraph)
library(stringr)
library (argparser)
source("/SGRNJ03/PiplineTest01/Software_test/liuchenglong/netfunc/netplot_function.R")
argv <- arg_parser('')
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--prefix",help='prefix')
argv <- add_argument(argv,"--highGene", help="interested Gene",default = "F")
argv <- add_argument(argv,"--highGeneColor", help="",default = "9999FF")
argv <- add_argument(argv,"--TFColor", help="",default = "1E90FF")
argv <- add_argument(argv,"--geneColor", help="",default = "E87D72")
argv <- add_argument(argv,"--edgeColor", help="",default = "00CED1")
argv <- add_argument(argv,"--showlabel", help='',default=0.3)
argv <- add_argument(argv,"--vertexLabelSize", help='',default=0.8)
argv <- add_argument(argv,"--legend", help="interested Gene",default = "T")
argv <- parse_args(argv)
setwd(argv$outdir)
prefix = argv$prefix
edgeColor = paste0('#',argv$edgeColor)
vertexLabelSize = argv$vertexLabelSize
showlabel = as.numeric(argv$showlabel)
edge_df = read.csv('edge.csv')
vertex_df = read.csv('vertex.csv')
layouts = unlist(str_split('layout_with_lgl,layout_with_kk,layout_with_fr,layout_with_dh,layout_nicely,layout_components',','))

if(argv$highGene != 'F'){
        for (layout in layouts){
        if (showlabel==0){
                if (argv$legend != 'T'){
        pdf(paste0(prefix,'/',layout,'.pdf'), width=10, height=10, useDingbats=FALSE)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout)
        dev.off()
        png(paste0(prefix,'/',layout,'.png'), width = 10,height = 10,units = 'in',res=300)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout)
        dev.off()
                }else{
        pdf(paste0(prefix,'/',layout,'.pdf'), width=10, height=10, useDingbats=FALSE)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout,
                legend.vertex.label = c('TF','Highlight Gene','Other Gene'),legend.vertex.color = c(paste0('#',argv$TFColor),paste0('#',argv$highGeneColor),"#FFFFFF")
                )
        dev.off()
        png(paste0(prefix,'/',layout,'.png'), width = 10,height = 10,units = 'in',res=300)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout,
                legend.vertex.label = c('TF','Highlight Gene','Other Gene'),legend.vertex.color = c(paste0('#',argv$TFColor),paste0('#',argv$highGeneColor),"#FFFFFF")
                )
        dev.off()
                }
        }else{
        if (argv$legend != 'T'){
        pdf(paste0(prefix,'/',layout,'.pdf'), width=10, height=10, useDingbats=FALSE)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout)
        dev.off()
        png(paste0(prefix,'/',layout,'.png'), width = 10,height = 10,units = 'in',res=300)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout)
        dev.off()
                }else{
        pdf(paste0(prefix,'/',layout,'.pdf'), width=10, height=10, useDingbats=FALSE)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout,
                legend.vertex.label = c('TF',paste0('Top ',showlabel*100,'% of Gene'),'Highlight Gene','Other Gene'),legend.vertex.color = c(paste0('#',argv$TFColor),paste0('#',argv$geneColor),paste0('#',argv$highGeneColor),"#FFFFFF")
                )
        dev.off()
        png(paste0(prefix,'/',layout,'.png'), width = 10,height = 10,units = 'in',res=300)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout,
                legend.vertex.label = c('TF',paste0('Top ',showlabel*100,'% of Gene'),'Highlight Gene','Other Gene'),legend.vertex.color = c(paste0('#',argv$TFColor),paste0('#',argv$geneColor),paste0('#',argv$highGeneColor),"#FFFFFF")
                )
        dev.off()
                }
        }
        }
}else{
        for (layout in layouts){
                if (argv$legend != 'T'){
        pdf(paste0(prefix,'/',layout,'.pdf'), width=10, height=10, useDingbats=FALSE)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout)
        dev.off()
        png(paste0(prefix,'/',layout,'.png'), width = 10,height = 10,units = 'in',res=300)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout)
        dev.off()
                }else{
        pdf(paste0(prefix,'/',layout,'.pdf'), width=10, height=10, useDingbats=FALSE)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout,
                legend.vertex.label = c('TF',paste0('Top ',showlabel*100,'% of Gene'),'Other Gene'),legend.vertex.color = c(paste0('#',argv$TFColor),paste0('#',argv$geneColor),"#FFFFFF")
                )
        dev.off()
        png(paste0(prefix,'/',layout,'.png'), width = 10,height = 10,units = 'in',res=300)
        net_plot(vertex.label.cex=vertexLabelSize,edge.color=edgeColor,edge = edge_df,vertex = vertex_df,vertex.size = vertex_df$size,vertex.color  = vertex_df$color,layout = layout,
                legend.vertex.label = c('TF',paste0('Top ',showlabel*100,'% of Gene'),'Other Gene'),legend.vertex.color = c(paste0('#',argv$TFColor),paste0('#',argv$geneColor),"#FFFFFF")
                )
        dev.off()
                }
        }
}
