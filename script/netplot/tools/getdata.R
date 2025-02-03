library(stringr)
library (argparser)
source("/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/scripts/script/netplot/utils/func.R")
argv <- arg_parser('')
argv <- add_argument(argv,"--regulons", help="scenic_regulons_importance.tsv file path")
argv <- add_argument(argv,"--TF", help="interested TF")
argv <- add_argument(argv,"--highGene", help="interested Gene",default = 'F')
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--top", help='',default=30)
argv <- add_argument(argv,"--showlabel", help='',default=0.3)
argv <- add_argument(argv,"--highGeneSize", help="",default = 3)
argv <- add_argument(argv,"--highGeneColor", help="",default = "00AA55")
argv <- add_argument(argv,"--TFColor", help="",default = "1E90FF")
argv <- add_argument(argv,"--geneColor", help="",default = "E87D72")
argv <- parse_args(argv)

metadata = read.csv(argv$regulons,sep = '\t')
intgene = unlist(str_split(argv$TF,','))
top = as.numeric(argv$top)
showlabel = as.numeric(argv$showlabel)
outdir = argv$outdir


dir.create(outdir)

setwd(outdir)

if (argv$highGene != 'F'){
    highGene = unlist(str_split(argv$highGene,','))
    metadata_tmp = list()
    outHighGene = list()
    topgene = c()
    for (i in intgene){
        metadata_tmp[[i]] = metadata[which(metadata['TF']==i),]
        metadata_tmp[[i]]$topNub = 1:nrow(metadata_tmp[[i]])
        row.names(metadata_tmp[[i]]) = 1:nrow(metadata_tmp[[i]])
        outHighGene[[i]] = findHighGene(metadata_tmp[[i]],top,highGene,i)
        metadata_tmp[[i]] = top_data(metadata_tmp[[i]],col = 'Weight',top = top)
        #metadata_tmp[[i]]$Weight = myscale(2,metadata_tmp[[i]]$Weight)
        topgene = append(topgene,top_data(metadata_tmp[[i]],'Weight',round(top*showlabel))$Gene)
    }
    edge_df = rbind(Reduce(rbind,metadata_tmp),Reduce(rbind,outHighGene))
    gene_label = c(intgene,topgene,highGene)
    vv = unique(c(unique(edge_df$Gene),unique(edge_df$TF)))
    vertex_df = data.frame(name=vv,label='')
    for(i in gene_label){
        vertex_df$label[which(vertex_df$name == i)] = i
    }
    vertex_df$size = ifelse(vertex_df$name %in% intgene,10,ifelse(vertex_df$name %in% highGene,as.numeric(argv$highGeneSize),ifelse(vertex_df$name %in% gene_label,3,1)))
    vertex_df$color = ifelse(vertex_df$name %in% intgene,paste0('#',argv$TFColor),ifelse(vertex_df$name %in% highGene,paste0('#',argv$highGeneColor),ifelse(vertex_df$name %in% gene_label,paste0('#',argv$geneColor),"#FFFFFF")))
}else{
    metadata_tmp = list()
    topgene = c()
    for (i in intgene){
        metadata_tmp[[i]] = metadata[which(metadata['TF']==i),]
        metadata_tmp[[i]]$topNub = 1:nrow(metadata_tmp[[i]])
        row.names(metadata_tmp[[i]]) = 1:nrow(metadata_tmp[[i]])
        metadata_tmp[[i]] = top_data(metadata_tmp[[i]],col = 'Weight',top = top)
        #metadata_tmp[[i]]$Weight = myscale(2,metadata_tmp[[i]]$Weight)
        topgene = append(topgene,top_data(metadata_tmp[[i]],'Weight',round(top*showlabel))$Gene)
    }
    edge_df = Reduce(rbind,metadata_tmp)
    gene_label = c(intgene,topgene)
    vv = unique(c(unique(edge_df$Gene),unique(edge_df$TF)))
    vertex_df = data.frame(name=vv,label='')
    for(i in gene_label){
        vertex_df$label[which(vertex_df$name == i)] = i
    }
    vertex_df$size = ifelse(vertex_df$name %in% intgene,10,ifelse(vertex_df$name %in% gene_label,3,1))
    vertex_df$color = ifelse(vertex_df$name %in% intgene,paste0('#',argv$TFColor),ifelse(vertex_df$name %in% gene_label,paste0('#',argv$geneColor),"#FFFFFF"))
}


write.csv(edge_df,'edge.csv',row.names=FALSE)
write.csv(vertex_df,'vertex.csv',row.names=FALSE)
