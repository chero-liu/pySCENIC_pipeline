findHighGene = function(mdata,top,highGene,TF){
    outHighGene=c()
        for (i in highGene){
        if (i %in% mdata$Gene){
            if (which(mdata['Gene'] == i)<=top){
                print(paste0(i,' top = ',which(mdata['Gene'] == i),' <= ',top,' in TF:',TF))
            }else{
                print(paste0(i,' top = ',which(mdata['Gene'] == i),' > ',top,' in TF:',TF))
                outHighGene = append(outHighGene,i)
            }
        }else{
            print(paste0(i,' not in TF:', TF))
        }
    }
    
    return(mdata[which(mdata['Gene'] == outHighGene),])
}
myscale = function(range,col){
    k = range/(max(col)-min(col))
    scale_value = c()
    for (i in col){
        scale_value = append(scale_value,k*(i-min(col)))
    }
    return(scale_value)
}
top_data = function(data,col,top){
    data = data[rev(order(data[,col])),]
    data = head(data,n = top)
    return(data)
}