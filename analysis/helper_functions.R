library(pROC)
library(plyr)
library(reshape2)
library(pheatmap)
library(gridExtra)
library(ggplot2)

# FILTER FUNCTIONS ######################################################

# Function to eliminate particular column name
remove_columns = function(df,column_names){
  idx = which(colnames(df)%in%column_names)
  df = df[,-idx]
  colnames(df)[1] = "UID"
  return(df)
}

# Function to parse compiled data input
parse_input = function(input){
  input$name = rownames(input)
  df = melt(input)
  thresholds = c()
  directions = c()
  for (name in df$name){
    thresholds = c(thresholds,strsplit(name,"_")[[1]][3])
    directions = c(directions,strsplit(name,"_")[[1]][4])
  }
  df$thresh = thresholds
  df$direction = directions
  return(df)
}

# Function to parse individual data input (also has thresh and direction column)
parse_single_input = function(input){
  input$name = rownames(input)
  df = melt(input,id.var=c("thresh","name","direction"))
  number = paste("^",iter,"_",sep="")
  df$name = gsub("_0.0|_1.0|_2.0|_3.0|_4.0|_5.0|_6.0|_7.0|_8.0|_9.0|_10.0|_11.0|_12.0|_13.0","",df$name)
  df$name = gsub("_pos|_posneg","",df$name)
  df$name = gsub(number,"",df$name)
  df$variable = gsub(paste("X",iter,"_",sep=""),"",df$variable)
  colnames(df) = c("thresh","comparisonimage","direction","queryimage","score")
  return(df)
}


# PLOTTING FUNCTIONS ####################################################
plot_distribution = function(df,title,ymax=2.6){

  gg = ggplot(df,aes(x=score,fill=strategy)) +
    geom_density(alpha=0.25) + 
    ylab("") + 
    xlab("") +
    ylim(0,ymax) +
    xlim(-1,1) +
    facet_wrap(~thresh,nrow=1) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.margin=unit(c(0,0,0,5),"mm")) +
    guides(fill=FALSE) + 
    ggtitle(title)
    
  return(gg)
}

plot_result = function(res,thresh,direction="posneg",outfile) {

  pdp = res$pdp[res$pdp$direction==direction,-which(colnames(res$pdp)=="direction")]
  pip = res$pip[res$pip$direction==direction,-which(colnames(res$pip)=="direction")]
  pds = res$pds[res$pds$direction==direction,-which(colnames(res$pds)=="direction")]
  pis = res$pis[res$pis$direction==direction,-which(colnames(res$pis)=="direction")]
  
  # Filter down to threshold of interest
  pdp = pdp$value[pdp$thresh==thresh] 
  pip = pip$value[pip$thresh==thresh] 
  pds = pds$value[pds$thresh==thresh] 
  pis = pis$value[pis$thresh==thresh] 
    
  pdp = cbind(pdp,"cca.pearson")
  pip = cbind(pip,"svi.pearson")
  pds = cbind(pds,"cca.spearman")
  pis = cbind(pis,"svi.spearman")
  
  flat = as.data.frame(rbind(pdp,pip,pds,pis),stringsAsFactors=FALSE)
  colnames(flat) = c("value","variable")
  flat$value = as.numeric(flat$value)

  # NAs mean the comparison was not possible - separate the data to plot
  number_na = length(which(is.na(flat$value)))
  if (number_na>0){
    nans = flat[which(is.na(flat$value)),]
    flat = flat[-which(is.na(flat$value)),]
    nans$value = 1
    values = c(0,0,0,0)
    variables = c("cca.pearson","svi.pearson","cca.spearman","svi.spearman")
    add = data.frame(value=values,variable=variables)
    nans = rbind(nans,add)
  } else {
    values = c(0,0,0,0)
    variables = c("cca.pearson","svi.pearson","cca.spearman","svi.spearman")
    # Make something tiny no one will see for empty plot
    nans = data.frame(value = values, variable=variables)
  }
  nans$variable = as.character(nans$variable)
  
  if (dim(flat)[1]!=0) {
    
    if (thresh!="9.0"){
      densityplot = ggplot(flat,aes(x=value, fill=variable)) + 
        geom_density(alpha=0.25) + 
        ylab("Density") + 
        xlab(paste("Threshold",thresh)) + 
        xlim(-1,1) +
        ylim(0,12) +
        theme(legend.position="none")
    } else {
      densityplot = ggplot(flat,aes(x=value, fill=variable)) + 
        geom_density(alpha=0.25) + 
        ylab("Density") + 
        xlab(paste("Threshold",thresh)) + 
        xlim(-1,1) +
        theme(legend.position="none")  
    }  
  nanplot = ggplot(nans,aes(value, fill=variable)) + 
    geom_bar(alpha=0.25,binwidth=12) + 
    xlab("NaN Count") + 
    ylab("") +
    ylim(0,2500) +
    facet_wrap(~variable,nrow=1) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    guides(fill=guide_legend(title=NULL)) +
    scale_x_discrete(breaks=NULL)
  
  g = arrangeGrob(densityplot, nanplot, ncol=2)
  ggsave(file=outfile,g)    
  }
}

# STATISTICS FUNCTIONS ####################################################

# Function to calculate confidence intervals
get_ci = function(dat,direction="upper"){
  error = qnorm(0.975)*sd(dat)/sqrt(length(dat))
  if (direction=="upper"){
    return(mean(dat)+error)
  } else {
    return(mean(dat)-error)    
  }
}

# Calculate accuracy across all contrasts
calculate_accuracy = function(input,queryimages,thresholds,directions,label) {

  res = c()
  
  for (thresh in thresholds){
    sub = input[which(input$thresh==thresh),]
    for (direction in directions){
      sub1 = sub[which(sub$direction==direction),]
      # For each query image, calculate accuracy
      acc = c()
      for (image in queryimages){
        sub2 = sub1[which(sub1$queryimage==image),]
        sub2 = sub2[with(sub2, order(-score)), ]
        if (sub2$comparisonimage[1] == image) {
          acc = c(acc,1)
        } else {
          acc = c(acc,0)
        }
      }
      acc = sum(acc) / length(acc)
      res = rbind(res,cbind(acc,thresh,direction))
    }
  }
  res = as.data.frame(res)
  res$acc = as.numeric(as.character(res$acc))
  res$thresh = as.numeric(as.character(res$thresh))
  res$strategy = label
  return(res)
}

# Calculate accuracy for single contrast
calculate_single_accuracy = function(input,image,thresholds,directions,label,iter) {
  
  res = c()
  
  for (thresh in thresholds){
    sub = input[which(input$thresh==thresh),]
    for (direction in directions){
      sub1 = sub[which(sub$direction==direction),]
      # For each query image, calculate accuracy
      sub2 = sub1[which(sub1$queryimage==image),]
      sub2 = sub2[with(sub2, order(-score)), ]
      if (sub2$comparisonimage[1] == image) {
          acc = 1
      } else {
          acc = 0
      }
      res = rbind(res,cbind(acc,thresh,direction,iter))
    }
  }
  res = as.data.frame(res)
  res$acc = as.numeric(as.character(res$acc))
  res$thresh = as.numeric(as.character(res$thresh))
  res$strategy = label
  return(res)
}


make_confmatrix = function(input,confmatrix,thresholds,directions,label) {
  
  for (thresh in thresholds){
    sub = input[which(input$thresh==thresh),]
    for (direction in directions){
      sub1 = sub[which(sub$direction==direction),]
      # For each query image, find actual and predicted
      for (actual in queryimages){
        sub2 = sub1[which(sub1$queryimage==actual),]
        sub2 = sub2[with(sub2, order(-score)), ]
        predicted = sub2$comparisonimage[1]
        confmatrix[actual,predicted] = confmatrix[actual,predicted] + 1        
      }     
    }
  }
  return(confmatrix)
}

# Get mean accuracy and sd across iterations for one set of variables
accuracy_table = function(threshold=1.0,direction="posneg",label="cca.pearson",queryimages) {
  
    accuracies = c()
    for (iter in iters){
        # Read in each of the input files
        cat(iter,"\n")
        pd = parse_single_input(read.csv(paste(iter,"_pearson_cca.tsv",sep=""),sep="\t",row.names=1)) 
        # Subset to threshold and direction
        pd = pd[pd$thresh==threshold,]
        pd = pd[pd$direction==direction,]
        for (actual in queryimages){
          sub = pd[which(pd$queryimage==actual),]
          sub2 = sub[with(sub, order(-score)), ]
          predicted = sub2$comparisonimage[1]
          if (predicted==actual){
            result = cbind(actual,1,iter)
            accuracies = rbind(accuracies,result)
          } else {
            result = cbind(actual,0,iter)
            accuracies = rbind(accuracies,result)
          }
        }
    }
    df = as.data.frame(accuracies,stringsAsFactors=FALSE)
    colnames(df) = c("contrast","classification","iteration")
    # Now we need a mean accuracy and standard deviation for each contrast, across iterations
    tableone = c()
    for (image in queryimages){
      sub = df[df$contrast==image,]
      accuracy_mean = mean(as.numeric(sub$classification))
      accuracy_sd = sd(as.numeric(sub$classification))
      result = cbind(image,accuracy_mean,accuracy_sd)
      tableone = rbind(tableone,result)
    }
    return(tableone)
}

# Get accuracy metrics for one task
accuracy_table = function(threshold=1.0,direction="posneg",label="cca.pearson",queryimages) {
  
  accuracies = c()
  for (iter in iters){
    # Read in each of the input files
    cat(iter,"\n")
    pd = parse_single_input(read.csv(paste(iter,"_pearson_cca.tsv",sep=""),sep="\t",row.names=1)) 
    # Subset to threshold and direction
    pd = pd[pd$thresh==threshold,]
    pd = pd[pd$direction==direction,]
    for (actual in queryimages){
      sub = pd[which(pd$queryimage==actual),]
      sub2 = sub[with(sub, order(-score)), ]
      predicted = sub2$comparisonimage[1]
      if (predicted==actual){
        result = cbind(actual,1,iter)
        accuracies = rbind(accuracies,result)
      } else {
        result = cbind(actual,0,iter)
        accuracies = rbind(accuracies,result)
      }
    }
  }
  df = as.data.frame(accuracies,stringsAsFactors=FALSE)
  colnames(df) = c("contrast","classification","iteration")
  # Now we need a mean accuracy and standard deviation for each contrast, across iterations
  tableone = c()
  for (image in queryimages){
    sub = df[df$contrast==image,]
    accuracy_mean = mean(as.numeric(sub$classification))
    accuracy_sd = sd(as.numeric(sub$classification))
    result = cbind(image,accuracy_mean,accuracy_sd)
    tableone = rbind(tableone,result)
  }
  return(tableone)
}
