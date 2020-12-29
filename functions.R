##displace the imzml xy when you cannot read through Cardinal
#does not work to remove duplicate coordinates
#requires MALDIquant and Cardinal
fixIMSCoord<-function(imzml){
  #read the data
  mq_data <- MALDIquantForeign::import(imzml, verbose=F)
  
  #extract the intensity data
  intmat <- lapply(mq_data, function(x) x@intensity)
  intmat <- do.call(rbind, intmat)
  
  #obtain the x and y coordinates
  xy <- lapply(mq_data, function(x) x@metaData$imaging$pos )
  xy <- data.frame(do.call(rbind, xy))
  
  # modify range so that the minimum value is positive
  if(min(range(xy$x))<=0 | min(range(xy$y))<=0){
    xy <- xy - min(range(xy$x), range(xy$y))+1
  }
  
  #convert back to MSImageSet format
  ds <- MSImageSet(spectra = t(intmat), mz=mq_data[[1]]@mass,
                   coord = data.frame(xy, sample = 'mysample'))

  return(ds)
}

##obtain peak stats, separates by sample
#epxort data matrix from Cardinal and calculates SNR with MALDIquant
imzml_SNR<-function(imzml, imzml_label=NA, peaklist, mz_range=5, 
                    peakwidth=1, output="height", 
                    noise="SuperSmoother", SNR=3, stats_fun=NULL, 
                    ppm=250){
  data=imzml
  resultlist=list()
  ##obtain raw data for desired peak
  #peaklist of compound
  if(length(peaklist)==1) {peaks <- c(peaklist, peaklist+5)} else {peaks <- peaklist}
  ##do peak picking
  if(output != "area" | output != "height") {pmethod="height"} else {pmethod=output}
  data_reduced <- reduceDimension(data, ref=peaks,type=pmethod, width=ppm, units="ppm")
  #reorder data
  order=order(data_reduced$y, data_reduced$x)
  #obtain sample information
  data_reduced_ordered=data_reduced[,order]
  sample=data_reduced_ordered$sample
  sample=as.character(sample)
  #obtain intensity
  intensity_reduced=iData(data_reduced_ordered)
  #give row names
  rownames(intensity_reduced)=data_reduced_ordered@featureData$mz
  ##Determine SNR for each species
  for(i in 1:length(peaklist)){
    species=peaklist[i]
    #truncate the spectra
    data_trunc=data[species-mz_range < Cardinal::mz(data) & Cardinal::mz(data) < species+mz_range,]
    #obtain intensity
    data_trunc_ordered=data_trunc[,order]
    #obtain samplelist
    sample=data_trunc_ordered$sample
    sample=as.character(sample)
    intensity = iData(data_trunc_ordered)
    #obtain mz_range
    mz_trunc=data_trunc_ordered@featureData$mz
    #name matrix
    rownames(intensity)=mz_trunc
    colnames(intensity)=sample
    #create coordinate list
    x=data_reduced_ordered$x
    y=data_reduced_ordered$y
    xy=paste("x",x,"y",y, sep="")
    #create intensity dataframe
    intensity_list=data.frame(xy,intensity_reduced[i,])
    #obtain coordinates for matching
    x=data_trunc_ordered$x
    y=data_trunc_ordered$y
    xy=paste("x",x,"y",y, sep="")
    #obtain SNR list
    sn=list()
    print("Calculating Noise")
    #Progress Bar
    pb   <- txtProgressBar(0, 100, initial=0, style=3)
    for(j in 1:length(sample)){
      spectre<-createMassSpectrum(mz_trunc,intensity[,j])
      nss <- estimateNoise(spectre, method=noise)
      ns_rd=nss[species-peakwidth/2 < nss[,1] & nss[,1]< species+peakwidth/2,2]
      sn[j]=ns_rd[which.max(ns_rd)]
      setTxtProgressBar(pb,round(100/length(sample)*j, digits=0))
    }
    sn=unlist(sn)
    print("Calculating Noise: DONE")
    #combine data
    data_list=data.frame(xy,sample,sn)
    #name the data table
    names(data_list) <- c("xy", "sample","sn")
    #combine data
    final_list <- merge(data_list,intensity_list,by="xy")
    names(final_list)[4]<- c("intensity")
    #create m/z list
    mz_list=rep(species, length(sample))
    #create list of labels
    if(is.na(imzml_label)==TRUE){label <- i} else {label<-imzml_label}
    ds_list=rep(label, length(sample))
    #combine to make final list
    final_list <- data.frame(final_list, final_list$intensity/final_list$sn, mz_list, ds_list)
    names(final_list)[5]<- c("SNR")
    names(final_list)[6]<- c("mz")
    names(final_list)[7]<- c("dataset")
    
    result <- final_list[which(final_list$intensity>(final_list$sn*SNR)),]
    #create resultlist
    resultlist[[i]]=result
  }
  print("Done")
  
  #return final table list
  df=do.call(rbind, resultlist)
  
  #save as dataframe
  return(df)
}

##calculates stats based on defined quantile
#uses tidyverse
imzml_stats<-function(df, quantile=c(0.1, 0.9)){
  #calculating stats
  print("Calculating stats")
  df=df %>% dplyr::group_by(sample, mz, dataset) 
  if (is.null(quantile)==FALSE){
    df_split=group_split(df)
    for (i in 1:length(df_split)){
      x=as.data.frame(df_split[[i]])
      quantiles <-quantile(x$intensity, quantile)
      x <-x[which(x$intensity > quantiles[1] &  x$intensity < quantiles[2]),]
      df_split[[i]]=x
    }
    df=do.call(rbind,df_split)
    df=df %>% dplyr::group_by(sample, mz, dataset) 
  }
  stats=summarise(df, mean(intensity), sd(intensity), median(intensity)) %>% ungroup
  #obtain piexel count
  pixel=count(df, dataset)
  #create final table
  stats_final=bind_cols(stats, pixel[,4])
  #convert to dataframe
  stats_final=as.data.frame(stats_final)
  #
  print("Done")
  return(stats_final)
}
