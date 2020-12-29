##Step 0: Export imzml and xml files##
#imzml export parameters from flexImaging: Data point: 15000; Min mass: 50; Max mass: 900
#xml export: export > Spectra List > Spectra List from Regions > Select All

##Step 1: load required packages and codes
source("strawberryIMS_functions.R")
require(Cardinal)
require(XML)
require(Biobase)
require(MALDIquant)
require(tidyverse)

#set working directory
setwd("directory")

##Step 2: Read imzml and save to RData format for later analysis##

#obtain list of imzml
data_list=as.matrix(read.table("imzml.txt"))

#read all the data files
for (i in 1:length(data_list)){
  
  #Use function to fix IMS coord to catch negative x, y coordinates
  setwd("directory for imzml files")
  data=fixIMSCoord(paste(data_list[[i]],".imzml",sep=""))
  
  #save data
  setwd("direcoty for rdata file")
  save(data, file=paste(data_list[[i]], ".Rdata", sep=""))
  
}

##Step 3: Obtain mean intensity across pixels of samples
##Removes pixels below a certain SNR so that only pixels with actual samples
##are included in the mean

#obtain imzml list from rdata
setwd("directory for rdata file")
data_list=as.matrix(read.table("rdata.txt"))

#obtain dataset/filenames, requires tidyverse
filename_list=str_remove(data_list, ".Rdata")

#obtain xml list
setwd("directory for xml file")
xml_list=as.matrix(read.table("xml.txt"))

#species list
species_list=c(219.0, 230.9, 381.1, 433.1, 
               449.1, 519.1, 535.1, 579.1)

stats_all=list()
table_all=list()

for (i in 1:length(data_list)){
  #read data
  #import and subset data
  setwd("directory for rdata file")
  load(data_list[[i]])
  
  #separate data from samples
  setwd("directory for xml file")
  data_split<-sampleROIs(data,XMLdata=xml_list[i])
  
  #visualize results
  par(mfrow=c(3,3), mai = c(0.4,0.4,0.4,0.4))
  image(data,mz=species_list[1], plusminus=0.2)
  image(data_split,mz=species_list[1], plusminus=0.2)
  
  #normalize data
  #data_split_norm=normalize(data_split,method="tic")
  #image(data_split_norm,mz=species_list[1], plusminus=0.1)
  
  #calculate SNR, remove those below certain threshold
  table=imzml_SNR(imzml=data_split, imzml_label=filename_list[[i]],
                  peaklist=species_list, ppm=250, 
                  noise="SuperSmoother",SNR=3,output="height", 
                  mz_range=5, peakwidth=1)
  
  #calculate stats with tidyverse
  #edit directly in the function if need adjustment
  stats=imzml_stats(table, quantile=c(0.1, 0.9))
  
  #save into list
  table_all[[i]]=table
  stats_all[[i]]=stats
}

##Step 4: Combine all calculations into one table and save##
#combine into one large table
stats_all=do.call(rbind, stats_all)
table_all=do.call(rbind, table_all)

#export the data
setwd("directory for export")
name=c("stats.csv")
write.csv(stats_all,file=name)

#export the data
setwd("directory for export")
name=c("table.Rdata")
save(table_all,file=name)
