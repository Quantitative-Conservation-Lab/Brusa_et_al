library(tidyverse)
library(here)

sb_df <- read.csv(here('Data', 'seaducks.csv'))

sb_df <- sb_df[which(sb_df$platform == 1),]

codes <- read.csv(here('Data', 'PSEMPSpeciesCodes.csv'))

Count.FF <- aggregate(sb_df$Count.FF,by = list(SPECIES = sb_df$SPECIES),FUN=sum)
Count.POV <- aggregate(sb_df$Count.POV,by = list(SPECIES = sb_df$SPECIES),FUN=sum)
Count.BM <- aggregate(sb_df$Count.BM,by = list(SPECIES = sb_df$SPECIES),FUN=sum)
Count.TC <- aggregate(sb_df$Count.TC,by = list(SPECIES = sb_df$SPECIES),FUN=sum)

sb_out <- data.frame(Count.FF$SPECIES,Count.FF$x,Count.POV$x,Count.BM$x,Count.TC$x)
colnames(sb_out) <- c("Code","Count.FF","Count.POV","Count.BM","Count.TC")

for(i in 1:nrow(sb_out)){
  
  if(length(codes$DESCRIPT[which(codes$CODE == sb_out$Code[i])]) > 0){
    sb_out$Name[i] <- codes$DESCRIPT[which(codes$CODE == sb_out$Code[i])]
  }else (sb_out$Name[i] <- NA)  
}

write.csv(sb_out,file = "Sum_data.csv")