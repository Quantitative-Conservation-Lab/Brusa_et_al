library(tidyverse)
library(here)

sd_df <- read.csv(here('Data', 'seaducks.csv'))


#Mis-Identification by group, then by individual, then integrate mis-ID info into detection model
#Remove uninformative records
sb_df <- sb_df %>% filter(sb_df$SPECIES != "PEFA" & sb_df$SPECIES != "UNSD" & sb_df$SPECIES != "BAEA" 
                          & sb_df$SPECIES != "BLOY" & sb_df$SPECIES != "GBHE" & sb_df$SPECIES != "NWCR"
                          & sb_df$SPECIES != "HAPO" & sb_df$SPECIES != "HASE" & sb_df$SPECIES != "BLLA"
                          & sb_df$SPECIES != "UNMA" & sb_df$SPECIES !=  "UNPD" & sb_df$SPECIES !=  "UMSD" 
                          & sb_df$SPECIES != "UNDD" & sb_df$SPECIES != "UNDU" & sb_df$SPECIES !=  "USSD"
                          & sb_df$SPECIES != "UNSD" & sb_df$SPECIES != "UNSB")


sb_df <- sb_df %>%
  mutate(case_when(str_detect(sb_df$SPECIES, "SC$") ~ "Scoter", str_detect(sb_df$SPECIES, "LO$") ~ "Loon", 
                   str_detect(sb_df$SPECIES, "GR$") ~ "Grebe", str_detect(sb_df$SPECIES, "ME$") ~ "Merganser", 
                   str_detect(sb_df$SPECIES, "UNGO") | str_detect(sb_df$SPECIES, "COGO") | 
                     str_detect(sb_df$SPECIES, "BAGO") ~ "Goldeneye", str_detect(sb_df$SPECIES, "FF$") ~ "Bufflehead",
                   str_detect(sb_df$SPECIES, "UNGU") | str_detect(sb_df$SPECIES, "GWGU") | 
                     str_detect(sb_df$SPECIES, "WG$") ~ "Gull",
                   str_detect(sb_df$SPECIES, "MU$") | str_detect(sb_df$SPECIES, "ML$") ~ "Murrelet",
                   str_detect(sb_df$SPECIES, "CO$") ~ "Cormorant", str_detect(sb_df$SPECIES, "LL$") ~ "Mallard",
                   str_detect(sb_df$SPECIES, "WI$") ~ "Wigeon", 
                   str_detect(sb_df$SPECIES, "PI$") ~ "Pintail", str_detect(sb_df$SPECIES, "BR$") ~ "Brant",
                   str_detect(sb_df$SPECIES, "DS$") ~ "Long-tailed Duck",
                   str_detect(sb_df$SPECIES, "RUDU") ~ "Ruddy Duck", str_detect(sb_df$SPECIES, "CAAU") | 
                     str_detect(sb_df$SPECIES, "PAAU") | str_detect(sb_df$SPECIES, "RHAU") ~ "Auklet", 
                   str_detect(sb_df$SPECIES, "PIGU") ~ "Guillemot", str_detect(sb_df$SPECIES, "AC$") ~ "Alcid", 
                   str_detect(sb_df$SPECIES, "HADU") ~ "Harelquin Duck", str_detect(sb_df$SPECIES, "DO$") ~ "Dove",
                   str_detect(sb_df$SPECIES, "CAGO") ~ "Canadian Goose", str_detect(sb_df$SPECIES, "SCAU") ~ "Scaup"))

colnames(sb_df)[30] <- "sp.group"

