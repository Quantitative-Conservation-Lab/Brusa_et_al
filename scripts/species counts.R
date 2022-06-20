library(tidyverse)
library(here)

sb_df <- read.csv(here('Data', 'seaducks.csv'))

#Remove non-birds
sb_df <- sb_df %>% filter(sb_df$SPECIES != "HAPO" & sb_df$SPECIES != "HASE" & sb_df$SPECIES != "BLLA"
                          & sb_df$SPECIES != "UNMA")

sb_df$Transect <- str_trunc(sb_df$TransectGroup, 12, "right", ellipsis = "")

sb_df <- sb_df %>% filter(platform == 1)

tran <- unique(sb_df$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))

sb_df$numeric_tran <- 0
for(i in dumb$tran){
  sb_df$numeric_tran[sb_df$Transect == i] <- match(i, dumb$tran)
}


SB.FF <- sb_df %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.FF,
              values_fn = sum, values_fill = 0)
SB.FF <- pivot_longer(data = SB.FF, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

SB.FF.total <- SB.FF %>% group_by(Species) %>%
  summarise(FF.Total = sum(Count))


SB.POV <- sb_df %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.POV,
              values_fn = sum, values_fill = 0)
SB.POV <- pivot_longer(data = SB.POV, cols = !numeric_tran, names_to = "Species", 
                        values_to = "Count")

SB.POV.total <- SB.POV %>% group_by(Species) %>%
  summarise(POV.Total = sum(Count))


SB.BM <- sb_df %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.BM,
              values_fn = sum, values_fill = 0)

SB.BM <- pivot_longer(data = SB.BM, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

SB.BM.total <- SB.BM %>% group_by(Species) %>%
  summarise(BM.Total = sum(Count))


SB.TC <- sb_df %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.TC,
              values_fn = sum, values_fill = 0)

SB.TC <- pivot_longer(data = SB.TC, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

SB.TC.total <- SB.TC %>% group_by(Species) %>%
  summarise(TC.Total = sum(Count))



#Classify to genus
small.groups <- sb_df %>%
  mutate(case_when(str_detect(sb_df$SPECIES, "SC$") ~ "Melanitta", str_detect(sb_df$SPECIES, "LO$") ~ "Gavia", 
                   str_detect(sb_df$SPECIES, "HOGR")  | str_detect(sb_df$SPECIES, "RNGR") ~ "Podiceps", 
                   str_detect(sb_df$SPECIES, "ME$") ~ "Mergus", str_detect(sb_df$SPECIES, "UNGO") | 
                   str_detect(sb_df$SPECIES, "COGO") | str_detect(sb_df$SPECIES, "BAGO") | 
                   str_detect(sb_df$SPECIES, "BUFF") ~ "Bucephala", str_detect(sb_df$SPECIES, "UNGU") |
                   str_detect(sb_df$SPECIES, "UBWG")  ~ "Gull", 
                   str_detect(sb_df$SPECIES, "GWGU") ~ "Larus", str_detect(sb_df$SPECIES, "MAMU") ~ "Brachyramphus", 
                   str_detect(sb_df$SPECIES, "WEGR") ~ "Aechmophorus",
                   str_detect(sb_df$SPECIES, "ANMU") ~ "Synthliboramphus", str_detect(sb_df$SPECIES, "ML$") ~ "Murrelet",
                   str_detect(sb_df$SPECIES, "DCCO") | str_detect(sb_df$SPECIES, "PECO") ~ "Phalacrocorax",
                   str_detect(sb_df$SPECIES, "UNCO") ~ "Cormorant", str_detect(sb_df$SPECIES, "UNAC") | 
                   str_detect(sb_df$SPECIES, "USAC") ~ "Alcid", str_detect(sb_df$SPECIES, "UNDD") ~ "Diving duck",
                   str_detect(sb_df$SPECIES, "WI$") ~ "Mareca", str_detect(sb_df$SPECIES, "BAEA") ~ "Haliaeetus",
                   str_detect(sb_df$SPECIES, "PI$") | str_detect(sb_df$SPECIES, "LL$") ~ "Anas", 
                   str_detect(sb_df$SPECIES, "BR$") | str_detect(sb_df$SPECIES, "CAGO") ~ "Branta", 
                   str_detect(sb_df$SPECIES, "NWCR") ~ "Corvus", str_detect(sb_df$SPECIES, "BLOY") ~ "Haematopus",
                   str_detect(sb_df$SPECIES, "DS$") ~ "Clangula", str_detect(sb_df$SPECIES, "GBHE") ~ "Ardea",
                   str_detect(sb_df$SPECIES, "RUDU") ~ "Oxyura", str_detect(sb_df$SPECIES, "COMU") | 
                   str_detect(sb_df$SPECIES, "UNMU") ~ "Uria", str_detect(sb_df$SPECIES, "PEFA") ~"Falco",
                   str_detect(sb_df$SPECIES, "RHAU") ~ "Cerorhinca", str_detect(sb_df$SPECIES, "UNDU") ~ "Duck",
                   str_detect(sb_df$SPECIES, "PIGU") ~ "Cepphus", str_detect(sb_df$SPECIES, "PIGU") ~ "Guillemot",  
                   str_detect(sb_df$SPECIES, "HADU") ~ "Histrionicus", str_detect(sb_df$SPECIES, "DO$") ~ "Columba",
                   str_detect(sb_df$SPECIES, "UNGR") ~ "Grebe", str_detect(sb_df$SPECIES, "UMSD") | 
                   str_detect(sb_df$SPECIES, "USSD") ~ "Shorebird", str_detect(sb_df$SPECIES, "UNML") ~ "Murrelet", 
                   str_detect(sb_df$SPECIES, "UNSB") ~ "Seabird", str_detect(sb_df$SPECIES, "SCAU") ~ "Aythya", 
                   str_detect(sb_df$SPECIES, "UNPD") ~ "Dabbling duck"))


#fitting one model for each observer...better than running two models to not have
colnames(small.groups)[32] <- "Genus"

SG.FF <- small.groups %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Genus, values_from = Count.FF,
              values_fn = sum, values_fill = 0)
SG.FF <- pivot_longer(data = SG.FF, cols = !numeric_tran, names_to = "Genus", 
                      values_to = "Count")

SG.FF.total <- SG.FF %>% group_by(Genus) %>%
  summarise(FF.Total = sum(Count))


SG.POV <- small.groups %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Genus, values_from = Count.POV,
              values_fn = sum, values_fill = 0)
SG.POV <- pivot_longer(data = SG.POV, cols = !numeric_tran, names_to = "Genus", 
                       values_to = "Count")

SG.POV.total <- SG.POV %>% group_by(Genus) %>%
  summarise(POV.Total = sum(Count))


SG.BM <- small.groups %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Genus, values_from = Count.BM,
              values_fn = sum, values_fill = 0)

SG.BM <- pivot_longer(data = SG.BM, cols = !numeric_tran, names_to = "Genus", 
                      values_to = "Count")

SG.BM.total <- SG.BM %>% group_by(Genus) %>%
  summarise(BM.Total = sum(Count))


SG.TC <- small.groups %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Genus, values_from = Count.TC,
              values_fn = sum, values_fill = 0)

SG.TC <- pivot_longer(data = SG.TC, cols = !numeric_tran, names_to = "Genus", 
                      values_to = "Count")

SG.TC.total <- SG.TC %>% group_by(Genus) %>%
  summarise(TC.Total = sum(Count))


large.groups <- sb_df %>%
  mutate(case_when(str_detect(sb_df$SPECIES, "SC$") | str_detect(sb_df$SPECIES, "BUFF") | 
                   str_detect(sb_df$SPECIES, "UNGO") | str_detect(sb_df$SPECIES, "COGO") | 
                   str_detect(sb_df$SPECIES, "BAGO") | str_detect(sb_df$SPECIES, "ME$") |
                   str_detect(sb_df$SPECIES, "HADU") | str_detect(sb_df$SPECIES, "DS$") |
                   str_detect(sb_df$SPECIES, "RUDU") | str_detect(sb_df$SPECIES, "SCAU") |
                   str_detect(sb_df$SPECIES, "UNDD") ~ "diving duck", str_detect(sb_df$SPECIES, "BAEA") |
                   str_detect(sb_df$SPECIES, "PEFA") ~ "hawks and falcons", 
                   str_detect(sb_df$SPECIES, "NWCR") ~ "corvids", str_detect(sb_df$SPECIES, "LO$") ~ "loon", 
                   str_detect(sb_df$SPECIES, "GR$") ~ "grebe", str_detect(sb_df$SPECIES, "UNGU") | 
                   str_detect(sb_df$SPECIES, "GWGU") | str_detect(sb_df$SPECIES, "WG$") ~ "gull", 
                   str_detect(sb_df$SPECIES, "UNSB") | str_detect(sb_df$SPECIES, "CO$") ~ "sea bird", 
                   str_detect(sb_df$SPECIES, "LL$") | str_detect(sb_df$SPECIES, "WI$") |
                   str_detect(sb_df$SPECIES, "NOPI") | str_detect(sb_df$SPECIES, "UNPD") ~ "dabbling duck", 
                   str_detect(sb_df$SPECIES, "GBHE") | str_detect(sb_df$SPECIES, "BLOY") |
                   str_detect(sb_df$SPECIES, "SD$") ~ "shorebird", str_detect(sb_df$SPECIES, "UNDU") ~ "duck",
                   str_detect(sb_df$SPECIES, "BR$") | str_detect(sb_df$SPECIES, "CAGO") ~ "waterfowl",
                   str_detect(sb_df$SPECIES, "CAAU") | str_detect(sb_df$SPECIES, "PAAU") | 
                   str_detect(sb_df$SPECIES, "RHAU") | str_detect(sb_df$SPECIES, "AC$") | 
                   str_detect(sb_df$SPECIES, "MU$") | str_detect(sb_df$SPECIES, "PIGU") | 
                   str_detect(sb_df$SPECIES, "ML$") ~ "alcid", str_detect(sb_df$SPECIES, "DO$") ~ "dove"))

#fitting one model for each observer...better than running two models to not have
colnames(large.groups)[32] <- "Label"

LG.FF <- large.groups %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Label, values_from = Count.FF,
              values_fn = sum, values_fill = 0)
LG.FF <- pivot_longer(data = LG.FF, cols = !numeric_tran, names_to = "Label", 
                      values_to = "Count")

LG.FF.total <- LG.FF %>% group_by(Label) %>%
  summarise(FF.Total = sum(Count))


LG.POV <- large.groups %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Label, values_from = Count.POV,
              values_fn = sum, values_fill = 0)
LG.POV <- pivot_longer(data = LG.POV, cols = !numeric_tran, names_to = "Label", 
                       values_to = "Count")

LG.POV.total <- LG.POV %>% group_by(Label) %>%
  summarise(POV.Total = sum(Count))


LG.BM <- large.groups %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Label, values_from = Count.BM,
              values_fn = sum, values_fill = 0)

LG.BM <- pivot_longer(data = LG.BM, cols = !numeric_tran, names_to = "Label", 
                      values_to = "Count")

LG.BM.total <- LG.BM %>% group_by(Label) %>%
  summarise(BM.Total = sum(Count))


LG.TC <- large.groups %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Label, values_from = Count.TC,
              values_fn = sum, values_fill = 0)

LG.TC <- pivot_longer(data = LG.TC, cols = !numeric_tran, names_to = "Label", 
                      values_to = "Count")

LG.TC.total <- LG.TC %>% group_by(Label) %>%
  summarise(TC.Total = sum(Count))




