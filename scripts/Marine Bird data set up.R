library(tidyverse)
library(nimble)
library(ggmcmc)
library(coda)
library(MCMCvis)
library(viridis)
library(here)

sb_df <- read.csv(here('Data', 'seaducks.csv'))

sb_df <- sb_df %>% filter(platform == 1)

#Mis-Identification by group, then by individual, then integrate mis-ID info into detection model
#Remove uninformative records; all shorebirds have been removed; CAGO removed (only one)
sb_df <- sb_df %>% filter(sb_df$SPECIES != "PEFA" & sb_df$SPECIES != "BAEA" & sb_df$SPECIES != "RODO"
                          & sb_df$SPECIES != "BLOY" & sb_df$SPECIES != "GBHE" & sb_df$SPECIES != "NWCR"
                          & sb_df$SPECIES != "HAPO" & sb_df$SPECIES != "HASE" & sb_df$SPECIES != "CAGO"
                          & sb_df$SPECIES != "UNMA" & sb_df$SPECIES !=  "UMSD" & sb_df$SPECIES !=  "USSD")

sb_df$SPECIES <- str_replace(sb_df$SPECIES, "OLDS", "LTDU")


#allocate unknowns for specific data

#Get temp species groups for sums
temp_df <- sb_df %>%
  mutate(case_when(str_detect(sb_df$SPECIES, "UNDD") ~ "UNDD", str_detect(sb_df$SPECIES, "UNPD") ~ "UNPD",
                   str_detect(sb_df$SPECIES, "UNDU") ~ "UNDU", str_detect(sb_df$SPECIES, "UNSB") ~ "UNSB",
                   str_detect(sb_df$SPECIES, "UNAC") | str_detect(sb_df$SPECIES, "USAC") | 
                     str_detect(sb_df$SPECIES, "PIGU") | str_detect(sb_df$SPECIES, "ML$") | 
                     str_detect(sb_df$SPECIES, "MU$") |  
                     str_detect(sb_df$SPECIES, "RHAU") ~ "alcid", str_detect(sb_df$SPECIES, "BUFF") ~ "bufflehead",
                   str_detect(sb_df$SPECIES, "DCCO") | str_detect(sb_df$SPECIES, "PECO") |
                     str_detect(sb_df$SPECIES, "UNCO") ~ "cormorant", str_detect(sb_df$SPECIES, "WI$") |
                     str_detect(sb_df$SPECIES, "MALL") | str_detect(sb_df$SPECIES, "NOPI") ~ "dabbling duck", 
                   str_detect(sb_df$SPECIES, "UNGO") | str_detect(sb_df$SPECIES, "COGO") | 
                     str_detect(sb_df$SPECIES, "BAGO") ~ "goldeneye", str_detect(sb_df$SPECIES, "BLBR") ~ "goose", 
                   str_detect(sb_df$SPECIES, "GR$") ~ "grebe", str_detect(sb_df$SPECIES, "UNGU") | 
                     str_detect(sb_df$SPECIES, "UBWG") | str_detect(sb_df$SPECIES, "GWGU") ~ "gull", 
                   str_detect(sb_df$SPECIES, "HADU") ~ "Harlequin duck",
                   str_detect(sb_df$SPECIES, "LTDU") ~ "long-tailed duck", str_detect(sb_df$SPECIES, "LO$") ~ "loon",
                   str_detect(sb_df$SPECIES, "ME$") ~ "merganser", 
                   str_detect(sb_df$SPECIES, "RUDU") ~ "ruddy duck", str_detect(sb_df$SPECIES, "SCAU") ~ "scaup",
                   str_detect(sb_df$SPECIES, "SC$") ~ "scoter"
  ))

colnames(temp_df)[34] <- "sp.group"

temp_df <- temp_df %>% filter(!str_sub(temp_df$SPECIES,start = 1, end = 1) == "U" | SPECIES == "UNGU" | SPECIES == "UBWG")

#Get group totals and species proportions
spp.sums <- temp_df %>% group_by(SPECIES) %>% summarise(Cam.Count.sp = sum(Count.FF), sp.group = sp.group, 
                                                        Count.FF = Count.FF)
spp.sums <- spp.sums %>% group_by(sp.group) %>% summarise(SPECIES = SPECIES, sp.group = sp.group, 
                                                          Cam.Count.sp = Cam.Count.sp, Cam.Count.grp = sum(Count.FF))
spp.sums <- spp.sums %>% group_by(SPECIES) %>% summarise(sp.group = paste(unique(sp.group)), 
                                                         Cam.Count.sp = round(mean(Cam.Count.sp), 0),
                                                         Cam.Count.grp = round(mean(Cam.Count.grp),0))

spp.sums$Proportions <- spp.sums$Cam.Count.sp/spp.sums$Cam.Count.grp

s.allo_df <- sb_df
s.allo_df <- s.allo_df %>% left_join(x = s.allo_df, y = spp.sums, by = "SPECIES")
s.allo_df <- s.allo_df[,-c(34:36)]

s.allo_df <- s.allo_df %>% rename(Cam.Count.FF = Count.FF)
s.allo_df <- s.allo_df %>% rename(Cam.Count.POV = Count.POV)


#Get 4-letter codes for unknowns
unkvec <- unique(sb_df$SPECIES[str_sub(sb_df$SPECIES,start = 1, end = 1) == "U"])
unkvec <- unkvec[!unkvec %in% c("UNGU", "UBWG", "UNDD", "UNPD")]

#Make a list with element names as 4-letter codes of unknowns
grouplist <- list()
for(i in 1:length(unkvec)){
  u <- unkvec[i]
  grouplist[[u]]$code <- NA
  grouplist[[u]]$prop <- NA
}

#Define the 4-letter codes that should be allocated over for each unknown code
grouplist[["UNDU"]]$code <- c("AMWI", "BAGO", "BLSC", "BUFF", "COGO", "COME", "EUWI", "HADU", "LTDU", "MALL", 
                              "NOPI","RBME", "RUDU", "SCAU", "SUSC", "WWSC")
grouplist[["UNSB"]]$code <- unique(sb_df$SPECIES[str_sub(sb_df$SPECIES, start = 1, end =1) != "U"])
grouplist[["UNAC"]]$code <- c("ANMU", "COMU", "MAMU", "PIGU", "RHAU")
grouplist[["UNML"]]$code <- c("ANMU", "MAMU")
grouplist[["USAC"]]$code <- c("ANMU", "MAMU")
grouplist[["UNCO"]]$code <- c("DCCO", "PECO")
grouplist[["UNGO"]]$code <- c("BAGO", "COGO")
grouplist[["UNGR"]]$code <- c("HOGR", "RNGR", "WEGR")
grouplist[["UNLO"]]$code <- c("COLO", "PALO", "RTLO")
grouplist[["UNME"]]$code <- c("COME", "RBME")
grouplist[["UNSC"]]$code <- c("BLSC", "SUSC", "WWSC")

for(ug in names(grouplist)){
  grouplist[[ug]]$prop <- as.double(unlist(s.allo_df %>%
                                             filter(SPECIES %in% grouplist[[ug]]$code) %>%
                                             filter(!duplicated(SPECIES)) %>%
                                             arrange(SPECIES) %>%
                                             select(Proportions)))
}


for(g in unique(s.allo_df$TransectGroup)){
  for(u in names(grouplist)){
    go <- u %in% s.allo_df$SPECIES[s.allo_df$TransectGroup == g]
    if(!go){next}
    camcols <- which(str_detect(colnames(s.allo_df),"Cam"))
    if(!any(s.allo_df[s.allo_df$TransectGroup == g & s.allo_df$SPECIES == u,camcols] > 0)){next}
    nspp <- length(grouplist[[u]]$code)
    #Need something to make a condition that only adds new lines if Cam.Count.FF > 0 for the UN spp assigned to u
    newline <- s.allo_df[s.allo_df$SPECIES == u & s.allo_df$TransectGroup == g,]
    newlines <- as.data.frame(t(matrix(ncol = nspp, nrow = ncol(newline), data = rep(unlist(newline),nspp))))
    specol <- which(colnames(s.allo_df) == "SPECIES")
    newlines[,specol] <- grouplist[[u]]$code
    newlines[,camcols] <- apply(newlines[,camcols], 2,FUN = function(x){
      as.numeric(x)*grouplist[[u]]$prop
    })
    colnames(newlines) <- colnames(s.allo_df)
    s.allo_df <- rbind(s.allo_df, newlines)
  }
}

s.allo_df <- s.allo_df %>% rename(Count.FF = Cam.Count.FF)
s.allo_df <- s.allo_df %>% rename(Count.POV = Cam.Count.POV)

s.allo_df$SPECIES[s.allo_df$SPECIES == "UNGU"] <- "GULL"


s.allo_df$Count.FF[str_sub(s.allo_df$SPECIES, start = 1, end = 1) == "U"] <- 0
s.allo_df$Count.POV[str_sub(s.allo_df$SPECIES, start = 1, end = 1) == "U"] <- 0

s.allo_df$SPECIES[s.allo_df$SPECIES == "GULL"] <- "UNGU"

#Rename into specific species/species groups
s.allo_df <- s.allo_df %>%
  mutate(case_when(str_detect(s.allo_df$SPECIES, "ANMU") ~ "ancient murrelet", 
                   str_detect(s.allo_df$SPECIES, "MAMU") ~ "marbled murrelet", 
                   str_detect(s.allo_df$SPECIES, "COMU") ~ "common murre",
                   str_detect(s.allo_df$SPECIES, "PIGU") ~ "pigeon guillemot",
                   str_detect(s.allo_df$SPECIES, "RHAU") ~ "rhinoceros auklet", 
                   str_detect(s.allo_df$SPECIES, "BUFF") ~ "bufflehead",
                   str_detect(s.allo_df$SPECIES, "DCCO") ~ "double-crested cormorant",
                   str_detect(s.allo_df$SPECIES, "PECO") ~ "pelagic cormorant",
                   str_detect(s.allo_df$SPECIES, "MALL") ~ "mallard duck",
                   str_detect(s.allo_df$SPECIES, "AMWI") ~ "American wigeon", 
                   str_detect(s.allo_df$SPECIES, "EUWI") ~ "Eurasian wigeon",
                   str_detect(s.allo_df$SPECIES, "NOPI") ~ "northern pintail", 
                   str_detect(s.allo_df$SPECIES, "COGO") ~ "common goldeneye", 
                   str_detect(s.allo_df$SPECIES, "BAGO") ~ "Barrow's goldeneye", 
                   str_detect(s.allo_df$SPECIES, "BLBR") ~ "goose",
                   str_detect(s.allo_df$SPECIES, "RNGR") ~ "red-necked grebe", 
                   str_detect(s.allo_df$SPECIES, "HOGR") ~ "horned grebe", 
                   str_detect(s.allo_df$SPECIES, "WEGR") ~ "western grebe", str_detect(s.allo_df$SPECIES, "UNGU") | 
                   str_detect(s.allo_df$SPECIES, "UBWG") | str_detect(s.allo_df$SPECIES, "GWGU") ~ "gull", 
                   str_detect(s.allo_df$SPECIES, "HADU") ~ "Harlequin duck", 
                   str_detect(s.allo_df$SPECIES, "LTDU") ~ "long-tailed duck", 
                   str_detect(s.allo_df$SPECIES, "COLO") ~ "common loon", 
                   str_detect(s.allo_df$SPECIES, "PALO") ~ "Pacific loon",
                   str_detect(s.allo_df$SPECIES, "RTLO") ~ "red-throated loon",
                   str_detect(s.allo_df$SPECIES, "COME") ~ "common merganser", 
                   str_detect(s.allo_df$SPECIES, "RBME") ~ "red-breasted merganser",
                   str_detect(s.allo_df$SPECIES, "COME") ~ "common merganser", 
                   str_detect(s.allo_df$SPECIES, "RUDU") ~ "ruddy duck", str_detect(s.allo_df$SPECIES, "SCAU") ~ "scaup",
                   str_detect(s.allo_df$SPECIES, "BLSC") ~ "black scoter", 
                   str_detect(s.allo_df$SPECIES, "SUSC") ~ "surf scoter",
                   str_detect(s.allo_df$SPECIES, "WWSC") ~ "white-winged scoter",
                   str_detect(s.allo_df$SPECIES, "UNDD") ~ "UNDD",
                   str_detect(s.allo_df$SPECIES, "UNPD") ~ "UNPD",
                   str_detect(s.allo_df$SPECIES, "UNDU") ~ "UNDU",
                   str_detect(s.allo_df$SPECIES, "UNSB") ~ "UNSB",
                   str_detect(s.allo_df$SPECIES, "UNAC") ~ "UNAC",
                   str_detect(s.allo_df$SPECIES, "UNML") ~ "UNML",
                   str_detect(s.allo_df$SPECIES, "USAC") ~ "USAC",
                   str_detect(s.allo_df$SPECIES, "UNCO") ~ "UNCO",
                   str_detect(s.allo_df$SPECIES, "UNGO") ~ "UNGO",
                   str_detect(s.allo_df$SPECIES, "UNGR") ~ "UNGR",
                   str_detect(s.allo_df$SPECIES, "UNLO") ~ "UNLO",
                   str_detect(s.allo_df$SPECIES, "UNME") ~ "UNME",
                   str_detect(s.allo_df$SPECIES, "UNSC") ~ "UNSC"
                   ))

s.allo_df <- s.allo_df[,-34]
colnames(s.allo_df)[34] <- "sp.group"

#To work in model, round counts to whole numbers
s.allo_df$Count.FF <- round(s.allo_df$Count.FF, 0)
s.allo_df$Count.POV <- round(s.allo_df$Count.POV, 0)



#allocate unknowns for general data
#Group into known general names
g.allo_df <- sb_df %>%
  mutate(case_when(str_detect(sb_df$SPECIES, "UNDD") ~ "UNDD", str_detect(sb_df$SPECIES, "UNPD") ~ "UNPD",
                   str_detect(sb_df$SPECIES, "UNDU") ~ "UNDU", str_detect(sb_df$SPECIES, "UNSB") ~ "UNSB",
                   str_detect(sb_df$SPECIES, "UNAC") | str_detect(sb_df$SPECIES, "USAC") | 
                   str_detect(sb_df$SPECIES, "PIGU") | str_detect(sb_df$SPECIES, "ML$") | 
                   str_detect(sb_df$SPECIES, "MU$") |  
                   str_detect(sb_df$SPECIES, "RHAU") ~ "alcid", str_detect(sb_df$SPECIES, "BUFF") ~ "bufflehead",
                   str_detect(sb_df$SPECIES, "DCCO") | str_detect(sb_df$SPECIES, "PECO") |
                   str_detect(sb_df$SPECIES, "UNCO") ~ "cormorant", str_detect(sb_df$SPECIES, "WI$") |
                   str_detect(sb_df$SPECIES, "MALL") | str_detect(sb_df$SPECIES, "NOPI") ~ "dabbling duck", 
                   str_detect(sb_df$SPECIES, "UNGO") | str_detect(sb_df$SPECIES, "COGO") | 
                   str_detect(sb_df$SPECIES, "BAGO") ~ "goldeneye", str_detect(sb_df$SPECIES, "BLBR") ~ "goose", 
                   str_detect(sb_df$SPECIES, "GR$") ~ "grebe", str_detect(sb_df$SPECIES, "UNGU") | 
                   str_detect(sb_df$SPECIES, "UBWG") | str_detect(sb_df$SPECIES, "GWGU") ~ "gull", 
                   str_detect(sb_df$SPECIES, "HADU") ~ "Harlequin duck",
                   str_detect(sb_df$SPECIES, "LTDU") ~ "long-tailed duck", str_detect(sb_df$SPECIES, "LO$") ~ "loon",
                   str_detect(sb_df$SPECIES, "ME$") ~ "merganser", 
                   str_detect(sb_df$SPECIES, "RUDU") ~ "ruddy duck", str_detect(sb_df$SPECIES, "SCAU") ~ "scaup",
                   str_detect(sb_df$SPECIES, "SC$") ~ "scoter"
                   ))

colnames(g.allo_df)[34] <- "group"

temp_df <- g.allo_df

temp_df <- temp_df %>% filter(SPECIES != "UNDU" | SPECIES != "UNSB")

temp_df <- temp_df %>%
  mutate(case_when(str_detect(temp_df$group, "scoter") |str_detect(temp_df$group, "scaup") |
                     str_detect(temp_df$group, "ruddy duck") | str_detect(temp_df$group, "merganser") |
                     str_detect(temp_df$group, "long-tailed duck") | str_detect(temp_df$group, "Harlequin duck") |
                     str_detect(temp_df$group, "goldeneye") | str_detect(temp_df$group, "dabbling duck") |
                     str_detect(temp_df$group, "bufflehead") ~ "ducks",
                     str_detect(temp_df$group, "loon") |
                     str_detect(temp_df$group, "gull") | str_detect(temp_df$group, "grebe") |
                     str_detect(temp_df$group, "alcid") | str_detect(temp_df$group, "cormorant") |
                     str_detect(temp_df$group, "goose") ~ "nonduck"
  ))

colnames(temp_df)[35] <- "duckgroup"

temp_df <- temp_df %>% filter(!str_sub(temp_df$group,start = 1, end = 1) == "U")

temp_df <- temp_df %>%
  mutate(case_when(str_detect(temp_df$group, "scoter") | str_detect(temp_df$group, "scaup") |
                     str_detect(temp_df$group, "ruddy duck") | str_detect(temp_df$group, "merganser") |
                     str_detect(temp_df$group, "long-tailed duck") | str_detect(temp_df$group, "Harlequin duck") |
                     str_detect(temp_df$group, "goldeneye") | str_detect(temp_df$group, "dabbling duck") |
                     str_detect(temp_df$group, "bufflehead") | str_detect(temp_df$group, "loon") |
                     str_detect(temp_df$group, "gull") | str_detect(temp_df$group, "grebe") |
                     str_detect(temp_df$group, "alcid") | str_detect(temp_df$group, "cormorant") |
                     str_detect(temp_df$group, "goose") ~ "birds"
  ))

colnames(temp_df)[36] <- "birds"

#Get group totals and species proportions
group.sums <- temp_df %>% group_by(group) %>% summarise(Cam.Count.grp = sum(Count.FF), duckgroup = duckgroup, 
                                                        Count.FF = Count.FF, birds = birds, group = group)
group.sums <- group.sums %>% group_by(duckgroup) %>% summarise(group = group, duckgroup = duckgroup,
                                                               birds = birds, Cam.Count.grp = Cam.Count.grp, 
                                                               Cam.Count.d = sum(Count.FF), Count.FF = Count.FF)
group.sums <- group.sums %>% group_by(birds) %>% summarise(group = group, duckgroup = duckgroup,
                                                           birds = birds, Cam.Count.d = Cam.Count.d,
                                                           Cam.Count.grp = Cam.Count.grp, 
                                                           Cam.Count.b = sum(Count.FF))
group.sums <- group.sums %>% group_by(group) %>% summarise(duckgroup = paste(unique(duckgroup)), 
                                                           birds = paste(unique(birds)),
                                                           Cam.Count.d = round(mean(Cam.Count.d), 0),
                                                           Cam.Count.grp = round(mean(Cam.Count.grp),0),
                                                           Cam.Count.b = round(mean(Cam.Count.b), 0))

group.sums$duckProportions <- group.sums$Cam.Count.grp/group.sums$Cam.Count.d
group.sums$birdProportions <- group.sums$Cam.Count.grp/group.sums$Cam.Count.b

g.allo_df <- g.allo_df %>% left_join(x = g.allo_df, y = group.sums, by = "group")
g.allo_df <- g.allo_df[,-c(35:39)]

g.allo_df <- g.allo_df %>% rename(Cam.Count.FF = Count.FF)
g.allo_df <- g.allo_df %>% rename(Cam.Count.POV = Count.POV)

#Get 4-letter codes for unknowns
unkvec <- c("UNDU", "UNSB")

#Make a list with element names as 4-letter codes of unknowns
grouplist <- list()
for(i in 1:length(unkvec)){
  u <- unkvec[i]
  grouplist[[u]]$grp <- NA
  grouplist[[u]]$dprop <- NA
  grouplist[[u]]$bprop <- NA
}

#Define the 4-letter codes that should be allocated over for each unknown code
grouplist[["UNDU"]]$grp <- c("bufflehead", "dabbling duck", "goldeneye", "Harlequin duck", "long-tailed duck", 
                         "merganser", "ruddy duck", "scaup", "scoter")
grouplist[["UNSB"]]$grp <- unique(g.allo_df$group[str_sub(g.allo_df$group, start = 1, end =1) != "U"])

for(ug in names(grouplist)){
  grouplist[[ug]]$dprop <- as.double(unlist(g.allo_df %>%
                                             filter(group %in% grouplist[[ug]]$grp) %>%
                                             filter(!duplicated(group)) %>%
                                             arrange(group) %>%
                                             select(duckProportions, birdProportions)))
}



for(g in unique(g.allo_df$TransectGroup)){
  for(u in names(grouplist)){
    go <- u %in% g.allo_df$group[g.allo_df$TransectGroup == g]
    if(!go){next}
    camcols <- which(str_detect(colnames(g.allo_df),"Cam"))
    if(!any(g.allo_df[g.allo_df$TransectGroup == g & g.allo_df$SPECIES == u,camcols] > 0)){next}
    nspp <- length(grouplist[[u]]$grp)
    newline <- g.allo_df[g.allo_df$group == u & g.allo_df$TransectGroup == g,]
    newlines1 <- as.data.frame(t(matrix(ncol = nspp, nrow = ncol(newline), data = rep(unlist(newline),nspp))))
    newlines2 <- as.data.frame(t(matrix(ncol = nspp, nrow = ncol(newline), data = rep(unlist(newline),nspp))))
    camcols <- which(str_detect(colnames(g.allo_df),"Cam"))
    newlines1[,camcols] <- apply(newlines1[,camcols], 2,FUN = function(x){
      as.numeric(x)*grouplist[[u]]$dprop
    })
    newlines2[,camcols] <- apply(newlines2[,camcols], 2,FUN = function(x){
      as.numeric(x)*grouplist[[u]]$bprop
    })
    colnames(newlines1) <- colnames(g.allo_df)
    colnames(newlines2) <- colnames(g.allo_df)
    g.allo_df <- rbind(g.allo_df, newlines1, newlines2)
  }
}


g.allo_df <- g.allo_df %>% rename(Count.FF = Cam.Count.FF)
g.allo_df <- g.allo_df %>% rename(Count.POV = Cam.Count.POV)

g.allo_df$SPECIES[g.allo_df$SPECIES == "UNGU"] <- "GULL"

g.allo_df$Count.FF[str_sub(g.allo_df$SPECIES, start = 1, end = 1) == "U"] <- 0
g.allo_df$Count.POV[str_sub(g.allo_df$SPECIES, start = 1, end = 1) == "U"] <- 0

g.allo_df$SPECIES[g.allo_df$SPECIES == "GULL"] <- "UNGU"

g.allo_df <- g.allo_df[,-c(35, 36)]

#To work in model, round counts to whole numbers
g.allo_df$Count.FF <- round(g.allo_df$Count.FF, 0)
g.allo_df$Count.POV <- round(g.allo_df$Count.POV, 0)





