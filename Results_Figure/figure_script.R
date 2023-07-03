library(here)
library(ggplot2)

results <- read.csv(here("Results_Figure","Results.csv"),header = TRUE)
results <- as.data.frame(results) %>%
  filter(Species != "Black Scoter")

ggplot() + 
  geom_pointrange(data = results %>% filter(Observer != "Camera"), aes(x = Species, y = Count, ymin=Obs.L, ymax=Obs.U, color = Observer, shape = Type), 
                  size = 0.5, alpha = 0.5, position=position_dodge(width=1))+
  geom_point(data = results %>% filter(Observer == "Camera"), 
                  aes(x = Species, y = Count, color = "Camera"), 
                  size = 2, alpha = 0.5, position=position_dodge(width=1))+
  theme_bw() +
  theme(text = element_text(size=14),axis.text.x = element_text(angle=60, hjust=1),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        plot.background = element_blank(),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  labs(title="",x ="", y = "Count or Abundance")+
  scale_y_log10()+
  scale_color_manual(values=c('Camera'='chartreuse2', 'Observer 1'='mediumorchid4', 'Observer 2' = 'darkorange1'))




#For averaging observers
ggplot() + 
  geom_point(data = results, aes(x = Species, y = FF.Count, color = "Camera"), shape=16, color="chartruese2", size=5) +
  geom_point(data = results, aes(x = Species, y = Obs.Ave.Count, color = "Observer"), color = "mediumorchid4", size = 3) +
  geom_pointrange(data = results, aes(x = Species, y = FF.Count, ymin=Obs.Ave.L, ymax=Obs.Ave.U, color = "Estimate"), size = 0.5)+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=60, hjust=1)) +
  labs(title="",x ="", y = "Count or Abundance", color = "Legend")+
  scale_color_manual(name='',
                     breaks=c('Camera', 'Observer', 'Estimate'),
                     values=c('Estimate'='goldenrod3', 'Camera'='chartreuse2', 'Observer'='mediumorchid4'))



