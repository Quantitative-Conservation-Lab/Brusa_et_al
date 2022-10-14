library(here)
library(ggplot2)

results <- read.csv(here("Results_Figure","Results.csv"),header = TRUE)
results <- as.data.frame(results)

ggplot() + 
geom_point(data = results, aes(x = Species, y = FF.Count, color = "Forward Facing"), shape=16, color="red", size=5) +
geom_point(data = results, aes(x = Species, y = Obs1.Count, color = "Observer"), color = "blue", size = 3) +
geom_pointrange(data = results, aes(x = Species, y = FF.Count, ymin=Obs1.L, ymax=Obs1.U, color = "Estimate"), size = 0.5)+
theme(text = element_text(size=20),axis.text.x = element_text(angle=60, hjust=1)) +
labs(title="",x ="Species", y = "Count or Abundance", color = "Legend")+
scale_color_manual(name='',
                  breaks=c('Forward Facing', 'Observer', 'Estimate'),
                  values=c('Estimate'='black', 'Forward Facing'='red', 'Observer'='blue'))
