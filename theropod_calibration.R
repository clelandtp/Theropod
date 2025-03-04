library(tidyverse)

#Set PATHNAME to folder with myo- and cah2- files
myo_files<- list.files(path = "PATHNAME",
                           pattern = "myo")

cah2_files<- list.files(path = "PATHNAME",
                       pattern = "cah2")


myo <- readr::read_tsv(myo_files, id = "file_name")
cah2 <- readr::read_tsv(cah2_files, id = "file_name")

cah2 <- read_csv('20240924_cah2_275nM_0-5MIT_0-3ul_0SID_full-filtered.csv')
myo <- read_csv('20240924_myo_225nM_0-5MIT_0-3ul_0SID_full-filtered.csv')


myo$setcharge <- ifelse(myo$`m/z`>848.1 & myo$`m/z`<852.5|myo$`m/z`>863 & myo$`m/z`<869, 20, 
                           ifelse(myo$`m/z`>892.6 & myo$`m/z`<898.7|myo$`m/z`>909 & myo$`m/z`<916.5, 19,
                                  ifelse(myo$`m/z`>942.2 & myo$`m/z`<949|myo$`m/z`>960.5 & myo$`m/z`<966, 18,
                                         ifelse(myo$`m/z`>997.5 & myo$`m/z`<1001, 17,
                                                ifelse(myo$`m/z`>807.7 & myo$`m/z`<812.5|myo$`m/z`>822.8 & myo$`m/z`<827.5, 21,
                                                       ifelse(myo$`m/z`>771.1 & myo$`m/z`<776|myo$`m/z`>784.3 & myo$`m/z`<788.5, 22,""))))))

cah2$setcharge <- ifelse(cah2$`m/z`>854.25 & cah2$`m/z`<856.82, 34, 
                           ifelse(cah2$`m/z`>880.1 & cah2$`m/z`<883.85, 33,
                                  ifelse(cah2$`m/z`>829.8 & cah2$`m/z`<832, 35,
                                         ifelse(cah2$`m/z`>907.5 & cah2$`m/z`<909.87, 32,
                                                ifelse(cah2$`m/z`>936.8 & cah2$`m/z`<938.79, 31,"")))))


cah2myo <- full_join(cah2,myo)
cah2myo <- cah2myo %>% 
  arrange(`m/z`)

cah2myo <- cah2myo %>% 
  filter(setcharge>0) %>% 
  mutate(slopecharge = SlopeComb2/as.numeric(setcharge))

density(cah2myo$slopecharge)

idx <- which.max(density(cah2myo$slopecharge)$y)

density(cah2myo$slopecharge)$x[idx]
density(cah2myo$slopecharge)$y[idx]



peak <- cah2myo %>% 
  summarize(densx = density(cah2myo$slopecharge)$x[idx],
            densy = density(cah2myo$slopecharge)$y[idx])

ggplot(cah2myo, aes(x = slopecharge))+geom_density()+theme_bw()+
  geom_point(data=peak, aes(x=densx, y = densy), color = 'purple', size = 3)+
  annotate('text', x = 2.7, y = 0.8, label = round(density(cah2myo$slopecharge)$x[idx],digits = 6))+ylab('Density')+xlab('Slope to Charge Ratio')

ggsave('calibrationfactor.png')
           