#########################################################
##### Modelling Mann-Kendall trend analysis #############
##########################################################

#### Calculating PIN and LMD prevalences by year and by region 
library(dplyr)

st<-read.table("FINAL_Nanopore_data_metadata.csv", header = T, sep = "\t")

st_deno<- st %>%
  filter(wsaf>=0.5, region%in%c("East","North"),px1_coverage=="High") %>%
  group_by(region,year) %>% 
  summarise(denominator=n())

st_deno<-data.frame(st_deno)  
st_deno<- st_deno %>%
  mutate(unique_name=paste(region,year)) 

st_count<- st %>%
  filter(wsaf>=0.5, region%in%c("East","North"),px1_coverage=="High") %>%
  group_by(region,year,px1_haplo) %>%
  summarise(count=n())
st_count<-data.frame(st_count)  
st_count<- st_count %>%
  mutate(unique_name=paste(region,year)) %>% 
  select(unique_name,count,px1_haplo)

st_sum<-merge(x=st_count,y=st_deno, by="unique_name")


st_sum<-st_sum %>%
  mutate(prevalence=100*round((count/denominator),3))

st_sum$year[st_sum$year==2008]<-2014
st_sum$year[st_sum$year==2012]<-2015
st_sum$year[st_sum$year==2004]<-2013


kl<-data.frame(unique_name="East_2004", count=0,
               px1_haplo="PIN",region="East",year=2013,denominator=88,prevalence=0.0) 
st_sum<-rbind(st_sum,kl)


et<- st_sum %>%
  filter(region=="East", year!="2013", px1_haplo=="PIN")

nt<- st_sum %>%
  filter(region=="North", px1_haplo=="PIN")


##################### Performing Kendall Test #######################

library(Kendall)
mk_nt<-MannKendall(nt$prevalence) 
mk_nt<-data.frame(pvalue=mk_nt$sl, tau=mk_nt$tau,region="North")
mk_et<-MannKendall(et$prevalence)
mk_et<-data.frame(pvalue=mk_et$sl, tau=mk_et$tau,region="East")
mk<-rbind(mk_et,mk_nt)

###
tau_nt <- cor.test(nt$prevalence, nt$year, method = "kendall")
tau_et <- cor.test(et$prevalence, et$year, method = "kendall")

trend_line_nt <- predict(loess(nt$prevalence ~ nt$year))
trend_line_et <- predict(loess(et$prevalence ~ et$year))

####

st_sum1<-st_sum %>%
  filter(px1_haplo%in%c("PIN"), year!="NA")

######################### PLotting Haplotype Frequency by year and by region with sample sizes and 
############################ Kendall Test parameters added #####
library(ggplot2)
library(ggpubr)

Fig4a<-ggplot() +
  geom_col(data=st_sum1,aes(x=year, y=prevalence, fill=px1_haplo),position = "identity")+facet_wrap(~region, scales = "free")+theme_pubclean()+
 
  theme(axis.line = element_line(size=1,colour = "black"),panel.border = element_blank(),
        axis.title = element_text(size=28, face = "bold"), axis.text.x = element_text(size = 22,angle = 0,vjust = .5),
        legend.title = element_text(size = 24, face = "bold"),legend.text = element_text(size = 22),strip.background = element_blank(),
        strip.text = element_text(size = 28, face = "bold"),axis.text.y = element_text(size = 22))+
  labs(fill="PfPX1 haplotype ")+ylab("Prevalence of PIN haplotype (%)")+xlab("Year")+
  scale_x_continuous(breaks = c(2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024),
                     labels = c("2004","2008","2012","2016","2017","2018","2019","2020","2021","2022","2023","2024"))+
  guides(color="none", fill="none")+scale_y_continuous(expand = c(0.025,0.025),limits = c(0,100))+
  scale_fill_manual(values=c("grey30","blue"))+
  geom_text(data=st_sum1,aes(x=year,y=0,label=paste("n=",denominator,sep = "")),size=6)+
  geom_text(data=mk,x=2018,y=85,aes(label=paste("P=",round(mk$pvalue,3),",",sep = "","  Tau=",round(mk$tau,2),",  Mann-Kendall Test")) ,size=8) #τ


Fig4a
ggsave("Figure4A_version3.pdf", width =60, height = 31.5, units = "cm", dpi =600, pointsize = 18,plot = Fig4a)

