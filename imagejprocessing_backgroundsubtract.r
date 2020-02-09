library("drc")
require("ggplot2")
require("reshape2")

require("scales")
require("tidyverse")

require("ggthemes")
require("ggdark")

direc <- "~/Dropbox/BondLab/Data/Typhoon Images/200205/rna/fijidata/"

fl<- list.files(direc, pattern = ".csv")

fl.cy3 <- list.files(direc, pattern = "Cy3")
fl.fret <- list.files(direc, pattern = "FRET")

strsplit(fl.fret[1], split="-", )


Lapply <- lapply(paste(direc, fl, sep = ""), read.csv, header=TRUE)
names(Lapply) <- fl
for(i in fl) {
  Lapply[[i]]$Source = i}
comb <- do.call(rbind, Lapply)

comb <- comb %>% separate(Source, c("plate", "channel"), sep = "-")
comb <- comb %>% separate(channel, c("filter", NA), sep = ".cs")

details <- as.data.frame(c("04", '05', '06','NA', '07', '08', '09', '10', '11', '12', '13', '03'))
details
details$max <- c(40,40,20,NA,20,10,2,2,10,100,20,40)

colnames(details) <- c("platerna", "max")
details
i <- 1
# newlist <- list()
for (i in 1:6){
  i <- i*2
  details <- details  
  
  rn1 <- details$platerna[i-1]
  rn2 <- details$platerna[i]
  pl <- paste(rn1, rn2, sep = "_")

  

  tmp1<-
  comb %>% filter(plate==pl) %>%
    filter(Row == "A" | Row =="B" | Row == "C") %>%
    mutate(sample = paste(rn1))
  
  tmp1$conc <- details$max[i-1]/ 2e-3 / (2^tmp1$Column)
  # tmp1%>%head()
  
  
  
  tmp2 <- 
    comb %>% filter(plate==pl) %>%
    filter(Row == "D" | Row =="E" | Row == "F") %>%
    mutate(sample = paste(rn2))
  
  tmp2$conc <- details$max[i]/ 2e-3 / (2^tmp1$Column)
  tmp2%>%head()
  
  
 background <- comb %>% filter(plate == pl) %>%
    filter(Row =="G" & Column == 1:3) %>%
   mutate(sample=paste("background"))
 background$conc <- 0
 
 negative <- comb %>% filter(plate == pl) %>%
   filter(Row =="G" & Column == 4:6) %>%
   mutate(sample=paste("negative"))
 negative$conc <- 0
 
 # for(k in 1:3) {
 # print(mean(negative$Mean[1*k:3*k]))
 # print(negative$)
 #   }
 
 positive <- comb %>% filter(plate == pl) %>%
   filter(Row =="G" & Column == 7:9) %>%
   mutate(sample=paste("positive"))
 positive$conc <- 0
  
 tmp3 <- 
   comb %>% filter(plate==pl) %>%
   filter(Row == "G" & Column ==10:12) %>%
   mutate(sample = "empty")
  
 tmp3$conc <- NA
 
 tmp4 <- comb %>% filter(plate == pl) %>%
   filter(Row =="H") %>%
   mutate(sample=paste("empty"))
 tmp4$conc <- 0
 
 
  

  tmp5 <- rbind(rbind(tmp1, tmp2), rbind(tmp3, rbind(negative, rbind(positive, rbind(background,tmp4)))))

  assign(paste("output",i,sep = "."), tmp5)
  }

dfl <- list(output.2, output.4, output.6, output.8, output.10, output.12)

newcomb <- do.call(rbind, dfl)

# head(newcomb)
# newcomb[grepl("03", newcomb$sample),]
# sel <- rownames(with(newcomb, newcomb[sample=="03" & Column =="1",]))



# newcomb%>% filter(sample =="03" & Column =="1") 

#clean up bad value
# ggplot(newcomb %>% filter(), aes(conc, Mean)) + geom_point() + scale_x_log10() + 
  # facet_wrap(~sample~filter, scales = "free_y")


fretdf <- newcomb %>% filter(filter=="FRET") 
cy3df <-  newcomb %>% filter(filter=="Cy3") 

calcfretdf <- fretdf
calcfretdf$calcfret <- fretdf$Mean / (fretdf$Mean + cy3df$Mean)

# head(calcfretdf)




#create plate lyout diagram
invert_geom_defaults()

ggplot(newcomb%>%filter(filter=="FRET"), aes(Column, fct_rev(Row), alpha=Mean)) +
  geom_point(size = 8) +
  # scale_x_discrete() +
  scale_x_continuous("",breaks = 1:12, position = "bottom", limits = c(0.75, 12.25)) +
  scale_y_discrete("") +
  theme_presentation(base_size = 15, base_family = "Courier") +
  theme(aspect.ratio = 8/12, legend.position = "") + facet_wrap(~plate)


 calcfretdf %>% filter(sample == "positive" | sample=="negative") %>% group_by(plate,sample) %>%
   summarise(meanval = mean(calcfret))
# corrections
# corrections[13:24,] <- corrections
# misclist <- c(04, 04, 06,06, 07, 07, 09, 09, 11, 11, 13, 13, 05,05, NA,  NA, 08, 08, 10, 10, 12, 12, 03,03)
# 
# corrections$rna <- misclist
# 

# calcfretdf$calcfret <- (calcfretdf$calcfret- 0.472)/(0.485-0.472)


testnew <- calcfretdf %>% filter(sample != "03" | Column != 01)
# unique(testnew$sample)



testnew %>% filter(sample=="03" & Column==01)
  # filter(sample == "04" & Column == "02")

bindings <- drm(calcfret~conc, curveid = sample, data=testnew%>%
                  # filter(sample != "03" & Column != 01) %>%
                  filter(
                           sample == "03" | 
                           sample == "04" | 
                           sample == "07" | 
                           sample == "08" | 
                           sample == "09" | 
                           sample == "10" | 
                           sample == "11" | 
                           sample == "13"
                           ),
                fct = LL.4(names = c("b", "min", "max", "Kd"), fixed = c(NA, NA, NA, NA)))
 pl <- plot(bindings)
 plot(bindings)
 
 bindings$curve
 
 head(pl)
  colnames(pl) <- c("conc", "04","07","08","09","10","11","13","03")
head(pl)
 
melt.pl <- melt(pl, id.vars = c("conc"))
colnames(melt.pl) <- c("conc", "sample", "value")
# bindings$coefficients[25:32]
kds <- as.data.frame(round(bindings$coefficients[25:32],0))
colnames(kds) <- c("Kd")
kds
kds$sample <- c("04","07","08","09","10","11","13","03")
head(kds)
head(melt.pl)
# 
plotdf <- calcfretdf %>%filter(sample != "06" & 
                                 sample != "background" &
                                 sample != "empty" & 
                                 sample != "NA" & 
                                 sample != "negative" & 
                                 sample != "positive" )
head(plotdf)
p3 <- ggplot(plotdf, aes(conc, calcfret, colour=sample)) + 
  scale_x_log10(limits=c(0.1,10000), breaks =c(10^(-1:5))) + 
  # geom_line(data=melt.pl, aes(conc, value)) +
  geom_label(data=kds , aes(x = 0.1, y= 0.48, label=paste("Kd =",Kd, "nM")), hjust = 0, colour="white") +
  geom_line(data=melt.pl, aes(conc, value, group=sample), colour="white", alpha=0.8) +
  geom_point(data=plotdf, stat="summary") + 
  geom_errorbar(data=plotdf, stat="summary", width=0.05) +
  labs(title = "FRET binding", y = "Arbitrary FRET Units", colour="RNA (BAJ_)") + 
  scale_y_continuous(limits = c(0.469,0.485), breaks=c(0.470, 0.485)) +
  # coord_cartesian(ylim = c(0,1)) + 
  # theme_classic(base_size = 18) +
  dark_theme_classic(base_size = 18, base_family = "Ubuntu")+
  theme(aspect.ratio = 0.5, legend.position = "") +
   facet_wrap(~sample, ncol = 1, strip.position = "right") 
  
newplate <- plotdf
newcurve <- melt.pl
newkd <- kds

p3 + ggsave("~/Dropbox/BondLab/conferences/Lorne2020/posterfigures/testcolourFRET.svg")

