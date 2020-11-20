rm(list=ls())

library(reshape2)
library(ggpubr)
library(tidyverse)

source("edag.r")

g.fun <- function(lt){

    lt$lx <- lt$lx/100000
    lt$dx <- lt$dx/100000

    lt <- ex.5(lt)

    g <- lt$lx*(edxs(lt)-lt$ex5*(1+log(lt$lx)))

    return(as.data.frame(cbind(Age=lt$Age, g)))
}

rho.fun <- function(lt1, lt2){
    y1 <- lt1$Year[1]
    y2 <- lt2$Year[1]
    rho <- -log(lt2$mx/lt1$mx)/(y2-y1)
    return(as.data.frame(cbind(Year=y1+5, Age=lt1$Age, rho)))
}


shift_axis <- function(p, y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(y=y)
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  p + annotation_custom(grid::grobTree(ax, vp = grid::viewport(y=1, height=sum(ax$height))), 
                        ymax=y, ymin=y) +
    geom_hline(aes(yintercept=y), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank())

}


## ========================================


library(HMDHFDplus)

## readHMD is a function of the package HMDHFDplus
## Please replace "local.filepath" with your own local filepath.
jpnf <- readHMD("local.filepath/hmd_countries/JPN/STATS/fltper_1x1.txt", fixup=T)%>%filter(Year%in%c(1990, 1995))
rusm <- readHMD("local.filepath/hmd_countries/RUS/STATS/mltper_1x1.txt", fixup=T)%>%filterfilter(Year%in%c(1990, 1995))
enwm <- readHMD("local.filepath/hmd_countries/GBRCENW/STATS/mltper_1x1.txt", fixup=T)%>%filter(Year%in%c(2000, 2005))
belm <- readHMD("local.filepath/hmd_countries/BEL/STATS/mltper_1x1.txt", fixup=T)%>%filter(Year%in%c(1960, 1965))

jpn.ed <- get.ed.ctry(jpnf)
rus.ed <- get.ed.ctry(rusm)

enw.ed <- get.ed.ctry(enwm)
bel.ed <- get.ed.ctry(belm)

## Table 1
tab1 <- as.data.frame(round(cbind(t(rus.ed), t(jpn.ed), t(bel.ed), t(enw.ed)),3))
write.csv(tab1, "table1.csv")

## Table 2
tab2 <- round(cbind(log(tab1[c(3,5:7),2]/tab1[c(3,5:7),1])/5,
                    log(tab1[c(3,5:7),4]/tab1[c(3,5:7),3])/5,
                    log(tab1[c(3,5:7),6]/tab1[c(3,5:7),5])/5,
                    log(tab1[c(3,5:7),8]/tab1[c(3,5:7),7])/5),4)

out.tab2 <- rbind(tab2[c(1,4,3,2),], tab2[3,]-tab2[2,])

write.csv(out.tab2, "table2.csv")


## Figure 1

## g(x,t)

jpn5 <- readHMD("local.filepath/hmd_countries/JPN/STATS/fltper_1x5.txt", fixup=T)
rus5 <- readHMD("local.filepath/hmd_countries/RUS/STATS/mltper_1x5.txt", fixup=T)
enw5 <- readHMD("local.filepath/hmd_countries/GBRCENW/STATS/mltper_1x5.txt", fixup=T)
bel5 <- readHMD("local.filepath/hmd_countries/BEL/STATS/mltper_1x1.tx5", fixup=T)


g.jpn <- g.fun(filter(jpn5, Year==1990))
g.rus <- g.fun(filter(rus5, Year==1990))
g.bel <- g.fun(filter(bel5, Year==1960))
g.enw <- g.fun(filter(enw5, Year==2000))

g.4ctry <- rbind(cbind(Country="Russia 1990-95", g.rus),
                 cbind(Country="Japan 1990-95", g.jpn),
                 cbind(Country="Belgium 1960-65", g.bel),
                 cbind(Country="England & Wales 2000-05", g.enw))

g.4ctry$Country  <-  factor(g.4ctry$Country, levels=c("Russia 1990-95",
                                                      "Belgium 1960-65",
                                                      "England & Wales 2000-05",
                                                      "Japan 1990-95"))

ad.jpn <- get.ed.lt(filter(jpn5, Year==1990))[3]
ad.rus <- get.ed.lt(filter(rus5, Year==1990))[3]
ad.bel <- get.ed.lt(filter(bel5, Year==1960))[3]
ad.enw <- get.ed.lt(filter(enw5, Year==2000))[3]

fig.g <- ggplot(g.4ctry, aes(Age, g, col=Country, linetype=Country))+geom_line()+
     geom_hline(yintercept=0, linetype=1, size=.3)+ylim(-77, 20)+
     scale_x_continuous(name="Age", breaks=seq(0, 110, by=10),
                        labels=seq(0, 110, by=10))+
     scale_y_continuous(breaks=seq(-80, 20, by=10),
                        labels = scales::number_format(accuracy = 0.01))+
     geom_vline(xintercept=c(ad.rus, ad.bel, ad.enw, ad.jpn),
                color=c("tomato", "darkgreen", "skyblue3", "violet"),
                linetype=c(1,3,5,2), size=.5, alpha=.8)+
     ylab(expression(g(x)))+
     theme(legend.position=c(.8, .8), legend.title=element_blank(), panel.grid.minor=element_blank())

shift.fig.g <- shift_axis(fig.g, y=0)


## --- rho(x,t) 
r.rus <- rho.fun(filter(rus5, Year==1990-5), filter(rus5, Year==1990+5))
r.jpn <- rho.fun(filter(jpn5, Year==1990-5), filter(jpn5, Year==1990+5))
r.bel <- rho.fun(filter(bel5, Year==1960-5), filter(bel5, Year==1965+5))
r.enw <- rho.fun(filter(enw5, Year==2000-5), filter(enw5, Year==2005+5))

r.4ctry <- rbind(cbind(Country="Russia 1990-95", r.rus),
                 cbind(Country="Belgium 1960-65", r.bel),
                 cbind(Country="England & Wales 2000-05", r.enw),
                 cbind(Country="Japan 1990-95", r.jpn))

r.4ctry$Country  <-  factor(r.4ctry$Country, levels=c("Russia 1990-95",
                                                      "Belgium 1960-65",
                                                      "England & Wales 2000-05",
                                                      "Japan 1990-95"))

fig.rho <- ggplot(r.4ctry, aes(Age, rho, linetype=Country, shape=Country, col=Country))+
     geom_point(size=1.2)+
     stat_smooth(span=.4, fill="gray80", size=.7)+
     ylab(expression(rho(x)))+
     ## ylab(expression(paste("Age-specific rates of mortality decline  ", rho(x))))+
     scale_x_continuous(name="Age", breaks=seq(0, 110, by=10), labels=seq(0, 110, by=10))+
     geom_hline(yintercept=0, linetype=1)+
     geom_vline(xintercept=c(ad.rus, ad.bel, ad.enw, ad.jpn),
                color=c("tomato", "darkgreen", "skyblue3", "violet"),
                linetype=c(1,3,5,2), size=.5, alpha=.8)+
     theme(legend.position=c(.8, .3), panel.grid.minor = element_blank())

ggarrange(shift.fig.g, fig.rho, nrow=2, ncol=1, labels=c("A", "B"), common.legend=T, legend="right")

ggsave("./Figs/Fig 1.jpg", height=8, width=8)


## Figure 2
rusm$Year <- as.factor(rusm$Year)
jpnf$Year <- as.factor(jpnf$Year)
belm$Year <- as.factor(belm$Year)
enwm$Year <- as.factor(enwm$Year)

fig.rus.dx <- ggplot(rusm, aes(Age, dx, col=Year, linetype=Year, group=Year))+ylim(0, 5000)+
    geom_line()+ggtitle("Male, Russia")+ylab("d(x)")+theme(legend.position="")+
    annotate(geom="text", 55, 4850, label=expression(paste(e^{'+'},"=14.70, REC=0.53")), size=3.2)+
    annotate(geom="text", 55, 4450, label=expression(paste(e^{'+'},"=15.86, REC=0.62")), size=3.2)+
    theme(legend.position=c(.2, .9), legend.background = element_blank(), legend.title=element_blank())
fig.jpn.dx <- ggplot(jpnf, aes(Age, dx, col=Year, linetype=Year, group=Year))+ylim(0, 5000)+
    geom_line()+ggtitle("Female, Japan")+ ylab("d(x)")+theme(legend.position="")+
    annotate(geom="text", 55, 4850, label=expression(paste(e^{'+'},"=9.29, REC=0.53")), size=3.2)+
    annotate(geom="text", 55, 4450, label=expression(paste(e^{'+'},"=9.42, REC=0.52")), size=3.2)+
    theme(legend.position=c(.2, .9), legend.background = element_blank(), legend.title=element_blank())
fig.bel.dx <- ggplot(belm, aes(Age, dx, col=Year, linetype=Year, group=Year))+ylim(0, 5000)+
    geom_line()+ggtitle("Male, Belgium")+ ylab("d(x)")+theme(legend.position="")+
    annotate(geom="text", 55, 4850, label=expression(paste(e^{'+'},"=13.56, REC=0.46")), size=3.2)+
    annotate(geom="text", 55, 4450, label=expression(paste(e^{'+'},"=13.12, REC=0.51")), size=3.2)+
    theme(legend.position=c(.2, .9), legend.background = element_blank(), legend.title=element_blank())
fig.enw.dx <- ggplot(enwm, aes(Age, dx, col=Year, linetype=Year, group=Year))+ylim(0, 5000)+
     geom_line()+ggtitle("Male, England & Wales")+ ylab("d(x)")+theme(legend.position="")+
     annotate(geom="text", 55, 4850, label=expression(paste(e^{'+'},"=10.80, REC=0.56")), size=3.2)+
    annotate(geom="text", 55, 4450, label=expression(paste(e^{'+'},"=10.59, REC=0.52")), size=3.2)+
    theme(legend.position=c(.2, .9), legend.background = element_blank(), legend.title=element_blank())

ggarrange(fig.rus.dx, fig.jpn.dx, fig.bel.dx, fig.enw.dx, labels=LETTERS[1:4], nrow=2, ncol=2, common.legend=T)

ggsave("./Figs/Fig 2.jpg", height=8, width=8)



## ++++++++++++++++++++++++++++++++++++++++++++++++++
## period vs. cohort


## selected countries

per.pk <- c("SWE", "FRACNP", "NLD", "GBRCENW", "NOR", "CHE")
coh.pk <- c("SWE", "FRACNP", "NLD", "GBRCENW", "NOR", "CHE")

## period, female

rec.per <- NULL
for(ct in per.pk){
    file.name <- paste("local.filepath/hmd_countries/",ct, "/STATS/fltper_1x1.txt", sep="")
    tmp <- readHMD(file.name, fixup=T)%>%get.ed.ctry()
    ct.name <- countries[which(ct==sh.ctry)]
    rec.per <- rbind(rec.per, cbind(Code=ct, Country=ct.name, tmp))
}

## period, male
rec.per.m <- NULL
for(ct in per.pk){
    file.name <- paste("local.filepath/hmd_countries/",ct, "/STATS/mltper_1x1.txt", sep="")
    tmp <- readHMD(file.name, fixup=T)%>%get.ed.ctry()
    ct.name <- countries[which(ct==sh.ctry)]
    rec.per.m <- rbind(rec.per.m, cbind(Code=ct, Country=ct.name, tmp))
}


## cohort, female
rec.coh <- NULL
for(ct in coh.pk){
    file.name <- paste("local.filepath/hmd_countries/",ct, "/STATS/fltcoh_1x1.txt", sep="")
    tmp <- readHMD(file.name, fixup=T)%>%get.ed.ctry()
    ct.name <- ctry.coh.name[which(ct==ctry.coh)]
    rec.coh <- rbind(rec.coh, cbind(Code=ct, Country=ct.name, tmp))
}

## cohort, male
rec.coh.m <- NULL
for(ct in coh.pk){
    file.name <- paste("local.filepath/hmd_countries/",ct, "/STATS/mltcoh_1x1.txt", sep="")
    tmp <- readHMD(file.name, fixup=T)%>%get.ed.ctry()
    ct.name <- ctry.coh.name[which(ct==ctry.coh)]
    rec.coh.m <- rbind(rec.coh.m, cbind(Code=ct, Country=ct.name, tmp))
}

## selected countries
per.pk <- c("SWE", "FRACNP", "NLD", "GBRCENW", "NOR", "CHE")
coh.pk <- c("SWE", "FRACNP", "NLD", "GBRCENW", "NOR", "CHE")


ggplot(filter(rec.per, Year>=1860), aes(Year, ed, facet=Country))+
    geom_point(shape=2, size=.8, col="blue")+ 
    geom_point(aes(Year, REC*40), shape=19, size=1, col="red")+
    ylab(expression(paste(e^{"+"})))+
    scale_y_continuous(sec.axis = sec_axis(~ ./40, name="REC"))+
    theme(axis.text.y.left = element_text(color = "blue"),
          axis.text.y.right = element_text(color = "red"))+
    facet_wrap(~Country)

ggsave("./Figs/Fig 4.jpg", height=5, width=8)

ggplot(filter(rec.coh, Year>=1860), aes(Year, ed, facet=Country))+
    geom_point(shape=2, size=.8, col="blue")+ 
    geom_point(aes(Year, REC*80), shape=19, size=1, col="red")+
    ylab(expression(paste(e^{"+"})))+xlab("Birth year")+
    scale_y_continuous(sec.axis = sec_axis(~ ./80, name="REC"))+
    theme(axis.text.y.left = element_text(color = "blue"),
          axis.text.y.right = element_text(color = "red"))+
    facet_wrap(~Country)

ggsave("./Fig/Figs/Fig 5.jpg", height=5, width=8)

