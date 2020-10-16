rm(list=ls())

library(reshape2)
library(ggpubr)
library(tidyverse)

source("edag.r")

simu2lt <- function(mx, x, N, r.sample=c("bi", "dr")){

    nx <- length(x)
    ax <- rep(0.5, nx)
    x <- x+ax

    qx <- 1-exp(-mx)

    lx <- dx <- rep(0, nx)
    lx[1] <- N
    for(i in 1:(nx-1)){
        if(r.sample=="bi") dx[i] <- rbinom(1, lx[i], prob=qx[i])
        if(r.sample=="dr") dx[i] <- round(lx[i]*qx[i])
        lx[i+1] <- lx[i]-dx[i]
    }
    
    Lx <- (lx-dx)+dx*ax
    mx <- dx/Lx

    Tx <- rev(cumsum(rev(Lx)))
    ex <- round(Tx/lx, 4)
    ## ex[n.last] <- ax[n.last] <- 1/mx[n.last]

    mx[is.nan(mx)] <- NA
    ex[is.nan(ex)] <- NA
    if(any(is.nan(ex))){
        posi <- which(is.na(ex))
        n.last <- posi[1]
    }else{
        n.last <- length(x)
    }

    mx <- round(mx, 5)
    Lx <- round(Lx)
    Tx <- round(Tx)
    ex[n.last-1] <- ax[n.last-1] <- 1/mx[n.last-1]
    
    ax[n.last:nx] <- 0
    out <- as.data.frame(cbind(x-.5, lx, dx, mx, Lx, Tx, ex, ax, n.last))
    out[is.na(out)] <- 0
    
    names(out) <- c("Age", "lx", "dx", "mx",  "Lx", "Tx", "ex", "ax", "n.last")
    return(out)
}

## ----------------------------------------

## Siler
siler <- function(a1, b1, r1, a2, r2, a3, b3, r3, t, x,
                  N=1e+6, r.sample="dr"){
    mx <- a1*exp(-r1*t)*exp(-b1*x)+a2*exp(-r2*t)+a3*exp(-r3*t)*exp(b3*x)
    out <- simu2lt(mx, x, N, r.sample)
    return(out)
}

siler.mx.y <- function(a1, b1, r1, a2, r2, a3, b3, r3, t, x){
    a1*exp(-r1*t)*exp(-b1*x)+a2*exp(-r2*t)+a3*exp(-r3*t)*exp(b3*x)
}


## ----------------------------------------


## initial values of parameters
## equvilantly, Japanese female life table 2017 

a1 <- 1e-9
a2 <- 3e-4
a3 <- 1e-6

b1 <- 1e4
b3 <- 0.13

r1 <- r2 <- r3 <- .01

r3.fast <- 0.013


ys <- 0:300
x <- 0:150


## ++++++++++++++++++++++++++++++++++++++++++++++++++
## three scenarios for the Siler model

lt.fast <- NULL
siler.shift <- siler.fast <- siler.slow <- NULL
ad.sl.fast <- ad.sl.slow <- 0
for(y in ys){
    ## shifting model
    siler.y <- siler(a1, b1, r1, a2, r2, a3, b3, r3, y, x)%>%get.ed.lt()

    ## slow progress at the  ages below a.dagger
    siler.y.slow <- siler(a1, b1, r1, a2, r2, a3, b3, ifelse(x<ad.sl.slow, r3.fast, r3), y, x) %>% get.ed.lt()
    ad.sl.slow <- siler.y.slow[3]

    ## faster progress at the ages above a.dagger
    lt.y.fast <- siler(a1, b1, r1, a2, r2, a3, b3, ifelse(x>=ad.sl.fast, r3.fast, r3), y, x)
    lt.fast <- rbind(lt.fast, lt.y.fast) # life tables
    siler.y.fast <- get.ed.lt(lt.y.fast)
    ad.sl.fast <- siler.y.fast[3]

    siler.shift <- rbind(siler.shift, siler.y)
    siler.slow <- rbind(siler.slow, siler.y.slow)
    siler.fast <- rbind(siler.fast, siler.y.fast)
}


simu <- as.data.frame(rbind(cbind(Year=0:1000, siler.fast, Scenario="Faster progress above a+"),
                            cbind(year=0:1000, siler.shift, Scenario="Shifting model"),
                            cbind(year=0:1000, siler.slow, Scenario="Faster progress below a+")))
names(simu) <- c("Year", "e0", "e.dag", "a.dag", "early", "late","REC", "Scenario")

write.csv(simu, "simu.tmp.csv", row.names=F)

simu <- read.csv("simu.tmp.csv")

names(simu)[8] <- "Scenario"
simu$Scenario <- recode(simu$Scenario, "Faster progress above a.dag"="Faster progress above a+",
                        "Shifting model"="Shifting model",
                        "Faster progress below a.dag"="Faster progress below a+")

simu$Scenario <- factor(simu$Scenario, levels=c("Faster progress above a+", "Shifting model", "Faster progress below a+"))


fig.rec <- ggplot(simu, aes(Year, REC, col=Scenario, linetype=Scenario))+xlab("Year (simulated)")+
     geom_smooth(size=.7, se=F)

fig.ed <- ggplot(simu, aes(Year, e.dag, col=Scenario, linetype=Scenario))+
     ylab(expression(e^"+"))+ xlab("Year (simulated)")+
     geom_smooth(size=.7, se=F)

fig.early <- ggplot(simu, aes(Year, early, col=Scenario, linetype=Scenario))+xlab("Year (simulated)")+
     ylab(expression(e[c]^"+"))+ geom_smooth(size=.7, se=F)

fig.late <- ggplot(simu, aes(Year, late, col=Scenario, linetype=Scenario))+xlab("Year (simulated)")+
     ylab(expression(e[e]^"+"))+ geom_smooth(size=.7, se=F)

fig.ad <- ggplot(simu, aes(Year, a.dag, col=Scenario, linetype=Scenario))+ylim(85, 115)+xlab("Year (simulated)")+
     ylab(expression(a^"+"))+ geom_smooth(size=.7, se=F)

fig.e0 <- ggplot(simu, aes(Year, e0, col=Scenario, linetype=Scenario))+ylim(85, 115)+xlab("Year (simulated)")+
     ylab(expression(e^"o"))+ geom_smooth(size=.7, se=F)

ggarrange(fig.e0, fig.ad, fig.ed, fig.rec, fig.early, fig.late, nrow=3, ncol=2, common.legend=T, labels=c("A", "B", "C", "D", "E", "F"), legend="right")

ggsave("./Figs/Fig 3.jpg", height=8, width=8)

