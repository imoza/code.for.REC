## 2013-04-15 16:18:14

## a collection of functions for analysis on e.dagger

## calculate a.dagger

edx <- function(lt, x0=0){
    ## calculate e.dagger at age x0
    lt <- ex.5(lt)
    lt.age <- subset(lt, Age>=x0)
    lt.age$f.x <- lt.age$dx/sum(lt.age$dx)
    ed.out <- round(sum(lt.age$f.x * lt.age$ex5),2)
    return(ed.out)
}

edxs <- function(lt){
    ## calculate e.dagger at age x's
    lt <- ex.5(lt)
    edxs <- apply(as.matrix(lt$Age),1,function(age) edx(lt,age))
    return(edxs)
}

ex.5 <- function(lt){
    ## calculate ex at the midpoint of interval (a, a+1/5/10)
    len <- nrow(lt)
    prop <- lt$ax[-len]/(lt$Age[-1]-lt$Age[-len])
    ex.out <- round(lt$ex[-len]+prop*(lt$ex[-1]-lt$ex[-len]), 2)
    ## the last age interval

    ex.last <- lt$ex[len]
    lt$ex5 <- c(ex.out, ex.last)
    return(lt)
}

## consider the early and late part of e.dagger
get.ad <- function(lt, x0=0){
    ## calculate adag for life table above age a (default is zero) 
    lt <- ex.5(lt)

    lt <- subset(lt, Age>=x0)
    lt$lx <- lt$lx/lt$lx[1]
    lt$dx <- lt$dx/sum(lt$dx)

    y <- edxs(lt)-lt$ex5*(1+log(lt$lx))
    x <- lt$Age

    fit <- approx(x=x,y=y,n=10000)
    xi <- fit$x; yi <- fit$y
    xn <- xi[which.min(abs(yi-0))]
    return(round(xn,3))
}

get.ed.ctry <- function(clt){
    clt[is.na(clt)] <- 0
    yrs <- unique(clt$Year)
    out <- matrix(unlist(apply(as.matrix(yrs), 1, function(i){
        filter(clt, Year==i)%>%get.ed.lt()})), nrow=length(yrs), ncol=6, byrow=T)
    out <- as.data.frame(cbind(Year=yrs, out))
    names(out) <- c("Year", "e0", "ed", "thr", "early", "late", "REC")
    return(out)
}


get.ed.lt <- function(lt){

    ## the old name is get.ed.all

    ad <- get.ad(lt)

    grp <- ifelse(max(lt$Age[-1]-lt$Age[-nrow(lt)])==1,1,5)

    ## calculate edag before (Compression) and after (Expansion) a.dagger
    lt <- lt[!is.na(lt$mx),]
    lt <- ex.5(lt)
    lt$dx <- lt$dx/sum(lt$dx)
    eps <- lt$dx*lt$ex5
    
    ed <- round(sum(eps), 3)

    if(grp==1){
        which.x <- which(lt$Age==floor(ad)); frac <- ad-floor(ad)
    }
    if(grp==5){
        which.x <- rev(which(lt$Age<=floor(ad)))[1]
        age0 <- lt$Age[which.x]
        frac <- (ad-age0)/5
    }

    dx.fr.early <- lt$dx[which.x]*frac
    ex.fr.early <- (lt$ex[which.x+1]-lt$ex[which.x])*frac
    ex.early <- lt$ex[which.x]+0.5*ex.fr.early
    eps.early <- dx.fr.early*ex.early

    ed.early <- sum(eps[1:which.x-1],eps.early)
    ed.late <- ed-ed.early
    
    e0 <- lt$ex[1]
    
    return(round(c(e0, ed, ad, ed.early, ed.late, ed.late/ed.early), 4))
}
