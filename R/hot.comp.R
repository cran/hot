"hot.comp" <-
function(expr.dat,status,set.size,nset)
  {
    p1 <- sum(status==0)
    p2 <- sum(status==1)
    p <- p1+p2
    ngene <- nrow(expr.dat)

    res <- matrix(nr=nset,nc=set.size+1,0)
    m1 <- m2 <- zz <- ynum <- rep(0,set.size)
    v1 <- v2 <- w <- bb <- matrix(nr=set.size,nc=set.size,0)
    index <- rep(-999,set.size)
    
    res <- .Fortran(name="hot",
                    PACKAGE="hot",
                     as.double(expr.dat),
                     as.integer(ngene),
                     as.integer(set.size),
                     as.integer(p),
                     as.integer(p1),
                     as.integer(p2),
                     as.integer(nset),
                     as.double(res),
                     as.double(m1),
                     as.double(m2),
                     as.double(v1),
                     as.double(v2),
                     as.double(w),
                     as.double(bb),
                     as.double(zz),
                     as.double(ynum),
                     as.integer(index))
    
    res <- matrix(nr=nset,nc=set.size+1,res[[8]])
    pval <- 1-pf(res[,set.size+1],df1=set.size,df2=p1+p2-set.size-1)
    res <- cbind(res,pval)
    
    return(res[order(res[,set.size+1],decreasing=TRUE),])
  }

