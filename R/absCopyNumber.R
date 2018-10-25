cluster.solution <- function(x, alpha.cut, tau.cut) {

    cur.alpha <- x[1,"alpha"]
    cur.tau <- x[1,"tau.def"]
    cur.cnt <- x[1,"count"]
    cur.mse.r <- 1.0/x[1,"mse"]
    clust <- data.frame(numeric(), numeric(), numeric(), numeric(), integer())

    for(i in 2:nrow(x)) {
        next.alpha <- x[i,"alpha"]
        next.tau <- x[i,"tau.def"]
        next.cnt <- x[i,"count"]
        next.mse.r <- 1.0/x[i,"mse"]
        if(abs(next.alpha-cur.alpha)<alpha.cut & abs(next.tau-cur.tau)<tau.cut) {  # same cluster
            # update: weighted average
            cur.wt <- cur.cnt/(cur.cnt+next.cnt)
            next.wt <- next.cnt/(cur.cnt+next.cnt)
            cur.alpha <- (cur.alpha*cur.wt*cur.mse.r + next.alpha*next.wt*next.mse.r) /(cur.wt*cur.mse.r + next.wt*next.mse.r)
            cur.tau <- (cur.tau*cur.wt*cur.mse.r + next.tau*next.wt*next.mse.r) /(cur.wt*cur.mse.r + next.wt*next.mse.r)
            cur.mse.r <- cur.mse.r*cur.wt + next.mse.r*next.wt
            cur.cnt <- cur.cnt + next.cnt
        } else {
            # save the current cluster
            clust <- rbind(clust, data.frame(cur.alpha, cur.tau, cur.tau, cur.mse.r, cur.cnt))

            cur.alpha <- next.alpha
            cur.tau <- next.tau
            cur.cnt <- next.cnt
            cur.mse.r <- next.mse.r
        }
    }
    clust <- rbind(clust, data.frame(cur.alpha, cur.tau, cur.tau, cur.mse.r, cur.cnt))
    colnames(clust) <- c("alpha", "tau", "tau.def", "mse", "count")
    clust[,1:3] <- round(clust[,1:3], 2)
    clust$mse <- 1.0/clust$mse

    #print(x)
    #print(clust)

    clust
}

grid.search.alpha <- function(seg.data, snv.data, alpha.min=0.2, alpha.max=1.0, tau.min=1.5, tau.max=5.0, min.sol.freq=0, min.seg.len=200, qmax=7, lamda=0.5, verbose=FALSE) {

    if(verbose) {
        cat("min.seg.len=",min.seg.len,"\n")
        cat("qmax=",qmax,"\n")
        cat("lamda=",lamda,"\n")			# 1.0: CN only
    }

    orig.seg.dat <- seg.data

    # filtering
    seg.data <- na.omit(seg.data)
    n.seg.1 <- nrow(seg.data)
    seg.data <- seg.data[seg.data[,"eff.seg.len"]>=min.seg.len,]
    n.seg.2 <- nrow(seg.data)
    if(verbose)
        cat("retained and all segments",c(n.seg.2, n.seg.1),"\n")

    gf <- seg.data[,"loc.end"]-seg.data[,"loc.start"]+1
    gf <- gf/sum(gf)
    eta <- 1.0

    r <- seg.data[,"normalized.ratio"]

    max.r.cutoff <- 3.0
    min.r.cutoff <- 1.0/3.0

    outlier.frac.1 <- length(which(r>max.r.cutoff))/length(r)
    outlier.gf.1 <- sum(gf[r>max.r.cutoff])
    if(verbose)
        cat(length(which(r>max.r.cutoff))/length(r),"segments with r>3.0 before rescaling\n")

    outlier2.frac.1 <- length(which(r<min.r.cutoff))/length(r)
    outlier2.gf.1 <- sum(gf[r<min.r.cutoff])
    if(verbose)
        cat(length(which(r<min.r.cutoff))/length(r),"segments with r<0.33 before rescaling\n")

    # assign SNP to each segment after filtering
    snp2seg <- NULL
    fb <- NULL
    het.ind <- NULL
    for(i in 1:nrow(snv.data)) {
        ind <- which(seg.data[,"chrom"]==snv.data[i,"chrom"] & seg.data[,"loc.start"]<=snv.data[i,"position"] & seg.data[,"loc.end"]>=snv.data[i,"position"])
        if(length(ind)==1) {
            if(r[ind] <= max.r.cutoff & r[ind] >= min.r.cutoff) {
                snp2seg <- c(snp2seg, ind)
                fb <- c(fb, snv.data[i,"tumor_var_freq"])
                het.ind <- c(het.ind,i)
            }
        }
    }
    n.snv <- length(het.ind)
    if(verbose)
        cat("# of SNVs used:",length(het.ind),"\n")


    snv.type <- "somatic"

    # weights
    wts <- rep(1, length(r))
    wts[r>max.r.cutoff] <- 0.0
    wts[r<min.r.cutoff] <- 0.0

    wts.cn <- wts/sum(wts)*lamda
    wts.het <- 1.0/length(snp2seg)*(1.0-lamda)

    alpha.grid <- seq(from=alpha.min,to=alpha.max,by=0.05)
    tau.grid <- seq(from=tau.min,to=tau.max,by=0.05)
    search.grid <- expand.grid(alpha.grid=alpha.grid, tau.grid=tau.grid)
    n.grid <- nrow(search.grid)

    search.res <- data.frame(alpha=numeric(), tau=numeric(), tau.def=numeric(), mse=numeric(), alpha0=numeric(), tau0=numeric())
    for(k in 1:n.grid) {
        if(verbose) {
            cat(k)
        }
        alpha0 <- search.grid[k,1]
        tau0 <- search.grid[k,2]

        a.res <- optimize.alpha(alpha0, alpha.min, alpha.max, tau0, wts.cn, r, gf, eta, wts.het, snp2seg, fb, qmax, snv.type)
        if(!is.null(a.res)) {
            search.res <- rbind(search.res, a.res)
        }
    }
    if(verbose) {
        cat("\n")
    }
    colnames(search.res) <- c("alpha", "tau", "tau.def", "mse", "alpha0", "tau0")

    search.res[,c("alpha", "tau", "tau.def")] <- round(search.res[,c("alpha", "tau", "tau.def")],2)
    tmp <- aggregate(search.res[,c("tau.def","mse")], search.res[,c("alpha", "tau")], mean)
    tmp2 <- aggregate(search.res[,1], search.res[,c("alpha", "tau")], length)
    if(verbose)
        print(c(nrow(search.res),sum(tmp2[,3])))

    search.res <- data.frame(tmp, count=tmp2[,3])
    if(verbose)
        cat(nrow(search.res), "unique results\n")

    oo <- order(search.res$mse)
    search.res <- search.res[oo,]
    rownames(search.res) <- NULL

    if(1) {  # do clustering
        search.res <- cluster.solution(search.res, alpha.cut=0.10, tau.cut=0.15)
    }

    # filtering out impossible solutions
    if(min.sol.freq==0)
        min.sol.freq <- 0.05*sum(search.res[,"count"])
    proper.ind <- which(search.res[,"alpha"]>=alpha.min & search.res[,"alpha"]<=alpha.max & search.res[,"tau.def"]>=tau.min & search.res[,"tau.def"]<=tau.max & search.res[,"count"]>=min.sol.freq)
    if(length(proper.ind)>0) {
        search.res <- search.res[proper.ind,,drop=FALSE]
    }


    if(1) {  # output the best result
        min.ind <- which(search.res$mse==search.res$mse[1])
        #cat("# of optimal parameters:",length(min.ind),"\n")

        alpha.opt <- search.res[1,"alpha"]
        tau.opt <- search.res[1,"tau.def"]

        tmp <- (r*(alpha.opt*tau.opt + (1.0-alpha.opt)*2) - (1.0-alpha.opt)*2)/alpha.opt
        qs <- round(tmp)
        qs[qs<0] <- 0
        qs[qs>qmax] <- qmax
        rhat <- (alpha.opt*qs + (1.0-alpha.opt)*2)/(alpha.opt*tau.opt + (1.0-alpha.opt)*2)
        abs.cn <- data.frame(seg.data[,1:4], r=r, rhat=rhat, CN=qs)

        q.het <- qs[snp2seg]
        if(snv.type=="somatic") {
            tmp <- ((alpha.opt*q.het+(1.0-alpha.opt)*2)*fb)/alpha.opt
        } else if(snv.type=="germline") {
            tmp <- ((alpha.opt*q.het+(1.0-alpha.opt)*2)*fb - (1.0-alpha.opt))/alpha.opt
        }
        mg <- round(tmp)
        mg[mg<0] <- 0
        tmp.ind <- which(mg>q.het)
        if(length(tmp.ind)>0) {
            mg[tmp.ind] <- q.het[tmp.ind]
        }
        abs.mg <- data.frame(snv.data[het.ind,],multiplicity=mg)

    }

    res.list <- list(searchRes=search.res, absCN=abs.cn, absSNV=abs.mg, orig.seg.dat=orig.seg.dat, seg.dat=seg.data, orig.snv.dat=snv.data, snv.dat=snv.data[het.ind,])

    res.list
}

optimize.alpha <- function(alpha0, alpha.min=0.2, alpha.max=1.0, tau0, wts.cn, r, w, eta, wts.het, snp2seg, fb, qmax, snv.type) {

    # intialization
    alpha <- alpha0
    tau <- tau0

    tmp <- (r*(alpha*tau + (1.0-alpha)*2) - (1.0-alpha)*2)/alpha
    qs <- round(tmp)
    qs[qs<0] <- 0
    qs[qs>qmax] <- qmax
    #tau <- sum(qs*w*eta)/sum(w*eta)
    tau <- sum(qs*w)/sum(w)
    qs.last <- qs

    q.het <- qs[snp2seg]
    if(snv.type=="somatic") {
        tmp <- ((alpha*q.het+(1.0-alpha)*2)*fb)/alpha
    } else if(snv.type=="germline") {
        tmp <- ((alpha*q.het+(1.0-alpha)*2)*fb - (1.0-alpha))/alpha
    }
    mg <- round(tmp)
    mg[mg<0] <- 0
    tmp.ind <- which(mg>q.het)
    if(length(tmp.ind)>0) {
        mg[tmp.ind] <- q.het[tmp.ind]
    }

    cn.df <- data.frame(x=qs, y=r, wt=wts.cn, lqs=2.0, ix=1)
    het.df <- data.frame(x=mg, y=fb, wt=wts.het, lqs=qs[snp2seg], ix=0)
    reg.df <- rbind(cn.df, het.df)

    # iteration
    optim.fail <- 'converged'
    opt.cycle <- 0
    fitEst <- NULL
    alphaEst <- NULL
    tauEst <- NULL
    mdq <- NULL

    repeat {

        opt.cycle <- opt.cycle+1
        if(opt.cycle>100) {
            res <- NULL
            break
        }

        reg.df[,"x"] <- c(qs, mg)
        reg.df[,"lqs"] <- c(rep(2.0,length(qs)),qs[snp2seg])

        if(snv.type=="somatic") {
            try.tmp <- try(nls.res <- nls(y ~ ix*(alpha*x+(1.0-alpha)*2)/(alpha*tau+(1.0-alpha)*2) + (1.0-ix)*(alpha*x)/(alpha*lqs+(1.0-alpha)*2) , data=reg.df, start=list(alpha=alpha), weights=wt, lower=alpha.min, upper=alpha.max, algorithm="port"), silent=TRUE)

        } else if(snv.type=="germline") {
            try.tmp <- try(nls.res <- nls(y ~ ix*(alpha*x+(1.0-alpha)*2)/(alpha*tau+(1.0-alpha)*2) + (1.0-ix)*(alpha*x+(1.0-alpha))/(alpha*lqs+(1.0-alpha)*2) , data=reg.df, start = list(alpha=alpha), weights=wt, lower=alpha.min, upper=alpha.max, algorithm="port"), silent=TRUE)
        }

        if(class(try.tmp)=="try-error") {
            optim.fail <- 'not_converged'
            res <- NULL
            break
        }

        if(nls.res$convInfo$isConv) {  # proceed

            alpha <- coef(nls.res)[1]
            alpha <- unname(alpha)
            alphaEst <- c(alphaEst, alpha)
            fitEst <- c(fitEst, summary(nls.res)$sigma)

            # update qs
            tmp <- (r*(alpha*tau + (1.0-alpha)*2) - (1.0-alpha)*2)/alpha
            qs <- round(tmp)
            qs[qs<0] <- 0
            qs[qs>qmax] <- qmax

            # update tau
            #tau <- sum(qs*w*eta)/sum(w*eta)
            tau <- sum(qs*w)/sum(w)
            tauEst <- c(tauEst, tau)

            # update mg
            q.het <- qs[snp2seg]
            if(snv.type=="somatic") {
                tmp <- ((alpha*q.het+(1.0-alpha)*2)*fb)/alpha
            } else if(snv.type=="germline") {
                tmp <- ((alpha*q.het+(1.0-alpha)*2)*fb - (1.0-alpha))/alpha
            }
            mg <- round(tmp)
            mg[mg<0] <- 0
            tmp.ind <- which(mg>q.het)
            if(length(tmp.ind)>0) {
                mg[tmp.ind] <- q.het[tmp.ind]
            }

            # compute mean q changes
            mdq <- c(mdq, mean(abs(qs-qs.last)))
            qs.last <- qs

            nL <- length(fitEst)
            if(nL<=1) {
            } else if(abs(fitEst[nL]-fitEst[(nL-1)])/abs(fitEst[(nL-1)])<1e-4) { # converged
                tau.def <- sum(qs*w)/sum(w)
                res <- c(alpha=alpha, tau=tau, tau.def=tau.def, mse=fitEst[nL], alpha0=alpha0, tau0=tau0)
                break
            }

        } else {
            optim.fail <- 'not_converged'
            res <- NULL
            break
        }
    }

    res
}


grid.search.alpha.simple <- function(seg.data, alpha.min=0.2, alpha.max=1.0, tau.min=1.5, tau.max=5.0, min.sol.freq=0, min.seg.len=200, qmax=7, verbose=FALSE) {

    if(verbose) {
        cat("min.seg.len=",min.seg.len,"\n")
        cat("qmax=",qmax,"\n")
    }

    orig.seg.dat <- seg.data

    # filtering
    seg.data <- na.omit(seg.data)
    n.seg.1 <- nrow(seg.data)
    seg.data <- seg.data[seg.data[,"eff.seg.len"]>=min.seg.len,]
    n.seg.2 <- nrow(seg.data)
    if(verbose)
        cat("retained and all segments",c(n.seg.2, n.seg.1),"\n")

    gf <- seg.data[,"loc.end"]-seg.data[,"loc.start"]+1  # vector of genome fraction: spaced length
    gf <- gf/sum(gf)
    eta <- 1.0

    r <- seg.data[,"normalized.ratio"]

    max.r.cutoff <- 3.0
    min.r.cutoff <- 1.0/3.0

    outlier.frac.1 <- length(which(r>max.r.cutoff))/length(r)
    outlier.gf.1 <- sum(gf[r>max.r.cutoff])
    if(verbose)
        cat(length(which(r>max.r.cutoff))/length(r),"segments with r>3.0 before rescaling\n")

    outlier2.frac.1 <- length(which(r<min.r.cutoff))/length(r)
    outlier2.gf.1 <- sum(gf[r<min.r.cutoff])
    if(verbose)
        cat(length(which(r<min.r.cutoff))/length(r),"segments with r<0.33 before rescaling\n")


    # weights
    wts <- rep(1, length(r))
    wts[r>max.r.cutoff] <- 0.0
    wts[r<min.r.cutoff] <- 0.0

    wts.cn <- wts/sum(wts)

    alpha.grid <- seq(from=alpha.min,to=alpha.max,by=0.05)
    tau.grid <- seq(from=tau.min,to=tau.max,by=0.05)
    search.grid <- expand.grid(alpha.grid=alpha.grid, tau.grid=tau.grid)
    n.grid <- nrow(search.grid)

    search.res <- data.frame(alpha=numeric(), tau=numeric(), tau.def=numeric(), mse=numeric(), alpha0=numeric(), tau0=numeric())
    for(k in 1:n.grid) {
        if(verbose) {
            cat(k)
        }
        alpha0 <- search.grid[k,1]
        tau0 <- search.grid[k,2]

        a.res <- optimize.alpha.simple(alpha0, alpha.min, alpha.max, tau0, wts.cn, r, gf, eta, qmax)
        if(!is.null(a.res)) {
            search.res <- rbind(search.res, a.res)
        }
    }
    if(verbose) {
        cat("\n")
    }
    colnames(search.res) <- c("alpha", "tau", "tau.def", "mse", "alpha0", "tau0")

    search.res[,c("alpha", "tau", "tau.def")] <- round(search.res[,c("alpha", "tau", "tau.def")],2)
    tmp <- aggregate(search.res[,c("tau.def","mse")], search.res[,c("alpha", "tau")], mean)
    tmp2 <- aggregate(search.res[,1], search.res[,c("alpha", "tau")], length)
    if(verbose)
        print(c(nrow(search.res),sum(tmp2[,3])))

    search.res <- data.frame(tmp, count=tmp2[,3])
    if(verbose)
        cat(nrow(search.res), "unique results\n")

    oo <- order(search.res$mse)
    search.res <- search.res[oo,]
    rownames(search.res) <- NULL

    if(1) {  # do clustering
        search.res <- cluster.solution(search.res, alpha.cut=0.10, tau.cut=0.15)
    }

    # filtering out impossible solutions
    if(min.sol.freq==0)
        min.sol.freq <- 0.05*sum(search.res[,"count"])
    proper.ind <- which(search.res[,"alpha"]>=alpha.min & search.res[,"alpha"]<=alpha.max & search.res[,"tau.def"]>=tau.min & search.res[,"tau.def"]<=tau.max & search.res[,"count"]>=min.sol.freq)
    if(length(proper.ind)>0) {
        search.res <- search.res[proper.ind,,drop=FALSE]
    }



    if(1) {  # output the best result
        min.ind <- which(search.res$mse==search.res$mse[1])
        #cat("# of optimal parameters:",length(min.ind),"\n")

        alpha.opt <- search.res[1,"alpha"]
        tau.opt <- search.res[1,"tau.def"]

        tmp <- (r*(alpha.opt*tau.opt + (1.0-alpha.opt)*2) - (1.0-alpha.opt)*2)/alpha.opt
        qs <- round(tmp)
        qs[qs<0] <- 0
        qs[qs>qmax] <- qmax
        rhat <- (alpha.opt*qs + (1.0-alpha.opt)*2)/(alpha.opt*tau.opt + (1.0-alpha.opt)*2)
        abs.cn <- data.frame(seg.data[,1:4], r=r, rhat=rhat, CN=qs)

    }

    res.list <- list(searchRes=search.res, absCN=abs.cn, orig.seg.dat=orig.seg.dat, seg.dat=seg.data)

    res.list
}

optimize.alpha.simple <- function(alpha0, alpha.min=0.2, alpha.max=1.0, tau0, wts.cn, r, w, eta, qmax) {

    # intialization
    alpha <- alpha0
    tau <- tau0

    tmp <- (r*(alpha*tau + (1.0-alpha)*2) - (1.0-alpha)*2)/alpha
    qs <- round(tmp)
    qs[qs<0] <- 0
    qs[qs>qmax] <- qmax
    #tau <- sum(qs*w*eta)/sum(w*eta)
    tau <- sum(qs*w)/sum(w)
    qs.last <- qs

    reg.df <- data.frame(x=qs, y=r, wt=wts.cn)

    # iteration
    optim.fail <- 'converged'
    opt.cycle <- 0
    fitEst <- NULL
    alphaEst <- NULL
    tauEst <- NULL
    mdq <- NULL

    repeat {

        opt.cycle <- opt.cycle+1
        if(opt.cycle>100) {
            res <- NULL
            break
        }

        reg.df[,"x"] <- qs

        try.tmp <- try(nls.res <- nls(y ~ (alpha*x+(1.0-alpha)*2)/(alpha*tau+(1.0-alpha)*2), data=reg.df, start=list(alpha=alpha), weights=wt, lower=alpha.min, upper=alpha.max, algorithm="port"), silent=TRUE)

        if(class(try.tmp)=="try-error") {
            optim.fail <- 'not_converged'
            res <- NULL
            break
        }

        if(nls.res$convInfo$isConv) {  # proceed

            alpha <- coef(nls.res)[1]
            alpha <- unname(alpha)
            alphaEst <- c(alphaEst, alpha)
            fitEst <- c(fitEst, summary(nls.res)$sigma)

            # update qs
            tmp <- (r*(alpha*tau + (1.0-alpha)*2) - (1.0-alpha)*2)/alpha
            qs <- round(tmp)
            qs[qs<0] <- 0
            qs[qs>qmax] <- qmax

            # update tau
            #tau <- sum(qs*w*eta)/sum(w*eta)
            tau <- sum(qs*w)/sum(w)
            tauEst <- c(tauEst, tau)

            # compute mean q changes
            mdq <- c(mdq, mean(abs(qs-qs.last)))
            qs.last <- qs

            nL <- length(fitEst)
            if(nL<=1) {
            } else if(abs(fitEst[nL]-fitEst[(nL-1)])/abs(fitEst[(nL-1)])<1e-4) { # converged
                tau.def <- sum(qs*w)/sum(w)
                res <- c(alpha=alpha, tau=tau, tau.def=tau.def, mse=fitEst[nL], alpha0=alpha0, tau0=tau0)
                break
            }

        } else {
            optim.fail <- 'not_converged'
            res <- NULL
            break
        }
    }

    res
}



#' %% ~~ function to do. ~~ Convert copy ratios into abosolute copy numbers
#' when given the purity and ploidy estimates
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ This
#' function converts copy ratios into abosolute copy numbers when the tumor
#' purity and ploidy are provided by the user. The tumor purity and ploidy
#' parameter should be taken from the result after running the "run.absCNSeq"
#' function.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param seg.data %% ~~Describe \code{seg.data} here~~ a data frame with five
#' columns: "chrom", "loc.start", "loc.end", "eff.seg.len", "normalized.ratio".
#' @param alpha.opt %% ~~Describe \code{alpha.opt} here~~ tumor purity estimate
#' @param tau.opt %% ~~Describe \code{tau.opt} here~~ tumor ploidy estimate
#' @param qmax %% ~~Describe \code{qmax} here~~ maximum allowed absolute copy
#' number for any segments
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... The input data frame is augmented with two additional
#' columns: rhat (expected copy ratio) and CN (absolute copy number)
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Lei Bao
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @examples
#' 
#' % ##---- Should be DIRECTLY executable !! ----
#' % ##-- ==>  Define data, use random,
#' % ##--	or do  help(data=index)  for the standard data sets.
#' %% seg.data <- read.table("example.cn.txt",header=TRUE,sep="\t")
#' %% data(absCNseq)
#' %% seg.CN <- compute.absCN(seg.data, 0.58, 2.18)
#' 
#' library(absCNseq)
#' my.res.list <- run.absCNSeq("example.cn.txt", "example.snv.txt", "myResult", "Sample1", seq.type="WES", min.seg.len=200)
#' seg.CN <- compute.absCN(my.res.list$seg.dat, my.res.list$searchRes[i,"alpha"], my.res.list$searchRes[i,"tau"])  # the i-th solution
#' 
#' % ## The function is currently defined as
#' 
compute.absCN <- function(seg.data, alpha.opt, tau.opt, qmax=7) {

    r <- seg.data[,"normalized.ratio"]
    tmp <- (r*(alpha.opt*tau.opt + (1.0-alpha.opt)*2) - (1.0-alpha.opt)*2)/alpha.opt
    qs <- round(tmp)
    qs[qs<0] <- 0
    qs[qs>qmax] <- qmax
    rhat <- (alpha.opt*qs + (1.0-alpha.opt)*2)/(alpha.opt*tau.opt + (1.0-alpha.opt)*2)
    abs.cn <- data.frame(seg.data[,1:4], r=r, rhat=rhat, CN=qs)

    abs.cn
}




#' %% ~~ function to do. ~~ main function to run absCNseq
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ This
#' wrapper function accepts data files and user specified parameters and runs
#' the absCNseq algorithm.
#' 
#' %% ~~ If necessary, more details than the description above ~~ Please refer
#' to the "example" folder to see the format of the input segmentation file or
#' SNV file.
#' 
#' @param seg.fn %% ~~Describe \code{seg.fn} here~~ The name of the file
#' containing the segmentation data. An example file called "example.cn.txt"
#' can be found under the "example" folder within the tarball.
#' @param snv.fn %% ~~Describe \code{snv.fn} here~~ The name of the file
#' containing a set of somatic single nucleotide variants (SNVs). An example
#' file called "example.snv.txt" can be found under the "example" folder within
#' the tarball.
#' @param res.dir %% ~~Describe \code{res.dir} here~~ The output directory
#' @param smp.name %% ~~Describe \code{smp.name} here~~ Sample name
#' @param seq.type %% ~~Describe \code{seq.type} here~~ Either "WES" (whole
#' exome sequencing) or "WGS" (whole genome sequencing)
#' @param alpha.min %% ~~Describe \code{alpha.min} here~~ The minimum allowed
#' value for tumor purity. Default is 0.20. If you do have the pathologist
#' estimate, set it as the lower bound of the pathologist estimate is usually
#' preferred.
#' @param alpha.max %% ~~Describe \code{alpha.max} here~~ The maximum allowed
#' value for tumor purity. Default is 1.0. If you do have the pathologist
#' estimate, set it as the upper bound of the pathologist estimate is usually
#' preferred.
#' @param tau.min %% ~~Describe \code{tau.min} here~~ The minimum allowed value
#' for tumor ploidy
#' @param tau.max %% ~~Describe \code{tau.max} here~~ The maximum allowed value
#' for tumor ploidy
#' @param min.sol.freq %% ~~Describe \code{min.sol.freq} here~~ A solution
#' should appear at least this many times to be kept. Singleton solutions are
#' usually not trustable. By default (min.sol.freq=0), the program will only
#' retain solutions that cover at least 1 percent of the search space.
#' @param min.seg.len %% ~~Describe \code{min.seg.len} here~~ The minimum
#' length of a segment to be included in computation. The default value is 200
#' bp for WES and 3000 bp for WGS.
#' @param qmax %% ~~Describe \code{qmax} here~~ Maximum allowed absolute copy
#' number for any segments.
#' @param lamda %% ~~Describe \code{lamda} here~~ The relative weight of the
#' segment copy ratio data over the SNV data. Must be a value in(0.0,1.0].
#' @param verbose %% ~~Describe \code{verbose} here~~
#' @return %% ~Describe the value returned a list is returned %% If it is a
#' LIST, use \item{searchRes }{a data frame giving the solution set (purity and
#' ploidy pairs) ranked by their fitting errors } \item{absCN }{a data frame
#' giving the absolute copy number estimates of tumor cells for the top first
#' solution (purity and ploidy pair)} \item{absSNV }{a data frame giving the
#' absolute multiplicity estimates of SNVs for the top first solution (purity
#' and ploidy pair) } \item{orig.seg.dat }{original copy ratio data read from
#' the input segmentation file } \item{seg.dat }{filtered copy ratio data }
#' \item{orig.snv.dat }{original SNV data read from the input SNV file }
#' \item{snv.dat }{filtered SNV data }
#' 
#' %% ...
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Lei Bao
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @examples
#' 
#' 
#' library(absCNseq)
#' my.res.list <- run.absCNSeq("example.cn.txt", "example.snv.txt", "myResult", "Sample1", seq.type="WES", min.seg.len=200)
#' seg.CN <- compute.absCN(my.res.list$seg.dat, my.res.list$searchRes[i,"alpha"], my.res.list$searchRes[i,"tau"])  # the i-th solution
#' plot.absCN(seg.CN, chromnum=1)  # plot chromosome 1
#' 
#' % Add one or more standard keywords, see file 'KEYWORDS' in the
#' % R documentation directory.
#' %\keyword{ main function }
#' %\keyword{ absCNSeq } % __ONLY ONE__ keyword per line
#' 
run.absCNSeq <- function(seg.fn, snv.fn=NULL, res.dir, smp.name, seq.type=c("WES","WGS"), alpha.min=0.2, alpha.max=1.0, tau.min=1.5, tau.max=5.0, min.sol.freq=0, min.seg.len=0, qmax=7, lamda=0.5, verbose=FALSE) {

    dir.create(res.dir, showWarnings = FALSE)

    if(!is.null(seg.fn) & file.exists(seg.fn)) {
        seg.data <- read.table(seg.fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    } else {
        stop("ERROR: The input segmentation file is not provided or does not exist.")
    }

    seq.type <- match.arg(seq.type, c("WES","WGS"))
    if(seq.type=="WES") {
        if(min.seg.len==0)
            min.seg.len <- 200
    } else if(seq.type=="WGS") {
        if(min.seg.len==0)
            min.seg.len <- 3000
    } else {
        stop("ERROR: you must specify the sequencing platform: WES or WGS.")
    }

    # check the format of segmentation file
    seg.in.col <- colnames(seg.data)
    seg.req.col <- c("chrom", "loc.start", "loc.end", "eff.seg.len", "normalized.ratio")
    tmp <- setdiff(seg.req.col, seg.in.col)
    if(length(tmp)>0) {
        stop("ERROR: The input segmentation file must include all the named columns in order: chrom, loc.start, loc.end, eff.seg.len, normalized.ratio.")
    }

    if(is.null(snv.fn)) {
        res.list <- grid.search.alpha.simple(seg.data, alpha.min, alpha.max, tau.min, tau.max, min.sol.freq, min.seg.len, qmax, verbose)
    } else {
        if(file.exists(snv.fn)) {
            snv.data <- read.table(snv.fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
        } else {
            stop("ERROR: The input SNV file does not exist.")
        }

        # check the format of SNV file
        snv.in.col <- colnames(snv.data)
        snv.req.col <- c("chrom",  "position", "tumor_var_freq")
        tmp <- setdiff(snv.req.col, snv.in.col)
        if(length(tmp)>0) {
            stop("ERROR: The input SNV file must include all the named columns in order: chrom, position, tumor_var_freq.")
        }

        res.list <- grid.search.alpha(seg.data, snv.data, alpha.min, alpha.max, tau.min, tau.max, min.sol.freq, min.seg.len, qmax, lamda, verbose)
    }

    res.file <- file.path(res.dir, paste(smp.name, ".RData",sep=""))
    save(res.list, file=res.file)

    res.list
}




#' %% ~~ function to do. ~~ plot estimated integer copy numbers from absCNseq
#' for each chromosome.
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ plot
#' estimated integer copy numbers from absCNseq for each chromosome.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param seg.data %% ~~Describe \code{seg.data} here~~ A result data frame
#' generated by the compute.absCN function. It has five mandatory columns:
#' "chrom", "loc.start", "loc.end", "r", "CN".
#' @param chromnum %% ~~Describe \code{chromnum} here~~ which chromosome to be
#' plotted. Must be a number from 1 to 22
#' @param rawdata %% ~~Describe \code{rawdata} here~~ (optional): raw copy
#' ratio data before segmentation. This data frame has four mandatory columns:
#' chromosome number, start position, end position, and normalized copy ratio
#' of tumor vs. normal sample on log2 scale (log2_ratio). The column names are
#' provided by the user. Default column names are taken from the VARSCAN
#' output.
#' @param chromvar %% ~~Describe \code{chromvar} here~~ The variable name of
#' the chromosome number column for the optional raw data.
#' @param loc.start %% ~~Describe \code{loc.start} here~~ The variable name of
#' the start position column for the optional raw data.
#' @param loc.end %% ~~Describe \code{loc.end} here~~ The variable name of the
#' end position column for the optional raw data.
#' @param normalized.ratio %% ~~Describe \code{normalized.ratio} here~~ The
#' variable name of the log2_ratio column for the optional raw data.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~ Lei Bao
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @examples
#' 
#' % ##---- Should be DIRECTLY executable !! ----
#' % ##-- ==>  Define data, use random,
#' % ##--	or do  help(data=index)  for the standard data sets.
#' 
#' library(absCNseq)
#' my.res.list <- run.absCNSeq("example.cn.txt", "example.snv.txt", "myResult", "Sample1", seq.type="WES", min.seg.len=200)
#' seg.CN <- compute.absCN(my.res.list$seg.dat, my.res.list$searchRes[1,"alpha"], my.res.list$searchRes[1,"tau"])  # the top first solution
#' plot.absCN(seg.CN, chromnum=1)  # plot chromosome 1
#' 
#' % ## The function is currently defined as
#' 
plot.absCN <- function(seg.data, chromnum=1, rawdata=NULL, chromvar='chrom', loc.start='chr_start', loc.end='chr_stop', normalized.ratio='log2_ratio'){

    # check the format of input data
    seg.in.col <- colnames(seg.data)
    seg.req.col <- c("chrom", "loc.start", "loc.end", "r", "CN")
    tmp <- setdiff(seg.req.col, seg.in.col)
    if(length(tmp)>0) {
        stop("ERROR: The input dataset must include all the named columns in order: chrom, loc.start, loc.end, r, CN.")
    }

    if(!(chromnum %in% 1:22)) {
        stop("ERROR: chromosome number must be from 1 to 22.")
    }

    segs=subset(seg.data, chrom==chromnum)
    par(xpd=NA,oma=c(2,1,3,1))

    lty=c(1,1)
    col=c('blue','red')
    lwd=c(4,4)
    ylim=c(0,8)
    legtxt=c('2Xraw copy ratio','estimated integer CN')

    if (!is.null(rawdata)){
        rawdata$chrom=rawdata[,chromvar]
        rawdata$chr_start=rawdata[,loc.start]
        rawdata$chr_stop=rawdata[,loc.end]
        rawdata$log2_ratio=rawdata[,normalized.ratio]

        raw=subset(rawdata, chrom==chromnum)
        plot(raw[,'chr_start'], 2*2^(raw$log2_ratio),  xlab='Position', sub=paste('chr',chromnum,sep=''), ylab='Copy Ratio/Number',col='green', cex=.1, axes=F,
             xlim=c(min(raw[,'chr_start'],na.rm=T),max(raw[,'chr_stop'],na.rm=T)), ylim=ylim)
    }  else {

        plot(segs[,"loc.start"], 2*segs[,"r"],  xlab='Position', sub=paste('chr',chromnum,sep=''), ylab='Copy Ratio/Number',col='green', cex=.1, axes=F,
             xlim=c(min(segs[,"loc.start"],na.rm=T),max(segs[,"loc.end"],na.rm=T)), ylim=ylim)
    }
    box()
    axis(1)
    axis(2,las=2)
    segments(segs[,"loc.start"], 2*segs[,"r"],segs[,"loc.end"], 2*segs[,"r"],col=col[1],lwd=lwd[1],lty=lty[1])
    segments(segs[,"loc.start"], segs[,"CN"],segs[,"loc.end"], segs[,"CN"],col=col[2],lwd=lwd[2],lty=lty[2])

    legend(min(segs[,"loc.start"],na.rm=T),10, legtxt, lty=lty, col=col,lwd=lwd, bty='n',ncol=2,xjust=0)
}
