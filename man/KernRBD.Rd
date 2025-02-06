\name{KernRBD}
\alias{KernRBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KernRBD(Ks, Block.IDs, Treatment, n.res = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ks}{
%%     ~~Describe \code{Ks} here~~
}
  \item{Block.IDs}{
%%     ~~Describe \code{Block.IDs} here~~
}
  \item{Treatment}{
%%     ~~Describe \code{Treatment} here~~
}
  \item{n.res}{
%%     ~~Describe \code{n.res} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Ks, Block.IDs, Treatment, n.res = 5000) 
{
    ind.na <- which(is.na(Block.IDs) | is.na(Treatment))
    if (length(ind.na) > 0) {
        Block.IDs <- Block.IDs[-ind.na]
        Treatment <- Treatment[-ind.na]
        for (i in 1:length(Ks)) {
            Ks[[i]] <- Ks[[i]][-ind.na, -ind.na]
        }
    }
    ind <- order(Block.IDs, Treatment)
    Block.IDs <- Block.IDs[ind]
    Treatment <- Treatment[ind]
    for (i in 1:length(Ks)) {
        Ks[[i]] <- Ks[[i]][ind, ind]
    }
    X <- as.matrix(as.data.frame(model.matrix(~-1 + factor(Treatment))))
    rownames(X) <- NULL
    colnames(X) <- NULL
    n <- length(unique(Block.IDs))
    m <- length(unique(Treatment))
    N <- nrow(X)
    for (i in 1:length(Ks)) {
        nK <- Ks[[i]] - (as.matrix(rep(1, N)) \%*\% t((1/N) * (Ks[[i]] \%*\% 
            as.matrix(rep(1, N))))) - (((1/N) * (Ks[[i]] \%*\% 
            as.matrix(rep(1, N)))) \%*\% rep(1, N)) + as.numeric((1/N^2) * 
            (rep(1, N) \%*\% Ks[[i]] \%*\% as.matrix(rep(1, N))))
        Ks[[i]] <- nK/sum(diag(nK))
    }
    if (m == 1) {
        stop("at least two treatment labels are needed.")
    }
    if (m == 2) {
        global.ran.ind <- list()
        for (i in 1:n.res) {
            global.ran.ind[[i]] <- as.numeric(unlist(tapply(1:N, 
                Block.IDs, function(x) x[sample(1:length(x))])))
        }
        T.obs.list <- list()
        T.null.list <- list()
        itembyitem.pvals <- numeric()
        for (i in 1:length(Ks)) {
            T.obs.list[[i]] <- T.obs <- sum(diag(t(X) \%*\% Ks[[i]] \%*\% 
                X))
            T.null.list[[i]] <- T.null <- sapply(global.ran.ind, 
                function(r) sum(diag(t(X[r, ]) \%*\% Ks[[i]] \%*\% 
                  X[r, ])))
            itembyitem.pvals[i] <- (sum(T.null > T.obs) + 1)/(n.res + 
                1)
        }
        T.obs.minP <- min(itembyitem.pvals)
        T.null.tab <- as.data.frame(T.null.list)
        null.itembyitem.pvals.list <- list()
        for (i in 1:length(Ks)) {
            null.itembyitem.pvals <- numeric()
            for (j in 1:n.res) {
                null.itembyitem.pvals[j] <- (sum(T.null.tab[-j, 
                  i] > T.null.tab[j, i]))/n.res
            }
            null.itembyitem.pvals.list[[i]] <- null.itembyitem.pvals
        }
        T.null.minP <- apply(as.data.frame(null.itembyitem.pvals.list), 
            1, min)
        minP.pval <- (sum(T.null.minP < T.obs.minP) + 1)/(n.res + 
            1)
        names(itembyitem.pvals) <- paste("Kern.", 1:length(Ks), 
            sep = "")
        global.out <- list(itembyitem.pvals = itembyitem.pvals, 
            minP.pval = minP.pval)
        output <- list(Ks = Ks, block.IDs = Block.IDs, treatment = Treatment, 
            global.out = global.out)
    }
    if (m > 2) {
        global.ran.ind <- list()
        for (i in 1:n.res) {
            global.ran.ind[[i]] <- as.numeric(unlist(tapply(1:N, 
                Block.IDs, function(x) x[sample(1:length(x))])))
        }
        T.obs.list <- list()
        T.null.list <- list()
        itembyitem.pvals <- numeric()
        for (i in 1:length(Ks)) {
            T.obs.list[[i]] <- T.obs <- sum(diag(t(X) \%*\% Ks[[i]] \%*\% 
                X))
            T.null.list[[i]] <- T.null <- sapply(global.ran.ind, 
                function(r) sum(diag(t(X[r, ]) \%*\% Ks[[i]] \%*\% 
                  X[r, ])))
            itembyitem.pvals[i] <- (sum(T.null > T.obs) + 1)/(n.res + 
                1)
        }
        T.obs.minP <- min(itembyitem.pvals)
        T.null.tab <- as.data.frame(T.null.list)
        null.itembyitem.pvals.list <- list()
        for (i in 1:length(Ks)) {
            null.itembyitem.pvals <- numeric()
            for (j in 1:n.res) {
                null.itembyitem.pvals[j] <- (sum(T.null.tab[-j, 
                  i] > T.null.tab[j, i]))/n.res
            }
            null.itembyitem.pvals.list[[i]] <- null.itembyitem.pvals
        }
        T.null.minP <- apply(as.data.frame(null.itembyitem.pvals.list), 
            1, min)
        minP.pval <- (sum(T.null.minP < T.obs.minP) + 1)/(n.res + 
            1)
        names(itembyitem.pvals) <- paste("Kern.", 1:length(Ks), 
            sep = "")
        global.out <- list(itembyitem.pvals = itembyitem.pvals, 
            minP.pval = minP.pval)
        pair.grps <- combinations(m, 2, unique(Treatment))
        n.pairs <- nrow(pair.grps)
        pairwise.ran.ind <- list()
        for (i in 1:n.pairs) {
            treat.ind <- is.element(Treatment, pair.grps[i, ])
            ran.arr <- list()
            for (j in 1:n.res) {
                int.arr <- 1:N
                for (k in 1:n) {
                  block.ind <- Block.IDs == unique(Block.IDs)[k]
                  treat.block.ind <- which(treat.ind & block.ind)
                  if (length(treat.block.ind) > 1) {
                    int.arr[treat.block.ind] <- sample(int.arr[treat.block.ind])
                  }
                }
                ran.arr[[j]] <- int.arr
            }
            pairwise.ran.ind[[i]] <- ran.arr
        }
        itembyitem.pvals.list <- list()
        itembyitem.adj.pvals.list <- list()
        T.null.list.list <- list()
        for (i in 1:length(Ks)) {
            T.obs.list <- list()
            T.null.list <- list()
            itembyitem.pvals <- numeric()
            for (j in 1:n.pairs) {
                T.obs.list[[j]] <- T.obs <- sum(diag(t(X) \%*\% 
                  Ks[[i]] \%*\% X))
                T.null.list[[j]] <- T.null <- sapply(pairwise.ran.ind[[j]], 
                  function(r) sum(diag(t(X[r, ]) \%*\% Ks[[i]] \%*\% 
                    X[r, ])))
                itembyitem.pvals[j] <- (sum(T.null > T.obs) + 
                  1)/(n.res + 1)
            }
            T.null.list.list[[i]] <- T.null.list
            T.adj.obs.list <- T.obs.list[order(itembyitem.pvals)]
            T.adj.null.list <- T.null.list[order(itembyitem.pvals)]
            P.ord <- itembyitem.pvals[order(itembyitem.pvals)]
            T.tab <- t(as.data.frame(T.adj.null.list))
            P.tab <- t(apply(T.tab, 1, function(x) rank(x, ties.method = "max"))/n.res)
            P.tab.ext <- rbind(P.tab, rep(1, n.res))
            for (j in n.pairs:1) {
                P.tab.ext[j, ] <- apply(P.tab.ext[j:(j + 1), 
                  ], 2, min)
            }
            Q.tab <- P.tab.ext[1:n.pairs, ]
            itembyitem.adj.pvals <- numeric()
            for (j in 1:n.pairs) {
                itembyitem.adj.pvals[j] <- (sum(Q.tab[j, ] < 
                  P.ord[j]) + 1)/(n.res + 1)
            }
            for (j in 2:n.pairs) {
                itembyitem.adj.pvals[j] <- max(itembyitem.adj.pvals[(j - 
                  1):j])
            }
            itembyitem.adj.pvals <- itembyitem.adj.pvals[order(order(itembyitem.pvals))]
            itembyitem.pvals.list[[i]] <- itembyitem.pvals
            itembyitem.adj.pvals.list[[i]] <- itembyitem.adj.pvals
        }
        itembyitem.pvals.tab <- as.data.frame(itembyitem.pvals.list)
        T.obs.minP <- apply(itembyitem.pvals.tab, 1, min)
        null.itembyitem.pvals.tab.list <- list()
        for (j in 1:n.pairs) {
            null.itembyitem.pvals.list <- list()
            for (i in 1:length(Ks)) {
                null.itembyitem.pvals <- numeric()
                for (k in 1:n.res) {
                  null.itembyitem.pvals[k] <- (sum(T.null.list.list[[i]][[j]][-k] > 
                    T.null.list.list[[i]][[j]][k]) + 1)/(n.res + 
                    1)
                }
                null.itembyitem.pvals.list[[i]] <- null.itembyitem.pvals
            }
            null.itembyitem.pvals.tab.list[[j]] <- as.data.frame(null.itembyitem.pvals.list)
        }
        T.null.minP <- list()
        for (j in 1:n.pairs) {
            T.null.minP[[j]] <- apply(null.itembyitem.pvals.tab.list[[j]], 
                1, min)
        }
        T.null.minP <- as.data.frame(T.null.minP)
        minP.pvals <- numeric()
        for (j in 1:n.pairs) {
            minP.pvals[j] <- (sum(T.null.minP[, j] < T.obs.minP[j]) + 
                1)/(n.res + 1)
        }
        T.adj.obs.list <- as.list(1 - T.obs.minP)
        T.adj.null.list <- as.list(1 - T.null.minP)
        T.adj.obs.list <- T.obs.list[order(minP.pvals)]
        T.adj.null.list <- T.null.list[order(minP.pvals)]
        P.ord <- minP.pvals[order(minP.pvals)]
        T.tab <- t(as.data.frame(T.adj.null.list))
        P.tab <- t(apply(T.tab, 1, function(x) rank(x, ties.method = "max"))/n.res)
        P.tab.ext <- rbind(P.tab, rep(1, n.res))
        for (j in n.pairs:1) {
            P.tab.ext[j, ] <- apply(P.tab.ext[j:(j + 1), ], 2, 
                min)
        }
        Q.tab <- P.tab.ext[1:n.pairs, ]
        minP.adj.pvals <- numeric()
        for (j in 1:n.pairs) {
            minP.adj.pvals[j] <- (sum(Q.tab[j, ] < P.ord[j]) + 
                1)/(n.res + 1)
        }
        for (j in 2:n.pairs) {
            minP.adj.pvals[j] <- max(minP.adj.pvals[(j - 1):j])
        }
        minP.adj.pvals <- minP.adj.pvals[order(order(minP.pvals))]
        conc.pairs <- character()
        for (i in 1:n.pairs) {
            conc.pairs[i] <- paste(pair.grps[i, 1], pair.grps[i, 
                2], sep = "-")
        }
        itembyitem.pvals <- t(as.data.frame(itembyitem.pvals.list))
        rownames(itembyitem.pvals) <- paste("Kern.", 1:length(Ks), 
            sep = "")
        colnames(itembyitem.pvals) <- conc.pairs
        itembyitem.adj.pvals <- t(as.data.frame(itembyitem.adj.pvals.list))
        rownames(itembyitem.adj.pvals) <- paste("Kern.", 1:length(Ks), 
            sep = "")
        colnames(itembyitem.adj.pvals) <- conc.pairs
        names(minP.pvals) <- conc.pairs
        names(minP.adj.pvals) <- conc.pairs
        pairwise.out <- list(itembyitem.pvals = itembyitem.pvals, 
            itembyitem.adj.pvals = itembyitem.adj.pvals, minP.pvals = minP.pvals, 
            minP.adj.pvals = minP.adj.pvals)
        output <- list(Ks = Ks, block.IDs = Block.IDs, treatment = Treatment, 
            global.out = global.out, pairwise.out = pairwise.out)
    }
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
