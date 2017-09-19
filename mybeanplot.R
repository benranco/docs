if(!require(beanplot))
{
  install.packages('beanplot', repos='http://cran.us.r-project.org')
}
 
library(beanplot)


beanplot <-
function (..., bw = "SJ-dpi", kernel = "gaussian", cut = 3, cutmin = -Inf, 
    cutmax = Inf, grownage = 10, what = c(TRUE, TRUE, TRUE, TRUE), 
    add = FALSE, col, axes = TRUE, log = "auto", handlelog = NA, 
    ll = 0.16, wd = NA, maxwidth = 0.8, maxstripline = 0.96, 
    method = "stack", names, overallline = "mean", beanlines = overallline, 
    horizontal = FALSE, side = "no", jitter = NULL, beanlinewd = 2, 
    frame.plot = axes, border = NULL, innerborder = NA, at = NULL, 
    boxwex = 1, ylim = NULL, xlim = NULL, show.names = NA) 
{
##################### New Code ######################
    write("Custom version of beanplot modified by Ben.", stdout())
################## End of New Code ##################

    #internal functions (later on, mlog and mexp will be defined)

##################### Old Code ######################
    # This is the original unmodified mdensityxy function supplied with beanplot.
    mdensityxy <- function(x) {
        if (length(x) > 0) {
            from <- max(cutmin, (min(mlog(x)) - cut * bw))
            to <- min(cutmax, max(mlog(x)) + cut * bw)
            density(mlog(x), bw = bw, kernel = kernel, from = from, 
                to = to)[c("x", "y")]
        }
        else list(x = numeric(), y = numeric())
    }
##################### New Code ######################
    # This is the function written by Ben, intended to replace the original mdensityxy 
    # to enable a unique bw, cut, cutmin and cutmax for each dataset.
    mdensityxy_ben <- function(datasets, bwVect, cutVect, cutminVect, cutmaxVect) {
        densities <- list()
        if (length(datasets) > 0) {
            if (length(bwVect) == 1) 
                bwVect <- rep(bwVect, length(datasets))
            if (length(cutVect) == 1) 
                cutVect <- rep(cutVect, length(datasets))
            if (length(cutminVect) == 1) 
                cutminVect <- rep(cutminVect, length(datasets))
            if (length(cutmaxVect) == 1) 
                cutmaxVect <- rep(cutmaxVect, length(datasets))

            if ( length(datasets) != length(bwVect) || 
                length(datasets) != length(cutVect) ||
                length(datasets) != length(cutminVect) ||
                length(datasets) != length(cutmaxVect) )
            {
                stop("The length of the bw vector, cut vector, cutMin vector and cutMax vector must all equal either 1 or the number of data sets to be plotted.")
            }

            for (i in 1:length(datasets)) {
                x <- datasets[[i]]
                from <- max(cutminVect[i], (min(mlog(x)) - cutVect[i] * bwVect[i]))
                to <- min(cutmaxVect[i], max(mlog(x)) + cutVect[i] * bwVect[i])
                densities[[i]] <- density(mlog(x), bw = bwVect[i], kernel = kernel, from = from, 
                    to = to)[c("x", "y")]
            }
            
        }
        else densities[[i]] <- list(x = numeric(), y = numeric())

        names(densities) <- names(datasets)
        simplify2array(densities)
    }
################## End of New Code ##################

    #get and store function arguments
    args <- match.call()
    mcall <- as.list(args)

    #settings with multiple options
    method <- pmatch(method, c("overplot", "stack", "jitter"))
    if (is.na(method) || method == 0) 
        stop("invalid plotting method")
    beanlines <- pmatch(beanlines, c("mean", "median", "quantiles"))
    if (is.na(beanlines) || beanlines == 0) 
        stop("invalid beanlines")
    overallline <- pmatch(overallline, c("mean", "median"))
    if (is.na(overallline) || overallline == 0) 
        stop("invalid overallline")
    side <- pmatch(side, c("no", "first", "second", "both"))
    if (is.na(side) || side == 0) 
        stop("invalid side")

    #get the groups dataset, we will generate one bean(-side) per group
    #and set the name and position settings
    groups <- beanplot:::getgroupsfromarguments(args)
    groups <- lapply(groups, na.omit)
    n <- length(groups)
    displayn <- if (side == 4) 
        ceiling(n/2)
    else n
    if (n == 0) 
        stop("no data found to beanplot")
    if (missing(names)) {
        if (is.null(base::names(groups))) 
            attr(groups, "names") = 1:displayn
        names <- base::names(groups)
    }
    else {
        attr(groups, "names") <- names
        if (is.na(show.names)) 
            show.names <- TRUE
    }
    if (is.null(at)) {
        at <- 1:displayn
    }
    if ((side == 4) && (length(names) > length(at))) {
        for (i in 1:length(at)) {
            names[i] <- makecombinedname(names[i * 2 - 1], names[i * 
                2])
        }
        length(names) <- length(at)
    }

    #color settings
    combinedpolygons <- ((side == 4) && (length(border) < 2) && 
        (n%%2 == 0))
    if (missing(col)) 
        col <- par("fg")
    if (!is.list(col)) 
        col <- list(col)
    else combinedpolygons <- FALSE
    col <- lapply(col, fixcolorvector)
    col <- rep(col, length.out = n)
	if (!is.null(border))
        border <- rep(border, length.out = n)

    #set the logarithm handling settings
    if (!add && log == "auto") {
        if (seemslog(groups)) {
            log <- "y"
            message("log=\"y\" selected")
        }
        else log <- ""
    }
    if (is.na(handlelog)) 
        if (add && ((horizontal & par()$xlog) || (!horizontal & 
            par()$ylog))) 
            handlelog <- TRUE
        else if (!add && (log != "")) 
            handlelog <- TRUE
        else handlelog <- FALSE
    if (handlelog) {
        mlog <- base::log
        mexp <- base::exp
    }
    else {
        mlog <- function(x) {
            x
        }
        mexp <- mlog
    }

    #generate the necessary data for the density shapes from the group data
    if (!is.numeric(bw)) {
##################### Old Code ######################
        #bw <- mean(sapply(groups, function(x) {
        #    ifelse(length(x) > 1, density(mlog(x), kernel = kernel, 
        #        bw = bw)$bw, NA)
        #}), na.rm = TRUE)
##################### New Code ######################
        bw <- sapply(groups, function(x) {
            ifelse(length(x) > 1, density(mlog(x), kernel = kernel, 
                bw = bw)$bw, NA)
        })
################## End of New Code ##################
        if (is.nan(bw)) 
            bw <- 0.5
    }
##################### Old Code ######################
    #dens <- sapply(groups, mdensityxy)
##################### New Code ######################
    dens <- mdensityxy_ben(groups, bw, cut, cutmin, cutmax)
################## End of New Code ##################
    for (i in 1:n) dens[["y", i]] <- dens[["y", i]] * min(1, 
        length(groups[[i]])/grownage)
    if (is.na(wd)) 
        wd <- maxwidth/max(unlist(dens["y", ]))
    wd2 <- wd * boxwex/2

    #plot windows and axes
    axespars <- lapply(mcall[base::names(mcall) %in% c("xaxt", 
        "yaxt", "las", "cex.axis", "col.axis", "format", "tick", 
        "xaxp", "yaxp")], eval, parent.frame())
    if (!add) {
        if (!is.numeric(xlim)) {
            if (side == 2) 
                xlim <- c(0, displayn)
            else if (side == 3) 
                xlim <- c(1, displayn + 1)
            else xlim <- c(0.5, displayn + 0.5)
        }
        if (!is.numeric(ylim)) 
            ylim <- range(groups, mexp(unlist(dens["x", ])))
        plot.new()
        windowpars <- lapply(mcall[base::names(mcall) %in% c("yaxs", 
            "xaxs")], eval)
        if (horizontal) {
            names(windowpars)[names(windowpars) %in% c("xaxs", 
                "yaxs")] <- rev(names(windowpars)[names(windowpars) %in% 
                c("xaxs", "yaxs")])
            if (log == "y") 
                log <- "x"
            do.call("plot.window", c(list(xlim = ylim, ylim = xlim, 
                log = log), windowpars))
        }
        else {
            do.call("plot.window", c(list(xlim = xlim, ylim = ylim, 
                log = log), windowpars))
        }
        if (frame.plot) 
            box()
        if (axes) 
            do.call("axis", c(list(side = 2 - horizontal), axespars))
    }
    if (axes) {
        if (is.na(show.names)) 
            show.names <- (n > 1)
        if (show.names) 
            do.call("axis", c(list(1 + horizontal, at = at, labels = names), 
                axespars))
    }

    #plot the window contents
    if (what[1]) {
        if (overallline == 2) 
            val <- mexp(median(mlog(unlist(groups))))
        else val <- mexp(mean(mlog(unlist(groups))))
        if (horizontal) 
            abline(v = val, lty = 3)
        else abline(h = val, lty = 3)
    } else {
		val = "value not calculated, overall line was omitted"
	}
    if (what[2]) {
        beanplotpolyshapes(side, dens, at, wd2, combinedpolygons, 
            displayn, n, col, border, horizontal, mlog, mexp)
    }
    if (what[3]) {
        stats = beanplotbeanlines(groups, side, beanlines, beanlinewd, 
            at, boxwex, n, col, horizontal, mlog, mexp)
    } else {
		stats = "not calculated, beanlines were omitted"
	}
    if (what[4]) {
        beanplotscatters(groups, side, method, jitter, dens, 
            at, wd2, boxwex, n, ll, maxstripline, col, horizontal, 
            mlog, mexp)
    }
    if (any(!is.na(innerborder))) {
        beanplotinnerborders(innerborder, at, dens, side, displayn, 
            n, horizontal, mexp)
    }

    #finally, prints labels
    titlepars <- lapply(mcall[base::names(mcall) %in% c("main", 
        "sub", "xlab", "ylab", "cex.main", "col.main", "cex.lab", 
        "col.lab", "cex.sub", "col.sub")], eval, parent.frame())
    do.call("title", titlepars)

    #return generated data that can be used for subsequent calls
    invisible(list(bw = bw, wd = wd, names = names, stats = stats, overall = val, log = log))
}


assignInNamespace('beanplot',beanplot,ns='beanplot')
environment(beanplot)<-asNamespace('beanplot');

