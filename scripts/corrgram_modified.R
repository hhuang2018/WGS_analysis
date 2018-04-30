original_corrgram <- function (x, type = NULL, order = FALSE, labels, panel = panel.shade, 
          lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
          text.panel = textPanel, label.pos = c(0.5, 0.5), label.srt = 0, 
          cex.labels = NULL, font.labels = 1, row1attop = TRUE, dir = "", 
          gap = 0, abs = FALSE, col.regions = colorRampPalette(c("red", 
                                                                 "salmon", "white", "royalblue", "navy")), cor.method = "pearson", 
          ...) 
{
  if (is.null(order)) 
    order <- FALSE
  if (length(label.pos) < 2) 
    stop("label.pos needs a vector of length 2")
  if (dir == "") {
    if (row1attop) 
      dir <- "left"
    else dir <- "right"
  }
  if (dir == "\\") 
    dir <- "left"
  if (dir == "/") 
    dir <- "right"
  if (ncol(x) < 2) 
    stop("Only one column in the argument to 'corrgram'")
  if (is.matrix(x) && isSymmetric(x) && min(x, na.rm = TRUE) >= 
      -1 - .Machine$double.eps && max(x, na.rm = TRUE) <= 1 + 
      .Machine$double.eps) 
    maybeCorr <- TRUE
  else maybeCorr <- FALSE
  if (is.null(type)) {
    if (maybeCorr) 
      type <- "corr"
    else type <- "data"
  }
  else if (type == "data") {
    if (maybeCorr) 
      warning("This looks like a correlation matrix.")
  }
  else if (type == "cor" | type == "corr") {
    type <- "corr"
    if (!maybeCorr) 
      stop("This is NOT a correlation matrix.")
  }
  else {
    stop("unknown data type in 'corrgram'")
  }
  if (type == "data" & !is.matrix(x)) 
    x <- x[, sapply(x, is.numeric)]
  if (type == "data") 
    cmat <- cor(x, use = "pairwise.complete.obs", method = cor.method)
  else cmat <- x
  cmat <- if (abs) 
    abs(cmat)
  else cmat
  if (order == TRUE | order == "PC" | order == "PCA") {
    x.eigen <- eigen(cmat)$vectors[, 1:2]
    e1 <- x.eigen[, 1]
    e2 <- x.eigen[, 2]
    alpha <- ifelse(e1 > 0, atan(e2/e1), atan(e2/e1) + pi)
    ord <- order(alpha)
    x <- if (type == "data") 
      x[, ord]
    else x[ord, ord]
  }
  else if (order == "OLO") {
    distx <- dist(cmat)
    ss <- seriate(distx, method = "OLO")
    ord <- get_order(ss)
    x <- if (type == "data") 
      x[, ord]
    else x[ord, ord]
  }
  else if (order == "GW") {
    distx <- dist(cmat)
    ss <- seriate(distx, method = "GW")
    ord <- get_order(ss)
    x <- if (type == "data") 
      x[, ord]
    else x[ord, ord]
  }
  else if (order == "HC") {
    distx <- dist(cmat)
    ss <- seriate(distx, method = "HC")
    ord <- get_order(ss)
    x <- if (type == "data") 
      x[, ord]
    else x[ord, ord]
  }
  else if (order != FALSE) {
    stop("Unknown order argument in 'corrgram'.")
  }
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font, srt) {
    text(x, y, txt, cex = cex, font = font, srt = srt)
  }
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if (side%%2 == 1) 
      Axis(x, side = side, xpd = NA, ...)
    else Axis(y, side = side, xpd = NA, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq(along = names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]])) 
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]]))) 
        stop("non-numeric argument to 'corrgram'")
    }
  }
  else if (!is.numeric(x)) 
    stop("non-numeric argument to 'corrgram'")
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  has.diag <- !is.null(diag.panel)
  if (has.diag && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (dir == "left") {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(x)
    if (is.null(labels)) 
      labels <- paste("var", 1:nc)
  }
  else if (is.null(labels)) 
    has.labs <- FALSE
  if (is.null(text.panel)) 
    has.labs <- FALSE
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  else NULL
  main <- if ("main" %in% nmdots) 
    dots$main
  else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) 
      oma[3] <- 6
  }
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))
  for (i in if (dir == "left") 
    1:nc
    else nc:1) for (j in 1:nc) {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ...)
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        if (i == j) {
          if (has.diag) {
            if (type == "data") 
              localDiagPanel(as.vector(x[, i]), NULL, ...)
            else localDiagPanel(NULL, x[i, i], ...)
          }
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            text.panel(label.pos[1], label.pos[2], labels[i], 
                       cex = cex.labels, font = font.labels, srt = label.srt)
          }
        }
        else if (i < j) {
          if (type == "data") 
            localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                           i]), NULL, col.regions, cor.method, ...)
          else localLowerPanel(NULL, NULL, x[j, i], col.regions, 
                               cor.method, ...)
        }
        else {
          if (type == "data") 
            localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                           i]), NULL, col.regions, cor.method, ...)
          else localUpperPanel(NULL, NULL, x[j, i], col.regions, 
                               cor.method, ...)
        }
      }
      else {
        par(new = FALSE)
      }
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}

#######################
modified_corrgram <- function (x, type = NULL, order = FALSE, labels, panel = panel.shade, 
                               lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
                               text.panel = textPanel, label.pos = c(0.5, 0.5), label.srt = 0, 
                               cex.labels = NULL, font.labels = 1, row1attop = TRUE, dir = "", 
                               gap = 0, abs = FALSE, col.regions = colorRampPalette(c("red", 
                                                                                      "salmon", "white", "royalblue", "navy")), cor.method = "pearson", 
                               ...) 
{
  #  whether sort the variables
  if (is.null(order)) 
    order <- FALSE
  
  ## label position sanity check
  if (length(label.pos) < 2) 
    stop("label.pos needs a vector of length 2")
  
  ## direction of matrix
  if (dir == "") {
    if (row1attop) 
      dir <- "left"
    else dir <- "right"
  }
  if (dir == "\\") 
    dir <- "left"
  if (dir == "/") 
    dir <- "right"
  
  ## data sanity check
  if (ncol(x) < 2) 
    stop("Only one column in the argument to 'corrgram'")
  
  ### Calculate correlation or not -- maybeCorrelation matrix --> symetric matrix
  if (is.matrix(x) && isSymmetric(x) && min(x, na.rm = TRUE) >= 
      -1 - .Machine$double.eps && max(x, na.rm = TRUE) <= 1 + 
      .Machine$double.eps) 
    maybeCorr <- TRUE   
  else maybeCorr <- FALSE
  
  ### input type --  data or correlation matrix
  if (is.null(type)) { # default: depending on maybeCorr; (input data)
    if (maybeCorr)     # maybeCorr - true
      type <- "corr"   # then it's a correlation matrix
    else type <- "data" # otherwise it's a data matirx
  }
  else if (type == "data") { # data type
    if (maybeCorr)           # give a warning
      warning("This looks like a correlation matrix.")
  }
  else if (type == "cor" | type == "corr") { # if it's a correlation matrix
    type <- "corr"   
    if (!maybeCorr)  # correlation matrix sanity check
      stop("This is NOT a correlation matrix.")
  }
  else { # unsupported data type
    stop("unknown data type in 'corrgram'")
  }
  
  ###### if it's a data, and is not in MATRIX format, then format into the MATRIX
  if (type == "data" & !is.matrix(x)) 
    x <- x[, sapply(x, is.numeric)]
  
  ##### if it's data matrix, calcualte the pairwise correlation using the desired method: "peason", "spearman", "kendall"
  if (type == "data") 
    cmat <- cor(x, use = "pairwise.complete.obs", method = cor.method)
  else cmat <- x # if it's already a correlation matrix, then define it as cmat
  
  ##### If need the absolute values of the correlations, 
  cmat <- if (abs) 
    abs(cmat)
  else cmat
  
  ##### If order the varialble: order by PCs from PCA
  if (order == TRUE | order == "PC" | order == "PCA") {
    x.eigen <- eigen(cmat)$vectors[, 1:2]
    e1 <- x.eigen[, 1]
    e2 <- x.eigen[, 2]
    alpha <- ifelse(e1 > 0, atan(e2/e1), atan(e2/e1) + pi)
    ord <- order(alpha)
    x <- if (type == "data") 
      x[, ord]
    else x[ord, ord]
  } ## Optimal leaf ordering (OLO) - hierarchical clustering
  else if (order == "OLO") {
    distx <- dist(cmat)
    ss <- seriate(distx, method = "OLO")
    ord <- get_order(ss)
    x <- if (type == "data") 
      x[, ord]
    else x[ord, ord]
  } ### an algorithm developed by Gruvaeus and Wainer (1972) 
    ###  -- the objects at the edge of each cluster are adjacent to that object outside the cluster to which it is nearest.
  else if (order == "GW") {
    distx <- dist(cmat)
    ss <- seriate(distx, method = "GW")
    ord <- get_order(ss)
    x <- if (type == "data") 
      x[, ord]
    else x[ord, ord]
  } ## Hierarchical clustering.
  else if (order == "HC") {
    distx <- dist(cmat)
    ss <- seriate(distx, method = "HC")
    ord <- get_order(ss)
    x <- if (type == "data") 
      x[, ord]
    else x[ord, ord]
  } ### Unsupported order method
  else if (order != FALSE) {
    stop("Unknown order argument in 'corrgram'.")
  }
  
  ######### subfunction for text panel in the plot 
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font, srt) {
    text(x, y, txt, cex = cex, font = font, srt = srt)
  }
  ######## subfucntion for axis position in the plot
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if (side%%2 == 1) 
      Axis(x, side = side, xpd = NA, ...)
    else Axis(y, side = side, xpd = NA, ...)
  }
  ####### Other plot paramters
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  
  ### Sanity check: make sure the input elements of the matrix is numeric
  ###
  ### if the input is not a matrix, 
  ### convert the input elements into the numeric type
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq(along = names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]])) 
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]]))) 
        stop("non-numeric argument to 'corrgram'")
    }
  }
  else if (!is.numeric(x)) 
    stop("non-numeric argument to 'corrgram'")
  
  #### plot panel parameters
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  has.diag <- !is.null(diag.panel)
  if (has.diag && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (dir == "left") {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(x)
    if (is.null(labels)) 
      labels <- paste("var", 1:nc)
  }
  else if (is.null(labels)) 
    has.labs <- FALSE
  if (is.null(text.panel)) 
    has.labs <- FALSE
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  else NULL
  main <- if ("main" %in% nmdots) 
    dots$main
  else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) 
      oma[3] <- 6
  }
  
  ##### Plot correlation matrix
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))
  for (i in if (dir == "left") 
    1:nc
    else nc:1) for (j in 1:nc) {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ...)
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        if (i == j) {
          if (has.diag) {
            if (type == "data") 
              localDiagPanel(as.vector(x[, i]), NULL, ...)
            else localDiagPanel(NULL, x[i, i], ...)
          }
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            text.panel(label.pos[1], label.pos[2], labels[i], 
                       cex = cex.labels, font = font.labels, srt = label.srt)
          }
        }
        else if (i < j) {
          if (type == "data") 
            localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                           i]), NULL, col.regions, cor.method, ...)
          else localLowerPanel(NULL, NULL, x[j, i], col.regions, 
                               cor.method, ...)
        }
        else {
          if (type == "data") 
            localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                           i]), NULL, col.regions, cor.method, ...)
          else localUpperPanel(NULL, NULL, x[j, i], col.regions, 
                               cor.method, ...)
        }
      }
      else {
        par(new = FALSE)
      }
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}
