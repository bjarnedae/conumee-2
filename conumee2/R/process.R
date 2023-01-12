##### PROCESSING methods #####



#' CNV.fit
#' @description Normalize query sample intensities by fitting intensities to reference set using a linear regression model.
#' @param query \code{CNV.data} object of query sample (multiple samples).
#' @param ref \code{CNV.data} object of reference set.
#' @param anno \code{CNV.anno} object. Use \code{CNV.create_anno} do create.
#' @param intercept logical. Should intercept be considered? Defaults to \code{TRUE}.
#' @param ... Additional parameters (\code{CNV.fit} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The log2 ratio of query intensities versus a linear combination of reference set intensities that best reflects query intensities is calculated (as determined by linear regression). Every query sample in the \code{CNV.analysis} object is compared to the linear combination fo control samples individually. The annotations provided to \code{CNV.fit} are saved within the returned \code{CNV.analysis} object and used for subsequent analysis steps.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d, ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' #x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients from linear regression
#' x@@fit$coef
#'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.fit", function(query, ref, anno, ...) {
    standardGeneric("CNV.fit")
})

#' @rdname CNV.fit
setMethod("CNV.fit", signature(query = "CNV.data", ref = "CNV.data", anno = "CNV.anno"),
          function(query, ref, anno, intercept = TRUE, reduce_noise = FALSE, perc_cpgs = 0.2) {
            if (ncol(query@intensity) == 0)
              stop("query intensities unavailable, run CNV.load")
            if (ncol(ref@intensity) == 0)
              stop("reference set intensities unavailable, run CNV.load")

            if (ncol(query@intensity) != 1)
              message("using multiple query samples")
            if (ncol(ref@intensity) == 1)
              warning("reference set contains only a single sample. use more samples for better results.")

            p <- unique(names(anno@probes))  # ordered by location
            if (!all(is.element(p, rownames(query@intensity))))
              stop("query intensities not given for all probes.")
            if (!all(is.element(p, rownames(ref@intensity))))
              stop("reference set intensities not given for all probes.")

            object <- new("CNV.analysis")
            object@date <- date()
            object@fit$args <- list(intercept = intercept)

            object@anno <- anno

            object@fit$coef <- data.frame(matrix(ncol = 0, nrow = ncol(ref@intensity)))
            object@fit$ratio <- data.frame(matrix(ncol = 0, nrow = length(p)))
            for (i in 1:ncol(query@intensity)) {

              message(paste(colnames(query@intensity)[i]), " (",round(i/ncol(query@intensity)*100, digits = 3), "%", ")", sep = "")
              r <- cor(query@intensity[p, ], ref@intensity[p, ])[i, ] < 0.99
              if (any(!r)) message("query sample seems to also be in the reference set. not used for fit.")
              if (intercept) {
                ref.fit <- lm(y ~ ., data = data.frame(y = log2(query@intensity[p,i]), X = log2(ref@intensity[p, r])))
              } else {
                ref.fit <- lm(y ~ . - 1, data = data.frame(y = log2(query@intensity[p,i]), X = log2(ref@intensity[p, r])))
              }
              object@fit$coef <- cbind(object@fit$coef,as.numeric(ref.fit$coefficients[-1]))

              ref.predict <- predict(ref.fit)
              ref.predict[ref.predict < 1] <- 1

              object@fit$ratio <- cbind(object@fit$ratio, log2(query@intensity[p,i]) - ref.predict[p])
            }


            colnames(object@fit$coef) <- colnames(query@intensity)
            rownames(object@fit$coef) <- colnames(ref@intensity)
            colnames(object@fit$ratio) <- colnames(query@intensity)
            rownames(object@fit$ratio) <- p

            object@fit$noise <- as.numeric()
            for (i in 1:ncol(query@intensity)) {
              object@fit$noise <- c(object@fit$noise, sqrt(mean((object@fit$ratio[-1,i] - object@fit$ratio[-nrow(object@fit$ratio),i])^2,na.rm = TRUE)))
            }

            names(object@fit$noise) <- colnames(query@intensity)
            return(object)
          })


#' CNV.bin
#' @description Combine single probe intensitiy values into predefined bins.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.bin} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per bin is calculated. Bins are defined using \code{CNV.create_anno}. A value by which all probe and bin intensity values are shifted in subsequent analysis steps is calculated by minimizing the median absolute deviation from all bins to zero (ideally shifting the copy-neutral state to 0).
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.bin", function(object, ...) {
    standardGeneric("CNV.bin")
})

#' @rdname CNV.bin
setMethod("CNV.bin", signature(object = "CNV.analysis"), function(object) {
  if (length(object@fit) == 0)
    stop("fit unavailable, run CNV.fit")

  o1 <- as.matrix(findOverlaps(query = object@anno@bins, subject = object@anno@probes))
  o2 <- data.frame(bin = names(object@anno@bins)[o1[, "queryHits"]],
                   probe = names(object@anno@probes)[o1[, "subjectHits"]], stringsAsFactors = FALSE)

  object@bin$ratio <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@bin$variance <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@bin$shift <- as.numeric()
  for (i in 1:ncol(object@fit$ratio)) {

    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    object@bin$ratio[[i]] <- sapply(split(object@fit$ratio[o2[, "probe"],i], o2[,"bin"]),
                                    median, na.rm = TRUE)[names(object@anno@bins)]

    object@bin$variance[[i]] <- sapply(split(object@fit$ratio[o2[, "probe"],i], o2[,"bin"]),
                                       var, na.rm = TRUE)[names(object@anno@bins)]

    object@bin$shift <- c(object@bin$shift, optim(0, function(s) median(abs(object@bin$ratio[[i]] -s),na.rm = TRUE),
                                                  method = "Brent", lower = -100, upper = 100)$par)
  }

  names(object@bin$shift) <- colnames(object@fit$ratio)
  names(object@bin$ratio) <- colnames(object@fit$ratio)
  names(object@bin$variance) <- colnames(object@fit$ratio)

  return(object)
})


#' CNV.detail
#' @description Combine single probe values within detail regions.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.detail} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per detail region is calculated. Detail regions are defined using \code{CNV.create_anno(detail_bed=)}
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detail", function(object, ...) {
    standardGeneric("CNV.detail")
})

#' @rdname CNV.detail
setMethod("CNV.detail", signature(object = "CNV.analysis"), function(object) {
  if (length(object@fit) == 0)
    stop("fit unavailable, run CNV.fit")
  # if(length(object@bin) == 0) stop('bin unavailable, run CNV.bin')

  if (length(object@anno@detail) == 0) {
    message("no detail regions provided, define using CNV.create_anno")
  } else {
    d1 <- as.matrix(findOverlaps(query = object@anno@detail, subject = object@anno@probes))
    d2 <- data.frame(detail = values(object@anno@detail)$name[d1[,"queryHits"]], probe = names(object@anno@probes[d1[, "subjectHits"]]),stringsAsFactors = FALSE)


    object@detail$ratio <- vector(mode = "list", length = ncol(object@fit$ratio))
    for (i in 1:ncol(object@fit$ratio)) {
      object@detail$ratio [[i]]<- sapply(split(object@fit$ratio[d2[, "probe"],i],
                                               d2[, "detail"]), median, na.rm = TRUE)[values(object@anno@detail)$name]
    }
    names(object@detail$ratio) <- colnames(object@fit$ratio)
    object@detail$probes <- table(d2[, 1])[values(object@anno@detail)$name]

  }
  return(object)
})


critical.value.ks.test <- function(n, conf, alternative = "two.sided") {

  if(alternative == "one-sided") conf <- 1- (1-conf)*2;

  # if the sample size is large(>50), under the null hypothesis, the absolute value of the difference
  # of the empirical cdf and the theoretical cdf should follow a kolmogorov distribution

  # pdf of the kolmogorov distribution minus the confidence level
  kolmogorov.pdf <- function(x) {
    i <- 1:10^4;
    sqrt(2*pi) / x * sum(exp(-(2*i - 1)^2*pi^2/(8*x^2))) - conf;
  }

  # the root of the function above
  # is the critical value for a specific confidence level multiplied by sqrt(n);
  critical.value <- uniroot(kolmogorov.pdf , lower = 10^(-6), upper = 3)$root / sqrt(n);


  return(critical.value);
}


create.qqplot.fit.confidence.interval <- function(x, distribution = qnorm, conf = 0.95, conf.method = "both", reference.line.method = "quartiles") {

  # remove the NA and sort the sample
  # the QQ plot is the plot of the sorted sample against the corresponding quantile from the theoretical distribution
  sorted.sample <- sort(x[!is.na(x)]);

  # the corresponding probabilities, should be (i - 0.5)/n, where i = 1,2,3,...,n
  probabilities <- ppoints(length(sorted.sample));

  # the corresponding quantile under the theoretical distribution
  theoretical.quantile <- distribution(probabilities);

  if(reference.line.method == "quartiles") {
    # get the numbers of the 1/4 and 3/4 quantile in order to draw a reference line
    quantile.x.axis <- quantile(sorted.sample, c(0.25, 0.75));
    quantile.y.axis <- distribution(c(0.25, 0.75));

    # the intercept and slope of the reference line
    b <- (quantile.x.axis[2] - quantile.x.axis[1]) / (quantile.y.axis[2] - quantile.y.axis[1]);
    a <- quantile.x.axis[1] - b * quantile.y.axis[1];
  }
  if(reference.line.method == "diagonal") {
    a <- 0;
    b <- 1;
  }
  if(reference.line.method == "robust") {
    coef.linear.model <- coef(lm(sorted.sample ~ theoretical.quantile));
    a <- coef.linear.model[1];
    b <- coef.linear.model[2];
  }


  # the reference line
  fit.value <- a + b * theoretical.quantile;

  # create some vectors to store the returned values
  upper.pw <- NULL;
  lower.pw <- NULL;
  upper.sim <- NULL;
  lower.sim <- NULL;
  u <- NULL;	# a vector of logical value of whether the probabilities are in the interval [0,1] for the upper band
  l <- NULL;	# a vector of logical value of whether the probabilities are in the interval [0,1] for the lower band

  ### pointwise method
  if (conf.method == "both" | conf.method == "pointwise") {

    # create the numeric derivative of the theoretical quantile distribution
    numeric.derivative <- function(p) {
      # set the change in independent variable
      h <- 10^(-6);
      if (h * length(sorted.sample) > 1) { h <- 1 / (length(sorted.sample) + 1); }
      # the function
      return((distribution(p + h/2) - distribution(p - h/2)) / h);
    }

    # the standard errors of pointwise method
    data.standard.error <- b * numeric.derivative(probabilities) * qnorm(1 - (1 - conf)/2) * sqrt(probabilities * (1 - probabilities) / length(sorted.sample));

    # confidence interval of pointwise method
    upper.pw <- fit.value + data.standard.error;
    lower.pw <- fit.value - data.standard.error;
  }

  ### simultaneous method
  if (conf.method == "both" | conf.method == "simultaneous") {

    # get the threshold value for the statistics---the absolute difference of the empirical cdf and the theoretical cdf
    # Note that this statistics should follow a kolmogorov distribution when the sample size is large

    # the critical value from the Kolmogorov-Smirnov Test
    critical.value <- critical.value.ks.test(length(sorted.sample), conf);

    # under the null hypothesis, get the CI for the probabilities
    # the probabilities of the fitted value under the empirical cdf
    expected.prob <- ecdf(sorted.sample)(fit.value);

    # the probability should be in the interval [0, 1]
    u <- (expected.prob + critical.value) >= 0 & (expected.prob + critical.value) <= 1;
    l <- (expected.prob - critical.value) >= 0 & (expected.prob - critical.value) <= 1;

    # get the corresponding quantiles from the theoretical distribution
    z.upper <- distribution((expected.prob + critical.value)[u]);
    z.lower <- distribution((expected.prob - critical.value)[l]);

    # confidence interval of simultaneous method
    upper.sim <- a + b * z.upper;
    lower.sim <- a + b * z.lower;
  }


  # return the values for constructing the Confidence Bands of one sample QQ plot
  # the list to store the returned values
  returned.values <- list(
    a = a,
    b = b,
    z = theoretical.quantile,
    upper.pw = upper.pw,
    lower.pw = lower.pw,
    u = u,
    l = l,
    upper.sim = upper.sim,
    lower.sim = lower.sim
  );
  return (returned.values);
}


#' CNV.focal
#' @description This optional function provides filtering for diagnostically relevant CNVs (high level amplification or homozygous deletion).
#' @param object \code{CNV.analysis} object.
#' @param conf numeric. This parameter affects the plotted confidence intervals. Which confidence level should be used? Default to \code{0.99}.
#' @param minoverlap integer. The function determines the bins that overlap with the genes of interest. Which minimum number of basepairs should be considered for an overlap? Defaul to \code{1000L}.
#' @param ... Additional parameters (\code{CNV.detailplot} generic, currently not used).
#' @return A \code{CNV.analysis} object with significantly altered bins and genes from the Cancer Gene Census (curated by the Sanger Institute).
#' @details This function should facilitate the detection of diagnostically relevant CNVs that affect single genes. Therefore, a qqplot illustrating the bins' log2-ratios is created for each chromosome arm. In the first step, bins that lie outside the confidence interval are identified and sorted based on their residuals to the confidence curves.
#' Secondly, these bins are overlapped with the Cancer Gene Census to identify diagnostically relevant genes. The resulting bins and genes are sorted in regards to their significance.
#' @examples
#'
#' x <- CNV.focal(object, conf = 0.99, minoverlap = 10000L)
#' x@@detail$cancer_genes
#' x@@detail$significant_bins
#'
#' @author Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.focal", function(object, ...) {
  standardGeneric("CNV.focal")
})

#' @rdname CNV.focal
setMethod("CNV.focal", signature(object = "CNV.analysis"), function(object, conf = 0.99, minoverlap = 1000L,...){

  if(ncol(object@anno@genome) == 2) {
    stop("CNV.focal is not compatible with mouse arrays.")
  }

  data("consensus_cancer_genes_hg19")

  significant.bins <- vector(mode = "list", length = ncol(object@fit$ratio))
  cancer.genes <- vector(mode = "list", length = ncol(object@fit$ratio))

  for (i in 1:ncol(object@fit$ratio)){

    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    bin.ratios <- object@bin$ratio[[i]] - object@bin$shift[i]

    residuals <- as.numeric() #inner loop
    for (j in 1:nrow(object@anno@genome)) {
      first <- IRanges(start = 1, end = object@anno@genome[j,3])
      first <- GRanges(seqnames = rownames(object@anno@genome)[j], first)
      second <- IRanges(start = object@anno@genome[j,3]+1, end = object@anno@genome[j,2])
      second <- GRanges(seqnames = rownames(object@anno@genome)[j], second)

      h.1 <- findOverlaps(query = object@anno@bins, subject = first, type = "within", ignore.strand = TRUE)
      ind.1 <- queryHits(h.1)

      h.2 <- findOverlaps(query = object@anno@bins, subject = second, type = "within", ignore.strand = TRUE)
      ind.2 <- queryHits(h.2)

      names.1 <- names(object@anno@bins[ind.1])
      names.2 <- names(object@anno@bins[ind.2])

      if (length(names.1)>0) {
        c.intervals <- create.qqplot.fit.confidence.interval(bin.ratios[names.1], distribution = qnorm, conf = conf, conf.method = "pointwise")
        qq.plot <- qqnorm(bin.ratios[names.1], plot.it = FALSE)
        y.c <- qq.plot$y[order(qq.plot$x)]
        upper.outliers <- which(y.c>c.intervals$upper.pw)
        upper.residuals <- y.c[upper.outliers] - c.intervals$upper.pw[upper.outliers]
        lower.outliers <- which(y.c<c.intervals$lower.pw)
        lower.residuals <- y.c[lower.outliers] - c.intervals$lower.pw[lower.outliers]
        residuals.1 <- c(abs(upper.residuals), abs(lower.residuals))
      }

      if (length(names.2)>0) {
        c.intervals <- create.qqplot.fit.confidence.interval(bin.ratios[names.2], distribution = qnorm, conf = conf, conf.method = "pointwise")
        qq.plot <- qqnorm(bin.ratios[names.2], plot.it = FALSE)
        y.c <- qq.plot$y[order(qq.plot$x)]
        upper.outliers <- which(y.c>c.intervals$upper.pw)
        upper.residuals <- y.c[upper.outliers] - c.intervals$upper.pw[upper.outliers]
        lower.outliers <- which(y.c<c.intervals$lower.pw)
        lower.residuals <- y.c[lower.outliers] - c.intervals$lower.pw[lower.outliers]
        residuals.2 <- c(abs(upper.residuals), abs(lower.residuals))
      }

      my_out <- c(residuals.1, residuals.2)
      residuals <- c(residuals, my_out)
      residuals.1 <- NULL
      residuals.2 <- NULL
    }

    outliers <- sort(residuals, decreasing = TRUE)
    significant.bins[[i]] <- object@anno@bins[names(outliers)]

    h.cancer_genes <- findOverlaps(query = object@anno@bins[names(outliers)], subject = consensus_cancer_genes_hg19, minoverlap = minoverlap)
    significant.genes <- consensus_cancer_genes_hg19$SYMBOL[unique(subjectHits(h.cancer_genes))]
    cancer.genes[[i]] <- significant.genes

  }

  names(significant.bins) <- colnames(object@fit$ratio)
  names(cancer.genes) <- colnames(object@fit$ratio)
  object@detail$significant_bins <- significant.bins
  object@detail$cancer_genes <- cancer.genes

  return(object)

})


#' @import DNAcopy
NULL

#' CNV.segment
#' @description Segment bin values (wrapper of \code{DNAcopy} package).
#' @param object \code{CNV.analysis} object.
#' @param alpha See details. Defaults to 0.001.
#' @param nperm See details. Defaults to 50000.
#' @param min.width See details. Defaults to 5.
#' @param undo.splits See details. Defaults to 'sdundo'.
#' @param undo.SD See details. Defaults to 2.2.
#' @param verbose See details. Defaults to 0.
#' @param ... Additional parameters supplied to the \code{segment} method of the \code{DNAcopy} package.
#' @return \code{CNV.analysis} object.
#' @details This method is a wrapper of the CNA, segment, segments.summary and segments.p methods of the DNAcopy package. Please refer to the respective man pages for more detailed information. The default parameters of \code{CNV.segment} override some of the default parameters of segment and are optimized for 450k data CNV analysis.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.segment", function(object, ...) {
    standardGeneric("CNV.segment")
})

#' @rdname CNV.segment
setMethod("CNV.segment", signature(object = "CNV.analysis"), function(object,
                                                                      alpha = 0.001, nperm = 50000, min.width = 5, undo.splits = "sdundo",
                                                                      undo.SD = 2.2, verbose = 0, ...) {
  if(length(object@fit) == 0){
    stop('fit unavailable, run CNV.fit')
  }
  if (length(object@bin) == 0){
    stop("bin unavailable, run CNV.bin")
  }

  a1 <- formals()
  a2 <- as.list(match.call())[-1]
  object@seg$args <- as.list(sapply(setdiff(unique(names(c(a1, a2))),
                                            c("object", "verbose")), function(an) if (is.element(an, names(a2)))
                                              a2[[an]] else a1[[an]], simplify = FALSE))

  object@seg$summary <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@seg$p <- vector(mode = "list", length = ncol(object@fit$ratio))

  for (i in 1:ncol(object@fit$ratio)) {

    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    x1 <- DNAcopy::CNA(genomdat = object@bin$ratio[[i]][names(object@anno@bins)],
                       chrom = as.vector(seqnames(object@anno@bins)), maploc = values(object@anno@bins)$midpoint,
                       data.type = "logratio", sampleid = "sampleid")

    x2 <- DNAcopy::segment(x = x1, weights = 1/object@bin$variance[[i]][names(object@anno@bins)], verbose = verbose, min.width = min.width,
                           nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD,
                           ...)

    object@seg$summary[[i]] <- DNAcopy::segments.summary(x2)
    object@seg$summary[[i]]$chrom <- as.vector(object@seg$summary[[i]]$chrom)
    object@seg$summary[[i]]$ID <- colnames(object@fit$ratio)[i]
    object@seg$p[[i]] <- DNAcopy::segments.p(x2)
    object@seg$p[[i]]$chrom <- as.vector(object@seg$p[[i]]$chrom)
  }
  names(object@seg$summary) <- colnames(object@fit$ratio)
  names(object@seg$p) <- colnames(object@fit$ratio)
  return(object)
})

