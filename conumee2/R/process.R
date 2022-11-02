##### PROCESSING methods #####

#' CNV.fit
#' @description Normalize query sample intensities by fitting intensities to reference set using a linear regression model.
#' @param query \code{CNV.data} object of query sample (multiple samples).
#' @param ref \code{CNV.data} object of reference set.
#' @param anno \code{CNV.anno} object. Use \code{CNV.create_anno} do create.
#' @param intercept logical. Should intercept be considered? Defaults to \code{TRUE}.
#' @param reduce_noise logical. Should the noise be reduced by excluding the probes with the lowest combined signal intensities among the set of control probes? Defaults to \code{FALSE}. For details see the publication.
#' @param perc_cpgs numeric. What percentage of CpGs should be excluded if \code{reduce_noise} it set to \code{TRUE}? Defaults to \code{0.2}.
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

            if (reduce_noise) {
              message("identifying and excluding most variable probes among the set of control samples") #reduce noise beginning

              c_samples <- ref@intensity
              means <- apply(c_samples, 1, mean)
              st_dev <- apply(c_samples, 1, sd)
              ratio <- st_dev/means
              c_samples$ratio <- ratio
              cpgs_exclude <- rownames(c_samples[order(c_samples$ratio, decreasing = TRUE),][1:(perc_cpgs*nrow(c_samples)),])


              qload <- query
              refload <- ref

              ind_del_query <- which(rownames(qload@intensity) %in% cpgs_exclude) #identifying CpGs that should be excluded

              if (ncol(qload@intensity) == 1) {
                qload@intensity <- as.data.frame(qload@intensity[-ind_del_query,], row.names = rownames(qload@intensity)[-ind_del_query])
              } else {
              qload@intensity <- qload@intensity[-ind_del_query,]
              }

              ind_del_ref <- which(rownames(refload@intensity) %in% cpgs_exclude)
              refload@intensity <- refload@intensity[-ind_del_ref,]

              message ("modifying annotation object for CVN.fit()")
              aobject <- anno

              ind_del_ao <- which(names(aobject@probes) %in% cpgs_exclude)
              aobject@probes <- anno@probes[-ind_del_ao,]
              #changing the annotation object

              message("creating new bins")
              anno.tile <- CNV.create_bins(hg19.anno = aobject@genome, bin_minsize = aobject@args$bin_minsize,
                                           hg19.gap = aobject@gap, hg19.exclude = aobject@exclude)
              message(" - ", length(anno.tile), " bins created")

              message("merging new bins")

              if (is.null(aobject@genome$pg)) {

              aobject@bins <- CNV.merge_bins_mice(hg19.anno = aobject@genome, hg19.tile = anno.tile,
                                                  bin_minprobes = aobject@args$bin_minprobes, hg19.probes = aobject@probes, bin_maxsize = aobject@args$bin_maxsize)
              } else {

              aobject@bins <- CNV.merge_bins(hg19.anno = aobject@genome, hg19.tile = anno.tile,
                                                    bin_minprobes = aobject@args$bin_minprobes, hg19.probes = aobject@probes, bin_maxsize = aobject@args$bin_maxsize)
              }
              message(" - ", length(aobject@bins), " bins remaining")
              message("annotation object finished")
              message(paste(nrow(qload@intensity)," CpGs preserved", sep = ""))

              if (is.null(aobject@genome$pg)) {
              message("getting the gene annotations for each bin")
              o <- findOverlaps(aobject@probes, aobject@bins)
              bin_genes <- sapply(lapply(lapply(split(aobject@probes$genes[queryHits(o)],
                                                      names(aobject@bins)[subjectHits(o)]), unique), sort), paste0, collapse=";")

              c_bins <- aobject@bins
              c_bins$genes <- ""
              c_bins[names(bin_genes)]$genes <- sub(";","", bin_genes)
              aobject@bins <- c_bins
              } else {
              message("getting the gene annotations for each bin")

              o <- findOverlaps(aobject@probes, aobject@bins)
              bin_genes <- sapply(lapply(lapply(split(aobject@probes$genes[queryHits(o)],
                                                      names(aobject@bins)[subjectHits(o)]), unique), sort), paste0, collapse=";")

              aobject@bins$genes <- bin_genes
              }

              object <- new("CNV.analysis")
              object@date <- date()
              object@fit$args <- list(intercept = intercept)

              object@anno <- aobject

              p <- unique(names(aobject@probes))
              object@fit$coef <- data.frame(matrix(ncol = 0, nrow = ncol(refload@intensity)))
              object@fit$ratio <- data.frame(matrix(ncol = 0, nrow = length(p)))
              for (i in 1:ncol(qload@intensity)) {
                message(paste(colnames(query@intensity)[i]), " (",i/ncol(query@intensity)*100, "%", ")", sep = "")
                r <- cor(qload@intensity[p, ], refload@intensity[p, ])[i, ] < 0.99
                if (any(!r)) message("query sample seems to also be in the reference set. not used for fit.")
                if (intercept) {
                  ref.fit <- lm(y ~ ., data = data.frame(y = qload@intensity[p,
                                                                             i], X = refload@intensity[p, r]))
                } else {
                  ref.fit <- lm(y ~ . - 1, data = data.frame(y = qload@intensity[p,
                                                                                 i], X = refload@intensity[p, r]))
                }
                object@fit$coef <- cbind(object@fit$coef,as.numeric(ref.fit$coefficients[-1]))

                ref.predict <- predict(ref.fit)
                ref.predict[ref.predict < 1] <- 1

                object@fit$ratio <- cbind(object@fit$ratio, as.numeric(log2(qload@intensity[p, i]/ref.predict[p])))
              }
              colnames(object@fit$coef) <- colnames(query@intensity)
              rownames(object@fit$coef) <- colnames(refload@intensity)
              colnames(object@fit$ratio) <- colnames(query@intensity)
              rownames(object@fit$ratio) <- p




              object@fit$noise <- as.numeric()
              for (i in 1:ncol(qload@intensity)) {
                object@fit$noise <- c(object@fit$noise, sqrt(mean((object@fit$ratio[-1,i] - object@fit$ratio[-nrow(object@fit$ratio),i])^2,
                                                                  na.rm = TRUE)))
              }
              names(object@fit$noise) <- colnames(query@intensity)
              return(object)          #reduce noise end
            } else {

            object <- new("CNV.analysis")
            object@date <- date()
            object@fit$args <- list(intercept = intercept)

            object@anno <- anno

            object@fit$coef <- data.frame(matrix(ncol = 0, nrow = ncol(ref@intensity)))
            object@fit$ratio <- data.frame(matrix(ncol = 0, nrow = length(p)))
            for (i in 1:ncol(query@intensity)) {
              message(paste(colnames(query@intensity)[i]), " (",i/ncol(query@intensity)*100, "%", ")", sep = "")
              r <- cor(query@intensity[p, ], ref@intensity[p, ])[i, ] < 0.99
              if (any(!r)) message("query sample seems to also be in the reference set. not used for fit.")
              if (intercept) {
                ref.fit <- lm(y ~ ., data = data.frame(y = query@intensity[p,
                                                                           i], X = ref@intensity[p, r]))
              } else {
                ref.fit <- lm(y ~ . - 1, data = data.frame(y = query@intensity[p,
                                                                               i], X = ref@intensity[p, r]))
              }
              object@fit$coef <- cbind(object@fit$coef,as.numeric(ref.fit$coefficients[-1]))

              ref.predict <- predict(ref.fit)
              ref.predict[ref.predict < 1] <- 1

              object@fit$ratio <- cbind(object@fit$ratio, as.numeric(log2(query@intensity[p, i]/ref.predict[p])))
            }


            colnames(object@fit$coef) <- colnames(query@intensity)
            rownames(object@fit$coef) <- colnames(ref@intensity)
            colnames(object@fit$ratio) <- colnames(query@intensity)
            rownames(object@fit$ratio) <- p

            object@fit$noise <- as.numeric()
            for (i in 1:ncol(query@intensity)) {
              object@fit$noise <- c(object@fit$noise, sqrt(mean((object@fit$ratio[-1,i] - object@fit$ratio[-nrow(object@fit$ratio),i])^2,
                                                                na.rm = TRUE)))
            }
            names(object@fit$noise) <- colnames(query@intensity)
            return(object)
          }})


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
  object@bin$shift <- as.numeric()
  for (i in 1:ncol(object@fit$ratio)) {
    object@bin$ratio[[i]] <- sapply(split(object@fit$ratio[o2[, "probe"],i], o2[,
                                                                                "bin"]), median, na.rm = TRUE)[names(object@anno@bins)]
    if (any(is.na(object@bin$ratio[[i]]))) {
      stop("not every bin contains at least one CpG, please reduce perc_cpgs in CNV.fit()")
   }

    object@bin$shift <- c(object@bin$shift, optim(0, function(s) median(abs(object@bin$ratio[[i]] -s),
                                                                        na.rm = TRUE), method = "Brent", lower = -100, upper = 100)$par)
  }
  names(object@bin$shift) <- colnames(object@fit$ratio)
  names(object@bin$ratio) <- colnames(object@fit$ratio)
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
    d2 <- data.frame(detail = values(object@anno@detail)$name[d1[,
                                                                 "queryHits"]], probe = names(object@anno@probes[d1[, "subjectHits"]]),
                     stringsAsFactors = FALSE)


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
  # if(length(object@fit) == 0) stop('fit unavailable, run CNV.fit')
  if (length(object@bin) == 0)
    stop("bin unavailable, run CNV.bin")
  # if(length(object@detail) == 0) stop('bin unavailable, run
  # CNV.detail')

  a1 <- formals()
  a2 <- as.list(match.call())[-1]
  object@seg$args <- as.list(sapply(setdiff(unique(names(c(a1, a2))),
                                            c("object", "verbose")), function(an) if (is.element(an, names(a2)))
                                              a2[[an]] else a1[[an]], simplify = FALSE))

  object@seg$summary <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@seg$p <- vector(mode = "list", length = ncol(object@fit$ratio))
  for (i in 1:ncol(object@fit$ratio)) {
    message(colnames(object@fit$ratio)[i])
    x1 <- DNAcopy::CNA(genomdat = object@bin$ratio[[i]][names(object@anno@bins)],
                       chrom = as.vector(seqnames(object@anno@bins)), maploc = values(object@anno@bins)$midpoint,
                       data.type = "logratio", sampleid = "sampleid")
    x2 <- DNAcopy::segment(x = x1, verbose = verbose, min.width = min.width,
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

