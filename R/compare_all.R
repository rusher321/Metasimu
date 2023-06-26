#' Compare_meta
#'
#' @param dat, matrix or data.frame, row is sample id
#' @param config , data.frame. group information
#' @param prep , default = 0.2,cutoff of filter
#' @param method , character, "all" or specific method
#' @param C , confounders only available for locom or ancombc
#' @param ...
#'
#' @return list
#' @export
#'
#' @examples
Compare_meta <- function(dat, config, prep = 0.2, method = "all", C = NULL, ...){

  id <- intersect(rownames(dat), rownames(config))
  data_match <- dat[id, ]
  Y <- config[id, 1]
  data_match <- data_match[, apply(data_match, 2, function(x){(sum(x!=0)/length(x)) >= prep})]

  outlist <- list()
  ############ locom compare
  if(method == "all" | method == "locom"){
  res.locom <- LOCOM::locom(otu.table = data_match, Y = Y, C = C,
                     ref.otu = NULL, fdr.nominal = 0.2,
                     n.perm.max = 50000, prev.cut = prep, n.rej.stop = 100, ...)

  out.locom <- data.frame(effect = t(res.locom$effect.size)[,1],
                          pvalue = t(res.locom$p.otu)[,1],
                          qvalue = t(res.locom$q.otu)[,1])

  rownames(out.locom) <- colnames(res.locom$effect.size)
  outlist[["locom"]] <- out.locom
  }

  ########### ancom-bc method
  sam.ID <- rownames(config)
  if (is.null(C)){
    meta.data <- data.frame(group = Y, row.names = id, stringsAsFactors = FALSE)
    formula <- "group"
  } else {
    meta.data <- data.frame(group = Y, C, row.names = id, stringsAsFactors = FALSE)
    colnames(meta.data)[2:ncol(meta.data)] <- colnames(C)
    formula <- paste("group + ",  paste(colnames(C), collapse = " + "), sep = "")
  }

  otu.table.ancombc <- t(data_match)
  otus.keep <- rownames(otu.table.ancombc)
  tax.table <- matrix(sample(letters, 7*length(otus.keep), replace = TRUE),
                      nrow = length(otus.keep), ncol = 7)
  rownames(tax.table) <- rownames(otu.table.ancombc)
  colnames(tax.table) <- c("Kingdom", "Phylum", "Class", "Order", "Family",
                           "Genus", "Species")

  otu.phylo <- phyloseq::otu_table(otu.table.ancombc, taxa_are_rows = TRUE)
  meta.phylo <- phyloseq::sample_data(meta.data)
  tax.phylo <- phyloseq::tax_table(tax.table)
  phylo.obj <- phyloseq::phyloseq(otu.phylo, meta.phylo, tax.phylo)

  if(method == "all" | method == "ancombc"){
      res.ancombc <- ANCOMBC::ancombc(phyloseq = phylo.obj, formula = formula,
                             p_adj_method = "BH",
                             zero_cut = 1-prep,
                             lib_cut = 0,
                             group = "group", struc_zero = TRUE, neg_lb = FALSE,
                             tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05,
                             global = TRUE)

      otu.ancombc <- data.frame(effect = res.ancombc$res$W[,1],
                                pvalue = res.ancombc$res$p_val[,1],
                                qvalue = res.ancombc$res$q_val[,1])

      rownames(otu.ancombc) <- rownames(res.ancombc$res$W)
      outlist[["ancombc"]] <- otu.ancombc
  }
  ##################### wilcox-raw
  if(method == "all" | method == "wilcox"){
    wilcox_res <- wilcox_all(datamatrix = data_match, configdata = meta.data)
    outlist[["wilcox_raw"]] <- wilcox_res
  }
  #################### wilcox-clr-psudo
  if(method == "all" | method == "clr"){
      data_match_p <- data_match
      data_match_p[data_match_p == 0] <- min(data_match[data_match!=0])
      data_match_clr <- robCompositions::cenLR(data_match_p)$x.clr
      data_match_clr_f <- data_match_clr[, colnames(data_match)]
      wilcox_res_clr <- wilcox_all(datamatrix = data_match_clr_f, configdata = meta.data)
      outlist[["wilcox_clr"]] <- wilcox_res_clr
  }

  ################### wilcox-clr-zcomposition
  if(method == "all" |method == "zclr"){
    data_match_z <- zCompositions::cmultRepl(data_match)
    data_match__zclr <- cenLR(data_match_z)$x.clr
    data_match_zclr_f <- data_match__zclr[, colnames(data_match)]
    wilcox_res_zclr <- wilcox_all(datamatrix = data_match_zclr_f, configdata = meta.data)
    outlist[["wilcox_res_zclr"]] <- wilcox_res_zclr
  }

  return(outlist)

}



wilcox_all <- function (datamatrix, configdata, ratio = "zero")
{

  out <- matrix(NA, nrow = ncol(datamatrix), ncol = 9)
  config <- as.factor(as.character(configdata[, 1]))
  nlevels = levels(droplevels(config))

  for (i in 1:ncol(datamatrix)) {
    tmp <- as.numeric(datamatrix[, i])
    g1 <- tmp[which(config == nlevels[1])]
    g2 <- tmp[which(config == nlevels[2])]
    wilcox_sign <- wilcox.test(g1, g2)$p.value
    data_tmp <- data.frame(tmp = tmp, config = config)
    effect <- as.numeric(rstatix::wilcox_effsize(tmp~config, data = data_tmp)[1,4])

    if (ratio == "zero") {
      or <- tapply(tmp, config, function(x) {
        sum(x != 0, na.rm = T)/length(x)
      })
    }
    else {
      or <- tapply(tmp, config, function(x) {
        sum(!is.na(x))/length(x)
      })
    }
    mean_abun <- tapply(tmp, config, mean, na.rm = T)
    median_abun <- tapply(tmp, config, median, na.rm = T)
    z_score <- statistic(coin::wilcox_test(tmp~config, data_tmp))
    out[i, 1:9] <- c(wilcox_sign, or, mean_abun, median_abun,
                     z_score, effect)
  }
  rownames(out) <- colnames(datamatrix)
  colnames(out) <- c("sign_p.value", paste0(rep(nlevels,
                                                3), rep(c("_ratio", "_mean", "_median"),
                                                        each = 2)), "z_score", "effect_size")
  out <- as.data.frame(out)
  out$p.adjust <- p.adjust(out$sign_p.value, method = "BH")
  out$enrich <- ifelse(out$p.adjust < 0.05, ifelse(out$z_score >
                                                     0, nlevels[1], nlevels[2]), "none")
  return(out)

}


