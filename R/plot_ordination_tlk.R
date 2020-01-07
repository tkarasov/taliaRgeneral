#' @title Plot ordination of phyloseq objects
#' @description Modified version of functions from phyloseq package
#' @param otu.file File containing the vcf
#'



################################################################################
# Define S3 generic extract_eigenvalue function; formerly S4 generic get_eigenvalue()
# Function is used by `plot_scree` to get the eigenvalue vector from different
# types of ordination objects. 
# Used S3 generic in this case because many ordination objects, the input, are
# not formally-defined S4 classes, but vaguely-/un-defined S3. This throws
# warnings during package build if extract_eigenvalue were S4 generic method,
# because the ordination classes don't appear to have any definition in phyloseq
# or dependencies.
#' @keywords internal
extract_eigenvalue = function(ordination) UseMethod("extract_eigenvalue", ordination)
# Default is to return NULL (e.g. for NMDS, or non-supported ordinations/classes).
extract_eigenvalue.default = function(ordination) NULL
# for pcoa objects
extract_eigenvalue.pcoa = function(ordination) ordination$values$Relative_eig
# for CCA objects
extract_eigenvalue.cca = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for RDA objects
extract_eigenvalue.rda = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for dpcoa objects
extract_eigenvalue.dpcoa = function(ordination) ordination$eig
# for decorana (dca) objects
extract_eigenvalue.decorana = function(ordination) ordination$evals

my.rarefy, ordinate(my.rarefy), "MDS", color = "Species"

plot_ordination_tlk <-function (physeq, ordination, type = "samples", axes = 1:2, color = NULL, 
          shape = NULL, label = NULL, title = NULL, justDF = FALSE) 
{
  if (length(type) > 1) {
    warning("`type` can only be a single option,\n            but more than one provided. Using only the first.")
    type <- type[[1]]
  }
  if (length(color) > 1) {
    warning("The `color` variable argument should have length equal to 1.", 
            "Taking first value.")
    color = color[[1]][1]
  }
  if (length(shape) > 1) {
    warning("The `shape` variable argument should have length equal to 1.", 
            "Taking first value.")
    shape = shape[[1]][1]
  }
  if (length(label) > 1) {
    warning("The `label` variable argument should have length equal to 1.", 
            "Taking first value.")
    label = label[[1]][1]
  }
  official_types = c("sites", "species", "biplot", "split", 
                     "scree")
  if (!inherits(physeq, "phyloseq")) {
    if (inherits(physeq, "character")) {
      if (physeq == "list") {
        return(official_types)
      }
    }
    warning("Full functionality requires `physeq` be phyloseq-class ", 
            "with multiple components.")
  }
  type = gsub("^.*site[s]*.*$", "sites", type, ignore.case = TRUE)
  type = gsub("^.*sample[s]*.*$", "sites", type, ignore.case = TRUE)
  type = gsub("^.*species.*$", "species", type, ignore.case = TRUE)
  type = gsub("^.*taxa.*$", "species", type, ignore.case = TRUE)
  type = gsub("^.*OTU[s]*.*$", "species", type, ignore.case = TRUE)
  type = gsub("^.*biplot[s]*.*$", "biplot", type, ignore.case = TRUE)
  type = gsub("^.*split[s]*.*$", "split", type, ignore.case = TRUE)
  type = gsub("^.*scree[s]*.*$", "scree", type, ignore.case = TRUE)
  if (!type %in% official_types) {
    warning("type argument not supported. `type` set to 'samples'.\n", 
            "See `plot_ordination('list')`")
    type <- "sites"
  }
  if (type %in% c("scree")) {
    return(plot_scree(ordination, title = title))
  }
  is_empty = function(x) {
    length(x) < 2 | suppressWarnings(all(is.na(x)))
  }
  specDF = siteDF = NULL
  trash1 = try({
    siteDF <- scores(ordination, choices = axes, display = "sites", 
                     physeq = physeq)
  }, silent = TRUE)
  trash2 = try({
    specDF <- scores(ordination, choices = axes, display = "species", 
                     physeq = physeq)
  }, silent = TRUE)
  siteSampIntx = length(intersect(rownames(siteDF), sample_names(physeq)))
  siteTaxaIntx = length(intersect(rownames(siteDF), taxa_names(physeq)))
  specSampIntx = length(intersect(rownames(specDF), sample_names(physeq)))
  specTaxaIntx = length(intersect(rownames(specDF), taxa_names(physeq)))
  if (siteSampIntx < specSampIntx & specTaxaIntx < siteTaxaIntx) {
    co = specDF
    specDF <- siteDF
    siteDF <- co
    rm(co)
  }else {
    if (siteSampIntx < specSampIntx) {
      siteDF <- specDF
      specDF <- NULL
    }
    if (specTaxaIntx < siteTaxaIntx) {
      specDF <- siteDF
      siteDF <- NULL
    }
  }
  if (is_empty(siteDF) & is_empty(specDF)) {
    warning("Could not obtain coordinates from the provided `ordination`. \n", 
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
    return(NULL)
  }
  if (is_empty(specDF) & type != "sites") {
    message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
    specDF <- data.frame(wascores(siteDF, w = veganifyOTU(physeq)), 
                         stringsAsFactors = FALSE)
  }
  if (is_empty(siteDF) & type != "species") {
    message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
    siteDF <- data.frame(wascores(specDF, w = t(veganifyOTU(physeq))), 
                         stringsAsFactors = FALSE)
  }
  specTaxaIntx <- siteSampIntx <- NULL
  siteSampIntx <- length(intersect(rownames(siteDF), sample_names(physeq)))
  specTaxaIntx <- length(intersect(rownames(specDF), taxa_names(physeq)))
  if (siteSampIntx < 1L & !is_empty(siteDF)) {
    warning("`Ordination site/sample coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
    siteDF <- NULL
  }
  if (specTaxaIntx < 1L & !is_empty(specDF)) {
    warning("`Ordination species/OTU/taxa coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
    specDF <- NULL
  }
  if (is_empty(siteDF) & is_empty(specDF)) {
    warning("Could not obtain coordinates from the provided `ordination`. \n", 
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
    return(NULL)
  }
  if (type %in% c("biplot", "split") & (is_empty(siteDF) | 
                                        is_empty(specDF))) {
    if (is_empty(siteDF)) {
      warning("Could not access/evaluate site/sample coordinates. Switching type to 'species'")
      type <- "species"
    }
    if (is_empty(specDF)) {
      warning("Could not access/evaluate species/taxa/OTU coordinates. Switching type to 'sites'")
      type <- "sites"
    }
  }
  if (type != "species") {
    sdf = NULL
    sdf = data.frame(access(physeq, slot = "sam_data"), stringsAsFactors = FALSE)
    if (!is_empty(sdf) & !is_empty(siteDF)) {
      siteDF <- cbind(siteDF, sdf[rownames(siteDF), ])
    }
  }
  if (type != "sites") {
    tdf = NULL
    tdf = data.frame(access(physeq, slot = "tax_table"), 
                     stringsAsFactors = FALSE)
    if (!is_empty(tdf) & !is_empty(specDF)) {
      specDF = cbind(specDF, tdf[rownames(specDF), ])
    }
  }
  if (!inherits(siteDF, "data.frame")) {
    siteDF <- as.data.frame(siteDF, stringsAsFactors = FALSE)
  }
  if (!inherits(specDF, "data.frame")) {
    specDF <- as.data.frame(specDF, stringsAsFactors = FALSE)
  }
  DF = NULL
  DF <- switch(EXPR = type, sites = siteDF, species = specDF, 
               {
                 specDF$id.type <- "Taxa"
                 siteDF$id.type <- "Samples"
                 colnames(specDF)[1:2] <- colnames(siteDF)[1:2]
                 DF = merge(specDF, siteDF, all = TRUE)
                 if (!is.null(shape)) {
                   DF <- rp.joint.fill(DF, shape, "Samples")
                 }
                 if (!is.null(shape)) {
                   DF <- rp.joint.fill(DF, shape, "Taxa")
                 }
                 if (!is.null(color)) {
                   DF <- rp.joint.fill(DF, color, "Samples")
                 }
                 if (!is.null(color)) {
                   DF <- rp.joint.fill(DF, color, "Taxa")
                 }
                 DF
               })
  if (justDF) {
    return(DF)
  }
  if (!is.null(color)) {
    if (!color %in% names(DF)) {
      warning("Color variable was not found in the available data you provided.", 
              "No color mapped.")
      color <- NULL
    }
  }
  if (!is.null(shape)) {
    if (!shape %in% names(DF)) {
      warning("Shape variable was not found in the available data you provided.", 
              "No shape mapped.")
      shape <- NULL
    }
  }
  if (!is.null(label)) {
    if (!label %in% names(DF)) {
      warning("Label variable was not found in the available data you provided.", 
              "No label mapped.")
      label <- NULL
    }
  }
  x = colnames(DF)[1]
  y = colnames(DF)[2]
  if (ncol(DF) <= 2) {
    message("No available covariate data to map on the points for this plot `type`")
    ord_map = aes_string(x = x, y = y)
  }else if (type %in% c("sites", "species", "split")) {
    ord_map = aes_string(x = x, y = y, color = color, shape = shape, 
                         na.rm = TRUE)
  }else if (type == "biplot") {
    if (is.null(color)) {
      ord_map = aes_string(x = x, y = y, size = "id.type", 
                           color = "id.type", shape = shape, na.rm = TRUE)
    }
    else {
      ord_map = aes_string(x = x, y = y, size = "id.type", 
                           color = color, shape = shape, na.rm = TRUE)
    }
  }
  p <- ggplot(DF, ord_map) + geom_point(na.rm = TRUE, fill = "id.type")
  if (type == "split") {
    p <- p + facet_wrap(~id.type, nrow = 1)
  }
  if (type == "biplot") {
    if (is.null(color)) {
      p <- update_labels(p, list(colour = "Ordination Type"))
    }
    p <- p + scale_size_manual("type", values = c(Samples = 5, 
                                                  Taxa = 2))
  }
  if (!is.null(label)) {
    label_map <- aes_string(x = x, y = y, label = label)
    p = p + geom_text(label_map, data = rm.na.phyloseq(DF, 
                                                       label), size = 2, vjust = 1.5, na.rm = TRUE)
  }
  if (!is.null(title)) {
    p = p + ggtitle(title)
  }
  if (length(extract_eigenvalue(ordination)[axes]) > 0) {
    eigvec = extract_eigenvalue(ordination)
    fracvar = eigvec[axes]/sum(eigvec)
    percvar = round(100 * fracvar, 1)
    strivar = as(c(p$label$x, p$label$y), "character")
    strivar = paste0(strivar, "   [", percvar, "%]")
    p = p + xlab(strivar[1]) + ylab(strivar[2])
  }
  p_fin <- p + theme_bw() + scale_fill_viridis_d()
  return(p_fin)
}
