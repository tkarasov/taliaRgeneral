write_phyloseq_tlk <- function(x, type="all", path=getwd()) {
  
  # TODO make read_phyloseq as well
  if (type == "OTU" || type == "all") {
    f <- paste(path, "otu_table.csv", sep="/")
    message("Writing OTU in the file ", f)
    # y <- as.data.frame(x@otu_table);
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    # Let us use abundances function here as it is guaranteed to be taxa x
    # samples always
    y <- abundances(x)
    write.csv(y, file=f)
  } else if (type == "TAXONOMY" || type == "all") {
    # Renamed from TAXA to TAXONOMY as the latter is used elsewhere
    f <- paste(path, "taxonomy_table.csv", sep="/")
    message("Writing TAXONOMY in the file ", f)
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    y <- as.data.frame(tax_table(x))
    write.csv(y, file=f)
  } else if (type == "METADATA" || type == "all") {
    f <- paste(path, "metadata_table.csv", sep="/")
    message("Writing METADATA in the file ", f)
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    y <- meta(x)
    write.csv(y, file=f)
  }
  
  return(path)
  
}
