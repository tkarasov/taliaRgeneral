#' @title Make Kinship Matrix
#' @description Reads in vcf file and creates kinship covariance matrix
#' @param otu.file File containing the vcf
#'



#tutorial on relatedness stuff here: https://bioconductor.statistik.tu-dortmund.de/packages/3.2/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#format-conversion-from-vcf-files
#vcf = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/poolsGVCF.filtered_snps_final.PASS.bi.vcf"
#genot_dir = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/"
#library(SNPRelate)


make_kinship<-function(vcf){
  vcf.fn <- vcf
  
  # Reformat
  genofile <- snpgdsVCF2GDS(vcf.fn, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/test.gds", method="biallelic.only")
  #genofile <- snpgdsOpen("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/test.gds")
  #snpgdsClose("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/test.gds")
  ibs <- snpgdsIBS(snpgdsOpen(genofile), num.thread = 4)
  
  set.seed(1000)
  
  # Try different LD thresholds for sensitivity analysis
  # snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
  return(ibs)
}

