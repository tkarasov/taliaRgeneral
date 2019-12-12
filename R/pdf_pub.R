# Function that will generate a pdf 
# @param title Title for the figure
# @param location Location for file to be saved
# @param width Width of figure
# @param height Height of figure

# 88mm = 3.46in
pdf_pub<-function(title, location = NULL, width = 3.46, height = NULL){
  pdf(paste(location, title, sep=""), width = width, height = height, family = "Arial", useDingbats=FALSE)
  }
