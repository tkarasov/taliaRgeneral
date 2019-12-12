# Function that will generate a pdf 
# @param title Title for the figure
# @param location Location for file to be saved
# @param width Width of figure
# @param height Height of figure


pdf_pub<-function(title, location = NULL, width = 88mm, height = NULL){
  pdf(paste(location, title, sep=""), width = 88mm, height = height, family = "Arial", useDingbats=FALSE))
  }
