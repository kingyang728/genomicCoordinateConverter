
#install.packages("R.utils")
#remotes::install_github("kingyang728/genomicCoordinateConverter",build = FALSE)
remove.packages("genomicCoordinateConverter")


devtools::install_github("kingyang728/genomicCoordinateConverter")
library(genomicCoordinateConverter)




DF1<-  Coordinate_Covnerter1("chrX",48823056,48823056,"G","C")
DF2 <- Coordinate_Covnerter2("chrX",48823056,48823056,"G","C")
DF3 <- Coordinate_Covnerter3("chrX",48823056,48823056)
DF4 <- Coordinate_Covnerter4("chrX",48823056,48823056)
DF5_1<-Coordinate_Covnerter5_1("HDAC6-201",2657,2657)
DF5_2<-Coordinate_Covnerter5_2("ENST00000334136",2657)
DF6 <- Coordinate_Covnerter6("chrX",48823056,48823056,"hg16","hg18")
