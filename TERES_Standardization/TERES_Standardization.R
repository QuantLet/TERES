#Load required packages
library("fGarch")

#Uncomment to change your working directory
#setwd()

Stock     = read.csv("DataIndices.csv")  
y         = Stock[,4]
y         = diff(log(y))
y         = na.omit(y)

#pre-white data with a GARCH model
GARCHvola= garchFit(~garch(1, 1), data =y)
ywhite=y / volatility(GARCHvola)

yclean=ywhite-mean(ywhite)

# uncomment to save the results
# write.table(yclean,"SP500.dat",row.names=F,col.names=F)
