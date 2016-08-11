
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **TERES_Standardization** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of QuantLet : TERES_Standardization

Published in : Tail Event Risk Expected Shortfall

Description : Uses a GARCH(1,1) model and a constant mean to standardze the data input.

Keywords : 'dynamic, garch, heavy-tailed, heteroskedasticity, nonstationary, preprocessing,
returns, standardization, variance'

See also : 'MSEconfexpectile0.95, SFSconfexpectile0.95, SFSconfexpectile0.95,
TERES_ExpectileQuantileDiffMulti'

Author : Philipp Gschöpf, Andrija Mihoci

Submitted : Tue, August 2 2016 by Roman Lykhnenko

```


### R Code:
```r
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

```
