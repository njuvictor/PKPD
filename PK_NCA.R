install.packages("PKNCA")

library(PKNCA)
library(digitize)
#import conc-t profile from literature 

cal <- ReadAndCal('/Users/siweili/Documents/Debio_200mgPK.jpg')

#extract color point
conc_t <- DigitData(col = 'red')

#Calibrate(data,calpoints,x1,x2,y1,y2)
data <- Calibrate(conc_t, cal, 1, 24, 0.01, 1)

MW = 561
data["conc_ng/mL"] = data[,2] * MW #ng/mL

plot(data$x, data$y, pch=20, col='grey',
     xlab = 'Nutrients concentration',
     ylab = 'Divisions per hour')

#function to digitize log-scale plot from imagine 
digitize.graph <- function(name,x1,x2,y1,y2,sets=1,setlabels=1:sets,log='',xlab='x axis',ylab='y axis')
{
  dataset <- data.frame(x=NULL,y=NULL,lab=NULL)
  cat('Mark axes min and max values \n')
  axes.points <- ReadAndCal(name)
  if(log=='x'){x1 <- log10(x1);x2 <- log10(x2)}
  if(log=='y'){y1 <- log10(y1);y2 <- log10(y2)}
  if(log=='xy'){x1 <- log10(x1);x2 <- log10(x2);y1 <- log10(y1);y2 <- log10(y2)}
  for(i in 1:sets)
  {
    cat(paste('Mark point set "',setlabels[i],'"\n',sep=''))
    data.points <- DigitData(col = 'red')
    dat <- Calibrate(data.points, axes.points, x1, x2, y1, y2)      
    dat$lab <- rep(setlabels[i],nrow(dat))
    dataset <- rbind(dat, dataset)
  }
  if(log=='x'){dataset$x <- 10^(dataset$x)}
  if(log=='y'){dataset$y <- 10^(dataset$y)}
  if(log=='xy'){dataset$x <- 10^(dataset$x);dataset$y <- 10^(dataset$y)}
  print(dataset)
  plot(dataset$x,dataset$y,log=log,pch=as.numeric(as.factor(dataset$lab)),col=as.numeric(as.factor(dataset$lab)),xlab=xlab,ylab=ylab)
  #legend('bottomright',setlabels, pch=1:sets,col=1:sets, bty='n')
  return(dataset)    
}

temp <- digitize.graph('/Users/siweili/Documents/Debio_200mgPK.jpg',
               x1=1, x2=24, y1=0.1, y2=10, 
               log = 'y', xlab = "Time(h)", ylab = "umol/L")

data <- temp
MW = 561
data["conc_ng/mL"] = data[,2] * MW #ng/mL

write.csv(data, "/Users/siweili/Downloads/debio_200mg_conc.csv")



###Back up
#convert to c.conc for NCA analysis 
d.conc= cbind(data[,4],data[,1],data[,3])
d.conc = data.frame(d.conc)
names(d.conc) <- c("conc","time", "subject")

temp <- cbind(300000, 0, 1)
temp <- data.frame(temp)
names(temp) <- c("dose", "time", "subject")
d.dose <- temp
my.dose <- PKNCAdose(d.dose, dose~time|subject)

# Load concentration-time data into a data.frame called d.conc
# with columns named "conc", "time", and "subject".
my.conc <- PKNCAconc(d.conc, conc~time|subject)
# Load dose-time data into a data.frame called d.dose
# with columns named "dose", "time", and "subject".
my.dose <- PKNCAdose(d.dose, dose~time|subject)
# Combine the concentration-time and dose-time data into an object
# ready for calculations.
my.data <- PKNCAdata(my.conc, my.dose)
# Perform the calculations
my.results <- pk.nca(my.data)
# Look at summary results
summary(my.results)


