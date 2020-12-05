setwd("/Users/siweili/Documents/")

library(digitize)

# 
# MonodGrowth <- function(params, M) {
#   with(params, rK*(M/(M1+M)))
# }
# 
# MonodError <- function(params, M, y) {
#   with(params,
#        sum((MonodGrowth(params, M)-y)^2))
# }


cal <- ReadAndCal('test.jpg')
#extract color point
growth <- DigitData(col = 'red')
#Calibrate(data,calpoints,x1,x2,y1,y2)
data <- Calibrate(growth, cal, 1, 8, 0.5, 1)

plot(data$x, data$y, pch=20, col='grey',
     xlab = 'Nutrients concentration',
     ylab = 'Divisions per hour')

# points(xcal, MonodGrowth(out$set, xcal),
#        type = 'l', lty = 1, lwd = 2)
# points(xcal, MonodGrowth(paper, xcal),
#        type = 'l', lty = 2)
# legend('bottomright',
#        legend = c('data', 'best fit', 'paper value'),
#        pch = c(20, NA, NA),
#        lty = c(NA, 1, 2),
#        lwd = c(NA, 2, 1),
#        col = c('grey', 'black', 'black'),
#        bty = 'n')
