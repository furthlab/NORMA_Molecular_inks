# Linear sequence from 0 to 1
x <- seq(0, 1, length.out = 7)

# Apply inverse power transform to cluster near the end
y <- round(1 + (56 - 1) * (1 - (1 - x)^2.5))
y

dat <- read.table('./data/csv/multiphoton/average_instensity_of_spot_over_timeResults.csv', sep=',', header=TRUE)
dat$Mean1<-dat$Mean1-min(dat$Mean1)

quartz(width = 134.6871/30, height = 108.783/30)
par(mar=c(4,5,0.25,1), xaxs='i')
tim <- seq(0, nrow(dat)-1)*30/60

plot(tim, dat$Mean1, las=1, type='l', lwd=2, xlab='Time (min)', ylab='Fluorescent \n intensity', axes=F, ylim=c(0,600))
axis(2, at=seq(0,600,by=200), las=1)
polygon(c(tim, rev(tim)), c(dat$Mean1+dat$StdDev1, rev(dat$Mean1-dat$StdDev1)), col='lightblue', lty=1)
lines(tim, dat$Mean1, lwd=2)
points(tim, dat$Mean1, pch=21, bg='white')
box()
axis(1)
par(xpd=T)
points( (seq(0, nrow(dat)-1)*30/60)[y], dat$Mean1[y], pch=21, bg='orange', cex=2)
par(xpd=F)

quartz.save(file = "./figures/pdf/kinetics_multiphoto.pdf", type="pdf")


