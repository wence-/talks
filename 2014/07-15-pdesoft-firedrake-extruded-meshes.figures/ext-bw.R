p0p0 <- read.table("ext-bw-data/P0-P0-mpi-vs-vbw.dat", header=T)
p0p1 <- read.table("ext-bw-data/P0-P1-mpi-vs-vbw.dat", header=T)
p0p1dg <- read.table("ext-bw-data/P0-P1dg-mpi-vs-vbw.dat", header=T)
p1dgp0 <- read.table("ext-bw-data/P1dg-P0-mpi-vs-vbw.dat", header=T)
p1dgp1 <- read.table("ext-bw-data/P1dg-P1-mpi-vs-vbw.dat", header=T)#
p1p0 <- read.table("ext-bw-data/P1-P0-mpi-vs-vbw.dat", header=T)
p1p1 <- read.table("ext-bw-data/P1-P1-mpi-vs-vbw.dat", header=T)
p1p1dg <- read.table("ext-bw-data/P1-P1dg-mpi-vs-vbw.dat", header=T)

source("tufteaxis.R")

xlim <- range(p0p0$nproc)
ylim <- range(p0p0$X75, p0p1$X75, p0p1dg$X75,
              p1dgp0$X75, p1dgp1$X75,
              p1p0$X75, p1p1$X75, na.rm=T)

pdf.options()
pdf(file="ext-vbw.pdf", width=12, height=8,
    pointsize=20)
par(las=1, mar=c(4, 5, 1, 0), mfrow=c(1, 2))

plot(X75 ~ nproc, p0p0, axes=F, ann=F,
     type="b", lwd=4, xlim=xlim, ylim=ylim,
     col="blue", pch=21, bg="blue")

points(X75 ~ nproc, p0p1, type="b", lwd=4, col="skyblue", lty="dotted",
       pch=22, bg="skyblue")

points(X75 ~ nproc, p0p1dg, type="b", lwd=4, col="steelblue", lty="dashed",
       pch=23, bg="steelblue")

points(X75 ~ nproc, p1p0, type="b", lwd=4, col="red",
       bg="red", pch=21)

points(X75 ~ nproc, p1p1, type="b", lwd=4, col="tomato", lty="dotted",
       bg="tomato", pch=22)

points(X75 ~ nproc, p1p1dg, type="b", lwd=4, col="red3", lty="dashed",
       bg="red3", pch=23)

points(X75 ~ nproc, p1dgp0, type="b", lwd=4, col="seagreen",
       bg="seagreen", pch=21)

points(X75 ~ nproc, p1dgp1, type="b", lwd=4, col="darkgreen", lty="dotted",
       bg="darkgreen", pch=22)


tufteaxis(1, at=c(1, 4, 12, 24))
tufteaxis(2, amin=min(ylim), amax=max(ylim), digits=0)
mtext("MPI processes", side=1, line=3)
par(las=0)
mtext("VBW [GB/s]", side=2, line=4)
title("75 cell layers")

## dev.off()
## pdf.options()
## pdf(file="ext-vbw-layers.pdf", width=6, height=8,
##     pointsize=16)
## par(las=1, mar=c(5, 5, 1, 1))

layers <- c(1, 4, 8, 16, 25, 50)
xlim <- range(layers)

par(mar=c(4, 1, 1, 1))
plot(layers, p0p0[3, 2:7], axes=F, ann=F,
     type="b", lwd=4, xlim=xlim, ylim=ylim,
     col="blue", pch=21, bg="blue")

points(layers, p0p1[3, 2:7], type="b", lwd=4, col="skyblue", lty="dotted",
       pch=22, bg="skyblue")

points(layers, p0p1dg[3, 2:7], type="b", lwd=4, col="steelblue", lty="dashed",
       pch=23, bg="steelblue")

points(layers, p1p0[3, 2:7], type="b", lwd=4, col="red",
       bg="red", pch=21)

points(layers, p1p1[3, 2:7], type="b", lwd=4, col="tomato", lty="dotted",
       bg="tomato", pch=22)

points(layers, p1p1dg[3, 2:7], type="b", lwd=4, col="red3", lty="dashed",
       bg="red3", pch=23)

points(layers, p1dgp0[3, 2:7], type="b", lwd=4, col="seagreen",
       bg="seagreen", pch=21)

points(layers, p1dgp1[3, 2:7], type="b", lwd=4, col="darkgreen", lty="dotted",
       bg="darkgreen", pch=22)


tufteaxis(1, at=layers)
mtext("Number of cell layers", side=1, line=3)
title("12 MPI processes")

par(mar=c(5,10,0,0), mfrow=c(1, 1), new=T)
plot.new()
legend("bottom",
       c("P0xP0", "P0xP1", "P0xP1dg",
         "P1xP0", "P1xP1", "P1xP1dg",
         "P1dgxP0", "P1dgxP1"),
       ncol=3,
       bty="n",
       col=c("blue", "skyblue", "steelblue",
           "red", "tomato", "red3",
           "seagreen", "darkgreen"),
       pt.bg=c("blue", "skyblue", "steelblue",
           "red", "tomato", "red3",
           "seagreen", "darkgreen"),
       pch=c(21, 22, 23),
       lty=c("solid", "dotted", "dashed"),
       lwd=4, seg.len=4, adj=0.5, x.intersp=3)
dev.off()
