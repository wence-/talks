setEPS()
L3bw <- read.table("L3-bandwidth.dat")

names(L3bw) <- c("layers", "T1", "T2", "T3", "T4", "STREAM")

## Fix layer numbering
L3bw$layers <- L3bw$layers - 1
source("tufteaxis.R")
xlim <- range(L3bw$layers)
ylim <- range(L3bw[2:5])

postscript(file="L3-bandwidth.eps", width=10, height=6.18, pointsize=20)
par(las=1, mar=c(4.5, 5.5, 0.5, 0.5))
plot(T1 ~ layers, L3bw, xlim=xlim, ylim=ylim,
       type="b", lwd=3, col="black", pch=21, bg="black",
     axes=F, ann=F)
points(T2 ~ layers, L3bw,
       type="b", lwd=3, col="blue", pch=22, bg="blue",
       lty="dashed")
points(T3 ~ layers, L3bw,
       type="b", lwd=3, col="purple", pch=23, bg="purple",
       lty="dotted")
points(T4 ~ layers, L3bw,
       type="b", lwd=3, col="brown", pch=24, bg="brown",
       lty="dotdash")
tufteaxis(1)
par(las=1)
tufteaxis(2, amin=min(ylim), amax=max(ylim), at=c(5000, 10000, 15000), mingap=0.2, digits=0)
mtext("Number of cell layers", side=1, line=3)
par(las=0)
mtext("L3 bandwidth/(MB/s)", side=2, line=4)
legend("topright", legend=c("1 Thread", "2 Threads", "3 Threads", "4 Threads"),
       bty="n", pch=c(21, 22, 23, 24), col=c("black", "blue", "purple", "brown"),
       pt.bg=c("black", "blue", "purple", "brown"), lwd=3)

dev.off()
Vbw <- read.table("valuable-bandwidth-by-layer.dat")
names(Vbw) <- c("layers", "T1", "T2", "T3", "T4")
Vbw$layers <- Vbw$layers - 1
STREAM <- 11341
xlim <- range(Vbw$layers)
ylim <- range(c(Vbw[2:5], STREAM))
postscript(file="valuable-bandwidth-by-layer.eps", width=10, height=6.18, pointsize=16)
par(las=1, mar=c(4.5, 5.5, 0.5, 0.5))
plot(T1 ~ layers, Vbw, xlim=xlim, ylim=ylim,
       type="b", lwd=3, col="black", pch=21, bg="black",
     axes=F, ann=F)
points(T2 ~ layers, Vbw,
       type="b", lwd=3, col="blue", pch=22, bg="blue",
       lty="dashed")
points(T3 ~ layers, Vbw,
       type="b", lwd=3, col="purple", pch=23, bg="purple",
       lty="dotted")
points(T4 ~ layers, Vbw,
       type="b", lwd=3, col="brown", pch=24, bg="brown",
       lty="dotdash")
tufteaxis(1)
par(las=1)
tufteaxis(2, amin=min(ylim), amax=STREAM, mingap=0.2, digits=0)
mtext("Number of cell layers", side=1, line=3)
par(las=0)
mtext("Valuable bandwidth/(MB/s)", side=2, line=4)
legend("topright", legend=c("1 Thread", "2 Threads", "3 Threads", "4 Threads"),
       bty="n", pch=c(21, 22, 23, 24), col=c("black", "blue", "purple", "brown"),
       pt.bg=c("black", "blue", "purple", "brown"), lwd=3)
dev.off()

Vbw8 <- read.table("valuable-bandwidth-by-thread.dat")
names(Vbw8) <- c("threads", "L1nsf", "L1", "L3", "L10", "L20", "L50", "L100", "L150", "L200")
xlim <- range(Vbw8$threads)
Vbw8 <- rbind(Vbw8[1:4,], Vbw8[8,])
ylim <- range(c(Vbw8[,2:9], STREAM), na.rm=T)
postscript(file="valuable-bandwidth-by-thread.eps", width=10, height=6.18, pointsize=16)
par(las=1, mar=c(4.5, 5.5, 0.5, 0.5))
plot(L1 ~ threads, Vbw8, xlim=xlim, ylim=ylim,
     type="b", lwd=3, col="black", pch=18, bg="black",
     axes=F, ann=F)
points(L3 ~ threads, Vbw8, type="b", lwd=3,
       col="blue", pch=19, bg="blue")
points(L10 ~ threads, Vbw8, type="b", lwd=3,
       col="purple", pch=20, bg="purple")
points(L20 ~ threads, Vbw8, type="b", lwd=3,
       col="gray", pch=21, bg="gray")
points(L50 ~ threads, Vbw8, type="b", lwd=3,
       col="brown", pch=22, bg="brown")
points(L100 ~ threads, Vbw8, type="b", lwd=3,
       col="cyan", pch=23, bg="cyan")
points(L150 ~ threads, Vbw8, type="b", lwd=3,
       col="magenta", pch=24, bg="magenta")
points(L200 ~ threads, Vbw8, type="b", lwd=3,
       col="pink", pch=25, bg="pink")
points(L1nsf ~ threads, Vbw8, type="b", lwd=3,
       col="seagreen", pch=17, bg="seagreen")

tufteaxis(1)
par(las=1)
tufteaxis(2, amin=min(ylim), amax=STREAM, mingap=0.2, digits=0)
mtext("Number of threads", side=1, line=3)
par(las=0)
mtext("Valuable bandwidth/(MB/s)", side=2, line=4)
legend("topleft", legend=c("1 Layer (bad base numbering)", "1 Layer", "3 Layers",
                      "10 Layers", "20 Layers",
                      "50 Layers", "100 Layers", "150 Layers", "200 Layers"),
       bty="n", pch=c(17:25),
       col=c("seagreen", "black", "blue", "purple", "gray", "brown", "cyan", "magenta", "pink"),
       pt.bg=c("seagreen", "black", "blue", "purple", "gray", "brown", "cyan", "magenta", "pink"), lwd=3)
dev.off()

xlim <- range(Vbw8$threads[1:4])
ylim <- range(Vbw8[,2:3], na.rm=T)
postscript(file="bad-numbering.eps", width=10, height=6.18, pointsize=16)
par(las=1, mar=c(4.5, 5.5, 0.5, 0.5))
plot(L1nsf ~ threads, Vbw8[1:4,], type="b", lwd=3,
     col="black", pch=21, bg="black", xlim=xlim, ylim=ylim, ann=F, axes=F)
points(L1 ~ threads, Vbw8[1:4,], type="b", lwd=3,
       col="blue", pch=22, bg="blue", lty="dashed")
tufteaxis(1, at=c(1,2,3,4))
par(las=1)
mtext("Number of threads", side=1, line=3)
tufteaxis(2, amin=min(ylim), amax=max(ylim), mingap=0.2, digits=0)
par(las=0)
mtext("Valuable bandwidth/(MB/s)", side=2, line=4)
legend("topleft", legend=c("Original numbering (1 layer)", "Good numbering (1 layer)"),
       bty="n", pch=21:22, col=c("black", "blue"),
       pt.bg=c("black", "blue"))
dev.off()
