newton <- read.table("newton-convergence.dat")
fas <- read.table("fas-convergence.dat")
newton.fas <- read.table("newton-fas-convergence.dat")

xlim <- range(newton$V1)

ylim <- range(c(newton$V2, fas$V2, newton.fas$V2))

par(mar=c(7, 7, 1, 1), las=0)
pdf("convergence.pdf", width=10, height=6.18, pointsize=16)
plot(V2 ~ V1, newton, xlim=xlim, ylim=ylim, pch=20, type="b", lwd=3, log="y", ann=F,
     axes=F)

points(V2 ~ V1, fas, pch=21, type="b", lwd=3, col="blue", lty="dashed")

points(V2 ~ V1, newton.fas, pch=22, type="b", lwd=3, col="brown", lty="dotted")


axis(1)
mtext("Iteration", side=1, line=3)

mtext("Residual", side=2, line=3)
par(las=1)
axis(2)

legend("bottomright", c("Newton", "FAS", "Newton-FAS"),
       lwd=3, pch=c(20, 21, 22), col=c("black", "blue", "brown"),
       bty="n")

dev.off()
