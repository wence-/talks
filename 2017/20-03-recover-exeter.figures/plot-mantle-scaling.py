import numpy
from matplotlib import pyplot
import seaborn

data = numpy.asarray([[16384, 1],
                      [32768, 1.03],
                      [65536, 1.00],
                      [131072, 0.88],
                      [262144, 0.76],
                      [524288, 0.61],
                      [1048576, 0.43],
                      [1572864, 0.32]])

odata = data.copy()
data[:, 1] = data[:, 1] * (data[:, 0] / data[0, 0])

seaborn.set(style="ticks")

fig = pyplot.figure(figsize=(9, 5), frameon=False)

ax = fig.add_subplot(111)

ax.set_xlabel("Number of MPI ranks\n(dofs/rank)")

ax.loglog()
ax.minorticks_off()

ax.set_ylabel("Speedup")


def doflabel(n):
    dofs = 8.3e9 / n
    if dofs >= 1e9:
        return '%.1fB' % (dofs/(1e9))
    if dofs >= 1e6:
        return '%.1fM' % (dofs/(1e6))
    elif dofs >= 1e3:
        return '%.1fk' % (dofs/(1e3))
    return '%d' % (dofs)


xticks = data[:, 0].astype(numpy.int32)

ax.set_xticks(xticks)

ax.set_xticklabels(["%s\n(%s)" % (n, doflabel(n))
                    for n in xticks])

ax.plot(data[:, 0], data[:, 1], marker="o", color="black",
        linewidth=2, clip_on=False)

yticks = [0.5, 1, 2, 4, 8, 16, 32, 64]
ax.set_yticks(yticks)
ax.set_yticklabels(map(str, yticks))

ax.set_ylim([0.5, 64])

ax2 = ax.twinx()

ax2.set_ylabel("Efficiency")

ax2.plot(data[:, 0], odata[:, 1],
         marker="s", color="black",
         linewidth=2, linestyle="dashed", clip_on=False)
ax2.yaxis.set_ticks_position("right")

ax2.set_ylim([0.3, 1.1])
seaborn.despine(fig=fig, right=False)

fig.tight_layout()

fig.savefig("mantle-strong-scale.pdf",
            orientation="landscape",
            format="pdf",
            transparent=True)
