import numpy
from matplotlib import pyplot
import seaborn

data = numpy.asarray([[5964, 0.084],
                      [8946, 0.120],
                      [13419, 0.173],
                      [20480, 0.25],
                      [27256, 0.341],
                      [35800, 0.423],
                      [47800, 0.513],
                      [72800, 0.662],
                      [109224, 0.843],
                      [163840, 1.004]])

data[:, 0] *= 64
seaborn.set(style="ticks")

fig = pyplot.figure(figsize=(9, 5), frameon=False)

ax = fig.add_subplot(111)

ax.set_xlabel("Number of cores\n(dofs/core)")

ax.loglog()
ax.minorticks_off()

ax.set_ylabel("Simulated years per day")

def doflabel(n):
    dofs = 31.8e9 / n
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

yticks = [0.06, 0.15, 0.3, 0.60, 1.2, 2.4]
ax.set_yticks(yticks)
ax.set_yticklabels(map(str, yticks))

ax.set_ylim([0.06, 2.4])

ax2 = ax.twinx()

ax2.set_ylabel("Efficiency")

ax2.plot(data[:, 0], (data[0, 0] / data[0, 1])/numpy.divide(data[:, 0], data[:, 1]),
         marker="s", color="black",
         linewidth=2, linestyle="dashed", clip_on=False)
ax2.yaxis.set_ticks_position("right")
ax2.set_ylim([0.3, 1.1])

seaborn.despine(fig=fig, right=False)

fig.tight_layout()

fig.savefig("sunway-strong-scale.pdf",
            orientation="landscape",
            format="pdf",
            transparent=True)
