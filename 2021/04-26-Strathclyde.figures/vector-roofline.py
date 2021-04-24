import numpy
import pandas as pd
from matplotlib import pyplot

skylake_naive_laplace = pd.read_csv("data/skylake_laplacian_hex_40_1_ve_gcc.csv")
skylake_vector_laplace = pd.read_csv("data/skylake_laplacian_hex_40_8_ve_gcc.csv")
skylake_naive_helmholtz = pd.read_csv("data/skylake_helmholtz_hex_40_1_ve_gcc.csv")
skylake_vector_helmholtz = pd.read_csv("data/skylake_helmholtz_hex_40_8_ve_gcc.csv")
skylake_naive_hyper = pd.read_csv("data/skylake_hyperelasticity_hex_40_1_ve_gcc.csv")
skylake_vector_hyper = pd.read_csv("data/skylake_hyperelasticity_hex_40_8_ve_gcc.csv")

FONTSIZE = 16
fig = pyplot.figure(figsize=(9, 5), frameon=False)
ax = fig.add_subplot(111)

# ax.set_title("Skylake: automated vectorisation roofline", fontsize=FONTSIZE)
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=10)
ax.set_xlim(2**-3, 2**8)
ax.set_ylim(20, 2000)
ax.set_xlabel("Arithmetic intensity [Flops/byte]")
ax.set_ylabel("Throughput [GFlops/s]")


def add_roofline(BW, FLOPS, label=None, **kwargs):
    xes = [2 ** n for n in range(-6, 9)]
    xes = numpy.insert(xes, numpy.searchsorted(xes, FLOPS / BW), FLOPS / BW)
    yes = [min(FLOPS, BW * x) for x in xes]
    ax.plot(xes, yes, linewidth=2, color="black", zorder=1,
            label=label, **kwargs)


def add_derived_metrics(df):
    df["flop"] = df["add"] + df["sub"] + df["mul"] + df["div"]
    df["ai"] = df["flop"] * df["cell"] / df["byte"]
    df["flop/s"] = df["flop"] * df["cell"] / df["time"]


SKYLAKE_BW = 81.04
SKYLAKE_FLOP = 1536
SKYLAKE_LINPACK = 977

add_roofline(SKYLAKE_BW, SKYLAKE_LINPACK, label="LINPACK peak", linestyle="dashed")
add_roofline(SKYLAKE_BW, SKYLAKE_FLOP, label="Roofline")

for df in [skylake_naive_laplace, skylake_vector_laplace,
           skylake_naive_helmholtz, skylake_vector_helmholtz,
           skylake_naive_hyper, skylake_vector_hyper]:
    add_derived_metrics(df)

# ax.plot(skylake_naive_laplace["ai"], skylake_naive_laplace["flop/s"]/1e9, "o", markersize=7, label="Laplace naive")
ax.plot(skylake_vector_laplace["ai"], skylake_vector_laplace["flop/s"]/1e9, "s", markersize=7, label="Laplacian")

# ax.plot(skylake_naive_helmholtz["ai"], skylake_naive_helmholtz["flop/s"]/1e9, "o", markersize=7, label="Helmholtz naive")
# ax.plot(skylake_vector_helmholtz["ai"], skylake_vector_helmholtz["flop/s"]/1e9, "s", markersize=7, label="Helmholtz vectorised")

# ax.plot(skylake_naive_hyper["ai"], skylake_naive_hyper["flop/s"]/1e9, "^", markersize=7, label="Hyperelasticity naive")
ax.plot(skylake_vector_hyper["ai"], skylake_vector_hyper["flop/s"] / 1e9, "v", markersize=7, label="Hyperelasticity")

legend = ax.legend(loc="best", frameon=False)

fig.savefig("vector-roofline.pdf",
            orientation="landscape",
            transparent=True,
            bbox_inches="tight",
            bbox_extra_artists=[legend])
