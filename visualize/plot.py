import proplot as plt
from numpy import array, arange, corrcoef
from matplotlib.lines import Line2D
import os

# Get global Min and Max of two lists


def getMinMax(ini, fin):
    minIni, maxIni = min(ini), max(ini)
    minFin, maxFin = min(fin), max(fin)
    Min, Max = min([minIni, minFin]), max([maxIni, maxFin])
    return Min, Max

# Plot hypocentral dislocation


def plotHypocenterDiff(ini, fin):
    print("+++ Plotting some statistics ...")
    xMin, xMax = getMinMax(ini["Lon"], fin["Lon"])
    yMin, yMax = getMinMax(ini["Lat"], fin["Lat"])
    zMin, zMax = getMinMax(ini["Dep"], fin["Dep"])
    gMin, gMax = getMinMax(ini["Gap"], fin["Gap"])
    axShape = [
        [1, 2, 3],
        [4, 5, 6]
    ]
    fig, axs = plt.subplots(axShape, share=False)
    axs.format(
        abc=True, abcloc="ul", suptitle="Dislocation plots")
    [ax.grid(ls=":") for ax in axs]
    for i, label in enumerate(["longitude (deg)", "latitude (deg)", "depth (km)", "gap (deg)"]):
        axs[i].set_xlabel("Raw - {label:s}".format(label=label))
        axs[i].set_ylabel("Relocated - {label:s}".format(label=label))
    for i, label in enumerate(["Epicentral", "Depth"], start=4):
        axs[i].set_xlabel("{label:s} error (km)".format(label=label))
        axs[i].set_ylabel("Number of event (#)")
    r = []
    for v in ["Lon", "Lat", "Dep", "Gap"]:
        r.append(corrcoef(ini[v], fin[v])[0][1])

    axs[0].plot(ini["Lon"], fin["Lon"], color="gray",
                mec="k", mew=0.5, marker="o", ls="")
    axs[0].plot([xMin, xMax], [xMin, xMax], color="k", ms=0.5)
    axs[0].set_xlim(xMin, xMax)
    axs[0].set_ylim(xMin, xMax)
    axs[0].format(lrtitle="r={r:f}".format(r=r[0]))

    axs[1].plot(ini["Lat"], fin["Lat"], color="gray",
                mec="k", mew=0.5, marker="o", ls="")
    axs[1].plot([yMin, yMax], [yMin, yMax], color="k", ms=0.5)
    axs[1].set_xlim(yMin, yMax)
    axs[1].set_ylim(yMin, yMax)
    axs[1].format(lrtitle="r={r:f}".format(r=r[1]))

    axs[2].plot(ini["Dep"], fin["Dep"], color="gray",
                mec="k", mew=0.5, marker="o", ls="")
    axs[2].plot([zMin, zMax], [zMin, zMax], color="k", ms=0.5)
    axs[2].set_xlim(zMin, zMax)
    axs[2].set_ylim(zMin, zMax)
    axs[2].format(lrtitle="r={r:f}".format(r=r[2]))

    axs[3].plot(ini["Gap"], fin["Gap"], color="gray",
                mec="k", mew=0.5, marker="o", ls="")
    axs[3].plot([gMin, gMax], [gMin, gMax], color="k", ms=0.5)
    axs[3].set_xlim(gMin, gMax)
    axs[3].set_ylim(gMin, gMax)
    axs[3].format(lrtitle="r={r:f}".format(r=r[3]))

    d = array([ini["ERH"], fin["ERH"]])
    axs[4].hist(d.T, arange(0, 8, 0.2), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur", legend_kw={"ncol": 1})

    d = array([ini["ERZ"], fin["ERZ"]])
    axs[5].hist(d.T, arange(0, 8, 0.2), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur", legend_kw={"ncol": 1})

    fig.save(os.path.join("relocation", "compareHyp.png"))

# Plot velocity models


def plotVelocityModels(velocityModels, maxDep):
    v1 = velocityModels["1"]["V"]
    v2 = velocityModels["2"]["V"]
    z1 = velocityModels["1"]["Z"]
    z2 = velocityModels["2"]["Z"]
    z1.append(maxDep)
    z2.append(maxDep)
    axShape = [
        [1]
    ]
    fig, axs = plt.subplots(axShape, share=False)
    axs.format(
        suptitle="Velocity models", xlabel="Velocity (km/s)", ylabel="Depth (km)")
    [ax.grid(ls=":") for ax in axs]
    for v, z, c in zip([v2, v1], [z2, z1], ["cyan7", "red7"]):
        for x, y1, y2 in zip(v, z[:-1], z[1:]):
            axs[0].vlines(x, y1, y2, color=c, lw=2)
        for x1, x2, y in zip(v[:-1], v[1:], z[1:]):
            axs[0].hlines(y, x1, x2, color=c, lw=2)
    custom_lines = [Line2D([0], [0], color="cyan7", lw=3),
                    Line2D([0], [0], color="red7", lw=3)]
    axs[0].invert_yaxis()
    axs[0].set_xlim(3, 7)
    axs[0].legend(custom_lines, ["Raw", "Relocation"], loc="ll", ncol=1)
    fig.save(os.path.join("relocation", "velocityModel.png"))
