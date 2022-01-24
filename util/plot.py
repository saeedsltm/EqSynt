import proplot as plt
from numpy import array, arange, corrcoef
from matplotlib.lines import Line2D
from pandas import read_csv, DataFrame
from obspy.geodetics.base import degrees2kilometers as d2k
import os

# Get global Min and Max of two lists


def getMinMax(*inpList):
    Min = min([min(x) for x in inpList])
    Max = max([max(x) for x in inpList])
    # minIni, maxIni = min(ini), max(ini)
    # minFin, maxFin = min(fin), max(fin)
    # Min, Max = min([minIni, minFin]), max([maxIni, maxFin])
    return Min, Max

# plot hypocentral map between initial and final event locations
def plotHypocentralMap(iniFile, finFile, stationFile, relocationPath):
    print("+++ Plotting seismicity map ...")
    # - Read input data
    ini = read_csv(iniFile, delim_whitespace=True)
    fin = read_csv(finFile, delim_whitespace=True)
    sta = DataFrame(stationFile)
    # - Setting global min, max of data
    xMin, xMax = getMinMax(ini["LON"], fin["LON"], sta["Lon"])
    yMin, yMax = getMinMax(ini["LAT"], fin["LAT"], sta["Lat"])
    zMin, zMax = getMinMax(ini["DEPTH"], fin["DEPTH"])
    # - Define shape of axis
    axShape = [
        [1]
    ]
    # - Somepreproccessing on figure
    fig, axs = plt.subplots(axShape, share=False)
    axs.format(
        abc=True, abcloc="ul", suptitle="Seismicity map")
    [ax.grid(ls=":") for ax in axs]

    axs[0].format(xlim=(xMin, xMax), ylim=(yMin, yMax), ylabel="Latitude (deg)", fontsize=7)
    axs[0].plot(ini["LON"].values, ini["LAT"].values, marker="o", ms=4, mec="k", mew=0.1, ls="", color="gray", alpha=0.2)
    axs[0].plot(fin["LON"].values, fin["LAT"].values, marker="o", ms=4, mec="k", mew=0.1, ls="", color="red")
    axs[0].plot(sta["Lon"].values, sta["Lat"].values, marker="^", ms=4, mec="k", mew=0.1, ls="", color="blue")
    
    px = axs[0].panel_axes(side="r", width="5em")
    px.grid(ls=":")
    px.format(xlim=(zMin, zMax), ylim=(yMin, yMax), xlabel="Depth (km)", fontsize=7)
    px.plot(ini["DEPTH"].values, ini["LAT"].values, marker="o", ms=4, mec="k", mew=0.1, ls="", color="gray", alpha=0.2)
    px.plot(fin["DEPTH"].values, fin["LAT"].values, marker="o", ms=4, mec="k", mew=0.1, ls="", color="red")

    px = axs[0].panel_axes(side="b", width="5em")
    px.grid(ls=":")
    px.format(xlim=(xMin, xMax), ylim=(zMax, zMin), xlabel="Longitude (deg)", ylabel="Depth (km)", fontsize=7)
    px.plot(ini["LON"].values, ini["DEPTH"].values, marker="o", ms=4, mec="k", mew=0.1, ls="", color="gray", alpha=0.2, autoreverse=False)
    px.plot(fin["LON"].values, fin["DEPTH"].values, marker="o", ms=4, mec="k", mew=0.1, ls="", color="red", autoreverse=False)
    # Saving figure
    fig.save(os.path.join(relocationPath, "seismicity.pdf"))
    

# Plot hypocentral dislocation


def plotHypocenterDiff(iniFile, finFile, relocationPath):
    print("+++ Plotting some statistics ...")
    # - Read input data
    ini = read_csv(iniFile, delim_whitespace=True)
    fin = read_csv(finFile, delim_whitespace=True)
    # - Setting global min, max of data
    xMin, xMax = getMinMax(ini["LON"], fin["LON"])
    yMin, yMax = getMinMax(ini["LAT"], fin["LAT"])
    zMin, zMax = getMinMax(ini["DEPTH"], fin["DEPTH"])
    gMin, gMax = getMinMax(ini["GAP"], fin["GAP"])
    gapMin, gapMax = 0, 360
    dMin, dMax = 0, 10
    # - Define shape of axis
    axShape = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
    ]
    # - Somepreproccessing on figure
    fig, axs = plt.subplots(axShape, share=False)
    axs.format(
        abc=True, abcloc="ul", suptitle="Dislocation plots")
    [ax.grid(ls=":") for ax in axs]
    for i, label in enumerate(["longitude (deg)", "latitude (deg)", "depth (km)", "longitude (deg)", "latitude (deg)", "depth (km)", "gap (deg)"]):
        axs[i].set_xlabel("Raw - {label:s}".format(label=label))
        axs[i].set_ylabel("Relocated - {label:s}".format(label=label))
    for i, label in enumerate(["Epicentral", "Depth"], start=7):
        axs[i].set_xlabel("{label:s} error (km)".format(label=label))
        axs[i].set_ylabel("Number of event (#)")
    r = []
    for v in ["LON", "LAT", "DEPTH", "GAP"]:
        r.append(corrcoef(ini[v], fin[v])[0][1])
    # - Visualizing data
    scr1 = axs[0].scatter(ini["LON"], fin["LON"], s=50, marker="o", c=fin["GAP"], lw=0.4, edgecolors="k",
                          cmap="lajolla", vmin=gapMin, vmax=gapMax)
    axs[0].plot([xMin, xMax], [xMin, xMax], color="k", ms=0.5)
    axs[0].format(lrtitle="r={r:f}".format(r=r[0]), xlim=(xMin, xMax), ylim=(xMin, xMax))
    ix = axs[0].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
    ix.format(title="raw-rel (km)", fontsize=8)
    ix.grid(ls=":")
    ix.spines["right"].set_visible(False)
    ix.spines["top"].set_visible(False)
    data = d2k(ini["LON"].values-fin["LON"].values)
    ix.hist(data, arange(-5, 6, 0.4), filled=True, alpha=0.7, edgecolor="k", color="gray")
    

    axs[1].scatter(ini["LAT"], fin["LAT"], s=50, marker="o", c=fin["GAP"], lw=0.4, edgecolors="k",
                   cmap="lajolla", vmin=gapMin, vmax=gapMax)
    axs[1].plot([yMin, yMax], [yMin, yMax], color="k", ms=0.5)
    axs[1].set_xlim(yMin, yMax)
    axs[1].format(lrtitle="r={r:f}".format(r=r[1]), xlim=(yMin, yMax), ylim=(yMin, yMax))
    ix = axs[1].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
    ix.format(title="raw-rel (km)", fontsize=8)
    ix.grid(ls=":")
    ix.spines["right"].set_visible(False)
    ix.spines["top"].set_visible(False)
    data = d2k(ini["LAT"].values-fin["LAT"].values)
    ix.hist(data, arange(-5, 6, 0.4), filled=True, alpha=0.7, edgecolor="k", color="gray")

    axs[2].scatter(ini["DEPTH"], fin["DEPTH"], s=50, marker="o", c=fin["GAP"], lw=0.4, edgecolors="k",
                   cmap="lajolla", vmin=gapMin, vmax=gapMax)
    axs[2].plot([zMin, zMax], [zMin, zMax], color="k", ms=0.5)
    axs[2].format(lrtitle="r={r:f}".format(r=r[2]), xlim=(zMin, zMax), ylim=(zMin, zMax))
    ix = axs[2].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
    ix.format(title="raw-rel (km)", fontsize=8)
    ix.grid(ls=":")
    ix.spines["right"].set_visible(False)
    ix.spines["top"].set_visible(False)
    data = ini["DEPTH"].values-fin["DEPTH"].values
    ix.hist(data, arange(-5, 6, 0.5), filled=True, alpha=0.7, edgecolor="k", color="gray")

    scr2 = axs[3].scatter(ini["LON"], fin["LON"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
                          cmap="lajolla", vmin=dMin, vmax=dMax)
    axs[3].plot([xMin, xMax], [xMin, xMax], color="k", ms=0.5)
    axs[3].set_xlim(xMin, xMax)
    axs[3].set_ylim(xMin, xMax)
    axs[3].format(lrtitle="r={r:f}".format(r=r[0]))

    axs[4].scatter(ini["LAT"], fin["LAT"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
                   cmap="lajolla", vmin=dMin, vmax=dMax)
    axs[4].plot([yMin, yMax], [yMin, yMax], color="k", ms=0.5)
    axs[4].set_xlim(yMin, yMax)
    axs[4].set_ylim(yMin, yMax)
    axs[4].format(lrtitle="r={r:f}".format(r=r[1]))

    axs[5].scatter(ini["DEPTH"], fin["DEPTH"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
                   cmap="lajolla", vmin=dMin, vmax=dMax)
    axs[5].plot([zMin, zMax], [zMin, zMax], color="k", ms=0.5)
    axs[5].set_xlim(zMin, zMax)
    axs[5].set_ylim(zMin, zMax)
    axs[5].format(lrtitle="r={r:f}".format(r=r[2]))

    scr3 = axs[6].scatter(ini["GAP"], fin["GAP"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
                          cmap="lajolla", vmin=dMin, vmax=dMax)
    axs[6].plot([gMin, gMax], [gMin, gMax], color="k", ms=0.5)
    axs[6].set_xlim(gMin, gMax)
    axs[6].set_ylim(gMin, gMax)
    axs[6].format(lrtitle="r={r:f}".format(r=r[3]))

    d = array([ini["ERH"], fin["ERH"]])
    axs[7].hist(d.T, arange(0, 8, 0.2), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur", legend_kw={"ncol": 1})

    d = array([ini["ERZ"], fin["ERZ"]])
    axs[8].hist(d.T, arange(0, 8, 0.2), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur", legend_kw={"ncol": 1})

    # - Colorbars
    fig.colorbar(
        scr1, row=1, loc="r", extend="both", label="Azimuthal gap ($\degree$)", shrink=0.9)
    fig.colorbar(
        scr2, row=2, loc="r", extend="both", label="Nearest station (km)", shrink=0.9)
    fig.colorbar(
        scr3, row=3, loc="r", extend="both", label="Nearest station (km)", shrink=0.9)
    # - Save figure
    fig.save(os.path.join(relocationPath, "compareHyp.pdf"))

# Plot velocity models


def plotVelocityModels(velocityModels, maxDep, relocationPath):
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
    fig.save(os.path.join(relocationPath, "velocityModel.pdf"))
