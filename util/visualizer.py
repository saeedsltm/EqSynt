import proplot as plt
from numpy import array, arange, corrcoef
from matplotlib.lines import Line2D
from pandas import read_csv, DataFrame
from obspy.geodetics.base import degrees2kilometers as d2k
import os

def getMinMax(*inpList):
    Min = min([min(x) for x in inpList])
    Max = max([max(x) for x in inpList])
    return Min, Max

def plotSeismicityMap(resultPath, stationsDict):
    print("+++ Plotting seismicity map ...")
    # - Read input data
    report_initial = read_csv(os.path.join(resultPath, "relocation", "report_initial.csv"))
    report_select_unweighted = read_csv(os.path.join(resultPath, "relocation", "report_select_unweighted.csv"))
    report_select_weighted = read_csv(os.path.join(resultPath, "relocation", "report_select_weighted.csv"))
    sta = DataFrame(stationsDict).T
    # - Setting global min, max of data
    xMin, xMax = getMinMax(report_initial["Lon"], report_select_unweighted["Lon"], report_select_weighted["Lon"], sta["Lon"])
    yMin, yMax = getMinMax(report_initial["Lat"], report_select_unweighted["Lat"], report_select_weighted["Lat"], sta["Lat"])
    zMin, zMax = getMinMax(report_initial["Dep"], report_select_unweighted["Dep"], report_select_weighted["Dep"])
    # - Define shape of axis
    axShape = [
        [1]
    ]
    # - Some preprocessing on figure
    fig, axs = plt.subplots(axShape, share=False)
    axs.format(
        abc=True, abcloc="ul", suptitle="Seismicity map")
    [ax.grid(ls=":") for ax in axs]

    axs[0].format(xlim=(xMin, xMax), ylim=(yMin, yMax),
                  ylabel="Latitude (deg)", fontsize=7)
    X = [report_initial["Lon"].values, report_select_unweighted["Lon"].values, report_select_weighted["Lon"].values, sta["Lon"].values]
    Y = [report_initial["Lat"].values, report_select_unweighted["Lat"].values, report_select_weighted["Lat"].values, sta["Lat"].values]
    C = ["green", "red", "blue", "black"]
    M = ["o", "o", "o", "^"]
    for x,y,c,m in zip(X, Y, C, M):
            axs[0].plot(x, y, marker=m, ms=4, mec="k", mew=0.1, ls="", color=c)

    px = axs[0].panel_axes(side="r", width="5em")
    px.grid(ls=":")
    px.format(xlim=(zMin, zMax), ylim=(yMin, yMax),
              xlabel="Depth (km)", fontsize=7)
    X = [report_initial["Dep"].values, report_select_unweighted["Dep"].values, report_select_weighted["Dep"].values]
    Y = [report_initial["Lat"].values, report_select_unweighted["Lat"].values, report_select_weighted["Lat"].values]
    C = ["green", "red", "blue"]
    M = ["o", "o", "o"]
    for x,y,c,m in zip(X, Y, C, M):
            px.plot(x, y, marker=m, ms=4, mec="k", mew=0.1, ls="", color=c)              

    px = axs[0].panel_axes(side="b", width="5em")
    px.grid(ls=":")
    px.format(xlim=(xMin, xMax), ylim=(zMax, zMin),
              xlabel="Longitude (deg)", ylabel="Depth (km)", fontsize=7)
    X = [report_initial["Lon"].values, report_select_unweighted["Lon"].values, report_select_weighted["Lon"].values, sta["Lon"].values]
    Y = [report_initial["Dep"].values, report_select_unweighted["Dep"].values, report_select_weighted["Dep"].values]
    C = ["green", "red", "blue"]
    M = ["o", "o", "o"]
    for x,y,c,m in zip(X, Y, C, M):
            px.plot(x, y, marker=m, ms=4, mec="k", mew=0.1, ls="", color=c, autoreverse=False)                 
                 
    # Saving figure
    fig.save(os.path.join(resultPath, "seismicity.pdf"))


def plotHypocenterDiff(resultPath, stationsDict, config):
    print("+++ Plotting some statistics ...")
    # - Read input data
    report_initial = read_csv(os.path.join(resultPath, "relocation", "report_initial.csv"))
    report_select_unweighted = read_csv(os.path.join(resultPath, "relocation", "report_select_unweighted.csv"))
    report_select_weighted = read_csv(os.path.join(resultPath, "relocation", "report_select_weighted.csv"))
    sta = DataFrame(stationsDict).T
    # - Setting global min, max of data
    xMin, xMax = getMinMax(report_initial["Lon"], report_select_unweighted["Lon"], report_select_weighted["Lon"], sta["Lon"])
    yMin, yMax = getMinMax(report_initial["Lat"], report_select_unweighted["Lat"], report_select_weighted["Lat"], sta["Lat"])
    zMin, zMax = getMinMax(report_initial["Dep"], report_select_unweighted["Dep"], report_select_weighted["Dep"])
    gMin, gMax = getMinMax(report_select_unweighted["GAP"], report_select_weighted["GAP"])
    gapMin, gapMax = 0, config["FGS"]["ColorbarGapMax"]
    dMin, dMax = 0, config["FGS"]["ColorbarNearestMax"]
    HistInsetERHMin, HistInsetERHMax = config["FGS"]["HistInsetERHMin"], config["FGS"]["HistInsetERHMax"]
    HistInsetERHInc, HistInsetERZInc = config["FGS"]["HistInsetERHInc"], config["FGS"]["HistInsetERZInc"]
    HistInsetERZMin, HistInsetERZMax = config["FGS"]["HistInsetERZMin"], config["FGS"]["HistInsetERZMax"]
    HistERHMax, HistERZMax = config["FGS"]["HistERHMax"], config["FGS"]["HistERZMax"]
    HistERHInc, HistERZInc = config["FGS"]["HistERHInc"], config["FGS"]["HistERZInc"]
    # - Define shape of axis
    axShape = [
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [9, 10, 11, 12],
    ]
    # - Some preprocessing on figure
    fig, axs = plt.subplots(axShape, share=False)
    axs.format(
        abc=True, abcloc="ul", suptitle="Dislocation plots")
    [ax.grid(ls=":") for ax in axs]
    axesLabels = ["longitude (deg)", "longitude (deg)", "latitude (deg)", "latitude (deg)", "depth (km)", "depth (km)", "Gap (deg)"]
    for i, label in enumerate(axesLabels):
        axs[i].set_xlabel("Raw - {label:s}".format(label=label))
        axs[i].set_ylabel("Relocated - {label:s}".format(label=label))
    for i, label in enumerate(["Epicentral", "Depth"], start=7):
        axs[i].set_xlabel("{label:s} error (km)".format(label=label))
        axs[i].set_ylabel("Number of event (#)")
    ru = [] # for unweighted
    rw = [] # for weighted
    for v in ["Lon", "Lat", "Dep", "GAP"]:
        ru.append(corrcoef(report_initial[v], report_select_unweighted[v])[0][1])
        rw.append(corrcoef(report_initial[v], report_select_weighted[v])[0][1])
    # - Visualizing data
    X = [report_initial["Lon"], report_initial["Lon"]]
    Y = [report_select_unweighted["Lon"], report_select_weighted["Lon"]]
    C = [report_select_unweighted["GAP"], report_select_weighted["GAP"]]
    T = ["$raw-rel{_u} (km)$", "$raw-rel{_w} (km)$"]
    R = [ru, rw]
    for i,(x,y,c,t,r) in enumerate(zip(X, Y, C, T, R)):
        scr1 = axs[i].scatter(x, y, s=50, marker="o", c=c, lw=0.4, edgecolors="k",
                            cmap="lajolla", vmin=gapMin, vmax=gapMax)
        axs[i].plot([xMin, xMax], [xMin, xMax], color="k", ms=0.5)
        axs[i].format(lrtitle="r={r:f}".format(r=r[0]),
                    xlim=(xMin, xMax), ylim=(xMin, xMax))
        ix = axs[i].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
        ix.format(title=t, fontsize=8)
        ix.grid(ls=":")
        ix.spines["right"].set_visible(False)
        ix.spines["top"].set_visible(False)
        data = d2k(x.values-y.values)  # type: ignore
        ix.hist(data, arange(HistInsetERHMin, HistInsetERHMax+1, HistInsetERHInc),
                filled=True, alpha=0.7, edgecolor="k", color="gray")

    # axs[1].scatter(ini["LAT"], fin["LAT"], s=50, marker="o", c=fin["GAP"], lw=0.4, edgecolors="k",
    #                cmap="lajolla", vmin=gapMin, vmax=gapMax)
    # axs[1].plot([yMin, yMax], [yMin, yMax], color="k", ms=0.5)
    # axs[1].set_xlim(yMin, yMax)
    # axs[1].format(lrtitle="r={r:f}".format(r=r[1]),
    #               xlim=(yMin, yMax), ylim=(yMin, yMax))
    # ix = axs[1].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
    # ix.format(title="raw-rel (km)", fontsize=8)
    # ix.grid(ls=":")
    # ix.spines["right"].set_visible(False)
    # ix.spines["top"].set_visible(False)
    # data = d2k(ini["LAT"].values-fin["LAT"].values)
    # ix.hist(data, arange(HistInsetERHMin, HistInsetERHMax+1, HistInsetERHInc),
    #         filled=True, alpha=0.7, edgecolor="k", color="gray")
    X = [report_initial["Lat"], report_initial["Lat"]]
    Y = [report_select_unweighted["Lat"], report_select_weighted["Lat"]]
    C = [report_select_unweighted["GAP"], report_select_weighted["GAP"]]
    T = ["$raw-rel{_u} (km)$", "$raw-rel{_w} (km)$"]
    R = [ru, rw]
    for i,(x,y,c,t,r) in enumerate(zip(X, Y, C, T, R), start=2):
        scr1 = axs[i].scatter(x, y, s=50, marker="o", c=c, lw=0.4, edgecolors="k",
                            cmap="lajolla", vmin=gapMin, vmax=gapMax)
        axs[i].plot([yMin, yMax], [yMin, yMax], color="k", ms=0.5)
        axs[i].format(lrtitle="r={r:f}".format(r=r[1]),
                    xlim=(yMin, yMax), ylim=(yMin, yMax))
        ix = axs[i].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
        ix.format(title=t, fontsize=8)
        ix.grid(ls=":")
        ix.spines["right"].set_visible(False)
        ix.spines["top"].set_visible(False)
        data = d2k(x.values-y.values)  # type: ignore
        ix.hist(data, arange(HistInsetERHMin, HistInsetERHMax+1, HistInsetERHInc),
                filled=True, alpha=0.7, edgecolor="k", color="gray")
    # axs[2].scatter(ini["DEPTH"], fin["DEPTH"], s=50, marker="o", c=fin["GAP"], lw=0.4, edgecolors="k",
    #                cmap="lajolla", vmin=gapMin, vmax=gapMax)
    # axs[2].plot([zMin, zMax], [zMin, zMax], color="k", ms=0.5)
    # axs[2].format(lrtitle="r={r:f}".format(r=r[2]),
    #               xlim=(zMin, zMax), ylim=(zMin, zMax))
    # ix = axs[2].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
    # ix.format(title="raw-rel (km)", fontsize=8)
    # ix.grid(ls=":")
    # ix.spines["right"].set_visible(False)
    # ix.spines["top"].set_visible(False)
    # data = ini["DEPTH"].values-fin["DEPTH"].values
    # ix.hist(data, arange(HistInsetERZMin, HistInsetERZMax+1, HistInsetERZInc),
    #         filled=True, alpha=0.7, edgecolor="k", color="gray")
    X = [report_initial["Dep"], report_initial["Dep"]]
    Y = [report_select_unweighted["Dep"], report_select_weighted["Dep"]]
    C = [report_select_unweighted["GAP"], report_select_weighted["GAP"]]
    T = ["$raw-rel{_u} (km)$", "$raw-rel{_w} (km)$"]
    R = [ru, rw]
    for i,(x,y,c,t,r) in enumerate(zip(X, Y, C, T, R), start=4):
        scr1 = axs[i].scatter(x, y, s=50, marker="o", c=c, lw=0.4, edgecolors="k",
                            cmap="lajolla", vmin=gapMin, vmax=gapMax)
        axs[i].plot([zMin, zMax], [zMin, zMax], color="k", ms=0.5)
        axs[i].format(lrtitle="r={r:f}".format(r=r[2]),
                    xlim=(zMin, zMax), ylim=(zMin, zMax))
        ix = axs[i].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
        ix.format(title=t, fontsize=8)
        ix.grid(ls=":")
        ix.spines["right"].set_visible(False)
        ix.spines["top"].set_visible(False)
        data = d2k(x.values-y.values)  # type: ignore
        ix.hist(data, arange(HistInsetERZMin, HistInsetERZMax+1, HistInsetERZInc),
                filled=True, alpha=0.7, edgecolor="k", color="gray")
    # scr2 = axs[3].scatter(ini["LON"], fin["LON"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
    #                       cmap="lajolla", vmin=dMin, vmax=dMax)
    # axs[3].plot([xMin, xMax], [xMin, xMax], color="k", ms=0.5)
    # axs[3].set_xlim(xMin, xMax)
    # axs[3].set_ylim(xMin, xMax)
    # axs[3].format(lrtitle="r={r:f}".format(r=r[0]))

    scr3 = axs[6].scatter(report_select_unweighted["GAP"], report_select_weighted["GAP"], s=50, marker="o", c=report_select_unweighted["GAP"], lw=0.4, edgecolors="k",
                          cmap="lajolla", vmin=dMin, vmax=dMax)
    axs[6].plot([gMin, gMax], [gMin, gMax], color="k", ms=0.5)
    axs[6].set_xlim(gMin, gMax)
    axs[6].set_ylim(gMin, gMax)
    axs[6].format(lrtitle="r={r:f}".format(r=R[0][3]),
                xlim=(gMin, gMax), ylim=(gMin, gMax))

        # ix = axs[i].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
        # ix.format(title=t, fontsize=8)
        # ix.grid(ls=":")
        # ix.spines["right"].set_visible(False)
        # ix.spines["top"].set_visible(False)
        # data = d2k(x.values-y.values)  # type: ignore
        # ix.hist(data, arange(HistInsetERZMin, HistInsetERZMax+1, HistInsetERZInc),
        #         filled=True, alpha=0.7, edgecolor="k", color="gray")
    # axs[4].scatter(ini["LAT"], fin["LAT"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
    #                cmap="lajolla", vmin=dMin, vmax=dMax)
    # axs[4].plot([yMin, yMax], [yMin, yMax], color="k", ms=0.5)
    # axs[4].set_xlim(yMin, yMax)
    # axs[4].set_ylim(yMin, yMax)
    # axs[4].format(lrtitle="r={r:f}".format(r=r[1]))

    # axs[5].scatter(ini["DEPTH"], fin["DEPTH"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
    #                cmap="lajolla", vmin=dMin, vmax=dMax)
    # axs[5].plot([zMin, zMax], [zMin, zMax], color="k", ms=0.5)
    # axs[5].set_xlim(zMin, zMax)
    # axs[5].set_ylim(zMin, zMax)
    # axs[5].format(lrtitle="r={r:f}".format(r=r[2]))

    # scr3 = axs[6].scatter(ini["GAP"], fin["GAP"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
    #                       cmap="lajolla", vmin=dMin, vmax=dMax)
    # axs[6].plot([gMin, gMax], [gMin, gMax], color="k", ms=0.5)
    # axs[6].set_xlim(gMin, gMax)
    # axs[6].set_ylim(gMin, gMax)
    # axs[6].format(lrtitle="r={r:f}".format(r=r[3]))

    # d = array([ini["ERH"], fin["ERH"]])
    # axs[7].hist(d.T, arange(0, HistERHMax+1, HistERHInc), filled=True, alpha=0.7, edgecolor="k",
    #             cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur", legend_kw={"ncol": 1})

    # d = array([ini["ERZ"], fin["ERZ"]])
    # axs[8].hist(d.T, arange(0, HistERZMax+1, HistERZInc), filled=True, alpha=0.7, edgecolor="k",
    #             cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur", legend_kw={"ncol": 1})

    # # - Colorbar
    # fig.colorbar(
    #     scr1, row=1, loc="r", extend="both", label="Azimuthal gap ($\degree$)", shrink=0.9)
    # fig.colorbar(
    #     scr2, row=2, loc="r", extend="both", label="Nearest station (km)", shrink=0.9)
    # fig.colorbar(
    #     scr3, row=3, loc="r", extend="both", label="Nearest station (km)", shrink=0.9)
    # - Save figure
    fig.save(os.path.join(resultPath, "compareHyp.pdf"))    