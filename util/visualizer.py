from matplotlib.pyplot import legend
import proplot as plt
from numpy import array, arange, corrcoef, sqrt, abs, loadtxt
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
    magnitudes = loadtxt(os.path.join(resultPath, "relocation", "magnitudes.dat"))
    report_initial["MAG"] = magnitudes
    report_select_unweighted["MAG"] = magnitudes
    report_select_weighted["MAG"] = magnitudes
    c = (report_select_unweighted["Lon"]>0)&(report_select_weighted["Lon"]>0)
    report_initial = report_initial[c]
    report_select_unweighted = report_select_unweighted[c]
    report_select_weighted = report_select_weighted[c]
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
    X = [report_initial["Lon"], report_select_unweighted["Lon"], report_select_weighted["Lon"]]
    Y = [report_initial["Lat"], report_select_unweighted["Lat"], report_select_weighted["Lat"]]
    M = [report_initial["MAG"], report_select_unweighted["MAG"], report_select_weighted["MAG"]]
    C = ["green", "red", "blue"]
    L = ["Raw", "Relocated$_u$", "Relocated$_w$"]
    for x,y,m,c,l in zip(X, Y, M, C, L):
            # axs[0].plot(x, y, marker=m, ms=4, mec="k", mew=0.1, ls="", color=c)
            axs[0].scatter(x, y, s=m, marker="o", c=c, lw=0.4, edgecolors="k", alpha=.5, label=l)
    axs[0].plot(sta["Lon"], sta["Lat"], marker="^", ms=4, mec="k", mew=0.1, ls="", color="k")

    px = axs[0].panel_axes(side="r", width="5em")
    px.grid(ls=":")
    px.format(xlim=(zMin, zMax), ylim=(yMin, yMax),
              xlabel="Depth (km)", fontsize=7)
    X = [report_initial["Dep"], report_select_unweighted["Dep"], report_select_weighted["Dep"]]
    Y = [report_initial["Lat"], report_select_unweighted["Lat"], report_select_weighted["Lat"]]
    M = [report_initial["MAG"], report_select_unweighted["MAG"], report_select_weighted["MAG"]]
    C = ["green", "red", "blue"]
    for x,y,c,m in zip(X, Y, C, M):
            # px.plot(x, y, marker=m, ms=4, mec="k", mew=0.1, ls="", color=c)              
            px.scatter(x, y, s=m, marker="o", c=c, lw=0.4, edgecolors="k", alpha=.5)

    px = axs[0].panel_axes(side="b", width="5em")
    px.grid(ls=":")
    px.format(xlim=(xMin, xMax), ylim=(zMax, zMin),
              xlabel="Longitude (deg)", ylabel="Depth (km)", fontsize=7)
    X = [report_initial["Lon"], report_select_unweighted["Lon"], report_select_weighted["Lon"]]
    Y = [report_initial["Dep"], report_select_unweighted["Dep"], report_select_weighted["Dep"]]
    M = [report_initial["MAG"], report_select_unweighted["MAG"], report_select_weighted["MAG"]]
    C = ["green", "red", "blue"]
    for x,y,c,m in zip(X, Y, C, M):
            # px.plot(x, y, marker=m, ms=4, mec="k", mew=0.1, ls="", color=c, autoreverse=False)  
            px.scatter(x, y, s=m, marker="o", c=c, lw=0.4, edgecolors="k", alpha=.5)
                 
    # Saving figure
    fig.save(os.path.join(resultPath, "seismicity.pdf"))


def plotHypocenterDiff(resultPath, stationsDict, config):
    print("+++ Plotting some statistics ...")
    # - Read input data
    report_initial = read_csv(os.path.join(resultPath, "relocation", "report_initial.csv"))
    report_select_unweighted = read_csv(os.path.join(resultPath, "relocation", "report_select_unweighted.csv"))
    report_select_weighted = read_csv(os.path.join(resultPath, "relocation", "report_select_weighted.csv"))
    c = (report_select_unweighted["Lon"]>0)&(report_select_weighted["Lon"]>0)
    report_initial = report_initial[c]
    report_select_unweighted = report_select_unweighted[c]
    report_select_weighted = report_select_weighted[c]
    sta = DataFrame(stationsDict).T
    # - Setting global min, max of data
    xMin, xMax = getMinMax(report_initial["Lon"], report_select_unweighted["Lon"], report_select_weighted["Lon"], sta["Lon"])
    yMin, yMax = getMinMax(report_initial["Lat"], report_select_unweighted["Lat"], report_select_weighted["Lat"], sta["Lat"])
    zMin, zMax = getMinMax(report_initial["Dep"], report_select_unweighted["Dep"], report_select_weighted["Dep"])
    gMin, gMax = getMinMax(report_select_unweighted["GAP"], report_select_weighted["GAP"])
    rMin, rMax = getMinMax(report_select_unweighted["RMS"], report_select_weighted["RMS"])
    gapMin, gapMax = 0, config["FGS"]["ColorbarGapMax"]
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
    axesLabels = ["longitude (deg)", "longitude (deg)", "latitude (deg)", "latitude (deg)", "depth (km)", "depth (km)"]
    W = ["u", "w", "u", "w", "u", "w"]
    for i, (label,w) in enumerate(zip(axesLabels, W)):
        axs[i].set_xlabel("Raw - {label:s}".format(label=label))
        axs[i].set_ylabel("Relocated$_{w:s}$ - {label:s}".format(label=label,w=w))
    for i, (label, unit) in enumerate(zip(["Gap", "RMS"], ["(deg)", "(s)"]), start=6):
        axs[i].set_xlabel("{label:s} Relocated$_u$ {unit:s}".format(label=label, unit=unit))
        axs[i].set_ylabel("{label:s} Relocated$_w$ {unit:s}".format(label=label, unit=unit))
    for i, label in enumerate(["Real horizontal", "Real depth", "Computed horizontal", "Computed depth"], start=8):
        axs[i].set_xlabel("{label:s} error (km)".format(label=label))
        axs[i].set_ylabel("Number of event (#)")

    ru = [] # for unweighted
    rw = [] # for weighted
    for v in ["Lon", "Lat", "Dep", "GAP", "RMS"]:
        ru.append(corrcoef(report_initial[v], report_select_unweighted[v])[0][1])
        rw.append(corrcoef(report_initial[v], report_select_weighted[v])[0][1])
    # - Visualizing data
    X = [report_initial["Lon"], report_initial["Lon"]]
    Y = [report_select_unweighted["Lon"], report_select_weighted["Lon"]]
    C = [report_select_unweighted["GAP"], report_select_weighted["GAP"]]
    T = ["raw-rel$_u$ (km)", "raw-rel$_w$ (km)"]
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
    T = ["raw-rel$_u$ (km)", "raw-rel$_w$ (km)"]
    R = [ru, rw]
    for i,(x,y,c,t,r) in enumerate(zip(X, Y, C, T, R), start=2):
        axs[i].scatter(x, y, s=50, marker="o", c=c, lw=0.4, edgecolors="k",
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
    T = ["raw-rel$_u$ (km)", "raw-rel$_w$ (km)"]
    R = [ru, rw]
    for i,(x,y,c,t,r) in enumerate(zip(X, Y, C, T, R), start=4):
        scr2 = axs[i].scatter(x, y, s=50, marker="o", c=c, lw=0.4, edgecolors="k",
                       cmap="lajolla", vmin=gapMin, vmax=gapMax)
        axs[i].plot([zMin, zMax], [zMin, zMax], color="k", ms=0.5)
        axs[i].format(lrtitle="r={r:f}".format(r=r[2]),
                    xlim=(zMin, zMax), ylim=(zMin, zMax))
        ix = axs[i].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
        ix.format(title=t, fontsize=8)
        ix.grid(ls=":")
        ix.spines["right"].set_visible(False)
        ix.spines["top"].set_visible(False)
        data = x.values-y.values  # type: ignore
        ix.hist(data, arange(HistInsetERZMin, HistInsetERZMax+1, HistInsetERZInc),
                filled=True, alpha=0.7, edgecolor="k", color="gray")
    # scr2 = axs[3].scatter(ini["LON"], fin["LON"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
    #                       cmap="lajolla", vmin=dMin, vmax=dMax)
    # axs[3].plot([xMin, xMax], [xMin, xMax], color="k", ms=0.5)
    # axs[3].set_xlim(xMin, xMax)
    # axs[3].set_ylim(xMin, xMax)
    # axs[3].format(lrtitle="r={r:f}".format(r=r[0]))

    r = []
    for v in ["GAP", "RMS"]:
        r.append(corrcoef(report_select_unweighted[v], report_select_weighted[v])[0][1])
    axs[6].scatter(report_select_unweighted["GAP"], report_select_weighted["GAP"], s=50, marker="o",
                   c=report_select_unweighted["GAP"], lw=0.4, edgecolors="k", cmap="lajolla", vmin=gapMin, vmax=gapMax)
    axs[6].plot([gMin, gMax], [gMin, gMax], color="k", ms=0.5)
    axs[6].set_xlim(gMin, gMax)
    axs[6].set_ylim(gMin, gMax)
    axs[6].format(lrtitle="r={r:f}".format(r=r[0]),
                  xlim=(gMin, gMax), ylim=(gMin, gMax))

    # axs[4].scatter(ini["LAT"], fin["LAT"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
    #                cmap="lajolla", vmin=dMin, vmax=dMax)
    # axs[4].plot([yMin, yMax], [yMin, yMax], color="k", ms=0.5)
    # axs[4].set_xlim(yMin, yMax)
    # axs[4].set_ylim(yMin, yMax)
    # axs[4].format(lrtitle="r={r:f}".format(r=r[1]))

    axs[7].scatter(report_select_unweighted["RMS"], report_select_weighted["RMS"], s=50, marker="o",
                   c=report_select_unweighted["GAP"], lw=0.4, edgecolors="k", cmap="lajolla", vmin=gapMin, vmax=gapMax)
    axs[7].plot([rMin, rMax], [rMin, rMax], color="k", ms=0.5)
    axs[7].set_xlim(rMin, rMax)
    axs[7].set_ylim(rMin, rMax)
    axs[7].format(lrtitle="r={r:f}".format(r=r[1]),
                xlim=(rMin, rMax), ylim=(rMin, rMax))

    # axs[5].scatter(ini["DEPTH"], fin["DEPTH"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
    #                cmap="lajolla", vmin=dMin, vmax=dMax)
    # axs[5].plot([zMin, zMax], [zMin, zMax], color="k", ms=0.5)
    # axs[5].set_xlim(zMin, zMax)
    # axs[5].set_ylim(zMin, zMax)
    # axs[5].format(lrtitle="r={r:f}".format(r=r[2]))
    
    d1 = distanceDiff(report_initial["Lon"], report_select_unweighted["Lon"], report_initial["Lat"], report_select_unweighted["Lat"])
    d2 = distanceDiff(report_initial["Lon"], report_select_weighted["Lon"], report_initial["Lat"], report_select_weighted["Lat"])
    d = d2k(array([d1, d2]))
    labels = ["$raw-rel{_u} (km)$", "$raw-rel{_w} (km)$"]
    axs[8].hist(d.T, arange(0, HistERHMax+1, HistERHInc), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=labels, legend="ur", legend_kw={"ncol": 1})

    # scr3 = axs[6].scatter(ini["GAP"], fin["GAP"], s=50, marker="o", c=fin["MIND"], lw=0.4, edgecolors="k",
    #                       cmap="lajolla", vmin=dMin, vmax=dMax)
    # axs[6].plot([gMin, gMax], [gMin, gMax], color="k", ms=0.5)
    # axs[6].set_xlim(gMin, gMax)
    # axs[6].set_ylim(gMin, gMax)
    # axs[6].format(lrtitle="r={r:f}".format(r=r[3]))

    d1 = abs(report_initial["Lon"] - report_select_unweighted["Lon"])
    d2 = abs(report_initial["Lon"] - report_select_weighted["Lon"])
    d = d2k(array([d1, d2]))
    labels = ["$raw-rel{_u} (km)$", "$raw-rel{_w} (km)$"]
    axs[9].hist(d.T, arange(0, HistERZMax+1, HistERZInc), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=labels, legend="ur", legend_kw={"ncol": 1})

    # d = array([ini["ERH"], fin["ERH"]])
    # axs[7].hist(d.T, arange(0, HistERHMax+1, HistERHInc), filled=True, alpha=0.7, edgecolor="k",
    #             cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur", legend_kw={"ncol": 1})

    d1 = sqrt(report_select_unweighted["LonErr"]**2 + report_select_unweighted["LatErr"]**2)
    d2 = sqrt(report_select_weighted["LonErr"]**2 + report_select_weighted["LatErr"]**2)
    d = array([d1, d2])
    labels = ["$HEr-rel{_u} (km)$", "$HEr-rel{_w} (km)$"]
    axs[10].hist(d.T, arange(0, HistERHMax+1, HistERHInc), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=labels, legend="ur", legend_kw={"ncol": 1})

    # d = array([ini["ERZ"], fin["ERZ"]])
    # axs[8].hist(d.T, arange(0, HistERZMax+1, HistERZInc), filled=True, alpha=0.7, edgecolor="k",
    #             cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur", legend_kw={"ncol": 1})

    d1 = report_select_unweighted["DepErr"]
    d2 = report_select_weighted["DepErr"]
    d = array([d1, d2])
    labels = ["$ZEr-rel{_u} (km)$", "$ZEr-rel{_w} (km)$"]
    axs[11].hist(d.T, arange(0, HistERZMax+1, HistERZInc), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=labels, legend="ur", legend_kw={"ncol": 1})

    # # - Colorbar
    fig.colorbar(
        scr1, row=1, loc="r", extend="both", label="Azimuthal gap ($\degree$)", shrink=0.9)
    fig.colorbar(
        scr2, row=2, loc="r", extend="both", label="Azimuthal gap ($\degree$)", shrink=0.9)
    # fig.colorbar(
    #     scr3, row=3, loc="r", extend="both", label="Nearest station (km)", shrink=0.9)
    # - Save figure
    fig.save(os.path.join(resultPath, "compareHyp.pdf"))    


def distanceDiff(xa, xb, ya, yb):
    return sqrt((xa-xb)**2+(ya-yb)**2)
