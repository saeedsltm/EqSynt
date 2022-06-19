import os

import proplot as plt
from numpy import abs, arange, array, corrcoef, loadtxt, sqrt
from obspy.geodetics.base import degrees2kilometers as d2k
from pandas import DataFrame, read_csv
from util.summarizer import summarize

def getMinMax(*inpList):
    Min = min([min(x) for x in inpList])
    Max = max([max(x) for x in inpList])
    return Min, Max


def loadData(resultPath, config):
    report_initial = read_csv(os.path.join(
        resultPath, "relocation", "xyzm_initial.dat"), delim_whitespace=True)
    report_select_unweighted = read_csv(os.path.join(
        resultPath, "relocation", "xyzm_select_unweighted.dat"), delim_whitespace=True)
    report_select_weighted = read_csv(os.path.join(
        resultPath, "relocation", "xyzm_select_weighted.dat"), delim_whitespace=True)
    magnitudes = loadtxt(os.path.join(
        resultPath, "relocation", "magnitudes.dat"))
    report_initial["MAG"] = magnitudes
    report_select_unweighted["MAG"] = magnitudes
    report_select_weighted["MAG"] = magnitudes
    summarize(report_select_unweighted, report_select_weighted, config, resultPath)
    return report_initial, report_select_unweighted, report_select_weighted


def plotSeismicityMap(resultPath, stationsDict, config):
    print("+++ Plotting seismicity map ...")
    # - Read input data
    report_initial, report_select_unweighted, report_select_weighted = loadData(
        resultPath, config)
    sta = DataFrame(stationsDict).T
    # - Setting global min, max of data
    xMin, xMax = getMinMax(
        report_initial["Lon"], report_select_unweighted["Lon"], report_select_weighted["Lon"], sta["Lon"])
    yMin, yMax = getMinMax(
        report_initial["Lat"], report_select_unweighted["Lat"], report_select_weighted["Lat"], sta["Lat"])
    zMin, zMax = getMinMax(
        report_initial["Dep"], report_select_unweighted["Dep"], report_select_weighted["Dep"])
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
    X = [report_initial["Lon"], report_select_unweighted["Lon"],
         report_select_weighted["Lon"]]
    Y = [report_initial["Lat"], report_select_unweighted["Lat"],
         report_select_weighted["Lat"]]
    M = [report_initial["MAG"], report_select_unweighted["MAG"],
         report_select_weighted["MAG"]]
    C = ["green", "red", "blue"]
    L = ["Raw", "Relocated$_u$", "Relocated$_w$"]
    for x, y, m, c, l in zip(X, Y, M, C, L):
        axs[0].scatter(x, y, s=m, marker="o", c=c, lw=0.4,
                       edgecolors="k", alpha=.5, label=l)
    axs[0].plot(sta["Lon"], sta["Lat"], marker="^",
                ms=4, mec="k", mew=0.1, ls="", color="k")

    px = axs[0].panel_axes(side="r", width="5em")
    px.grid(ls=":")
    px.format(xlim=(zMin, zMax), ylim=(yMin, yMax),
              xlabel="Depth (km)", fontsize=7)
    X = [report_initial["Dep"], report_select_unweighted["Dep"],
         report_select_weighted["Dep"]]
    Y = [report_initial["Lat"], report_select_unweighted["Lat"],
         report_select_weighted["Lat"]]
    M = [report_initial["MAG"], report_select_unweighted["MAG"],
         report_select_weighted["MAG"]]
    C = ["green", "red", "blue"]
    for x, y, c, m in zip(X, Y, C, M):
        px.scatter(x, y, s=m, marker="o", c=c,
                   lw=0.4, edgecolors="k", alpha=.5)

    px = axs[0].panel_axes(side="b", width="5em")
    px.grid(ls=":")
    px.format(xlim=(xMin, xMax), ylim=(zMax, zMin),
              xlabel="Longitude (deg)", ylabel="Depth (km)", fontsize=7)
    X = [report_initial["Lon"], report_select_unweighted["Lon"],
         report_select_weighted["Lon"]]
    Y = [report_initial["Dep"], report_select_unweighted["Dep"],
         report_select_weighted["Dep"]]
    M = [report_initial["MAG"], report_select_unweighted["MAG"],
         report_select_weighted["MAG"]]
    C = ["green", "red", "blue"]
    for x, y, c, m in zip(X, Y, C, M):
        px.scatter(x, y, s=m, marker="o", c=c,
                   lw=0.4, edgecolors="k", alpha=.5)

    # Saving figure
    fig.save(os.path.join(resultPath, "seismicity.pdf"))


def plotHypocenterDiff(resultPath, stationsDict, config):
    print("+++ Plotting some statistics ...")
    # - Read input data
    report_initial, report_select_unweighted, report_select_weighted = loadData(
        resultPath, config)
    sta = DataFrame(stationsDict).T
    # - Setting global min, max of data
    xMin, xMax = getMinMax(
        report_initial["Lon"], report_select_unweighted["Lon"], report_select_weighted["Lon"], sta["Lon"])
    yMin, yMax = getMinMax(
        report_initial["Lat"], report_select_unweighted["Lat"], report_select_weighted["Lat"], sta["Lat"])
    zMin, zMax = getMinMax(
        report_initial["Dep"], report_select_unweighted["Dep"], report_select_weighted["Dep"])
    gMin, gMax = getMinMax(
        report_select_unweighted["GAP"], report_select_weighted["GAP"])
    rMin, rMax = getMinMax(
        report_select_unweighted["RMS"], report_select_weighted["RMS"])
    mdsMin, mdsMax = getMinMax(report_initial["MDS"], report_initial["MDS"])
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
    axesLabels = ["longitude (deg)", "longitude (deg)", "latitude (deg)",
                  "latitude (deg)", "depth (km)", "depth (km)"]
    W = ["u", "w", "u", "w", "u", "w"]
    for i, (label, w) in enumerate(zip(axesLabels, W)):
        axs[i].set_xlabel("Raw - {label:s}".format(label=label))
        axs[i].set_ylabel(
            "Relocated$_{w:s}$ - {label:s}".format(label=label, w=w))
    for i, (label, unit) in enumerate(zip(["Gap", "RMS"], ["(deg)", "(s)"]), start=6):
        axs[i].set_xlabel(
            "{label:s} Relocated$_u$ {unit:s}".format(label=label, unit=unit))
        axs[i].set_ylabel(
            "{label:s} Relocated$_w$ {unit:s}".format(label=label, unit=unit))
    for i, label in enumerate(["Real horizontal", "Real depth", "Computed horizontal", "Computed depth"], start=8):
        axs[i].set_xlabel("{label:s} error (km)".format(label=label))
        axs[i].set_ylabel("Number of event (#)")

    ru = []  # for unweighted
    rw = []  # for weighted
    for v in ["Lon", "Lat", "Dep"]:
        ru.append(corrcoef(report_initial[v],
                  report_select_unweighted[v])[0][1])
        rw.append(corrcoef(report_initial[v], report_select_weighted[v])[0][1])
    # - Visualizing data
    X = [report_initial["Lon"], report_initial["Lon"]]
    Y = [report_select_unweighted["Lon"], report_select_weighted["Lon"]]
    C = [report_select_unweighted["GAP"], report_select_weighted["GAP"]]
    T = ["raw-rel$_u$ (km)", "raw-rel$_w$ (km)"]
    R = [ru, rw]
    for i, (x, y, c, t, r) in enumerate(zip(X, Y, C, T, R)):
        r = x.corr(y)
        scr1 = axs[i].scatter(x, y, s=50, marker="o", c=c, lw=0.4, edgecolors="k",
                              cmap="lajolla", vmin=gapMin, vmax=gapMax)
        axs[i].plot([xMin, xMax], [xMin, xMax], color="k", ms=0.5)
        axs[i].format(lrtitle="r={r:f}".format(r=r),
                      xlim=(xMin, xMax), ylim=(xMin, xMax))
        ix = axs[i].inset([0.1, 0.6, 0.3, 0.3], transform="axes", zoom=False)
        ix.format(title=t, fontsize=8)
        ix.grid(ls=":")
        ix.spines["right"].set_visible(False)
        ix.spines["top"].set_visible(False)
        data = d2k(x.values-y.values)  # type: ignore
        ix.hist(data, arange(HistInsetERHMin, HistInsetERHMax+1, HistInsetERHInc),
                filled=True, alpha=0.7, edgecolor="k", color="gray")

    X = [report_initial["Lat"], report_initial["Lat"]]
    Y = [report_select_unweighted["Lat"], report_select_weighted["Lat"]]
    C = [report_select_unweighted["GAP"], report_select_weighted["GAP"]]
    T = ["raw-rel$_u$ (km)", "raw-rel$_w$ (km)"]
    R = [ru, rw]
    for i, (x, y, c, t, r) in enumerate(zip(X, Y, C, T, R), start=2):
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

    X = [report_initial["Dep"], report_initial["Dep"]]
    Y = [report_select_unweighted["Dep"], report_select_weighted["Dep"]]
    C = [report_select_unweighted["MDS"], report_select_weighted["MDS"]]
    T = ["raw-rel$_u$ (km)", "raw-rel$_w$ (km)"]
    R = [ru, rw]
    for i, (x, y, c, t, r) in enumerate(zip(X, Y, C, T, R), start=4):
        scr2 = axs[i].scatter(x, y, s=50, marker="o", c=c, lw=0.4, edgecolors="k",
                              cmap="lajolla", vmin=mdsMin, vmax=mdsMax)
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

    r = []
    for v in ["GAP", "RMS"]:
        r.append(
            corrcoef(report_select_unweighted[v], report_select_weighted[v])[0][1])
    axs[6].scatter(report_select_unweighted["GAP"], report_select_weighted["GAP"], s=50, marker="o",
                   c=report_select_unweighted["MDS"], lw=0.4, edgecolors="k", cmap="lajolla", vmin=mdsMin, vmax=mdsMax)
    axs[6].plot([gMin, gMax], [gMin, gMax], color="k", ms=0.5)
    axs[6].set_xlim(gMin, gMax)
    axs[6].set_ylim(gMin, gMax)
    axs[6].format(lrtitle="r={r:f}".format(r=r[0]),
                  xlim=(gMin, gMax), ylim=(gMin, gMax))

    axs[7].scatter(report_select_unweighted["RMS"], report_select_weighted["RMS"], s=50, marker="o",
                   c=report_select_unweighted["MDS"], lw=0.4, edgecolors="k", cmap="lajolla", vmin=mdsMin, vmax=mdsMax)
    axs[7].plot([rMin, rMax], [rMin, rMax], color="k", ms=0.5)
    axs[7].set_xlim(rMin, rMax)
    axs[7].set_ylim(rMin, rMax)
    axs[7].format(lrtitle="r={r:f}".format(r=r[1]),
                  xlim=(rMin, rMax), ylim=(rMin, rMax))

    d1 = distanceDiff(report_initial["Lon"], report_select_unweighted["Lon"],
                      report_initial["Lat"], report_select_unweighted["Lat"])
    d2 = distanceDiff(report_initial["Lon"], report_select_weighted["Lon"],
                      report_initial["Lat"], report_select_weighted["Lat"])
    d = d2k(array([d1, d2]))
    labels = ["$raw-rel{_u} (km)$", "$raw-rel{_w} (km)$"]
    axs[8].hist(d.T, arange(0, HistERHMax+1, HistERHInc), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=labels, legend="ur", legend_kw={"ncol": 1})

    d1 = abs(report_initial["Lon"] - report_select_unweighted["Lon"])
    d2 = abs(report_initial["Lon"] - report_select_weighted["Lon"])
    d = d2k(array([d1, d2]))
    labels = ["$raw-rel{_u} (km)$", "$raw-rel{_w} (km)$"]
    axs[9].hist(d.T, arange(0, HistERZMax+1, HistERZInc), filled=True, alpha=0.7, edgecolor="k",
                cycle=("cyan7", "red7"), labels=labels, legend="ur", legend_kw={"ncol": 1})

    d1 = report_select_unweighted["ERH"]
    d2 = report_select_weighted["ERH"]
    d = array([d1, d2])
    labels = ["$erh-rel{_u} (km)$", "$erh-rel{_w} (km)$"]
    axs[10].hist(d.T, arange(0, HistERHMax+1, HistERHInc), filled=True, alpha=0.7, edgecolor="k",
                 cycle=("cyan7", "red7"), labels=labels, legend="ur", legend_kw={"ncol": 1})

    d1 = report_select_unweighted["ERZ"]
    d2 = report_select_weighted["ERZ"]
    d = array([d1, d2])
    labels = ["$erz-rel{_u} (km)$", "$erz-rel{_w} (km)$"]
    axs[11].hist(d.T, arange(0, HistERZMax+1, HistERZInc), filled=True, alpha=0.7, edgecolor="k",
                 cycle=("cyan7", "red7"), labels=labels, legend="ur", legend_kw={"ncol": 1})

    # - Colorbar
    fig.colorbar(
        scr1, row=1, loc="r", extend="both", label="Azimuthal gap ($\degree$)", shrink=0.9)
    fig.colorbar(
        scr2, row=2, loc="r", extend="both", label="Min distance to station ($km$)", shrink=0.9)

    # - Save figure
    fig.save(os.path.join(resultPath, "compareHyp.pdf"))


def distanceDiff(xa, xb, ya, yb):
    return sqrt((xa-xb)**2+(ya-yb)**2)
