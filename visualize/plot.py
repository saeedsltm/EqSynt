import proplot as plt
from numpy import array, arange
import os

def getMinMax(ini, fin):
    minIni, maxIni = min(ini), max(ini)
    minFin, maxFin = min(fin), max(fin)
    Min, Max = min([minIni, minFin]), max([maxIni, maxFin])
    return Min, Max

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
    for i,label in enumerate(["longitude (deg)", "latitude (deg)", "depth (km)", "gap (deg)"]):
        axs[i].set_xlabel("Raw - {label:s}".format(label=label))
        axs[i].set_ylabel("Relocated - {label:s}".format(label=label))
    for i,label in enumerate(["Epicentral", "Depth"], start=4):
        axs[i].set_xlabel("{label:s} error (km)".format(label=label))
        axs[i].set_ylabel("Number of event (#)")
    
    axs[0].plot(ini["Lon"], fin["Lon"], color="gray", mec="k", mew=0.5, marker="o", ls="")
    axs[0].set_xlim(xMin, xMax)
    axs[0].set_ylim(xMin, xMax)
    
    axs[1].plot(ini["Lat"], fin["Lat"], color="gray", mec="k", mew=0.5, marker="o", ls="")
    axs[1].set_xlim(yMin, yMax)
    axs[1].set_ylim(yMin, yMax)    
    
    axs[2].plot(ini["Dep"], fin["Dep"], color="gray", mec="k", mew=0.5, marker="o", ls="")
    axs[2].set_xlim(zMin, zMax)
    axs[2].set_ylim(zMin, zMax)    
    
    axs[3].plot(ini["Gap"], fin["Gap"], color="gray", mec="k", mew=0.5, marker="o", ls="")
    axs[3].set_xlim(gMin, gMax)
    axs[3].set_ylim(gMin, gMax)
    
    d = array([ini["ERH"], fin["ERH"]])
    axs[4].hist(d.T, arange(0, 8, 0.2), filled=True, alpha=0.7, edgecolor="k",
    cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur")
    
    d = array([ini["ERZ"], fin["ERZ"]])
    axs[5].hist(d.T, arange(0, 8, 0.2), filled=True, alpha=0.7, edgecolor="k",
    cycle=("cyan7", "red7"), labels=["raw", "relocated"], legend="ur")         

    fig.save(os.path.join("relocation", "compareHyp.png"))