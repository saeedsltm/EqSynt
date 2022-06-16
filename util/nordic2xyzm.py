from obspy import read_events
from obspy.geodetics.base import degrees2kilometers as d2k
from math import sqrt
from numpy import loadtxt
from pandas import DataFrame
import os
import warnings
warnings.filterwarnings("ignore")

def handleNone(value, degree=False):
    if value == None:
        return None
    else:
        if degree:
            return d2k(value)
        return value

def getHer(event):
    if event.origins[0].latitude_errors.uncertainty:
        return round(d2k(sqrt(event.origins[0].latitude_errors.uncertainty**2 + event.origins[0].longitude_errors.uncertainty**2)), 1)
    else:
        return None

def getZer(event):
    if event.origins[0].depth_errors.uncertainty:
        return event.origins[0].depth_errors.uncertainty*0.001
    else:
        return None


def catalog2xyzm(catalogFileName, resultPath):
    cat = read_events(os.path.join(resultPath, "relocation", catalogFileName))
    magnitudes = loadtxt(os.path.join(resultPath, "relocation", "magnitudes.dat"))
    outputFile = os.path.join(resultPath, "relocation", "xyzm.dat")
    catDict = {}
    for i,event in enumerate(cat):
        picks = event.picks
        arrivals = event.origins[0].arrivals
        ort = event.origins[0].time
        lat = event.origins[0].latitude
        lon = event.origins[0].longitude
        try:
            dep = event.origins[0].depth*0.001
        except TypeError:
            dep = None
        # mag = event.magnitudes[0].mag
        try:
            nus = event.origins[0].quality.used_station_count
        except AttributeError:
            nus = None
        nuP = len([arrival.phase for arrival in arrivals if "P" in arrival.phase.upper()])
        nuS = len([arrival.phase for arrival in arrivals if "S" in arrival.phase.upper()])
        # mds = handleNone(event.origins[0].quality.minimum_distance, degree=True)
        mds = handleNone(min([arrival.distance for arrival in event.origins[0].arrivals]), degree=True)
        try:
            gap = handleNone(event.origins[0].quality.azimuthal_gap)
        except AttributeError:
            gap = None
        try:
            rms = handleNone(event.origins[0].quality.standard_error)
        except AttributeError:
            rms = None
        erh = getHer(event)
        erz = getZer(event)
        catDict[i] = {
            "ort":ort,
            "lat":lat,
            "lon":lon,
            "dep":dep,
            # "mag":mag,
            "nus":nus,
            "nuP":nuP,
            "nuS":nuS,
            "mds":mds,
            "gap":gap,
            "rms":rms,
            "erh":erh,
            "erz":erz,
        }
    df = DataFrame(catDict).T
    df["mag"] = magnitudes
    with open(outputFile, "w") as f:
        df.to_string(f, index=False)


catalog2xyzm("hyp.out", "results/test/")



# # Write to "xyzm" format
# def write2xyzm(event, xyzmFile):
#     if len(event["Line1"][0]):
#         ot, lat, lon, dep, mag = event["Line1"][0]
#         ot = ot.strftime("%Y %m %d %H %M %S.%f")[:21]
#         if None in [lat, lon, dep, mag]:
#             formatter = "   nan    nan   nan  nan     nan     nan     nan   nan nan  nan    nan    nan {ot:s}\n"
#                         xyzmFile.write(formatter.format(ot=ot))
#             return
#         gap, erH, erZ = event["LineE"][0]
#         nStUsed = len(
#             list(set([line4[0] for line4 in event["Line4"] if line4[2] != 4])))
#         nPhaseP = len([line4[1] for line4 in event["Line4"]
#                     if (line4[1] == "P" and line4[2] != 4)])
#         nPhaseS = len([line4[1] for line4 in event["Line4"]
#                     if (line4[1] == "S" and line4[2] != 4)])
#         minD = min([line4[6] for line4 in event["Line4"] if line4[2] != 4])
#         rms = getRMS([line4[5] for line4 in event["Line4"] if line4[2] != 4])
#         formatter = "{lon:6.3f} {lat:6.3f} {dep:5.1f} {mag:4.1f} {nStUsed:7d} {nPhaseP:7d} {nPhaseS:7d} {minD:5.1f} {gap:3.0f} {rms:4.2f} {erH:6.2f} {erZ:6.2f} {ot:s}\n"
#         xyzmFile.write(formatter.format(
#             lon=lon, lat=lat, dep=dep, mag=mag, nStUsed=nStUsed, nPhaseP=nPhaseP, nPhaseS=nPhaseS, minD=minD, gap=gap, rms=rms, erH=erH, erZ=erZ, ot=ot
#         ))














