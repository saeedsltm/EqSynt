import warnings
from math import sqrt

from numpy import loadtxt, nan, mean, round_
from obspy import read_events
from obspy.geodetics.base import degrees2kilometers as d2k
from pandas import DataFrame

warnings.filterwarnings("ignore")


def handleNone(value, degree=False):
    """Handle missing values

    Args:
        value (float): a float value
        degree (bool, optional): whether convert to degree or not. Defaults to False.

    Returns:
        float: handled value
    """
    if value == None:
        return None
    else:
        if degree:
            return d2k(value)
        return value


def getHer(event):
    """Get horizontal error of event

    Args:
        event (obspy.event): an obspy event

    Returns:
        float: event horizontal error
    """
    if event.origins[0].latitude_errors.uncertainty:
        return round(d2k(sqrt(event.origins[0].latitude_errors.uncertainty**2 + event.origins[0].longitude_errors.uncertainty**2)), 1)
    else:
        return None


def getZer(event):
    """Get depth error of event

    Args:
        event (obspy.event): an obspy event

    Returns:
        float: event depth error
    """
    if event.origins[0].depth_errors.uncertainty:
        return event.origins[0].depth_errors.uncertainty*0.001
    else:
        return None


def catalog2xyzm(hypInp, catalogFileName):
    """Convert catalog to xyzm file format

    Args:
        hypInp (str): file name of NORDIC file
        catalogFileName (str): file name of xyzm.dat file
    """
    cat = read_events(hypInp)
    magnitudes = loadtxt("magnitudes.dat")
    outputFile = "xyzm_{catalogFileName:s}.dat".format(
        catalogFileName=catalogFileName.split(".")[0])
    catDict = {}
    for i, event in enumerate(cat):
        arrivals = event.origins[0].arrivals
        ort = event.origins[0].time
        lat = event.origins[0].latitude
        lon = event.origins[0].longitude
        try:
            dep = event.origins[0].depth*0.001
        except TypeError:
            dep = None
        try:
            nus = event.origins[0].quality.used_station_count
        except AttributeError:
            nus = None
        nuP = len(
            [arrival.phase for arrival in arrivals if "P" in arrival.phase.upper()])
        nuS = len(
            [arrival.phase for arrival in arrivals if "S" in arrival.phase.upper()])
        mds = handleNone(
            min([arrival.distance for arrival in event.origins[0].arrivals]), degree=True)
        ads = round_(handleNone(
            mean([arrival.distance for arrival in event.origins[0].arrivals]), degree=True), 2)
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
            "ORT": ort,
            "Lon": lon,
            "Lat": lat,
            "Dep": dep,
            # "mag":mag,
            "Nus": nus,
            "NuP": nuP,
            "NuS": nuS,
            "ADS": ads,
            "MDS": mds,
            "GAP": gap,
            "RMS": rms,
            "ERH": erh,
            "ERZ": erz,
        }
    df = DataFrame(catDict).T
    df["MAG"] = magnitudes
    df = df.replace({"None": nan})
    with open(outputFile, "w") as f:
        df.to_string(f, index=False)
