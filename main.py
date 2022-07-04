import os
import random
import sys
from glob import glob
from json import load
from pathlib import Path
from shutil import copy
from string import ascii_uppercase

import pykonal
from LatLon import lat_lon as ll
from numpy import array
from numpy.random import multivariate_normal, seed
from obspy import UTCDateTime as utc
from obspy.geodetics.base import gps2dist_azimuth as gps
from tqdm.contrib import tzip

from util import feedEvents2Catalog
from util.logger import myLogger
from util.nordic2xyzm import catalog2xyzm
from util.synthTime import extractTT, generateTTT
from util.visualizer import plotHypocenterDiff, plotSeismicityMap

# Define main class


class Main():
    # Init
    def __init__(self):
        self.resultsPath = os.path.join(
            "results", utc.now().strftime("%Y_%j_%H%M%S"))
        self.resultsPath = "results/2022_185_134736"
        Path(self.resultsPath).mkdir(parents=True, exist_ok=True)
        self.report = myLogger(os.path.join(
            self.resultsPath, "report"), mode="w")
        self.readConfiguration()
        self.config2Log()
        self.checkRequiredFiles()
        self.clearExistingFiles()
        self.node_intervals = [self.config["TravelTimeGenerator"]["XGridIntervals"],
                               self.config["TravelTimeGenerator"]["YGridIntervals"],
                               self.config["TravelTimeGenerator"]["ZGridIntervals"]]
        self.travelTimeDict = {}

    # Read configuration parameters
    def readConfiguration(self):
        # - check if configuration file exists
        if not os.path.exists("config.json"):
            msg = "+++ Could not find configuration file! Aborting ..."
            print(msg)
            self.report.info(msg)
            sys.exit()
        with open("config.json") as f:
            self.config = load(f)
        msg = "+++ Configuration file was loaded successfully."
        print(msg)
        self.report.info(msg)

    # Write configuration to log file
    def config2Log(self):
        for key in self.config:
            if isinstance(self.config[key], dict):
                for k, v in self.config[key].items():
                    msg = "+++ Configuration on {key} > {k}: {v}".format(
                        key=key, k=k, v=v)
                    self.report.info(msg)
            else:
                msg = "+++ {v} :".format(v=self.config[key])
                self.report.info(msg)

    # Check required input files
    def checkRequiredFiles(self):
        # - check earthquake arrivals and station files
        if self.config["FSS"]["Flag"] and not self.config["RSS"]["Flag"]:
            self.createNewStationFile()
        elif not self.config["FSS"]["Flag"] and self.config["RSS"]["Flag"]:
            self.eqFile = os.path.join(
                "EqInput", self.config["RSS"]["Inputs"]["EqFile"])
            self.stationFile = os.path.join(
                "EqInput", self.config["RSS"]["Inputs"]["StationFile"])
            if not os.path.exists(self.eqFile) or not os.path.exists(self.stationFile):
                msg = "+++ Could not find earthquake or station files! Aborting ..."
                print(msg)
                self.report.info(msg)
                sys.exit()        # - Create NLloc required directories
        Path("time").mkdir(parents=True, exist_ok=True)
        Path("model").mkdir(parents=True, exist_ok=True)
        Path("obs").mkdir(parents=True, exist_ok=True)
        Path(os.path.join(self.resultsPath, "relocation")).mkdir(
            parents=True, exist_ok=True)
        Path("tmp").mkdir(parents=True, exist_ok=True)

    # Clear pre-existing files
    def clearExistingFiles(self):
        for d in ["model", "obs", "time", "tmp"]:
            for f in glob(os.path.join(d, "*")):
                os.remove(f)
        msg = "+++ Removing preexisting files in 'model', 'obs', 'time' and 'tmp' directories ..."
        print(msg)
        self.report.info(msg)

    def readVelocityFile(self, stationFile):
        """reading velocity model from "STATION0.HYP" file 

        Args:
            stationFile (str): station file in "NORDIC" format

        Returns:
            (dict): a dictionary contains velocity model
        """
        emptyLines = 0
        velocityModelDict = {"Vp": [], "Z": [], "VpVs": 1.73, "Moho": 46.0}
        with open(stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 2 and l.strip():
                    Vp, Z = [float(x) for x in l.split()[:2]]
                    velocityModelDict["Vp"].append(Vp)
                    velocityModelDict["Z"].append(Z)
                if emptyLines == 2 and len(l) > 20 and l[21] == "N":
                    _, Z = [float(x) for x in l.split()[:2]]
                    velocityModelDict["Moho"] = Z
                if emptyLines == 3 and l.strip():
                    VpVs = float(l[16:20])
                    velocityModelDict["VpVs"] = VpVs
                    break
        return velocityModelDict

    def readStationFile(self, stationFile):
        """read station information from "STATION0.HYP" file

        Args:
            stationFile (str): station file in "NORDIC" format

        Returns:
            dict: a dictionary contains stations information
        """
        emptyLines = 0
        stations = {}
        with open(stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 1 and l.strip():
                    code = l[:6].strip()
                    lat = ll.Latitude(degree=int(
                        l[6:8]), minute=float(l[8:13])).decimal_degree
                    lon = ll.Longitude(degree=int(
                        l[15:17]), minute=float(l[17:22])).decimal_degree
                    elv = float(l[23:27])
                    stations[code] = {"Lat": lat, "Lon": lon, "Elv": elv}
        return stations

    def generateTTTable(self):
        """generate travel time tables and store a bank of files
        """
        # reading velocity file
        velocityModelDict = self.readVelocityFile(
            os.path.join("EqInput", "STATION0.HYP"))
        if self.config["TravelTimeGenerator"]["Flag"]:
            xMaxDist = self.config["TravelTimeGenerator"]["NumXGrid"] * \
                self.config["TravelTimeGenerator"]["XGridIntervals"]
            zMaxDist = self.config["TravelTimeGenerator"]["NumZGrid"] * \
                self.config["TravelTimeGenerator"]["ZGridIntervals"]
            print("+++ Generating travel time tables for X range (0, {xMax:4.0f})km, and Z range (0, {zMax:4.0f})km".format(
                xMax=xMaxDist, zMax=zMaxDist))
            print("+++ Generating travel time tables for P phase ...")
            generateTTT(
                velocityModelDict, "P",
                self.config["TravelTimeGenerator"]["NumXGrid"], self.config["TravelTimeGenerator"]["NumZGrid"],
                self.node_intervals, self.config["TravelTimeGenerator"]["DecimationFactor"]
            )
            print("+++ Generating travel time tables for S phase ...")
            generateTTT(
                velocityModelDict, "S",
                self.config["TravelTimeGenerator"]["NumXGrid"], self.config["TravelTimeGenerator"]["NumZGrid"],
                self.node_intervals, self.config["TravelTimeGenerator"]["DecimationFactor"]
            )

    # Create new station file
    def createNewStationFile(self):
        """Generating new STATION0.HYP file. It will produce the followings:
        - Inheriting RESET defaults from SEISAN,
        - New random station coordinates,
        - New velocity model.
        """
        stationCodes, stationLats, stationLons = [], [], []
        nStations = self.config["FSS"]["Station"]["NumberOfStations"]
        stationCodeLength = 4
        random.seed(self.config["FSS"]["Station"]["SeedRandomID"])
        stationCodes = ["".join(random.sample(
            ascii_uppercase, k=stationCodeLength)) for _ in range(nStations)]
        StationBoundCenter = self.config["FSS"]["Station"]["StationBoundCenter"]
        StationBoundCov = self.config["FSS"]["Station"]["StationBoundCov"]
        seed(self.config["FSS"]["Station"]["SeedRandomID"])
        stationLons, stationLats = multivariate_normal(
            StationBoundCenter, StationBoundCov, nStations).T
        Vp = self.config["FSS"]["Model"]["Vp"]
        Z = self.config["FSS"]["Model"]["Z"]
        Interfaces = self.config["FSS"]["Model"]["Interfaces"]
        VpVs = self.config["FSS"]["Model"]["VpVs"]
        xNear = self.config["FSS"]["Model"]["XNear"]
        xFar = self.config["FSS"]["Model"]["XFar"]
        trialDepth = self.config["FSS"]["Model"]["TrialDepth"]
        with open(os.path.join("files", "resets.dat")) as f, open(os.path.join("EqInput", "STATION0.HYP"), "w") as g:
            for l in f:
                g.write(l)
            g.write("\n")
            g.write("\n")
            for code, lat, lon in zip(stationCodes, stationLats, stationLons):
                lat = ll.Latitude(lat)
                lon = ll.Longitude(lon)
                elv = 0.0
                g.write("  {code:4s}{latDeg:2.0f}{latMin:5.2f}N {lonDeg:2.0f}{lonMin:5.2f}E{elv:4.0f}\n".format(
                    code=code,
                    latDeg=lat.decimal_degree, latMin=lat.decimal_minute,
                    lonDeg=lon.decimal_degree, lonMin=lon.decimal_minute,
                    elv=elv
                ))
            g.write("\n")
            for v, z, i in zip(Vp, Z, Interfaces):
                g.write(" {v:5.2f}  {z:6.3f}       {i:1s}     \n".format(
                    v=v, z=z, i=i
                ))
            g.write("\n")
            g.write("{trialDepth:4.0f}.{xNear:4.0f}.{xFar:4.0f}. {VpVs:4.2f}".format(
                trialDepth=trialDepth, xNear=xNear, xFar=xFar, VpVs=VpVs
            ))
            g.write("\nNew")

    def createNewCatalogFile(self, computeWeight=True):
        """Generate new catalog file. It will produce the followings:
        - New catalog based on user defined normal distribution,
        - Write new catalog in NORDIC format.
        """
        stationsDict = self.readStationFile(
            os.path.join("EqInput", "STATION0.HYP"))
        eventLons, eventLats, eventDeps = [], [], []
        eventsArrivals = []
        eventsInfo = []
        NumberOfEvents = self.config["FSS"]["Catalog"]["NumberOfEvents"]
        CatalogBoundCenter = self.config["FSS"]["Catalog"]["CatalogBoundCenter"]
        CatalogBoundCov = self.config["FSS"]["Catalog"]["CatalogBoundCov"]
        seed(self.config["FSS"]["Catalog"]["SeedRandomID"])
        eventLons, eventLats = multivariate_normal(
            CatalogBoundCenter, CatalogBoundCov, NumberOfEvents).T
        CatalogBoundDepMin = self.config["FSS"]["Catalog"]["CatalogBoundDepMin"]
        CatalogBoundDepMax = self.config["FSS"]["Catalog"]["CatalogBoundDepMax"]
        random.seed(self.config["FSS"]["Catalog"]["SeedRandomID"])
        eventDeps = [random.uniform(
            CatalogBoundDepMin, CatalogBoundDepMax) for _ in range(NumberOfEvents)]
        print(
            "+++ Generating new catalog, weighting={computeWeight} ...".format(computeWeight=computeWeight))
        eventOt = utc("2020-01-01")
        eventsTimeShift = 100.0
        for i, (eventLat, eventLon, eventDep) in enumerate(tzip(eventLats, eventLons, eventDeps)):  # type: ignore
            newEventArrivals = {}
            eventOt += eventsTimeShift
            newEventInfo = {"OriginTime": eventOt, "Longitude": eventLon,
                            "Latitude": eventLat, "Depth": eventDep}
            eventsInfo.append(newEventInfo)
            selectedStations, azimuths, receivers = self.computeNewDistances(
                i, eventLat, eventLon, stationsDict)
            hdf5FileP = os.path.join(
                "ttt", "dep{z:003d}{vt:s}.hdf5".format(z=int(eventDep), vt="P"))
            if hdf5FileP not in self.travelTimeDict.keys():
                traveltime = pykonal.fields.read_hdf(  # type: ignore
                    hdf5FileP)
                self.travelTimeDict[hdf5FileP] = traveltime
            extractedTTP = extractTT(
                self.travelTimeDict[hdf5FileP], receivers)
            hdf5FileS = os.path.join(
                "ttt", "dep{z:003d}{vt:s}.hdf5".format(z=int(eventDep), vt="S"))
            if hdf5FileS not in self.travelTimeDict.keys():
                traveltime = pykonal.fields.read_hdf(  # type: ignore
                    hdf5FileS)
                self.travelTimeDict[hdf5FileS] = traveltime
            extractedTTS = extractTT(
                self.travelTimeDict[hdf5FileS], receivers)
            for station, ttP, ttS, distance, azimuth in zip(selectedStations, extractedTTP, extractedTTS, receivers, azimuths):
                travelTimeP, weightP = self.setTTError(
                    eventOt+ttP, m=0.0, std=self.config["FPS"]["ErrorOnArrivals"]["P"], addWeight=computeWeight)
                travelTimeS, weightS = self.setTTError(
                    eventOt+ttS, m=0.0, std=self.config["FPS"]["ErrorOnArrivals"]["S"], addWeight=computeWeight)
                newEventArrivals[station] = {
                    "P": ["P", "emergent", weightP, travelTimeP, distance[0], azimuth],
                    "S": ["S", "emergent", weightS, travelTimeS, distance[0], azimuth]
                }
            eventsArrivals.append(newEventArrivals)
        feedCatalog = feedEvents2Catalog.feedCatalog()
        newCatalog = feedCatalog.setCatalog(eventsInfo, eventsArrivals)
        self.writeCatalog(newCatalog, computeWeight)

    def writeCatalog(self, catalog, computeWeight):
        """Write catalog of earthquakes in NORDIC format.

        Args:
            catalog (obspy.catalog): an obspy catalog of events,
            computeWeight (bool): whether or not computing weights.
        """
        if computeWeight:
            catalog.write(os.path.join(
                "EqInput", "select_weighted.out"), format="NORDIC")
        else:
            catalog.write(os.path.join(
                "EqInput", "select_unweighted.out"), format="NORDIC")
        with open(os.path.join("EqInput", "magnitudes.dat"), "w") as f:
            for event in catalog:
                mag = event.magnitudes[0].mag
                f.write("{mag:3.1f}\n".format(mag=mag))

    def setTTError(self, tt, m=0.0, std=1.0, addWeight=True):
        """Setting error with travel times.

        Args:
            tt (float): travel time of a desired phase (s),
            m (float, optional): mean of Gaussian distribution for generating error. Defaults to 0.0.
            std (float, optional): standard deviation of Gaussian distribution for generating error. Defaults to 1.0.
            addWeight (bool, optional): whether or not considering weights.. Defaults to True.

        Returns:
            tuple: a tuple contains of travel time and associated weight.
        """
        err = random.gauss(m, std)
        tt = tt + err
        w = self.mapError2Weight(err, m, std, addWeight)
        return tt, w

    def mapError2Weight(self, err, m, std, addWeight=True):
        """Map error magnitude to weight.

        Args:
            err (float): error in seconds,
            m (float, optional): mean of Gaussian distribution for generating error. Defaults to 0.0.
            std (float, optional): standard deviation of Gaussian distribution for generating error. Defaults to 1.0.
            addWeight (bool, optional): whether or not considering weights.. Defaults to True.

        Returns:
            _type_: _description_
        """
        if not addWeight:
            return 0
        absError = abs(abs(std) - abs(m))
        if 0.0 <= abs(err) <= 0.25*absError:
            return 0
        elif 0.25*absError <= abs(err) <= 0.5*absError:
            return 1
        elif 0.5*absError <= abs(err) <= 0.75*absError:
            return 2
        elif 0.75*absError <= abs(err) <= 1.0*absError:
            return 3
        elif 1.0*absError <= abs(err):
            return 4

    def computeNewDistances(self, eventID, eventLat, eventLon, stations):
        """Compute new distances from epicenter to stations

        Args:
            eventID (int): ID of an event,
            eventLat (float): latitude of an event,
            eventLon (float): longitude of an event,
            stations (dict): a dictionary contains of station information.

        Returns:
            tuple: a tuple contains a list of stations, associated azimuths and distances.
        """
        random.seed(eventID)
        MinNumberOfUsedStation = 5
        randomStations = random.sample(
            list(stations.keys()), k=random.randint(MinNumberOfUsedStation, len(stations)))
        stationLats = [stations[station]["Lat"] for station in randomStations]
        stationLons = [stations[station]["Lon"] for station in randomStations]
        distances = [gps(eventLat, eventLon, stationLat, stationLon)[
            0]*1e-3 for (stationLat, stationLon) in zip(stationLats, stationLons)]
        azimuths = [gps(eventLat, eventLon, stationLat, stationLon)[
            1] for (stationLat, stationLon) in zip(stationLats, stationLons)]
        return randomStations, azimuths, [array([distance, 1.0, 0.0]) for distance in distances]

    def makeReportFile(self):
        """Make report in xyzm.dat file format.
        """
        for selectFile in glob(os.path.join("EqInput", "select_*.out")):
            copy(selectFile, os.path.join(self.resultsPath, "relocation"))
        copy(os.path.join("EqInput", "STATION0.HYP"),
             os.path.join(self.resultsPath, "relocation"))
        copy(os.path.join("EqInput", "magnitudes.dat"),
             os.path.join(self.resultsPath, "relocation"))
        copy(os.path.join("files", "report.inp"),
             os.path.join(self.resultsPath, "relocation"))
        root = os.getcwd()
        os.chdir(os.path.join(self.resultsPath, "relocation"))
        for inpFile in glob("select*.out"):
            with open("hyp.inp", "w") as f:
                f.write("{inpFile:s}\nn\n".format(inpFile=inpFile))
            cmd = "hyp < hyp.inp > /dev/null"
            os.system(cmd)
            catalog2xyzm("hyp.out", inpFile)
        initial = "select_unweighted.out"
        catalog2xyzm(initial, "initial.out")
        os.chdir(root)

    def visualizeResult(self):
        """Visualize results by generating some maps.
        """
        stationsDict = self.readStationFile(
            os.path.join("EqInput", "STATION0.HYP"))
        plotSeismicityMap(self.resultsPath, stationsDict,
                          self.config)  # type: ignore
        plotHypocenterDiff(self.resultsPath, stationsDict,
                           self.config)  # type: ignore


# Run application
if __name__ == "__main__":
    app = Main()
    # app.generateTTTable()
    # app.createNewStationFile()
    # app.createNewCatalogFile(computeWeight=False)
    # app.createNewCatalogFile(computeWeight=True)
    app.makeReportFile()
    app.visualizeResult()
