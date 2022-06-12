from LatLon import lat_lon as ll
from datetime import datetime as dt
from datetime import timedelta as td
from statistics import mean, stdev
from string import ascii_uppercase
from json import load
from pathlib import Path
from shutil import copy, move
from glob import glob
from pytz import utc
import pykonal
from tqdm.contrib import tzip
from math import ceil
from util.plot import plotSeismicityMap, plotHypocenterDiff, plotVelocityModels
from util.synthTime import extractTT, generateTTT
from obspy.geodetics.base import gps2dist_azimuth as gps
from obspy import UTCDateTime as utc
from util import feedEvents2Catalog
from numpy import array
from util.nordic2xyzm import nordic2xyzm
from util.logger import myLogger
import random
import os
import sys

# Define main class


class Main():
    # Init
    def __init__(self):
        # self.resultsPath = os.path.join("results", dt.now().strftime("%Y_%j_%H%M%S"))
        self.resultsPath = os.path.join("results", "test")
        Path(self.resultsPath).mkdir(parents=True, exist_ok=True)
        self.report = myLogger(os.path.join(self.resultsPath, "report"), mode="w")
        self.readConfiguration()
        self.config2Log()
        self.checkRequiredFiles()
        self.clearExistingFiles()        
        self.node_intervals = [self.config["TravelTimeGenerator"]["velocityXGridIntervals"],
                               self.config["TravelTimeGenerator"]["velocityYGridIntervals"],
                               self.config["TravelTimeGenerator"]["velocityZGridIntervals"]]   
        self.travelTimeDict = {}     
    
    # Read configuration parameters
    def readConfiguration(self):
        # - check if configuration file exists
        if not os.path.exists("config.json"):
            msg = "+++ Could not find configuration file! Aborting ..."
            print(msg); self.report.info(msg)
            sys.exit()
        with open("config.json") as f:
            self.config = load(f)
        msg = "+++ Configuration file was loaded successfully."
        print(msg); self.report.info(msg)
    
    # Write configuration to log file
    def config2Log(self):
        for key in self.config:
            if isinstance(self.config[key], dict):
                for k,v in self.config[key].items():
                    msg = "+++ Configuration on {key} > {k}: {v}".format(key=key, k=k, v=v)
                    self.report.info(msg)
            else:
                msg = "+++ {v} :".format(v=self.config[key])
                self.report.info(msg)

    # Check required input files
    def checkRequiredFiles(self):
        # - check earthquake arrivals and station files
        if self.config["FSS"]["flag"] and not self.config["RSS"]["flag"]:
            self.createNewStationFile()
        elif not self.config["FSS"]["flag"] and self.config["RSS"]["flag"]:
            self.eqFile = os.path.join("EqInput", self.config["RSS"]["Inputs"]["EqFile"])
            self.stationFile = os.path.join(
                "EqInput", self.config["RSS"]["Inputs"]["StationFile"])
            if not os.path.exists(self.eqFile) or not os.path.exists(self.stationFile):
                msg = "+++ Could not find earthquake or station files! Aborting ..."
                print(msg); self.report.info(msg)
                sys.exit()        # - Create NLloc required directories
        Path("time").mkdir(parents=True, exist_ok=True)
        Path("model").mkdir(parents=True, exist_ok=True)
        Path("obs").mkdir(parents=True, exist_ok=True)
        Path(os.path.join(self.resultsPath, "relocation")).mkdir(parents=True, exist_ok=True)
        Path("tmp").mkdir(parents=True, exist_ok=True)
        
    # Clear pre-existing files
    def clearExistingFiles(self):
        for d in ["model", "obs", "time", "tmp"]:
            for f in glob(os.path.join(d, "*")):
                os.remove(f)
        msg = "+++ Removing preexisting files in 'model', 'obs', 'time' and 'tmp' directories ..."
        print(msg); self.report.info(msg)

    def generateTTTable(self):
        """generate travel time tables and store a bank of files
        """
        # reading velocity file
        velocityModelDict = self.readVelocityFile(os.path.join("EqInput", "STATION0.HYP"))
        if self.config["TravelTimeGenerator"]["flag"]:
            xMaxDist = self.config["TravelTimeGenerator"]["numXGrid"] * \
                self.config["TravelTimeGenerator"]["velocityXGridIntervals"]
            zMaxDist = self.config["TravelTimeGenerator"]["numZGrid"] * \
                self.config["TravelTimeGenerator"]["velocityZGridIntervals"]
            print("\n+++ Generating travel time tables for X range (0, {xMax:4.0f})km, and Z range (0, {zMax:4.0f})km\n".format(
                xMax=xMaxDist, zMax=zMaxDist))
            print("+++ Generating travel time tables for P phase ")
            generateTTT(
                velocityModelDict, "P",
                self.config["TravelTimeGenerator"]["numXGrid"], self.config["TravelTimeGenerator"]["numZGrid"],
                self.node_intervals, self.config["TravelTimeGenerator"]["decimationFactor"]
            )
            print("+++ Generating travel time tables for S phase ")
            generateTTT(
                velocityModelDict, "S",
                self.config["TravelTimeGenerator"]["numXGrid"], self.config["TravelTimeGenerator"]["numZGrid"],
                self.node_intervals, self.config["TravelTimeGenerator"]["decimationFactor"]
            )

    # Create new station file
    def createNewStationFile(self):
        """_summary_
        """
        stationCodes, stationLats, stationLons = [], [], []
        nStations = self.config["FSS"]["Station"]["NumberOfStations"]
        stationCodeLength = 4
        random.seed(self.config["FSS"]["Station"]["SeedRandomID"])
        stationCodes = ["".join(random.sample(ascii_uppercase, k=stationCodeLength)) for _ in range(nStations)]
        StationBoundLonMin = self.config["FSS"]["Station"]["StationBoundLonMin"]
        StationBoundLonMax = self.config["FSS"]["Station"]["StationBoundLonMax"]
        StationBoundLatMin = self.config["FSS"]["Station"]["StationBoundLatMin"]
        StationBoundLatMax = self.config["FSS"]["Station"]["StationBoundLatMax"]            
        Vp = self.config["FSS"]["Model"]["Vp"]
        Z = self.config["FSS"]["Model"]["Z"]
        Interfaces = self.config["FSS"]["Model"]["Interfaces"]
        VpVs = self.config["FSS"]["Model"]["VpVs"]
        xNear = self.config["FSS"]["Model"]["xNear"]
        xFar = self.config["FSS"]["Model"]["xFar"]
        trialDepth = self.config["FSS"]["Model"]["trialDepth"]
        if self.config["FSS"]["Station"]["GaussianDistribution"]:
            meanLon, stdLon = mean([StationBoundLonMin, StationBoundLonMax]), stdev([StationBoundLonMin, StationBoundLonMax])
            meanLat, stdLat = mean([StationBoundLatMin, StationBoundLatMax]), stdev([StationBoundLatMin, StationBoundLatMax])
            random.seed(self.config["FSS"]["Station"]["SeedRandomID"])
            stationLons = [random.gauss(meanLon, stdLon) for _ in range(nStations)]
            random.seed(self.config["FSS"]["Station"]["SeedRandomID"])
            stationLats = [random.gauss(meanLat, stdLat) for _ in range(nStations)]
        elif self.config["FSS"]["Station"]["UniformDistribution"]:
            random.seed(self.config["FSS"]["Station"]["SeedRandomID"])
            stationLons = [random.uniform(StationBoundLonMin, StationBoundLonMax) for _ in range(nStations)]
            random.seed(self.config["FSS"]["Station"]["SeedRandomID"])
            stationLats = [random.uniform(StationBoundLatMin, StationBoundLatMax) for _ in range(nStations)]
        elif self.config["FSS"]["Station"]["FaultDistribution"]:
            pass
        with open(os.path.join("files", "resets.dat")) as f, open(os.path.join("EqInput", "STATION0.HYP"), "w") as g:
            for l in f:
                g.write(l)
            g.write("\n"); g.write("\n")
            for code, lat, lon in zip(stationCodes, stationLats, stationLons):
                lat = ll.Latitude(lat)
                lon = ll.Longitude(lon)
                elv = 0.0
                g.write("  {code:4s}{latDeg:2.0f}{latMin:5.2f}N {lonDeg:2.0f}{lonMin:5.2f}E{elv:4.0f}\n".format(
                    code=code, latDeg=lat.decimal_degree, latMin=lat.decimal_minute, lonDeg=lon.decimal_degree, lonMin=lon.decimal_minute, elv=elv
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
    
    def createNewCatalogFile(self):
        """_summary_
        """
        stationsDict = self.readStationFile(os.path.join("EqInput", "STATION0.HYP"))
        eventLons, eventLats, eventDeps = [], [], []
        eventsArrivals = []
        eventsInfo = []
        NumberOfEvents = self.config["FSS"]["Catalog"]["NumberOfEvents"]
        CatalogBoundLonMin = self.config["FSS"]["Catalog"]["CatalogBoundLonMin"]
        CatalogBoundLonMax = self.config["FSS"]["Catalog"]["CatalogBoundLonMax"]
        CatalogBoundLatMin = self.config["FSS"]["Catalog"]["CatalogBoundLatMin"]
        CatalogBoundLatMax = self.config["FSS"]["Catalog"]["CatalogBoundLatMax"]
        CatalogBoundDepMin = self.config["FSS"]["Catalog"]["CatalogBoundDepMin"]
        CatalogBoundDepMax = self.config["FSS"]["Catalog"]["CatalogBoundDepMax"]
        if self.config["FSS"]["Catalog"]["GaussianDistribution"]:
            meanLon, stdLon = mean([CatalogBoundLonMin, CatalogBoundLonMax]), stdev([CatalogBoundLonMin, CatalogBoundLonMax])
            meanLat, stdLat = mean([CatalogBoundLatMin, CatalogBoundLatMax]), stdev([CatalogBoundLatMin, CatalogBoundLatMax])
            meanDep, stdDep = mean([CatalogBoundDepMin, CatalogBoundDepMax]), stdev([CatalogBoundDepMin, CatalogBoundDepMax])
            random.seed(self.config["FSS"]["Catalog"]["SeedRandomID"])
            eventLons = [random.gauss(meanLon, stdLon) for _ in range(NumberOfEvents)]
            random.seed(self.config["FSS"]["Catalog"]["SeedRandomID"])
            eventLats = [random.gauss(meanLat, stdLat) for _ in range(NumberOfEvents)]
            random.seed(self.config["FSS"]["Catalog"]["SeedRandomID"])
            eventDeps = [random.gauss(meanDep, stdDep) for _ in range(NumberOfEvents)]
        elif self.config["FSS"]["Catalog"]["UniformDistribution"]:
            random.seed(self.config["FSS"]["Catalog"]["SeedRandomID"])
            eventLons = [random.uniform(CatalogBoundLonMin, CatalogBoundLonMax) for _ in range(NumberOfEvents)]
            random.seed(self.config["FSS"]["Catalog"]["SeedRandomID"])
            eventLats = [random.uniform(CatalogBoundLatMin, CatalogBoundLatMax) for _ in range(NumberOfEvents)]
            random.seed(self.config["FSS"]["Catalog"]["SeedRandomID"])
            eventDeps = [random.uniform(CatalogBoundDepMin, CatalogBoundDepMax) for _ in range(NumberOfEvents)]
        elif self.config["FSS"]["Catalog"]["FaultDistribution"]:
            pass
        print("+++ Generating new catalog")
        eventOt = utc.now()  # type: ignore
        eventsTimeShift = 60.0
        for (eventLat, eventLon, eventDep) in tzip(eventLats, eventLons, eventDeps):  # type: ignore
            newEventArrivals = {}
            eventOt += eventsTimeShift
            newEventInfo = {"OriginTime":eventOt, "Longitude":eventLon, "Latitude":eventLat, "Depth":eventDep}
            eventsInfo.append(newEventInfo)
            selectedStations, receivers = self.computeNewDistances(eventLat, eventLon, stationsDict)
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
            for station, ttP, ttS, distance in zip(selectedStations, extractedTTP, extractedTTS, receivers):
                travelTimeP, weightP = self.addTTError(eventOt+ttP, m=0.0, std=self.config["FPS"]["ErrorOnArrivals"]["P"])
                travelTimeS, weightS = self.addTTError(eventOt+ttP, m=0.0, std=self.config["FPS"]["ErrorOnArrivals"]["S"])
                newEventArrivals[station] = {
                    "P":["P", "emergent", weightP, travelTimeP, distance[0]],
                    "S":["S", "emergent", weightS, travelTimeS, distance[0]]
                    }
            eventsArrivals.append(newEventArrivals)
        feedCatalog = feedEvents2Catalog.feedCatalog()
        newCatalog = feedCatalog.setCatalog(eventsInfo, eventsArrivals)
        newCatalog.write(os.path.join("EqInput", "select.out"), format="NORDIC")

    def addTTError(self, tt, m=0.0, std=1.0):
        """_summary_

        Args:
            tt (_type_): _description_
            m (float, optional): _description_. Defaults to 0.0.
            std (float, optional): _description_. Defaults to 1.0.

        Returns:
            _type_: _description_
        """
        err = random.gauss(m, std)
        tt = tt + err
        w = self.mapError2Weight(err, m, std)
        return tt, w
    
    def mapError2Weight(self, err, m, std):
        """_summary_

        Args:
            err (_type_): _description_
            m (_type_): _description_
            std (_type_): _description_

        Returns:
            _type_: _description_
        """
        if not self.config["FPS"]["ErrorOnArrivals"]["ActivateWeightingBasedOnError"]:
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


    def computeNewDistances(self, eventLat, eventLon, stations):
        """given event coordinates it computes new station's distance

        Args:
            eventLat (float): event latitude in degree
            eventLon (float): event longitude in degree
            picks (obspy.picks): obspy picks for an event
            stations (dict): a dictionary contains stations information

        Returns:
            numpy.array: an array contains stations distance
        """
        randomStations = random.sample(list(stations.keys()), k=random.randint(3, len(stations)))
        stationLats = [stations[station]["Lat"] for station in randomStations]
        stationLons = [stations[station]["Lon"] for station in randomStations]
        distances = [gps(eventLat, eventLon, stationLat, stationLon)[
            0]*1e-3 for (stationLat, stationLon) in zip(stationLats, stationLons)]
        return randomStations, [array([distance, 1.0, 0.0]) for distance in distances]
                
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

    # # Write NLLOC velocity model
    # def writeNLlocVelocityModel(self):
    #     msg = "+++ Writing velocity model in NLLOC format ..."
    #     self.report.info(msg)
    #     velocityModelDict = self.readVelocityFile(os.path.join("EqInput", "STATION0.HYP"))
    #     VPs, Zs, VpVs = velocityModelDict["Vp"], velocityModelDict["Z"], velocityModelDict["VpVs"]
    #     if self.config["FPS"]["ErrorOnVelocityModel"]["Flag"]:
    #         velocityModels = {"1": {"V": None, "Z": None},
    #                           "2": {"V": None, "Z": None}}
    #         velocityModels["1"]["V"] = VPs
    #         velocityModels["1"]["Z"] = Zs
    #         VelocityerrorPercentage = self.config["FPS"]["ErrorOnVelocityModel"]["VelocityerrorPercentage"]
    #         ThiknesserrorPercentage = self.config["FPS"]["ErrorOnVelocityModel"]["ThiknesserrorPercentage"]
    #         VPs = self.AddError(VPs, VelocityerrorPercentage)
    #         Zs = self.AddError(Zs, ThiknesserrorPercentage)
    #         velocityModels["2"]["V"] = VPs  # type: ignore
    #         velocityModels["2"]["Z"] = Zs  # type: ignore
    #         maxDepthToPlot = self.config["FPS"]["ErrorOnVelocityModel"]["MaxDepthToPlot"]
    #         plotVelocityModels(velocityModels, maxDepthToPlot, self.resultsPath, self.config)
    #     with open(os.path.join("EqInput", "model.dat"), "w") as f:
    #         f.write("#LAYER depth, Vp_top, Vp_grad, Vs_top, Vs_grad, p_top, p_grad\n")
    #         for Vp, Z in zip(VPs, Zs):
    #             f.write("LAYER {Z:6.2f} {Vp:5.2f} 0.00 {Vs:5.2f} 0.00 2.70 0.00\n".format(
    #                 Z=Z, Vp=Vp, Vs=Vp/VpVs
    #             ))

    # # Add random error to list items
    # def AddError(self, x, errorPercentage):
    #     seed = self.config["FPS"]["ErrorOnVelocityModel"]["Seed"]
    #     if seed:
    #         random.seed(seed)
    #     x = [v + random.gauss(0, v*errorPercentage*1e-2) for v in x]
    #     return x

    # # Parse station information
    # def parseStationInfo(self):
    #     msg = "+++ Parsing station information ..."
    #     self.report.info(msg)
    #     emptyLines = 0
    #     self.stations = {"Code": [], "Lat": [], "Lon": [], "Elv": []}
    #     with open(self.stationFile) as f:
    #         for l in f:
    #             if not l.strip():
    #                 emptyLines += 1
    #             if emptyLines == 1 and l.strip():
    #                 code = l[:6].strip()
    #                 lat = ll.Latitude(degree=int(
    #                     l[6:8]), minute=float(l[8:13])).decimal_degree
    #                 lon = ll.Longitude(degree=int(
    #                     l[15:17]), minute=float(l[17:22])).decimal_degree
    #                 elv = float(l[23:27])
    #                 for key, value in zip(["Code", "Lat", "Lon", "Elv"], [code, lat, lon, elv]):
    #                     self.stations[key].append(value)

    # # Write NLLOC station file
    # def writeNLlocStationFile(self):
    #     msg = "+++ Writing station information in NLLOC format ..."
    #     codes, lats, lons, elvs = self.stations.items()
    #     with open(os.path.join("EqInput", "stations.dat"), "w") as f:
    #         for code, lat, lon, elv in zip(codes[1], lats[1], lons[1], elvs[1]):
    #             if elv >= 0:
    #                 f.write("GTSRCE {code:4s} LATLON {lat:7.3f} {lon:7.3f}  0.000 {elv:6.3f}\n".format(
    #                     code=code, lat=lat, lon=lon, elv=elv*1e-3
    #                 ))
    #             else:
    #                 f.write("GTSRCE {code:4s} LATLON {lat:7.3f} {lon:7.3f} {elv:6.3f}  0.000\n".format(
    #                     code=code, lat=lat, lon=lon, elv=-elv*1e-3
    #                 ))

    # # Parse earthquake info line
    # def parseEventLine(self, l):
    #     for i in [6, 8, 11, 12, 13, 14, 16, 17, 19]:
    #         if l[i] == " ":
    #             l = l[:i]+"0"+l[i+1:]
    #     try:
    #         ot = dt.strptime(l[:20], " %Y %m%d %H%M %S.%f")
    #         lat = float(l[24:30])
    #         lon = float(l[32:38])
    #         dep = float(l[39:43])
    #         return ot, lat, lon, dep
    #     except ValueError:
    #         return None, None, None, None

    # # Parse earthquake file
    # def parseEarthquakeFile(self):
    #     msg = "+++ Parsing earthquake input file ..."
    #     self.report.info(msg)
    #     self.earthquakesInfo = {
    #         "OT": [],
    #         "LAT": [],
    #         "LON": [],
    #         "DEP": [],
    #         "PArrivals": [],
    #         "SArrivals": []}
    #     with open(self.eqFile) as f:
    #         for l in f:
    #             if l.strip() and l[79] == "1":
    #                 ot, lat, lon, dep = self.parseEventLine(l)
    #                 if ot == None:
    #                     continue
    #                 self.earthquakesInfo["OT"].append(ot)
    #                 self.earthquakesInfo["LAT"].append(lat)
    #                 self.earthquakesInfo["LON"].append(lon)
    #                 self.earthquakesInfo["DEP"].append(dep)
    #                 self.earthquakesInfo["PArrivals"].append([])
    #                 self.earthquakesInfo["SArrivals"].append([])
    #             if l.strip() and l[10] == "P" and l[79] in [" ", "4"]:
    #                 staCode = l[:5].strip()
    #                 self.earthquakesInfo["PArrivals"][-1].append(staCode)
    #             if l.strip() and l[10] == "S" and l[79] in [" ", "4"]:
    #                 staCode = l[:5].strip()
    #                 self.earthquakesInfo["SArrivals"][-1].append(staCode)

    # # Write NLLOC control file and generate synthetics
    # def writeNLlocControlFile(self):
    #     msg = "+++ Writing NLLOC control file ..."
    #     self.report.info(msg)
    #     lonMin = self.config["Region"]["LonMin"]
    #     lonMax = self.config["Region"]["LonMax"]
    #     latMin = self.config["Region"]["LatMin"]
    #     latMax = self.config["Region"]["LatMax"]
    #     yNum = self.config["FPS"]["VGGRID"]["yNum"]
    #     zNum = self.config["FPS"]["VGGRID"]["zNum"]
    #     lonOrig = mean([lonMin, lonMax])
    #     latOrig = mean([latMin, latMax])
    #     OTs, LATs, LONs, Deps, PArrivalsList, SArrivalsList = self.earthquakesInfo.items()
    #     mechType = self.config["FSS"]["Model"]["Fault"]["Type"]
    #     strike = self.config["FSS"]["Model"]["Fault"]["Strike"]
    #     dip = self.config["FSS"]["Model"]["Fault"]["Dip"]
    #     rake = self.config["FSS"]["Model"]["Fault"]["Rake"]
    #     velocityModelDict = self.readVelocityFile(os.path.join("EqInput", "STATION0.HYP"))
    #     with open("nlloc.conf", "w") as f:
    #         # - GENERIC CONTROL STATEMENT
    #         f.write("# +++ GENERIC CONTROL STATEMENT\n")
    #         f.write("CONTROL -1 54321\n")
    #         f.write("TRANS LAMBERT WGS-84 {latOrig:.2f} {lonOrig:.2f} {latMin:.2f} {latMax:.2f} 0.00\n".format(
    #             latOrig=latOrig, lonOrig=lonOrig, latMin=latMin, latMax=latMax
    #         ))
    #         f.write("\n")
    #         # - VELOCITY GRID STATEMENT
    #         f.write("# +++ VELOCITY GRID STATEMENT\n")
    #         f.write("VGOUT ./model/layer\n")
    #         f.write("VGTYPE P\nVGTYPE S\n")
    #         f.write("VGGRID 2 {yNum:d} {zNum:d} 0.0 0.0 -3.0 1.0 1.0 1.0 SLOW_LEN\n".format(
    #             yNum=yNum, zNum=zNum
    #         ))
    #         f.write("INCLUDE ./EqInput/model.dat\n")
    #         f.write("\n")
    #         # - GRID2TIME STATEMENTS
    #         f.write("# +++ GRID2TIME STATEMENTS\n")
    #         f.write("GTFILES ./model/layer ./time/layer P\n")
    #         f.write("GTMODE GRID2D ANGLES_YES\n")
    #         f.write("INCLUDE ./EqInput/stations.dat\n")
    #         f.write("GT_PLFD  1.0e-9  0\n")
    #         f.write("\n")
    #         # - Quality to Error Mapping
    #         f.write("# +++ Quality to Error Mapping\n")
    #         f.write("EQQUAL2ERR 0.1 0.2 0.4 0.8 99999.9\n")
    #         # - VpVs
    #         f.write("# +++ P Velocity to S Velocity Ratio\n")
    #         f.write("EQVPVS {VpVs:5.2f}\n".format(
    #             VpVs=velocityModelDict["VpVs"]))
    #         f.write("\n")
    #     # Run NLLOC for initiation
    #     self.runNLlocForward()
    #     # Now generating synthetics
    #     msg = "+++ Generating synthetic earthquakes ..."
    #     print(msg); self.report.info(msg)
    #     PArrivalsList, SArrivalsList = self.editUsedPhase(list(PArrivalsList), list(SArrivalsList))
    #     for i, (ot, lat, lon, dep, pArrivals, sarrivals) in enumerate(tzip(OTs[1], LATs[1], LONs[1], Deps[1], PArrivalsList[1], SArrivalsList[1])):  # type: ignore
    #         outConfig = os.path.join("tmp", "nlloc_{i:d}.conf".format(i=i+1))
    #         copy("nlloc.conf", outConfig)
    #         # - Shift each earthquake if requested by user
    #         lat += self.config["FPS"]["LocationShift"]["ShiftInLatitude"]
    #         lon += self.config["FPS"]["LocationShift"]["ShiftInLongitude"]
    #         with open(outConfig, "a") as f:
    #             # - TIME2EQ STATEMENTS
    #             f.write("# +++ TIME2EQ STATEMENTS\n")
    #             f.write(
    #                 "EQFILES time/layer obs/SYNEQ_{i:d}.obs\n".format(i=i+1))
    #             f.write("EQMECH {mechType:s} {strike:5.1f} {dip:5.1f} {rake:5.1f}\n".format(
    #                 mechType=mechType, strike=strike, dip=dip, rake=rake
    #             ))
    #             f.write("EQMODE SRCE_TO_STA\n")
    #             f.write("\n")
    #             f.write("EQSRCE SYNEQ_{i:d} LATLON {eventLat:7.3f} {eventLon:7.3f} {eventDep:5.1f} 0.0\n".format(
    #                 i=i, eventLat=lat, eventLon=lon, eventDep=dep
    #             ))
    #             for pArrival in pArrivals:
    #                 errorP = self.config["FPS"]["ErrorOnArrivals"]["P"]
    #                 f.write("EQSTA {staCode:4s} P GAU {errorP:6.2f} GAU {errorP:6.2f}\n".format(
    #                     staCode=pArrival, errorP=errorP))
    #             for sarrival in sarrivals:
    #                 errorS = self.config["FPS"]["ErrorOnArrivals"]["S"]
    #                 f.write("EQSTA {staCode:4s} S GAU {errorS:6.2f} GAU {errorS:6.2f}\n".format(
    #                     staCode=sarrival, errorS=errorS))
    #             f.write("\n")
    #         # - Synthetic arrival times
    #         cmd = "Time2EQ {outConfig:s} > /dev/null".format(
    #             outConfig=outConfig)
    #         os.system(cmd)
    #         self.correctOT(ot, os.path.join(
    #             "obs", "SYNEQ_{i:d}.obs".format(i=i+1)))
    #     if os.path.exists("nlloc.conf"): os.remove("nlloc.conf")

    # # Setting up how many phases (P, S) will be used 
    # def editUsedPhase(self, PArrivalsList, SArrivalsList):
    #     PArrivals = PArrivalsList[1]
    #     SArrivals = SArrivalsList[1]
    #     percentageP = self.config["FPS"]["UsageOfPhase"]["PercentageUseOfP"]
    #     percentageS = self.config["FPS"]["UsageOfPhase"]["PercentageUseOfS"]
    #     PArrivals = [random.sample(PArrival, ceil(len(PArrival)*percentageP*1e-2)) for PArrival in PArrivals]
    #     SArrivals = [random.sample(SArrival, ceil(len(SArrival)*percentageS*1e-2)) for SArrival in SArrivals]
    #     PArrivalsList[1] = PArrivals
    #     SArrivalsList[1] = SArrivals
    #     return PArrivalsList, SArrivalsList

    # # Correct origin time
    # def correctOT(self, ot, syntheticEqFile):
    #     with open(syntheticEqFile) as f, open("tmpFile", "w") as g:
    #         for i, l in enumerate(f):
    #             if i == 0:
    #                 line = l.split()
    #                 line.append("{0:f}".format(ot.timestamp()))
    #                 line.pop(-2)
    #                 l = " ".join(line)
    #                 g.write(l+"\n")
    #             elif l.strip() and "#" not in l and "PUBLIC_ID" not in l:
    #                 line = l.split()
    #                 arrivalTime = " ".join(line[6:9])
    #                 arrivalTime = dt.strptime(arrivalTime, "%Y%m%d %H%M %S.%f")
    #                 arrivalTime = arrivalTime - dt(1900, 1, 1, 0, 0, 0)
    #                 arrivalTime = ot + td(seconds=arrivalTime.total_seconds())
    #                 date = arrivalTime.strftime("%Y%m%m")
    #                 hourMin = arrivalTime.strftime("%H%M")
    #                 sec = arrivalTime.strftime("%S.%f")
    #                 g.write(
    #                     "{stationName:9s} {instrument:4s} {component:4s} {pPhaseOnset:4s} {PhaseDes:6s} {firstMotion:1s} {date:8s} {hourMin:6s} {sec:7s} {err:3s} {errMag:9s} {codaDur:9s} {amplitude:9s} {period:9s} {priorWt:9s}\n".format(
    #                         stationName=line[0], instrument=line[1], component=line[2], pPhaseOnset=line[3],
    #                         PhaseDes=line[4], firstMotion=line[5], date=date, hourMin=hourMin, sec=sec, err=line[9],
    #                         errMag=line[10], codaDur=line[11], amplitude=line[12], period=line[13], priorWt=line[14]
    #                     ))
    #     move("tmpFile", syntheticEqFile)

    # # Run NLLOC for generating velocity grid and travel time tables
    # def runNLlocForward(self):
    #     # - Make velocity grid files for P and S
    #     cmd = "Vel2Grid nlloc.conf"
    #     os.system(cmd)
    #     msg = "+++ P and S velocity grids have been created successfull ..."
    #     print(msg); self.report.info(msg)
    #     # - Generate P travel time grid files for stations
    #     cmd = "Grid2Time nlloc.conf"
    #     os.system(cmd)
    #     # - Now for S
    #     cmd = "sed -i 's/layer P/layer S/g' nlloc.conf"
    #     os.system(cmd)
    #     cmd = "Grid2Time nlloc.conf"
    #     os.system(cmd)
    #     msg = "+++ P and S travel time grids have been created successfully ..."
    #     print(msg); self.report.info(msg)

    # # Merge synthetic files
    # def mergeSyntheticFiles(self):
    #     msg = "+++ Merging synthetic files ..."
    #     self.report.info(msg)
    #     syntheticFiles = sorted(glob(os.path.join(
    #         "obs", "SYNEQ_*.obs")), key=lambda x: int(x.split("_")[1].split(".")[0]))
    #     for syntheticFile in syntheticFiles:
    #         cmd = "cat {syntheticFile:s} >> obs/SYNEQ.obs".format(
    #             syntheticFile=syntheticFile)
    #         os.system(cmd)
    #     for _ in glob(os.path.join("obs", "SYNEQ_*")):
    #         os.remove(_)

    # # Convert NLLOC to NORDIC
    # def convertNLLOC2NORDIC(self):
    #     msg = "+++ Converting NLLOC to NORDIC format ..."
    #     self.report.info(msg)
    #     nllocObsFile = os.path.join("obs/SYNEQ.obs")
    #     if not nllocObsFile:
    #         msg = "+++ Could not find synthetic earthquake file! Aborting ..."
    #         print(msg); self.report.info(msg)
    #         sys.exit(0)
    #     hypo_file = open(os.path.join(self.resultsPath, "relocation", "synthetic.out"), "w")
    #     with open(nllocObsFile) as f:
    #         line1_flag = False
    #         line7_flag = False
    #         start_flag = False
    #         for l in f:
    #             if "# EQEVENT:" in l:
    #                 line1_flag = True
    #                 line7_flag = True
    #                 if start_flag:
    #                     hypo_file.write("\n")
    #                     start_flag = False
    #             elif "#" not in l:
    #                 start_flag = True
    #                 if line1_flag and line7_flag:
    #                     y = int(l[34:38])
    #                     m = int(l[38:40])
    #                     d = int(l[40:42])
    #                     hh = int(l[43:45])
    #                     mm = int(l[45:47])
    #                     ss = 0
    #                     ms = 0
    #                     ot = dt(y, m, d, hh, mm, ss, ms)
    #                     line_1 = " {year:4d} {month:02d}{day:02d} {hour:02d}{minute:02d} {second:4.1f} L                                                         1\n".format(
    #                         year=ot.year, month=ot.month, day=ot.day, hour=ot.hour, minute=ot.minute, second=ot.second
    #                     )
    #                     line_7 = " STAT SP IPHASW D HRMN SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W  DIS CAZ7\n"
    #                     hypo_file.write(line_1)
    #                     hypo_file.write(line_7)
    #                     line1_flag = False
    #                     line7_flag = False
    #                 line = l.strip().split()
    #                 stationCode = line[0]
    #                 instrument = " "
    #                 comp = " "
    #                 qualityIndicator = " "
    #                 phaseID = line[4]+"   "
    #                 weightingIndicator = 0
    #                 if line[4] == "P":
    #                     firstMotion = line[5]
    #                 else:
    #                     firstMotion = " "
    #                 arrivalTime = dt(int(line[6][0:4]), int(line[6][4:6]), int(line[6][6:8]),
    #                                  int(line[7][0:2]), int(line[7][2:4]), int(
    #                     line[8].split(".")[0]),
    #                     int(line[8].split(".")[1]))
    #                 line_4 = " {stationCode:5s}{instrument:1s}{comp:1s} {qualityIndicator:1s}{phaseID:4s}{weightingIndicator:1d} {firstMotion:1s} {hour:02d}{minute:02d} {second:5.2f}                                                   4\n".format(
    #                     stationCode=stationCode, instrument=instrument, comp=comp, qualityIndicator=qualityIndicator, phaseID=phaseID, weightingIndicator=weightingIndicator, firstMotion=firstMotion, hour=arrivalTime.hour, minute=arrivalTime.minute, second=arrivalTime.second+arrivalTime.microsecond*1e-6
    #                 )
    #                 hypo_file.write(line_4)
    #     hypo_file.write("\n")
    #     hypo_file.close()
    #     msg = "+++ Output file 'nlloc.hyp' has been created.\n"
    #     print(msg)
    #     # - Copy nordic station file into 'relocation' directory
    #     copy(self.stationFile, os.path.join(self.resultsPath, "relocation", "STATION0.HYP"))

    # # Relocate data using "hypocenter" program
    # def relocateWithHyp(self, targetDir, nordicFileName):
    #     root = os.getcwd()
    #     os.chdir(targetDir)
    #     with open("hyp.inp", "w") as f:
    #         f.write("{nordicFileName:s}\n".format(
    #             nordicFileName=nordicFileName))
    #         f.write("n\n")
    #     cmd = "hyp < hyp.inp > /dev/null"
    #     os.system(cmd)
    #     os.chdir(root)

    # # Make summary on initial and final locations
    # def makeSummaryFile(self, targetDir, nordicFileName):
    #     root = os.getcwd()
    #     os.chdir(targetDir)
    #     nordic2xyzm(nordicFileName)
    #     os.chdir(root)

    # # Compare raw and relocated data
    # def compare(self):
    #     msg = "+++ Making some statistical figures ..."
    #     self.report.info(msg)
    #     relocationPath = os.path.join(self.resultsPath, "relocation")
    #     self.makeSummaryFile(targetDir="EqInput", nordicFileName="select.out")
    #     self.relocateWithHyp(targetDir=relocationPath,
    #                          nordicFileName="synthetic.out")
    #     self.makeSummaryFile(targetDir=relocationPath, nordicFileName="hyp.out")
    #     iniFile = os.path.join("EqInput", "xyzm.dat")
    #     finFile = os.path.join(relocationPath, "xyzm.dat")
    #     plotSeismicityMap(iniFile, finFile, self.stations, self.resultsPath)
    #     plotHypocenterDiff(iniFile, finFile, self.resultsPath, self.config)


# Run application
if __name__ == "__main__":
    app = Main()
    app.generateTTTable()
    app.createNewStationFile()
    app.createNewCatalogFile()
    # app.parseStationInfo()
    # app.parseVelocityModel()
    # app.writeNLlocVelocityModel()
    # app.writeNLlocStationFile()
    # app.parseEarthquakeFile()
    # app.writeNLlocControlFile()
    # app.mergeSyntheticFiles()
    # app.convertNLLOC2NORDIC()
    # app.compare()
