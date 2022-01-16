from LatLon import lat_lon as ll
from json import load
import os, sys

class Main():

    # Init
    def __init__(self):
        pass

    # Read configuration parameters
    def readConfiguration(self):
        # - check if configuration file exists
        if not os.path.exists("config.json"):
            print("+++ Could not find configuration file! Aborting ...")
            sys.exit()
        with open("config.json") as f:
            self.config = load(f)
        print("+++ Configuration file was loaded successfully.")

    # Check required input files
    def checkRequiredFiles(self):
        # - check earthquake arrivals and station files
        self.eqFile = os.path.join("EqInput", self.config["Inputs"]["EqFile"])
        self.stationFile = os.path.join("EqInput", self.config["Inputs"]["StationFile"])
        if not os.path.exists(self.eqFile) or not os.path.exists(self.stationFile):
            print("+++ Could not find earthquake or station files! Aborting ...")
            sys.exit()

    # Parse velocity model
    def parseVelocityModel(self):
        emptyLines = 0
        self.velocityModel = {"Vp":[], "Z":[]}
        with open(self.stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 2 and l.strip():
                    Vp, Z = [float(x) for x in l.split()[:2]]
                    self.velocityModel["Vp"].append(Vp)
                    self.velocityModel["Z"].append(Z)
                
    # Parse station information
    def parseStationInfo(self):
        emptyLines = 0
        self.stations = {"Code":[], "Lat":[], "Lon":[], "Elv":[]}
        with open(self.stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 2 and l.strip():
                    code = l[:6].strip()
                    lat = ll.Latitude()


    # Write NLLOC control file
    def writeNLlocControlFile(self):
        pass
    # Convert NLLOC to NORDIC
    def convertNLLOC2NORDIC(self):
        pass

if __name__ == "__main__":
    app = Main()
    app.readConfiguration()
    app.checkRequiredFiles()
    app.parseVelocityModel()