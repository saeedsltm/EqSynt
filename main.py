from LatLon import lat_lon as ll
from datetime import datetime as dt
from statistics import mean
from json import load
from pathlib import Path
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
        # - Create NLloc required directories
        Path("time").mkdir(parents=True, exist_ok=True)
        Path("model").mkdir(parents=True, exist_ok=True)
        Path("obs").mkdir(parents=True, exist_ok=True)

    # Parse velocity model
    def parseVelocityModel(self):
        emptyLines = 0
        self.velocityModel = {"Vp":[], "Z":[], "VpVs":1.73}
        with open(self.stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 2 and l.strip():
                    Vp, Z = [float(x) for x in l.split()[:2]]
                    self.velocityModel["Vp"].append(Vp)
                    self.velocityModel["Z"].append(Z)
                if emptyLines == 3 and l.strip():
                    VpVs = float(l[16:20])
                    self.velocityModel["VpVs"] = VpVs
                    break
    
    # Write NLLOC velocity model
    def writeNLlocVelocityModel(self):
        VPs, Zs, VpVs = self.velocityModel.items()
        with open(os.path.join("EqInput", "model.dat"), "w") as f:
            f.write("#LAYER depth, Vp_top, Vp_grad, Vs_top, Vs_grad, p_top, p_grad\n")
            for Vp, Z in zip(VPs[1], Zs[1]):
                f.write("LAYER {Z:6.2f} {Vp:5.2f} 0.00 {Vs:5.2f} 0.00 2.70 0.00\n".format(
                    Z=Z, Vp=Vp,Vs=Vp/VpVs[1]
                ))
                
    # Parse station information
    def parseStationInfo(self):
        emptyLines = 0
        self.stations = {"Code":[], "Lat":[], "Lon":[], "Elv":[]}
        with open(self.stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 1 and l.strip():
                    code = l[:6].strip()
                    lat = ll.Latitude(degree=int(l[6:8]), minute=float(l[8:13])).decimal_degree
                    lon = ll.Longitude(degree=int(l[15:17]), minute=float(l[17:22])).decimal_degree
                    elv = float(l[23:27])
                    for key,value in zip(["Code", "Lat", "Lon", "Elv"],[code, lat, lon, elv]):
                        self.stations[key].append(value)
    
    # Write NLLOC station file
    def writeNLlocStationFile(self):
        codes, lats, lons, elvs = self.stations.items()
        with open(os.path.join("EqInput", "stations.dat"), "w") as f:
            for code, lat, lon, elv in zip(codes[1], lats[1], lons[1], elvs[1]):
                if elv >= 0:
                    f.write("GTSRCE {code:4s} LATLON {lat:7.3f} {lon:7.3f}  0.000 {elv:6.3f}\n".format(
                        code=code, lat=lat, lon=lon, elv=elv*1e-3
                    ))
                else:
                    f.write("GTSRCE {code:4s} LATLON {lat:7.3f} {lon:7.3f} {elv:6.3f}  0.000\n".format(
                        code=code, lat=lat, lon=lon, elv=-elv*1e-3
                    ))

    # Parse earthquake info line
    def parseEventLine(self, l):
        for i in [6, 8, 11, 12, 13 ,14, 16, 17, 19]:
            if l[i] == " ":
                l = l[:i]+"0"+l[i+1:]
        ot = dt.strptime(l[:20], " %Y %m%d %H%M %S.%f")
        lat = float(l[24:30])
        lon = float(l[32:38])
        dep = float(l[39:43])
        return ot, lat, lon , dep

    # Parse earthquake file
    def parseEarthquakeFile(self):
        self.earthquakesInfo = {
            "OT":[],
            "LAT":[],
            "LON":[],
            "DEP":[],
            "PArrivals":[],
            "SArrivals":[]}
        with open(self.eqFile) as f:
            for l in f:
                if l.strip() and l[79] == "1":
                    ot, lat, lon, dep = self.parseEventLine(l)
                    self.earthquakesInfo["OT"].append(ot)
                    self.earthquakesInfo["LAT"].append(lat)
                    self.earthquakesInfo["LON"].append(lon)
                    self.earthquakesInfo["DEP"].append(dep)
                    self.earthquakesInfo["PArrivals"].append([])
                    self.earthquakesInfo["SArrivals"].append([])

                if l.strip() and l[10] == "P" and l[79] in [" ", "4"]:
                    staCode = l[:5].strip()
                    self.earthquakesInfo["PArrivals"][-1].append(staCode)
                if l.strip() and l[10] == "S"  and l[79] in [" ", "4"]:
                    staCode = l[:5].strip()
                    self.earthquakesInfo["SArrivals"][-1].append(staCode)

    # Write NLLOC control file
    def writeNLlocControlFile(self):
        lonMin = self.config["Region"]["LonMin"]
        lonMax = self.config["Region"]["LonMax"]
        latMin = self.config["Region"]["LatMin"]
        latMax = self.config["Region"]["LatMax"]
        yNum = self.config["VGGRID"]["yNum"]
        zNum = self.config["VGGRID"]["zNum"]        
        lonOrig=mean([lonMin, lonMax])
        latOrig=mean([latMin, latMax])
        OTs, LATs, LONs, Deps, PArrivalsList, SArrivalsList = self.earthquakesInfo.items()
        mechType = self.config["FirstMotionCalculation"]["Type"]
        strike = self.config["FirstMotionCalculation"]["Strike"]
        dip = self.config["FirstMotionCalculation"]["Dip"]
        rake = self.config["FirstMotionCalculation"]["Rake"]
        with open("nlloc.conf", "w") as f:
            # - GENERIC CONTROL STATEMENT
            f.write("# +++ GENERIC CONTROL STATEMENT\n")
            f.write("CONTROL -1 54321\n")
            f.write("TRANS LAMBERT WGS-84 {latOrig:.2f} {lonOrig:.2f} {latMin:.2f} {latMax:.2f} 0.00\n".format(
                latOrig=latOrig, lonOrig=lonOrig, latMin=latMin, latMax=latMax
            ))
            f.write("\n")
            # - VELOCITY GRID STATEMENT
            f.write("# +++ VELOCITY GRID STATEMENT\n")
            f.write("VGOUT ./model/layer\n")
            f.write("VGTYPE P\nVGTYPE S\n")
            f.write("VGGRID 2 {yNum:d} {zNum:d} 0.0 0.0 -3.0 1.0 1.0 1.0 SLOW_LEN\n".format(
                yNum=yNum, zNum=zNum
            ))
            f.write("INCLUDE ./EqInput/model.dat\n")
            f.write("\n")
            # - GRID2TIME STATEMENTS
            f.write("# +++ GRID2TIME STATEMENTS\n")
            f.write("GTFILES ./model/layer ./time/layer P\n")
            f.write("GTMODE GRID2D ANGLES_YES\n")
            f.write("INCLUDE ./EqInput/stations.dat\n")
            f.write("GT_PLFD  1.0e-9  0\n")
            f.write("\n")
            # - TIME2EQ STATEMENTS
            f.write("# +++ TIME2EQ STATEMENTS\n")
            f.write("EQFILES time/layer obs/SYNEQ.obs\n")
            f.write("EQMECH {mechType:s} {strike:5.1f} {dip:5.1f} {rake:5.1f}\n".format(
                mechType=mechType, strike=strike, dip=dip, rake=rake
            ))
            f.write("EQMODE SRCE_TO_STA\n")
            f.write("\n")
            for i, (ot, lat, lon, dep, parrivals, sarrivals) in enumerate(zip(OTs[1], LATs[1], LONs[1], Deps[1], PArrivalsList[1], SArrivalsList[1])):
                f.write("EQSRCE SYNEQ_{i:d} LATLON {eventLat:7.3f} {eventLon:7.3f} {eventDep:5.1f} 0.0\n".format(
                    i=i, eventLat=lat, eventLon=lon, eventDep=dep
                ))
                for parrival in parrivals:
                    errorP = self.config["ErrorOnArrivals"]["P"]
                    f.write("EQSTA {staCode:4s} P GAU {errorP:6.2f} GAU  0.00\n".format(staCode=parrival, errorP=errorP))
                for sarrival in sarrivals:
                    errorS = self.config["ErrorOnArrivals"]["S"]
                    f.write("EQSTA {staCode:4s} S GAU {errorS:6.2f} GAU  0.00\n".format(staCode=sarrival, errorS=errorS))
                f.write("\n")
            # - Quality to Error Mapping
            f.write("# +++ Quality to Error Mapping\n")
            f.write("EQQUAL2ERR 0.1 0.2 0.4 0.8 99999.9\n")
            # - VpVs
            f.write("# +++ P Velocity to S Velocity Ratio\n")
            f.write("EQVPVS {VpVs:5.2f}\n".format(VpVs=self.velocityModel["VpVs"]))
    
    # Run NLLOC
    def runNLloc(self):
        # - Make velocity grid files for P and S 
        cmd = "Vel2Grid nlloc.conf"
        os.system(cmd)
        print("+++ P and S velocity grids have been created successfull.")
        # - Generate P travel time grid files for stations
        cmd = "Grid2Time nlloc.conf"
        os.system(cmd)
        # - Now for S
        cmd = "sed -i 's/layer P/layer S/g' nlloc.conf"
        os.system(cmd)
        cmd = "Grid2Time nlloc.conf"
        os.system(cmd)
        print("+++ P and S travel time grids have been created successfull.")
        # - Now generate synthetics
        cmd = "Time2EQ nlloc.conf > /dev/null"
        os.system(cmd)
        nEQ = len(self.earthquakesInfo["OT"])
        print("+++ {nEQ:d} synthetic earthquakes have been generated successfull.".format(nEQ=nEQ))


    # Convert NLLOC to NORDIC
    def convertNLLOC2NORDIC(self):
        pass

if __name__ == "__main__":
    app = Main()
    app.readConfiguration()
    app.checkRequiredFiles()
    app.parseVelocityModel()
    app.writeNLlocVelocityModel()
    app.parseStationInfo()
    app.writeNLlocStationFile()
    app.parseEarthquakeFile()
    app.writeNLlocControlFile()
    app.runNLloc()
