from LatLon import lat_lon as ll
from datetime import datetime as dt
from datetime import timedelta as td
from statistics import mean
from json import load
from pathlib import Path
from shutil import copy, move
from glob import glob
from tqdm.contrib import tzip
import os
import sys

# Define main class


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
        self.stationFile = os.path.join(
            "EqInput", self.config["Inputs"]["StationFile"])
        if not os.path.exists(self.eqFile) or not os.path.exists(self.stationFile):
            print("+++ Could not find earthquake or station files! Aborting ...")
            sys.exit()
        # - Create NLloc required directories
        Path("time").mkdir(parents=True, exist_ok=True)
        Path("model").mkdir(parents=True, exist_ok=True)
        Path("obs").mkdir(parents=True, exist_ok=True)
        Path("relocation").mkdir(parents=True, exist_ok=True)
        Path("tmp").mkdir(parents=True, exist_ok=True)

    # Clear existing files
    def clearExistingFiles(self):
        for d in ["model", "obs", "time", "tmp", "relocation"]:
            for f in glob(os.path.join(d, "*")):
                os.remove(f)
        print("+++ Removing prexisting files in 'model', 'obs', 'time' and 'tmp' directories ...")

    # Parse velocity model
    def parseVelocityModel(self):
        emptyLines = 0
        self.velocityModel = {"Vp": [], "Z": [], "VpVs": 1.73}
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
                    Z=Z, Vp=Vp, Vs=Vp/VpVs[1]
                ))

    # Parse station information
    def parseStationInfo(self):
        emptyLines = 0
        self.stations = {"Code": [], "Lat": [], "Lon": [], "Elv": []}
        with open(self.stationFile) as f:
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
                    for key, value in zip(["Code", "Lat", "Lon", "Elv"], [code, lat, lon, elv]):
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
        for i in [6, 8, 11, 12, 13, 14, 16, 17, 19]:
            if l[i] == " ":
                l = l[:i]+"0"+l[i+1:]
        try:
            ot = dt.strptime(l[:20], " %Y %m%d %H%M %S.%f")
            lat = float(l[24:30])
            lon = float(l[32:38])
            dep = float(l[39:43])
            return ot, lat, lon, dep
        except ValueError:
            return None, None, None, None

    # Parse earthquake file
    def parseEarthquakeFile(self):
        self.earthquakesInfo = {
            "OT": [],
            "LAT": [],
            "LON": [],
            "DEP": [],
            "PArrivals": [],
            "SArrivals": []}
        with open(self.eqFile) as f:
            for l in f:
                if l.strip() and l[79] == "1":
                    ot, lat, lon, dep = self.parseEventLine(l)
                    if ot == None:
                        continue
                    self.earthquakesInfo["OT"].append(ot)
                    self.earthquakesInfo["LAT"].append(lat)
                    self.earthquakesInfo["LON"].append(lon)
                    self.earthquakesInfo["DEP"].append(dep)
                    self.earthquakesInfo["PArrivals"].append([])
                    self.earthquakesInfo["SArrivals"].append([])
                if l.strip() and l[10] == "P" and l[79] in [" ", "4"]:
                    staCode = l[:5].strip()
                    self.earthquakesInfo["PArrivals"][-1].append(staCode)
                if l.strip() and l[10] == "S" and l[79] in [" ", "4"]:
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
        lonOrig = mean([lonMin, lonMax])
        latOrig = mean([latMin, latMax])
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
            # - Quality to Error Mapping
            f.write("# +++ Quality to Error Mapping\n")
            f.write("EQQUAL2ERR 0.1 0.2 0.4 0.8 99999.9\n")
            # - VpVs
            f.write("# +++ P Velocity to S Velocity Ratio\n")
            f.write("EQVPVS {VpVs:5.2f}\n".format(
                VpVs=self.velocityModel["VpVs"]))
            f.write("\n")
        # Run NLLOC for initiation
        self.runNLlocForward()
        # Noe generating synthetics
        print("+++ Generating synthetic earthquakes ...")
        for i, (ot, lat, lon, dep, parrivals, sarrivals) in enumerate(tzip(OTs[1], LATs[1], LONs[1], Deps[1], PArrivalsList[1], SArrivalsList[1])):
            outConfig = os.path.join("tmp", "nlloc_{i:d}.conf".format(i=i+1))
            copy("nlloc.conf", outConfig)
            with open(outConfig, "a") as f:
                # - TIME2EQ STATEMENTS
                f.write("# +++ TIME2EQ STATEMENTS\n")
                f.write(
                    "EQFILES time/layer obs/SYNEQ_{i:d}.obs\n".format(i=i+1))
                f.write("EQMECH {mechType:s} {strike:5.1f} {dip:5.1f} {rake:5.1f}\n".format(
                    mechType=mechType, strike=strike, dip=dip, rake=rake
                ))
                f.write("EQMODE SRCE_TO_STA\n")
                f.write("\n")
                f.write("EQSRCE SYNEQ_{i:d} LATLON {eventLat:7.3f} {eventLon:7.3f} {eventDep:5.1f} 0.0\n".format(
                    i=i, eventLat=lat, eventLon=lon, eventDep=dep
                ))
                for parrival in parrivals:
                    errorP = self.config["ErrorOnArrivals"]["P"]
                    f.write("EQSTA {staCode:4s} P GAU {errorP:6.2f} GAU {errorP:6.2f}\n".format(
                        staCode=parrival, errorP=errorP))
                for sarrival in sarrivals:
                    errorS = self.config["ErrorOnArrivals"]["S"]
                    f.write("EQSTA {staCode:4s} S GAU {errorS:6.2f} GAU {errorS:6.2f}\n".format(
                        staCode=sarrival, errorS=errorS))
                f.write("\n")
            # - Synthetic arrival times
            cmd = "Time2EQ {outConfig:s} > /dev/null".format(
                outConfig=outConfig)
            os.system(cmd)
            self.correctOT(ot, os.path.join(
                "obs", "SYNEQ_{i:d}.obs".format(i=i+1)))

    # Correct origin time
    def correctOT(self, ot, syntheticEqFile):
        with open(syntheticEqFile) as f, open("tmpFile", "w") as g:
            for i, l in enumerate(f):
                if i == 0:
                    line = l.split()
                    line.append("{0:f}".format(ot.timestamp()))
                    line.pop(-2)
                    l = " ".join(line)
                    g.write(l+"\n")
                elif l.strip() and "#" not in l and "PUBLIC_ID" not in l:
                    line = l.split()
                    arrivalTime = " ".join(line[6:9])
                    arrivalTime = dt.strptime(arrivalTime, "%Y%m%d %H%M %S.%f")
                    arrivalTime = arrivalTime - dt(1900, 1, 1, 0, 0, 0)
                    arrivalTime = ot + td(seconds=arrivalTime.total_seconds())
                    date = arrivalTime.strftime("%Y%m%m")
                    hourMin = arrivalTime.strftime("%H%M")
                    sec = arrivalTime.strftime("%S.%f")
                    g.write(
                        "{stationName:9s} {instrument:4s} {component:4s} {pPhaseOnset:4s} {PhaseDes:6s} {firstMotion:1s} {date:8s} {hourMin:6s} {sec:7s} {err:3s} {errMag:9s} {codaDur:9s} {amplitude:9s} {period:9s} {priorWt:9s}\n".format(
                            stationName=line[0], instrument=line[1], component=line[2], pPhaseOnset=line[3],
                            PhaseDes=line[4], firstMotion=line[5], date=date, hourMin=hourMin, sec=sec, err=line[9],
                            errMag=line[10], codaDur=line[11], amplitude=line[12], period=line[13], priorWt=line[14]
                        ))
        move("tmpFile", syntheticEqFile)

    # Run NLLOC for generating velocity grid and travel time tables
    def runNLlocForward(self):
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

    # Merge synthetic files
    def mergeSyntheticFiles(self):
        syntheticFiles = sorted(glob(os.path.join(
            "obs", "SYNEQ_*.obs")), key=lambda x: int(x.split("_")[1].split(".")[0]))
        for syntheticFile in syntheticFiles:
            cmd = "cat {syntheticFile:s} >> obs/SYNEQ.obs".format(
                syntheticFile=syntheticFile)
            os.system(cmd)
        for _ in glob(os.path.join("obs", "SYNEQ_*")):
            os.remove(_)

    # Convert NLLOC to NORDIC
    def convertNLLOC2NORDIC(self):
        nllocObsFile = os.path.join("obs/SYNEQ.obs")
        if not nllocObsFile:
            print("+++ Could not find synthetic earthquake file! Aborting ...")
            sys.exit(0)
        hypo_file = open(os.path.join("relocation", "synthetic.out"), "w")
        with open(nllocObsFile) as f:
            line1_flag = False
            line7_flag = False
            start_flag = False
            for l in f:
                if "# EQEVENT:" in l:
                    line1_flag = True
                    line7_flag = True
                    if start_flag:
                        hypo_file.write("\n")
                        start_flag = False
                elif "#" not in l:
                    start_flag = True
                    if line1_flag and line7_flag:
                        y = int(l[34:38])
                        m = int(l[38:40])
                        d = int(l[40:42])
                        H = int(l[43:45])
                        M = int(l[45:47])
                        S = 0
                        MS = 0
                        ot = dt(y, m, d, H, M, S, MS)
                        line_1 = " %4d %02d%02d %02d%02d %4.1f L                                                         1\n"\
                                 % (ot.year, ot.month, ot.day, ot.hour, ot.minute, ot.second)
                        line_7 = " STAT SP IPHASW D HRMN SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W  DIS CAZ7\n"
                        hypo_file.write(line_1)
                        hypo_file.write(line_7)
                        line1_flag = False
                        line7_flag = False
                    line = l.strip().split()
                    sta_name = line[0]
                    Instrument = " "
                    comp = " "
                    Quality_Indicator = " "
                    Phase_ID = line[4]+"   "
                    Weighting_Indicator = 0
                    if line[4] == "P":
                        First_Motion = line[5]
                    else:
                        First_Motion = " "
                    arival_t = dt(int(line[6][0:4]), int(line[6][4:6]), int(line[6][6:8]),
                                  int(line[7][0:2]), int(line[7][2:4]), int(
                                      line[8].split(".")[0]),
                                  int(line[8].split(".")[1]))
                    line_4 = " %-5s%1s%1s %1s%-4s%1d %1s %02d%02d %5.2f                                                   4\n"\
                             % (sta_name, Instrument, comp, Quality_Indicator, Phase_ID, Weighting_Indicator, First_Motion, arival_t.hour, arival_t.minute, arival_t.second+arival_t.microsecond*1e-6)
                    hypo_file.write(line_4)
        hypo_file.write("\n")
        hypo_file.close()
        print("+++ Output file 'nlloc.hyp' has been created.\n")
        # - Copy nordic station file into 'relocation' directory
        copy(self.stationFile, os.path.join("relocation", "STATION0.HYP"))

    # Relocate data using "hypocenter" program
    def relocateWithHyp(self, targetDir, nordicFileName):
        root = os.getcwd()
        os.chdir(targetDir)
        with open("hyp.inp", "w") as f:
            f.write("{nordicFileName:s}\n".format(
                nordicFileName=nordicFileName))
            f.write("n\n")
        cmd = "hyp < hyp.inp > /dev/null"
        os.system(cmd)
        os.chdir(root)

    # Make report
    def makeReport(self, targetDir, nordicFileName):
        root = os.getcwd()
        os.chdir(targetDir)
        with open("report.inp", "w") as f:
            f.write("1    2         3  6 4  7 5   8       9    10  11\n")
        cmd = "report {nordicFileName:s} < report.inp > /dev/null".format(
            nordicFileName=nordicFileName)
        os.system(cmd)
        os.chdir(root)

    # Compare raw and relocated data
    def compare(self):
        self.makeReport(targetDir="EqInput", nordicFileName="select.out")
        self.relocateWithHyp(targetDir="relocation", nordicFileName="synthetic.out")
        self.makeReport(targetDir="relocation", nordicFileName="hyp.out")


# Run application
if __name__ == "__main__":
    app = Main()
    app.readConfiguration()
    app.checkRequiredFiles()
    app.clearExistingFiles()
    app.parseVelocityModel()
    app.writeNLlocVelocityModel()
    app.parseStationInfo()
    app.writeNLlocStationFile()
    app.parseEarthquakeFile()
    app.writeNLlocControlFile()
    app.mergeSyntheticFiles()
    app.convertNLLOC2NORDIC()
    app.compare()
