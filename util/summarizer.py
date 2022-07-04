import os
from statistics import mean
from pandas import DataFrame
from util.extra import distanceDiff
from obspy.geodetics.base import degrees2kilometers as d2k


def countNP(df, df_metrics, key, column):
    """Count number and percentage of specific event

    Args:
        df (DataFrame): data frame of catalog
        df_metrics (dict): a dictionary contains metrics
        key (str): desired key of data frame
        column (str): desired column of data frame

    Returns:
        tuple: number and percentage  of event
    """
    n = df[df_metrics[key]][column].size
    p = n/df[column].size*100.0
    return n, p


def computeClass(df_un, df_w, c, config):
    """Compute class statistics

    Args:
        df_un (data frame): unweighted data frame of events
        df_w (data frame): weighted data frame of events
        c (str): defined class
        config (dict): a dictionary contains main configuration

    Returns:
        tuple: number and percentage  of classes in unweighted and weighted catalog
    """
    Neq_un = df_un["ORT"].size
    Neq_w = df_w["ORT"].size
    ERH = config["RPS"]["Classes"][c]["ERH"]
    ERZ = config["RPS"]["Classes"][c]["ERZ"]
    GAP = config["RPS"]["Classes"][c]["GAP"]
    RMS = config["RPS"]["Classes"][c]["RMS"]
    MDS = config["RPS"]["Classes"][c]["MDS"]
    NuP = config["RPS"]["Classes"][c]["NuP"]
    NuS = config["RPS"]["Classes"][c]["NuS"]
    df_un = df_un[(df_un["ERH"] <= ERH) & (df_un["ERZ"] <= ERZ) & (df_un["GAP"] <= GAP) & (
        df_un["RMS"] <= RMS) & (df_un["MDS"] <= MDS) & (df_un["NuP"] >= NuP) & (df_un["NuS"] >= NuS)]
    df_w = df_w[(df_w["ERH"] <= ERH) & (df_w["ERZ"] <= ERZ) & (df_w["GAP"] <= GAP) & (
        df_w["RMS"] <= RMS) & (df_w["MDS"] <= MDS) & (df_w["NuP"] >= NuP) & (df_w["NuS"] >= NuS)]
    return df_un["ORT"].size, df_un["ORT"].size/Neq_un*100.0, df_w["ORT"].size, df_w["ORT"].size/Neq_w*100.0


def summarize(df_un, df_w, stationsDict, config, resultPath):
    """Make a statistical summary between unweighted and weighted catalogs

    Args:
        df_un (data frame): unweighted data frame of events
        df_w (data frame): weighted data frame of events
        config (dict): a dictionary contains main configuration
        resultPath (str): path to result directory
    """
    metrics_df_un = {
        "ERH_2km": 'df_un["ERH"]<=2.0',
        "ERH_5km": 'df_un["ERH"]<=5.0',
        "ERZ_2km": 'df_un["ERZ"]<=2.0',
        "ERZ_5km": 'df_un["ERZ"]<=5.0',
        "RMS_0_1s": 'df_un["RMS"]<=0.1',
        "RMS_0_3s": 'df_un["RMS"]<=0.3',
        "GAP_150d": 'df_un["GAP"]<=150',
        "GAP_200d": 'df_un["GAP"]<=200',
        "GAP_250d": 'df_un["GAP"]<=250',
        "MDS_5km": 'df_un["MDS"]<=5',
        "MDS_10km": 'df_un["MDS"]<=10',
        "MDS_15km": 'df_un["MDS"]<=15',
        "Nus_n": 'df_un["Nus"]>=5',
        "NuP_n": 'df_un["NuP"]>=5',
        "NuS_n": 'df_un["NuS"]>=5'
    }
    metrics_df_w = {
        "ERH_2km": 'df_w["ERH"]<=2.0',
        "ERH_5km": 'df_w["ERH"]<=5.0',
        "ERZ_2km": 'df_w["ERZ"]<=2.0',
        "ERZ_5km": 'df_w["ERZ"]<=5.0',
        "RMS_0_1s": 'df_w["RMS"]<=0.1',
        "RMS_0_3s": 'df_w["RMS"]<=0.3',
        "GAP_150d": 'df_w["GAP"]<=150',
        "GAP_200d": 'df_w["GAP"]<=200',
        "GAP_250d": 'df_w["GAP"]<=250',
        "MDS_5km": 'df_w["MDS"]<=5',
        "MDS_10km": 'df_w["MDS"]<=10',
        "MDS_15km": 'df_w["MDS"]<=15',
        "Nus_n": 'df_w["Nus"]>=5',
        "NuP_n": 'df_w["NuP"]>=5',
        "NuS_n": 'df_w["NuS"]>=5'
    }
    for metric in metrics_df_un.keys():
        metrics_df_un[metric] = eval(metrics_df_un[metric])
    for metric in metrics_df_w.keys():
        metrics_df_w[metric] = eval(metrics_df_w[metric])

    station_df = DataFrame(stationsDict).T

    nSta = station_df["Lat"].size
    nEvt = df_un["Lat"].size
    mDistEvt = d2k(mean([mean(distanceDiff(x, df_un["Lon"], y, df_un["Lat"]))
                         for x, y in zip(df_un["Lon"], df_un["Lat"])]))
    mDistSta = d2k(mean([mean(distanceDiff(x, station_df["Lon"], y, station_df["Lat"]))
                         for x, y in zip(station_df["Lon"], station_df["Lat"])]))
    mDistEvtSta = mean(df_un["ADS"])

    with open(os.path.join(resultPath, "summary.dat"), "w") as f:
        f.write("----Problem status:\n")
        f.write("Number of events: {nEvt:.0f}\n".format(nEvt=nEvt))
        f.write("Number of stations: {nSta:.0f}\n".format(nSta=nSta))
        f.write("Average distance between events (km): {mDistEvt:.2f}\n".format(
            mDistEvt=mDistEvt))
        f.write("Average distance between stations (km): {mDistSta:.2f}\n".format(
            mDistSta=mDistSta))
        f.write(
            "Average distance between events-stations pair (km): {mDistEvtSta:.1f}\n".format(mDistEvtSta=mDistEvtSta))

        f.write("----Statistics:\n")
        f.write("".center(50, "="))
        f.write("| Events relocated without weighting scheme |".center(40, "="))
        f.write("| Events relocated using weighting scheme |".center(40, "="))
        f.write("\n")

        f.write("+++ Number(%) of event with ERH<2.0km:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "ERH_2km", "ERH")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "ERH_2km", "ERH")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with ERH<5.0km:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "ERH_5km", "ERH")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "ERH_5km", "ERH")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with ERZ<2.0km:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "ERZ_2km", "ERZ")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "ERZ_2km", "ERZ")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with ERZ<5.0km:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "ERZ_5km", "ERZ")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "ERZ_5km", "ERZ")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with RMS<0.1s:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "RMS_0_1s", "RMS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "RMS_0_1s", "RMS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with RMS<0.3s:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "RMS_0_3s", "RMS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "RMS_0_3s", "RMS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with GAP<150:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "GAP_150d", "GAP")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "GAP_150d", "GAP")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with GAP<200:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "GAP_200d", "GAP")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "GAP_200d", "GAP")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with GAP<250:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "GAP_250d", "GAP")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "GAP_250d", "GAP")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with MDS<5:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "MDS_5km", "MDS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "MDS_5km", "MDS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with MDS<10:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "MDS_10km", "MDS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "MDS_10km", "MDS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with MDS<15:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "MDS_15km", "MDS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "MDS_15km", "MDS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with Nus>5:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "Nus_n", "Nus")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "Nus_n", "Nus")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with NuP>5:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "NuP_n", "NuP")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "NuP_n", "NuP")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        f.write("+++ Number(%) of event with NuS>5:".ljust(50, " "))
        n, p = countNP(df_un, metrics_df_un, "NuS_n", "NuS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        n, p = countNP(df_w, metrics_df_w, "NuS_n", "NuS")
        f.write("{n:.0f}({p:.1f})".format(n=n, p=p).rjust(40, " "))
        f.write("\n")

        # classes
        for c in config["RPS"]["Classes"].keys():
            f.write(" Class {c:s} ".format(c=c).center(50, "="))
            f.write("| Events relocated without weighting scheme |".center(40, "="))
            f.write("| Events relocated using weighting scheme |".center(40, "="))
            f.write("\n")
            NoEqInClass_un, PerEqInClass_un, NoEqInClass_w, PerEqInClass_w = computeClass(
                df_un, df_w, c, config)
            f.write("+++ Number(%) of event:".ljust(50, " "))
            f.write("{n:.0f}({p:.1f})".format(
                n=NoEqInClass_un, p=PerEqInClass_un).rjust(40, " "))
            f.write("{n:.0f}({p:.1f})".format(
                n=NoEqInClass_w, p=PerEqInClass_w).rjust(40, " "))
            f.write("\n")
