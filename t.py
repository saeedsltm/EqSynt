import obspy

event = obspy.core.event.Event()


def getPick(staCode, phaseHint):
    pick = obspy.core.event.Pick()
    pick.phase_hint = phaseHint
    pick.time = obspy.UTCDateTime(2020,1,1)
    pick.waveform_id = obspy.core.event.WaveformStreamID("BI", staCode)
    return pick
    
    
def getArival(staCode, phase, pick_id):
    arrival = obspy.core.event.Arrival()
    arrival.phase = phase
    arrival.time = obspy.UTCDateTime(2020,1,1)
    arrival.pick_id = pick_id
    return arrival

def getMagnitude(mag, magType, origin):
    magnitude = obspy.core.event.Magnitude()
    magnitude.mag = mag
    magnitude.magnitude_type = magType
    magnitude.origin_id = origin.resource_id
    return magnitude
    
picks, arrivals = [], []
arrivalsDict = {"AAA":{"P":["Pg"], "S":["Sg"]},
                "BBB":{"P":["Pn"], "S":[]},
                "CCC":{"P":["P"], "S":["Sg"]},
                "DDD":{"P":["Pn"], "S":["Sg"]},
                "EEE":{"P":["Pg"], "S":["Sg"]},
                "FFF":{"P":["Pg"], "S":["Sg"]},
                }
for staCode in arrivalsDict:
    phaseP = arrivalsDict[staCode]["P"]
    phaseS = arrivalsDict[staCode]["S"]
    for p in phaseP:
        pick = getPick(staCode, p)
        arrival = getArival(staCode, p, pick.resource_id)
        picks.append(pick)
        arrivals.append(arrival)

origin = obspy.core.event.Origin()
origin.time = obspy.UTCDateTime(2020,1,1)
origin.latitude = 32.0
origin.longitude = 52.0
origin.depth = 10000.0
origin.arrivals = arrivals

magnitude = getMagnitude(4.3, "Ml", origin)
event.origins.append(origin)
event.picks = picks
event.magnitudes.append(magnitude)
catalog = obspy.core.event.Catalog()
catalog.append(event)
