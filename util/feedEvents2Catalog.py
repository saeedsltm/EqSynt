from obspy.core import event
from random import gauss
from obspy.geodetics.base import kilometer2degrees as k2d

class feedCatalog():
    def __init__(self):
        pass

    def setPick(self, staCode, onset, phaseHint, weight, time):
        pick = event.Pick()
        pick.onset = onset
        pick.phase_hint = phaseHint
        pick.update({"extra":{"nordic_pick_weight":{"value":0}}})
        pick.extra.nordic_pick_weight.value = weight
        pick.time = time
        pick.waveform_id = event.WaveformStreamID("BI", staCode)
        return pick
        
    def setArrival(self, phase, time, distance, azimuth, pick_id):
        arrival = event.Arrival()
        arrival.phase = phase
        arrival.time = time
        arrival.distance = k2d(distance)
        arrival.azimuth = azimuth
        arrival.pick_id = pick_id
        return arrival

    def setMagnitude(self, mag, magType, origin):
        magnitude = event.Magnitude()
        magnitude.mag = mag
        magnitude.magnitude_type = magType
        magnitude.origin_id = origin.resource_id
        return magnitude
    
    def getPicksArrivals(self, arrivalsDict):   
        picks, arrivals = [], []
        for staCode in arrivalsDict:
            phaseP = arrivalsDict[staCode]["P"]
            phaseS = arrivalsDict[staCode]["S"]
            for Phase in [phaseP, phaseS]:
                phase, onset, weight, time, distance, azimuth = Phase
                pick = self.setPick(staCode, onset, phase, weight, time)
                arrival = self.setArrival(phase, time, distance, azimuth, pick.resource_id)
                picks.append(pick)
                arrivals.append(arrival)
        return picks, arrivals
    
    def setOrigin(self, eventInfo, arrivals):
        origin = event.Origin()
        origin.time = eventInfo["OriginTime"]
        origin.latitude = eventInfo["Latitude"]
        origin.longitude = eventInfo["Longitude"]
        origin.depth = eventInfo["Depth"]*1e3
        origin.arrivals = arrivals
        return origin
    
    def setEvent(self, eventInfo, arrivalsDict):
        picks, arrivals = self.getPicksArrivals(arrivalsDict)                
        Event = event.Event()
        origin = self.setOrigin(eventInfo, arrivals)
        magnitude = self.setMagnitude(0.5*len(arrivalsDict), "Ml", origin)
        Event.origins.append(origin)
        Event.picks = picks
        Event.magnitudes.append(magnitude)
        return Event
    
    def setCatalog(self, eventsInfo, eventArrivals):
        catalog = event.Catalog()
        for eventInfo, arrivalsDict in zip(eventsInfo, eventArrivals):
            Event = self.setEvent(eventInfo, arrivalsDict)
            catalog.append(Event)
        return catalog

# if __name__ == "__main__":
#     eventArrivals = [{"AAA":{"P":["Pg", utc(2020,1,1)], "S":["Sg", utc(2020,1,1)]},
#                  "BBB":{"P":["Pn", utc(2020,1,1)], "S":[]},
#                  "CCC":{"P":["P", utc(2020,1,1)], "S":["Sg", utc(2020,1,1)]},
#                  "DDD":{"P":["Pn", utc(2020,1,1)], "S":["Sg", utc(2020,1,1)]},
#                  "EEE":{"P":["Pg", utc(2020,1,1)], "S":["Sg", utc(2020,1,1)]},
#                  "FFF":{"P":["Pg", utc(2020,1,1)], "S":["Sg", utc(2020,1,1)]},
#                 }]
#     app = feedCatalog()
#     catalog = app.setCatalog(eventArrivals)
