from obspy.core import event
from obspy.geodetics.base import kilometer2degrees as k2d


class feedCatalog():
    """An obspy Catalog Contractor
    """

    def __init__(self):
        pass

    def setPick(self, staCode, onset, phaseHint, weight, time):
        """Fill the Obspy pick object

        Args:
            staCode (str): station code,
            onset (str): onset of phase
            phaseHint (str): define P or S phase
            weight (int): dedicated weight to the phase
            time (time): arrival time of phase

        Returns:
            obspy.pick: an obspy pick object
        """
        pick = event.Pick()
        pick.onset = onset
        pick.phase_hint = phaseHint
        pick.update({"extra": {"nordic_pick_weight": {"value": 0}}})
        pick.extra.nordic_pick_weight.value = weight
        pick.time = time
        pick.waveform_id = event.WaveformStreamID("BI", staCode)
        return pick

    def setArrival(self, phase, time, distance, azimuth, pick_id):
        """Fill Obspy arrival object

        Args:
            phase (str): phase code
            time (time): phase arrival time
            distance (float): epicentral distance
            azimuth (float): azimuth between event and station
            pick_id (str): string representing associated pick

        Returns:
            obspy.event.arrival: an obspy arrival object
        """
        arrival = event.Arrival()
        arrival.phase = phase
        arrival.time = time
        arrival.distance = k2d(distance)
        arrival.azimuth = azimuth
        arrival.pick_id = pick_id
        return arrival

    def setMagnitude(self, mag, magType, origin):
        """Fill Obspy magnitude object

        Args:
            mag (float): magnitude of event
            magType (str): type of magnitude
            origin (obspy.event.origin): an obspy event origin object

        Returns:
            obspy.event.magnitude: an obspy magnitude object
        """
        magnitude = event.Magnitude()
        magnitude.mag = mag
        magnitude.magnitude_type = magType
        magnitude.origin_id = origin.resource_id
        return magnitude

    def getPicksArrivals(self, arrivalsDict):
        """Make pick and arrival from input dictionary

        Args:
            arrivalsDict (dict): a dictionary contains arrival and pick

        Returns:
            tuple: a tuple contains obspy pick and arrival objects
        """
        picks, arrivals = [], []
        for staCode in arrivalsDict:
            phaseP = arrivalsDict[staCode]["P"]
            phaseS = arrivalsDict[staCode]["S"]
            for Phase in [phaseP, phaseS]:
                phase, onset, weight, time, distance, azimuth = Phase
                pick = self.setPick(staCode, onset, phase, weight, time)
                arrival = self.setArrival(
                    phase, time, distance, azimuth, pick.resource_id)
                picks.append(pick)
                arrivals.append(arrival)
        return picks, arrivals

    def setOrigin(self, eventInfoDict, arrivals):
        """Fill Obspy origin object

        Args:
            eventInfoDict (dict): a dictionary contains event information
            arrivals (obspy.arrival): an obspy arrivals

        Returns:
            obspy.origin: an obspy origin object
        """
        origin = event.Origin()
        origin.time = eventInfoDict["OriginTime"]
        origin.latitude = eventInfoDict["Latitude"]
        origin.longitude = eventInfoDict["Longitude"]
        origin.depth = eventInfoDict["Depth"]*1e3
        origin.arrivals = arrivals
        return origin

    def setEvent(self, eventInfoDict, arrivalsDict):
        """Fill Obspy event object

        Args:
            eventInfoDict (dict): a dictionary contains event information
            arrivalsDict (dict): a dictionary contains event arrivals

        Returns:
            obspy.event: an obspy event
        """
        picks, arrivals = self.getPicksArrivals(arrivalsDict)
        Event = event.Event()
        origin = self.setOrigin(eventInfoDict, arrivals)
        magnitude = self.setMagnitude(0.5*len(arrivalsDict), "Ml", origin)
        Event.origins.append(origin)
        Event.picks = picks
        Event.magnitudes.append(magnitude)
        return Event

    def setCatalog(self, eventsInfoList, eventArrivals):
        """Fill Obspy catalog object

        Args:
            eventsInfoList (list): a list contains event information
            eventArrivals (list): a list contains event information

        Returns:
            _type_: _description_
        """
        catalog = event.Catalog()
        for eventInfo, arrivalsDict in zip(eventsInfoList, eventArrivals):
            Event = self.setEvent(eventInfo, arrivalsDict)
            catalog.append(Event)
        return catalog
