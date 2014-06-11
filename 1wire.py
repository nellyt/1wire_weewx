#
#    Copyright (c) 2013 Neil Trimboy
#
#    See the file LICENSE.txt for your full rights.
#
#    $Revision: $
#    $Author: $
#    $Date: $
#
"""1wire for the weewx weather system"""

from __future__ import with_statement
import math
import time
import syslog

import weedb
import weeutil.weeutil
import weewx.abstractstation
import weewx.wxformulas
import Pysolar
import datetime

import ow

from collections import deque # For wind pipe

import os

def logmsg(dst, msg):
    syslog.syslog(dst, 'owfs: %s' % msg)

def logdbg(msg):
    logmsg(syslog.LOG_DEBUG, msg)

def loginf(msg):
    logmsg(syslog.LOG_INFO, msg)

def logcrt(msg):
    logmsg(syslog.LOG_CRIT, msg)

def logerr(msg):
    logmsg(syslog.LOG_ERR, msg)
    
def loader(config_dict, engine):

    # TODO see other drivers for altitude to see (maybe ) how to get lat/long into driver
    start_ts = resume_ts = None

    station = onewire(start_time=start_ts, resume_time=resume_ts, **config_dict['1wire'])
     # TODO check filesystem exists here
    return station


class onewire(weewx.abstractstation.AbstractStation):
    """Station 1Wire"""

    def __init__(self, **stn_dict):
        """Initialize the simulator

        # TODO Some of this next stuff is redundant as it come from Simulator
        NAMED ARGUMENTS:
        
        interface: Where to find the one-wire sensors.  Options include
        u, /dev/ttyS0
        [Required. Default is u (usb)]
        

        loop_interval: The time (in seconds) between emitting LOOP packets. [Optional. Default is 2.5]

        start_time: The start (seed) time for the generator in unix epoch time [Optional. If 'None',
                    or not present, then the present time will be used.]

        resume_time: The start time for the loop. [Optional. If 'None',
                     or not present, then start_time will be used].


        """

        logdbg("Station 1WIRE created")
        self.interface         = stn_dict.get('interface', 'u')
        self.LOOP_INTERVAL = float(stn_dict.get('loop_interval', 1))
        self.RAIN_INTERVAL = float(stn_dict.get('rain_interval', 30))
        self.OTHERS_INTERVAL = float(stn_dict.get('rain_interval', 30))
        
        # No start time specified. We are in realtime mode.
        self.real_time = True
        start_ts = self.the_time = time.time()

        #self.mode = stn_dict['mode']
        
        self.wind = WindSpeed()
        self.rain = RainCount()

        self.last_rain_time = time.time()
        self.last_other_time = time.time()

        # The following doesn't make much meteorological sense, but it is easy to program!
        self.observations = {'outTemp'    : inTemp(),
                             'inTemp'     : inTemp(),
                             'barometer'  : Barometer(),
                             'outHumidity': inHumid(),
                             'inHumidity' : inHumid(),
                             'radiation'  : Solar(),
                             'maxInstantRadiation' : MaxSolar()}
                             
        ow.init(self.interface)

    def genLoopPackets(self):

        while True:

            # Determine how long to sleep
            if self.real_time:
                # We are in real time mode. Try to keep synched up with the
                # wall clock
                sleep_time = self.the_time + self.LOOP_INTERVAL - time.time()
                if sleep_time > 0:
                    time.sleep(sleep_time)
            else:
                # A start time was specified, so we are not in real time. Just sleep
                # the appropriate interval
                time.sleep(self.LOOP_INTERVAL)

            # Update the simulator clock:
            self.the_time += self.LOOP_INTERVAL

            # Because a packet represents the measurements observed over the
            # time interval, we want the measurement values at the middle
            # of the interval.
            avg_time = self.the_time - self.LOOP_INTERVAL/2.0

            if (time.time() - self.last_rain_time > self.RAIN_INTERVAL) :  
                # TODO put assign last time here
                _rainpacket = {'dateTime': int(self.the_time+0.5),
                           'usUnits' : weewx.US }
                # calculate the rain increment from the rain total
                raindata = self.rain.value_at(avg_time)
                if raindata is None:
                    _rainpacket['rain'] = None
                    _rainpacket['rainRate'] = None
                    _rainpacket['rainRateA'] = None
                    _rainpacket['rainRateB'] = None
                else:
                    _rainpacket['rain'] = raindata[0]
                    _rainpacket['rainRate'] = raindata[1]
                    _rainpacket['rainRateA'] = raindata[2]
                    _rainpacket['rainRateB'] = raindata[3]
                yield _rainpacket
          self.last_rain_time = time.time()
            
            _windpacket = {'dateTime': int(self.the_time+0.5),
                       'usUnits' : weewx.METRIC }
            wspeeds = self.wind.value_at(avg_time)
            if wspeeds is None:
                _windpacket['windSpeed'] = None
                _windpacket['windDir'] = None
                _windpacket['windGust'] = None
                _windpacket['windGustDir'] = None
            else:
                _windpacket['windSpeed'] = wspeeds[0]
                _windpacket['windDir'] = wspeeds[2]
                _windpacket['windGust'] = wspeeds[1]
                _windpacket['windGustDir'] = wspeeds[3]
            yield _windpacket

            if (time.time() - self.last_other_time > self.OTHERS_INTERVAL) :
                #TODO put assign last time here
    	        _packet = {'dateTime': int(self.the_time+0.5),
                           'usUnits' : weewx.METRIC }
                for obs_type in self.observations:
                    _packet[obs_type] = self.observations[obs_type].value_at(avg_time)

                _packet['dewpoint'] = weewx.wxformulas.dewpointC(_packet['outTemp'], _packet['outHumidity'])
                _packet['heatindex'] = weewx.wxformulas.heatindexC(_packet['outTemp'], _packet['outHumidity'])
                _packet['windchill'] = weewx.wxformulas.windchillC(_packet['outTemp'], _windpacket['windSpeed'])
                yield _packet
                self.last_other_time = time.time()

    def getTime(self):
        return self.the_time

    @property
    def hardware_name(self):
        return "1Wire"

class inTemp(object):
    def __init__(self):
        """Some text"""

    def value_at(self, time_ts):
        try:
          temp = float(ow.owfs_get('10.25A5A0020800/temperature'))
          return temp
        except ow.exError, e:
          logerr("Failed attempt to get inTemp data: %s" % ( e))
          return None
          
        

class inHumid(object):
    def __init__(self):
        """Some text"""

    def value_at(self, time_ts):
        try:
          humidity = float(ow.owfs_get('26.378C21010000/HIH4000/humidity'))
          return humidity
        except ow.exError, e:
          logerr("Failed attempt to get inHumid data: %s" % ( e))
          return None    

class Barometer(object):
    def __init__(self):
        """Some text"""

    def value_at(self, time_ts):
        try:
          altimeter = float(ow.owfs_get('26.C9AABC000000/B1-R1-A/pressure'))
          return altimeter
        except ow.exError, e:
          logerr("Failed attempt to get Barometer data: %s" % ( e))
          return None        

class RainCount(object):
    def __init__(self):
        """Some text"""
        logdbg("1WIRE: RainCount created")
        self._last_rain_cnt = None
        self._last_rain_ts = None
        self._last_detected_tip_cnt = None
        self._last_detected_tip_ts = None
        self._last_detected_tip_cntb = None
        self._last_detected_tip_tsb = None
        self._this_tip_ts = None
        self._this_tip_cnt = None
        self.rainrateduration = 300 #NOTE TODO in seconds. This must be bigger than polling rate!
        self.bucket_size = 0.01 # Rainwise gauge is 0.01 inches per tip
        self.rainpipe = deque()
        self.TIME_INDEX = 0
        self.CNT_INDEX = 1

    def value_at(self, time_ts):
        """Returns rain and rain rate
        Various rain rate algorithms are implemented and returned.
        User can decide which they want and assign it to 
        """
        try:
          rainCnt = int(ow.owfs_get('1D.22AA0D000000/counters.B'))
        except ow.exError, e:
          logerr("Failed attempt to get rain data: %s" % ( e))
          return None
        # possible catch needed here
        #logdbg("1WIRE: RainCounter Count is %s %s" % (rainCnt, self._last_rain_cnt))
        if self._last_rain_cnt is None:
          # First call from power up
          self._last_rain_cnt = rainCnt
          self._last_rain_ts = time_ts
          return None
          
        if self._last_rain_cnt is not None:
            if rainCnt < self._last_rain_cnt:
                # This could be a wrap round or the DS2423 has lost power and reset to 0
                # The DS2423 counters are 32 bit
                # TODO still to verify this lot
                logdbg("1WIRE: RainCount counter wraparound detected: curr: %s last: %s (mm)" % (rainCnt, self._last_rain_cnt))
                rain = 0 #TODO
            else :
                rain = (rainCnt - self._last_rain_cnt) * self.bucket_size
        else:
            rain = None
        logdbg("1WIRE: RainCount %s %s" % (rainCnt, self._last_rain_cnt))
        
        # This returns the rain rate (weewx wants in/hr)
        #
        # rain rate may be under-reported.
        # There is a lot on the web about rainrate and the various ways to
        # calculate it, this load of code has not been sanity checked
        # http://wiki.sandaysoft.com/a/FAQ#How_is_my_rain_rate_calculated.3F
        # http://www.localweather.net.nz/smf/hardware-software-and-technology/how-does-your-station-calculate-your-rain-rate/
        # http://www.wxforum.net/index.php?topic=12188.0
        # http://www2.buoy.com/pipermail/weather/2006-December/007261.html
        
        # *** METHOD 1 ***
        # The Cumuls way
        # As described in http://wiki.sandaysoft.com/a/FAQ#How_is_my_rain_rate_calculated.3F
        # Cumulus simply takes the rain total from the last five minutes and calculates a rate based
        # on that; e.g. a single tip of 0.3mm in 5 minutes is a rate of 3.6mm/hr. 
        # TODO if count is null dont do this
        while self.rainpipe and self.rainpipe[0][self.TIME_INDEX] < time_ts - self.rainrateduration:
          # If the pipe is longer than self.rainrateduration, then truncate it
          _unused = self.rainpipe.popleft()
        self.rainpipe.append((time_ts,rainCnt)) # Put the latest count and time in the pipe
        count = len(self.rainpipe)
        if count < 2:
            logdbg("1WIRE: RainCounter Pipe is too short for Rainrate calculation")
            self._last_rain_cnt = rainCnt
            return (rain, None, None, None)
        time_delta_s = self.rainpipe[count-1][self.TIME_INDEX] - self.rainpipe[0][self.TIME_INDEX]
        cnt_delta = self.rainpipe[count-1][self.CNT_INDEX] - self.rainpipe[0][self.CNT_INDEX]
        rainrateA = (60*60) * self.bucket_size * cnt_delta/time_delta_s
        logdbg("1WIRE: Rainrate #1 Queue %s" % (self.rainpipe))
        logdbg("1WIRE: Rainrate #1 %s %s %s" % (rainrateA, time_delta_s, cnt_delta))
        

        
        # *** METHOD 2 ***
        # The Davis way
        # NOTE For this method we really want the polling rate to be quite high to avoid 
        # quantisation of time ?
        # As described in
        # http://www.davisnet.com/product_documents/weather/app_notes/AN_28-derived-weather-variables.pdf
        # http://www.wxforum.net/index.php?topic=12188.0
        #
        # TODO still need to cope with rain counter wrap round and Nones
        logdbg("1WIRE: Rainrate #2 XXX %s %s %s" % (rainCnt, self._last_detected_tip_cnt, self._last_detected_tip_ts))

        if rainCnt != self._last_rain_cnt:
          # Bucket tipped this iteration
          self._last_detected_tip_ts = self._this_tip_ts
          self._last_detected_tip_cnt = self._this_tip_cnt
          self._this_tip_ts = time_ts
          self._this_tip_cnt = rainCnt
          if (self._last_detected_tip_cnt is None) :
            # This is the first tip
            logdbg("1WIRE: Rainrate #2 Ignore First tip %s %s" % (rainCnt, self._last_detected_tip_cnt))
            self._last_detected_tip_cnt = self._this_tip_cnt
            self._last_detected_tip_ts = self._this_tip_ts
            rainrateB = None # We cant calculate a rate on first tip from power up so return
          else :
            time_delta_s = self._this_tip_ts - self._last_detected_tip_ts
            cnt_delta = self._this_tip_cnt - self._last_detected_tip_cnt
            rainrateB = (60*60)/time_delta_s * self.bucket_size * cnt_delta # same as wxformulas.calculate_rain_rate()
            # Tidy up ready for next tip
            logdbg("1WIRE: Rainrate #2 Tip %s %s %s" % (cnt_delta, time_delta_s, rainrateB))
        else :
          # Bucket did not tip this iteration
          if (self._last_detected_tip_ts is None) :
            logdbg("1WIRE: Rainrate #2 Silence %s %s" % (rainCnt, self._last_detected_tip_cnt))
            rainrateB = None # We cant calculate a rate on first tip from power up so return
          else :
            # This implements the decay refered to by Davis
            time_delta_s = time_ts - self._last_detected_tip_ts
            cnt_delta = rainCnt - self._last_detected_tip_cnt
            rainrateB = (60*60)/time_delta_s * self.bucket_size * cnt_delta # same as wxformulas.calculate_rain_rate()
            logdbg("1WIRE: Rainrate #2 Decay %s %s %s" % (cnt_delta, time_delta_s, rainrateB))
            if time_ts - self._this_tip_ts > 15 * 60 :
              # If we didnt tip this time and its is at least 15 minutes since this tip force the rate to zero
              rainrateB = 0
              logdbg("1WIRE: Rainrate #2 Truncated %s %s %s" % (cnt_delta, time_delta_s, rainrateB))
        
        # By default set rainrate to algorithm #2
        rainrate = rainrateB   
        
        self._last_rain_cnt = rainCnt
        return (rain, rainrate, rainrateA, rainrateB)

class WindSpeed(object):
    
    def __init__(self):
        """Some text"""
        logdbg("1WIRE: WindSpeed created")
        self._wind_ts = None
        self._windCnt = None
        self.GUSTDURATION = 3 # 3 second running mean
        self.MEANDURATION = (10*60) # 10 minute mean
        self.gustpipe =deque()
        self.speedpipe =deque()
        self.TIME_INDEX = 0

    def value_at(self, time_ts):
        self.MPH_PER_HZ = 2.5 
        kph_per_hz = self.MPH_PER_HZ * 1.60934400
        self.MIN_GUST_SPEED_KPH = 29.63200 # 16 knots
        self.GUST_VARIATION_SPEED_KPH = 16.68 # 9 knots
        self.MIN_GUST_SPEED_KPH = 10 # 16 knots
        self.GUST_VARIATION_SPEED_KPH = 5 # 9 knots
        
        return None
        # NOTE need to read non cached versions for wind
        # Are the counters shorted, is this for noise, could do clever stuff in here if so
        wspeedCnt = float(ow.owfs_get('uncached/1D.A4C80D000000/counters.A'))
        # possible catch needed here
        
        # TODO, wee need to read to registers but must ensure we do not get a race condition
        # OWFS should do this for us if it is configured correctly - check this
        # Direction is returned on VAD in the range of 5% to 95% of VDD
        wdir_vdd = float(ow.owfs_get('uncached/26.BE50E7000000/VDD'))
        # possible catch needed here
        wdir_vad = float(ow.owfs_get('uncached/26.BE50E7000000/VAD'))
        # possible catch needed here
        
            
        #logdbg("1WIRE: Inspeed Vortex, Count=%d" % (wspeedCnt))
        # TODO Catch wrap round
        if self._wind_ts is None:
            self._wind_ts = time_ts
            logdbg("1WIRE: WindSpeed %s" % self._wind_ts)
            return None
        delta = time_ts - self._wind_ts
        if self._windCnt is None:
            self._windCnt = wspeedCnt
            return None
        cntdelta = wspeedCnt - self._windCnt
        speed = kph_per_hz * cntdelta/delta
        logdbg("1WIRE: WindSpeed tdelta %s Cntdelta %s Speed %s" % (delta, cntdelta, speed))
        
   
        # We need to generate moving averages here but rather than storing speeds for each time instant, if we store the count
        # then we can do the averaging with just 2 subtractions and a divide rather than a full blown sum all speeds and divide
        while self.gustpipe and self.gustpipe[0][self.TIME_INDEX] < time_ts - self.GUSTDURATION:
          _unused = self.gustpipe.popleft()
        self.gustpipe.append((time_ts, wspeedCnt))
        count = len(self.gustpipe)
        if count < 2:
            return None
        delta = self.gustpipe[count-1][self.TIME_INDEX] - self.gustpipe[0][self.TIME_INDEX]
        cntdelta = self.gustpipe[count-1][1] - self.gustpipe[0][1]
        gustspeed = kph_per_hz * cntdelta/delta
        #logdbg("1WIRE: WindSpeed GustSpeed Queue %s" % (self.gustpipe))
        logdbg("1WIRE: WindSpeed GustSpeed %s %s %s" % (gustspeed, sum(item[1] for item in self.gustpipe) / count, count))
        # NOTE if pipe is less than GUSTDURATION then out speeds are correct but avaraged over the wrong period
        
        # Average Speed
        while self.speedpipe and self.speedpipe[0][self.TIME_INDEX] < time_ts - self.MEANDURATION:
          _unused = self.speedpipe.popleft()
        self.speedpipe.append((time_ts, wspeedCnt))
        count = len(self.speedpipe)
        if count < 2:
            logdbg("1WIRE: WindSpeed Average Pipe is too short for calculation")
            return None
        delta = self.speedpipe[count-1][self.TIME_INDEX] - self.speedpipe[0][self.TIME_INDEX]
        cntdelta = self.speedpipe[count-1][1] - self.speedpipe[0][1]
        avspeed = kph_per_hz * cntdelta/delta
        # TODO possible wrap round to protect fore
        #logdbg("1WIRE: WindSpeed AvSpeed Queue %s" % (self.speedpipe))
        logdbg("1WIRE: WindSpeed AvSpeed %s %s %s %s" % (avspeed, self.speedpipe[0][1], self.speedpipe[count-1][1], count))
        # NOTE if pipe is less than GUSTDURATION then out speeds are correct but avaraged over the wrong period
        
        self._windCnt = wspeedCnt
        self._wind_ts = time_ts
        
        normalised = wdir_vad - 0.05*wdir_vdd # Normalise returned value from 5-95% to 0-90% of VDD
        if (normalised < 0):
            logerr("1WIRE: Inspeed e-Vane error, Normalised < 0, VAD=%d, VDD=%d" % (normalise, wdir_vad, wdir_vdd))
            direction = None
        direction = (400 * (normalised))/wdir_vdd
        if (direction > 360):
            logerr("1WIRE: Inspeed e-Vane error, Direction > 360, normalised=%d, VDD=%d" % (normalised, wdir_vdd))
            direction = None
            
        roundeddirection = int(direction + 0.5)
        gustdirection = roundeddirection
        
        logdbg("1WIRE: Rounding WindDir %s %s" % (roundeddirection, direction))
        # TODO Next bit is optional
        if (roundeddirection == 0):
            roundeddirection = 360
        if (avspeed ==0):
            roundeddirection = 0
        if (avspeed ==0):
            roundeddirection = None
        
        
        if (gustspeed  < self.MIN_GUST_SPEED_KPH):
            gustspeed = None
            gustdirection = None
        elif (gustspeed - avspeed < self.GUST_VARIATION_SPEED_KPH): #TODO Careful here what units are 5 in knot mph/kph ??
            gustspeed = None
            gustdirection = None
        
        #gustdirection = None
        #gustspeed = None
        # Put 16.093440006146921597227828997904 here for 10 mph
        #avspeed = 16.093440006146921597227828997904
        #avspeed = avspeed
        #roundeddirection = 47
        
        return (avspeed, gustspeed, roundeddirection, gustdirection) 

class Solar(object):
    def __init__(self):
        """Some text"""

    def value_at(self, time_ts):
        light = float(ow.owfs_get('26.C29621010000/S3-R1-A/illuminance'))
        # possible catch needed here
        return light * 0.4 * 0.7
        


class MaxSolar(object):
    def __init__(self):
        """Some text"""

    def value_at(self, time_ts):
        d = datetime.datetime.utcnow()
        # TODO need correct lat lon in decimal
        # they are in lat_f = stn_info.latitude_f
        latitude_deg = -43.2698333333 # positive in the northern hemisphere
        longitude_deg = 172.37608 # negative reckoning west from prime meridian in Greenwich, England
        altitude_deg = Pysolar.GetAltitude(latitude_deg, longitude_deg, d)
        azimuth_deg = Pysolar.GetAzimuth(latitude_deg, longitude_deg, d)
        radiation = Pysolar.radiation.GetRadiationDirect(d, altitude_deg)
        if (radiation < 0):
          logerr("1WIRE: Max Solar returned negative value")
          radiation = 0
        return radiation



class Observation(object):

    def __init__(self, magnitude=1.0, average=0.0, period=96.0, phase_lag=0.0, start=None):
        """Initialize an observation function.

        magnitude: The value at max. The range will be twice this value
        average: The average value, averaged over a full cycle.
        period: The cycle period in hours.
        phase_lag: The number of hours after the start time when the observation hits its max
        start: Time zero for the observation in unix epoch time."""

        if not start:
            raise ValueError("No start time specified")
        self.magnitude = magnitude
        self.average   = average
        self.period    = period * 3600.0
        self.phase_lag = phase_lag * 3600.0
        self.start     = start

    def value_at(self, time_ts):
        """Return the observation value at the given time.

        time_ts: The time in unix epoch time."""

        phase = 2.0 * math.pi * (time_ts - self.start - self.phase_lag) / self.period
        return self.magnitude * math.cos(phase) + self.average

if __name__ == "__main__":

    station = onewire(mode='simulator',loop_interval=5.0)
    for packet in station.genLoopPackets():
        print weeutil.weeutil.timestamp_to_string(packet['dateTime']), packet


