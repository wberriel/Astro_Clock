#
#   The library uses sunrise-sunset.org for all time calculations.

# library that replaces argparse. Supposed to be easy to use.
import begin

# http and rest library.
import requests

# library to handle date and time
import arrow

import math

import sun_data


def lat_validate(input):
    # this will throw a type error if input is not correct.
    try:
        lat = float(input)

    # If it's not a valid floating point.
    except:
        return 0

    if abs(lat) <= 90:
        return 1
    else:
        return 0


def lng_validate(input):
    # this will throw a type error if input is not correct.
    try:
        lng = float(input)
    # if it's not a valid floating point number.
    except:
        return 0

    if abs(lng) <= 180:
        return 1
    else:
        return 0


def cosDegree(x):
    return math.cos(math.radians(x))


def acosDegree(x):
    # acos expects a value between -1 and 1 and returns it in radians
    return math.degrees(math.acos(x))


def sinDegree(x):
    return math.sin(math.radians(x))


def asinDegree(x):
    # asin expects a value between -1 and 1 and returns it in radians
    return math.degrees(math.asin(x))


def tanDegree(x):
    return math.tan(math.radians(x))


def atanDegree(x):
    # arc tan takes -infinity to infinity and returns -pi to pi
    return math.degrees(math.atan(x))


def atan2Degree(x, y):
    # atan2 uses a vector, thus 2 inputs
    return math.degrees(math.atan2(x, y))


def json_time(lat, lng, date):

    payload = {'lat': lat,
               'lng': lng,
               'date': date.format('YYYY-MM-DD'),
               'formatted': '0'}

    resp = requests.get('http://api.sunrise-sunset.org/json?', params=payload)
    if resp.status_code == requests.codes.ok:
        parsed_JSON = resp.json()

        if parsed_JSON['status'] == 'OK':
            results = parsed_JSON['results']
            sunrise_time = arrow.get(results['sunrise'],)

            sunset_time = arrow.get(results['sunset'],)

            return (sunrise_time, sunset_time)
    else:
        resp.raise_for_status()


def earth_periodic_term_sum(TERMS, jme):
    # This will sum the terms arrays using the following formula:
    # For = sum(map(A*cos(B+C*jme))) where Lx are primary rows and A,B,C are
    # parts of the tuple in each array.
    res_array = []
    for term_array in TERMS:
        # term_array should be an array of tuples with 3 values.
        # Values in L_TERMS work with radians.
        res_array.append(sum(map(lambda x: x[0] * math.cos(x[1] + (x[2] * jme)),
                                 term_array)))

    result = 0.0
    # we should now have an array of the L_Terms summed. Now we need the formula
    # L = (L0 + L1*JME + L2 *JME^2 + L3 *JME^3 + L4*JME^4 + L5*JME^5) / 10^8
    for i in range(len(res_array)):
        result += (res_array[i] * math.pow(jme, i))/math.pow(10.0, 8)

    return result


def xy_term_summation(X, terms):
    result = 0
    for j in range(len(X)):
        result += X[j] * terms[j]

    return result


def nutation_longitude_and_obliquity(jce):
    X = [0.0]*5
    # Mean elongation of the moon from the sun in degrees
    X[0] = 297.85036 + (445267.111480 * jce) - \
        (0.0019142 * math.pow(jce, 2)) + (math.pow(jce, 3) / 189474.0)

    # Mean anomaly of the sun(earth) in degrees
    X[1] = 367.52772 + (35999.050340 * jce) - \
        (0.001603 * math.pow(jce, 2)) - (math.pow(jce, 3) / 300000.0)

    # mean anomaly of the moon in degrees
    X[2] = 134.96298 + (477198.867398 * jce) - \
        (0.0086972 * math.pow(jce, 2)) - (math.pow(jce, 3) / 56250.0)

    # moon's argument of latitude in degrees
    X[3] = 93.27191 + (483202.017538 * jce) - \
        (0.0036825 * math.pow(jce, 2)) - (math.pow(jce, 3) / 327270.0)

    X[4] = 125.04452 - (1934.136261 * jce) - \
        (0.0020708 * math.pow(jce, 2)) - (math.pow(jce, 3) / 450000.0)

    nut_lng = 0
    obliquity = 0

    for i in range(len(sun_data.Y_TERMS)):
        xy_term_sum = xy_term_summation(X, sun_data.Y_TERMS[i])
        pe_term = sun_data.PE_TERMS[i]
        nut_lng += (pe_term[0] + (pe_term[1] * jce)) * sinDegree(xy_term_sum)
        obliquity += (pe_term[2] + (pe_term[3] * jce)) * cosDegree(xy_term_sum)

    # the values are in arc seconds, need to convert to degrees
    return (nut_lng/360000000, obliquity/36000000)


def ecliptic_mean_obliquity(jme):
    u = jme/10.0

    terms = (84381.448, - 4580.93, -1.55, 1999.25, -51.38, -249.67, -39.05,
             7.12, 27.87, 5.79, 2.45)
    result = 0
    for i in range(len(terms)):
        result += terms[i] * math.pow(u, i)
    return result


def mean_sidreal_time(jd, jc):
    return 280.46061837 + (360.98564736629 * (jd - 2451545)) + \
        (0.000387933 * math.pow(jc, 2)) - (math.pow(jc, 3) / 38710000)


# This returns the value in radians. perhaps it can be in degrees?
def sun_right_ascension(sun_lng, obliquity, geoB):
    return math.atan2((sinDegree(sun_lng) * cosDegree(obliquity) -
                       (tanDegree(geoB) * sinDegree(obliquity))),
                      cosDegree(sun_lng))


# The declination will be returned in radians. Can I refactor to be degrees?
def sun_declination(sun_lng, obliquity, geoB):
    return math.asin((sinDegree(geoB) * cosDegree(obliquity)) +
                     (sinDegree(obliquity) * sinDegree(sun_lng)))


# Returns the Julian Day, which is day from the start of the Julian Calendar.
# The formula (based on UTC) is:
# Floor(365.25 *(Year + 4716)) + Floor(30.6001 * (Month + 1)) + Day + B -
# 1524.5
def getJD(date):
    # The date needs to be in UTC
    utcDate = date.to('UTC')

    print "date: %s and utcDate: %s" % (date, utcDate)
    year = utcDate.year
    month = utcDate.month
    day = float(utcDate.day)
    hour = utcDate.hour
    minute = utcDate.minute

    # if the current month is less than or equal to 2, we need to decrement the
    # year and add 12 to the month.
    if month <= 2:
        year -= 1
        month += 12

    # adjust the day by the hour and minute as fractions of a day.
    day = day + hour/24.0 + minute/1440.0

    # leap year adjustments:
    B = 2 - int(year/100) + int(year/400)

    JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + \
        B - 1524.5

    return JD


def JDtoDate(JD):
    JD = JD+.5
    Z = int(JD)
    F = JD-Z

    alpha = int((Z-1867216.25)/36524.25)
    A = Z+1+alpha - int(alpha/4)

    B = A+1524
    C = int((B-122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)

    day = B - D - int(30.6001 * E)
    hour = int(F * 24)
    minute = int((F - (hour / 24.0)) * 1440)
    second = int((F - (hour / 24.0) - (minute / 1440.0)) * 86400)

    if(E > 13.5):
        month = int(E - 13)
    else:
        month = int(E - 1)

    if(month > 2):
        year = int(C-4716)
    else:
        year = int(C-4715)

    dateString = "%04d-%02d-%02dT%02d:%02d:%02d-00:00" % (year, month, day,
                                                          hour, minute, second)
    return arrow.get(dateString)


# http://www.nrel.gov/docs/fy08osti/34302.pdf
def astro_time(lat, lng, date):
    # Convert date to Julian Day
    # JDE is ephemeris, apparantly Terrestrial time differs from UTC by about 67
    # seconds.
    jd = getJD(date)
    jde = jd + 67/86400.0
    # Julian Century and millenium too:
    jc = (jd - 2451545)/36525.0
    jce = (jde - 2451545)/36525.0
    jme = jce/10.0

    # convert latitude and longitude to floats
    geo_lat = float(lat)
    geo_lng = float(lng)

    # L and B are in radians, and need to be converted to degrees
    L = earth_periodic_term_sum(sun_data.L_TERMS, jme)
    L = math.degrees(L) % 360
    B = earth_periodic_term_sum(sun_data.B_TERMS, jme)
    B = math.degrees(B) % 360
    R = earth_periodic_term_sum(sun_data.R_TERMS, jme)

    print "L = %f B = %f R = %f" % (L, B, R)

    # need to calculate the geocentric Longitude (geoL) and latitude (geoB)
    geoL = (L + 180) % 360
    geoB = -B
    print "geo L: %f geo B: %f" % (geoL, geoB)

    # now we need the nuttation in longitude and obliquity
    (nut_lng, nut_obliquity) = nutation_longitude_and_obliquity(jce)
    print "Nuttation of Longitude: %f, Obliquity %f" % (nut_lng, nut_obliquity)

    # the ecliptic obliquity is the mean obliquity/3600 + the nutation of the
    # obliquity. in degrees
    obliquity = ecliptic_mean_obliquity(jme)/3600 + nut_obliquity
    print "Ecliptic True Obliquity: %f" % obliquity

    # The aberration Correction = -20.4898/(3600*R)
    aberration_correction = -20.4898/(3600*R)

    sun_lng = geoL + nut_lng + aberration_correction

    print "Apparent Sun Longitude: %f" % sun_lng

    # Apparent sidreal time at greenwich is = mean sidreal time + nuttation of
    # longitude * cos(Ecliptic True Obliquity) limited to 360
    # Mean sidreal time is in degrees.
    sidreal_time = mean_sidreal_time(jd, jc) + nut_lng * cosDegree(obliquity)
    sidreal_time = sidreal_time % 360
    print "Apparent Sidereal Time: %f" % sidreal_time

    # sun right ascension will be in radians, convert to degrees and ensure it's
    # within 0-360
    r_asc = sun_right_ascension(sun_lng, obliquity, geoB)
    r_asc = math.degrees(r_asc) % 360

    print "Right Ascension: %f" % r_asc

    # sun declination will be in radians, convert to degrees
    sun_dec = sun_declination(sun_lng, obliquity, geoB)
    sun_dec = math.degrees(sun_dec) % 360
    print "Sun Declination: %f" % sun_dec

    # Local hour angle shoudl be in degrees, limited to 0-360
    H = (sidreal_time + geo_lng - r_asc) % 360

    print "Local Hour Angle: %f" % H


@begin.start(auto_convert=True)
def main(lat='40.7128',
         lng='-74.0059',
         date=None,
         printall=False,
         timezone='local'):
    "Returns the sunrise and/or sunset time for a given day."

    if date:
        try:
            a_date = arrow.get(date)
        except arrow.parser.ParserError:
            print "%s is not a valid date." % (date)
            print "Please provide a valid date of the form YYYY-MM-DD"
            return

    else:
        a_date = arrow.utcnow()
    if not lat_validate(lat):
        print "%s is not a valid latitude." % (lat)
        print "Please only use a floating point number betwen -90 and 90."
        return

    if not lng_validate(lng):
        print "%s is not a valid longitude." % (lng)
        print "Please only use a floating point number betwen -180 and 180."
        return
    try:
        arrow.now(timezone)
    except arrow.parser.ParserError:
        print "%s is not a valid timezone." % (timezone)
        print "Please use a valid timezone similar to \'US/Pacific\' or omit", \
            " for local timezone."
        return

    print "Astro Time"
    (sunrise_time, sunset_time) = astro_time(lat, lng, a_date)
    print "Sunrise: %s" % sunrise_time.to('local').format("HH:mm:SS")
    print "Sunset: %s" % sunset_time.to('local').format("HH:mm:SS")

    print "JSON Time"
    (sunrise_time, sunset_time) = json_time(lat, lng, a_date)
    sunrise_time.to('local')
    sunset_time.to('local')
    print "Sunrise: %s" % sunrise_time.to('local').format("HH:mm:SS")
    print "Sunset: %s" % sunset_time.to('local').format("HH:mm:SS")
    print "date: %s" % a_date
    print "JD : %f" % getJD(a_date)
    print JDtoDate(2453097)
