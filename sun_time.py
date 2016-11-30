#
#   The library uses sunrise-sunset.org for all time calculations.

# library that replaces argparse. Supposed to be easy to use.
import begin

# http and rest library.
import requests

# library to handle date and time
import arrow

import math


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


# Returns the Julian Day, which is day from the start of the Julian Calendar.
# The formula (based on UTC) is:
# Floor(365.25 *(Year + 4716)) + Floor(30.6001 * (Month + 1)) + Day + B -
# 1524.5
def getJD(date):
    # The date needs to be in UTC
    utcDate = date.to('UTC')

    print "date: %s and utcDate: %s" %(date, utcDate)
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


# astro_time calculates the solar position astronomically. It is taken from an
# algorithm first published here:
# http://www.nrel.gov/docs/fy08osti/34302.pdf
def astro_time(latS, lngS, date):

    # First we need the Julian Date.
    JD = getJD(date)

    lat = float(latS)
    lng = float(lngS)

    # We need to calculate the Mean Anomaly (difference in sun's motion from
    # a circle. We're going to use the one from AA:
    # http://aa.quae.nl/en/reken/zonpositie.html#mjx-eqn-eqm-aarde
    # It should be less than 1% off, so not a huge deal.
    M = (357.52 + 0.9856 * (JD - 2451545)) % 360

    # Need to find the true center vs. the mean center. According to
    # http://aa.quae.nl/en/reken/zonpositie.html#mjx-eqn-eqm-aarde again, we
    # need to do the following: Center = C1sinM + C2sin(2M) + C3sin(3M)
    # with C1 C2 and C3 being constants from a table. Also of note, python uses
    # radians for trigonometry functions.
    C = (1.9148 * sinDegree(M)) + (0.02 * sinDegree(2*M)) + \
        (0.0003 * sinDegree(3*M))

    # From the astronomy answers page above, the ecliptic longitude of the sun
    # is L = (M + 102.9373 + C + 180)%360. The modulus 360 is because it's in
    # degrees.
    L = (M + 282.9373 + C) % 360

    # Next, we need to find the right ascension and the declination.
    # Per AA, Ascension is approximately L + A2sin(2L) + A4sin(4L) + A6sin(6L)
    # Per AA, Declination is approximately D1sin(L) + D3sin^3(L) + D5sin^5(L)
    # A2, A4, A6 and D1, D3, D5 are constants.
    RA = L - (2.4657*sinDegree(2*L)) + (0.0529*sinDegree(4*L)) -\
        (0.0014*sinDegree(6*L))

    D = (22.7908*sinDegree(L)) + (0.5991 * math.pow(sinDegree(L), 3)) + \
        (0.0492*math.pow(sinDegree(L), 5))

    # Need to calculate the Mean Sidereal Time for the location, using the
    # following approximation from AA:
    # ST = (280.14 + 360.98*(int(JD)-2451545) - lat)mod 360
    ST = (280.14 + 360.98 * (int(JD)-2451545) - lat) % 360

    # The Hour Angle is the Sidereal Time - Right Ascension
    H = ST-RA

    # The Azimuth is arctan2(sin(H), cos(H)*sin(lng) - tan(D)*cos(lng))
    A = atan2Degree(sinDegree(H),
                    ((cosDegree(H)*sinDegree(lng)) -
                     (tanDegree(D) * cosDegree(lng))))

    # The altitude of the sun is:
    # h = arcsin(sin(lng)*sin(D)+cos(lng)*cos(D)*cos(H))
    h = asinDegree((sinDegree(lng) * sinDegree(D)) +
                   (cosDegree(lng) * cosDegree(D) * cosDegree(H)))

    # The sun transit time = (RA - D - ST)/360.0
    m = (RA-D-ST)/360.0

    # When from the AA, when the sun is at -.83 degrees, that is sunrise/set
    # The hour/angle is :
    # Ht = acos((sin(-.83) - sin(latitude)*sin(D)) / (cos(lat) * cos(D)))
    Ht = acosDegree((sinDegree(-.83) - (sinDegree(lat) * sinDegree(D))) /
                    (cosDegree(lat) * cosDegree(D)))

    # We need to calculate the Julian time for the transit near date.
    # This is done by determining nx:
    # nx = (JD-2451545-0.0009) - lng/360
    # Transit Time = JD+0.0009x(int(nx) - nx)
    nx = (JD - 2451545-0.0009) - lng/360.0
    print "JD = %f" % JD
    print  "nx = %f" % nx
    Jt = JD+(int(nx) - nx)

    # You can be more accurate by finding Ht' 2 more times by pluggin Ht into
    # this, but we don't need to.
    # Convert Ht to hours. For sunrise Ht = 360 - Ht, for sunset, Ht
    Jrise = Jt - (Ht/360)
    sunrise_time = JDtoDate(Jrise)
    Jset = Jt + (Ht/360)
    sunset_time = JDtoDate(Jset)
    print "Jtransit: %f" % Jt
    print sunrise_time
    print sunset_time
    return(sunrise_time, sunset_time)


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

    a_date = a_date.replace(hour=12, minute=0, second=0)

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
