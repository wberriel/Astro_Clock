# Script Name   : sun_time
# Author        : William Berriel
# Created       : Nov 6 2016
# Last Modified :
# Version       : 0.01

# Modifications :

# Description   : This script will tell you the time for sunrise/sunset based
#   based on a logitude/latitude and a day. Defaults to today. It also is a
#   test platform for the begins library.
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


def sinDegree(x):
    return math.sin(math.radians(x))


def tanDegree(x):
    return math.tan(math.radians(x))


def atanDegree(x):
    # arc tan takes -infinity to infinity and returns -pi to pi
    return math.degrees(math.atan(x))


def json_time(lat, lng, date, timezone):

    payload = {'lat': lat,
               'lng': lng,
               'date': date.format('YYYY-MM-DD'),
               'formatted': '0'}

    resp = requests.get('http://api.sunrise-sunset.org/json?', params=payload)
    if resp.status_code == requests.codes.ok:
        parsed_JSON = resp.json()

        if parsed_JSON['status'] == 'OK':
            results = parsed_JSON['results']
            sunrise_time = arrow.get(results['sunrise']).to(timezone)

            sunset_time = arrow.get(results['sunset']).to(timezone)

            return (sunrise_time, sunset_time)
    else:
        resp.raise_for_status()


# astro_time calculates the solar position astronomically. It is taken from an
# algorithm first published in Almanac for Computers, 1990. It's found here:
# http://williams.best.vwh.net/sunrise_sunset_algorithm.htm.
def astro_time(lat, lng, date):

    # The first step is convert date to day of the year.
    # In arrow format "DDDD" is the day of the year.
    # Convert this to an integer and we have the numerical day of the year.
    day = int(date.format("DDDD"))

    # The second step is to convert longitude to an hour value.
    # This is longitude /15 = longitude hours. This is because 360/24 = 15.
    lng_Hour = lng / 15

    t_sunRise = day + ((6-lng_Hour)/24)
    t_sunset = day + ((18-lng_Hour)/24)

    # We need to calculate the Mean Anomaly (difference in sun's motion from
    # a circle. We're going to use the one from AA:
    # http://aa.quae.nl/en/reken/zonpositie.html#mjx-eqn-eqm-aarde
    # It should be less than 1% off, so not a huge deal.
    M = (0.9856 * day - 3.59) % 360

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

    D = (22.7908*sinDegree(x)) + (0.5991 * math.pow(sinDegree(x), 3)) + \
        (0.0492*math.pow(sinDegree(x), 5))


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
        a_date = arrow.now()
    if not lat_validate(lat):
        print "%s is not a valid latitude." % (lat)
        print "Please only use a floating point number betwen -90 and 90."
        return

    if not lng_validate(lng):
        print "%s is not a valid longitude." % (lng)
        print "Please only use a floating point number betwen -180 and 180."
        return
    try:
        arrow.get(timezone)
    except arrow.parser.ParserError:
        print "%s is not a valid timezone." % (timezone)
        print "Please use a valid timezone similar to \'US/Pacific\' or omit", \
            " for local timezone."
        return


