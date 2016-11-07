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


@begin.start(auto_convert=True)
def main(lat='40.7128',
         lng='-74.0059',
         date=None):
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
    payload = {'lat': lat,
               'lng': lng,
               'date': a_date.format('YYYY-MM-DD'),
               'formatted': '0'}

    print payload
    resp = requests.get('http://api.sunrise-sunset.org/json?', params=payload)
    if resp.status_code == requests.codes.ok:
        parsed_JSON = resp.json()
        if parsed_JSON['status'] == 'OK':
            results = parsed_JSON['results']
            sunrise_time = arrow.get(results['sunrise']).to('local')
            sunset_time = arrow.get(results['sunset']).to('local')
            print "Sunrise: %s" % (sunrise_time.time())
            print "Sunset: %s" % (sunset_time.time())
            for key in results.keys():
                print "%s: %s" % (key, results[key])

    else:
        resp.raise_for_status()
