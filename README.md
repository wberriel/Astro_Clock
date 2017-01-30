# sun_time
An implentation of the NASA Solar Positioning Algorithm written in Python. It is based on the algorithm from  http://www.nrel.gov/docs/fy08osti/34302.pdf. 

Currently it calculates the time of sunrise and sunset. 

Astro_Clock requires the following libraries: Begin, Requests, Arrow.

The system also pulls data through JSON for verification. It uses http://api.sunrise-sunset.org for verification. Currently it's within 0.1% of that reference. 

Please run sun_time --help to see command line usage.
