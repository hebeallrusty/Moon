# Moon Calculator

This project calculates the Moon Rise and Set times using the method in Meeus Astronomical Algorithms.

The time is iterated 3 times to increase the accuracy of the result

Times are returned in respect of UT

Usage:
MoonTime(Year, Month, Day, Latitude, Longitude, Event)

Year should be the full year (i.e. 2023 not 23)
Event should be one of "Rise", "Set" or "Transit", anything else will halt execution.

Note that both Latitude and Longitude are to be entered as a decimal.
Longitude should be entered as positive (in the West) and negative (in the East)

Output is a tuple containing the Hour, Minute and Second of the Event, or False if it doesn't occur'

