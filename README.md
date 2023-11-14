# Moon Calculator

This project calculates the Moon Rise and Set times and Moon Phase using the method in Meeus Astronomical Algorithms.

The Moon Times are iterated 5 times to increase the accuracy of the result, although this could probably be reduced to 3

Times are returned in respect of UT

## Moon Rise, Set and Transit
Usage:
MoonTime(Year, Month, Day, Latitude, Longitude, Event)

Year should be the full year (i.e. 2023 not 23)
Event should be one of "Rise", "Set" or "Transit", anything else will halt execution.

Note that both Latitude and Longitude are to be entered as a decimal.
Longitude should be entered as **positive in the West and negative in the East**

Output is a tuple containing the Hour, Minute and Second of the Event, or False if it doesn't occur'

## Moon Phase
Phase(Year)

Year should be the full year + a fractional part for the month.

Note that the fractional part is centric around a New Moon, but, the fraction may produce an event that happens in the previous month, or following month. Care should be taken when an event happens twice in a month (such as a Blue Moon)

Output is a dictionary containing all events for that fractional part of the month, with keys *Full Moon*, *New Moon*, *First Quarter* and *Last Quarter*

