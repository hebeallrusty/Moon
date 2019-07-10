#import math
import datetime

# calculate if it's a leapyear
# if a year is div by 4 then check if it is also div by 100 and 400

def LeapYear (oDate):
	#print("--------- start ---------")
	# define some constants
	if (oDate.year % 4) == 0: # if no remainder after div by 4
		#print("div by 4")
		# check if the year is a century if it is then it also needs to be divisble by 400. If it isn't a century then it's a leap year
		if (oDate.year % 100) == 0: # if no remainder after div by 100
			#print("div by 100")
			if (oDate.year % 400) == 0: # if no remainder after div by 400
				#print("div by 400")
				# year is century and divisible by 400 and therefore a leapyear
				return 1
			else:
				# year is a century and not div by 400
				return 0
		else: # year is not a century therefore a leapyear
			return 1
			
	#print("---------- end -----------")
	# year isn't div by 4 and not a leapyear
	return 0

