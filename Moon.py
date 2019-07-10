##
## Find the percent of the moon that is illuminated based on the date.
##
## Input:
##	oDate - datetime object with date for illumination
##
## Output:
## 	single output which is a float containing the percentage of the moon illuminated
##
## Comments:
## 	Only accurate for the date (and not time) and not useful to find the exact time the moon is full, or new. Also only give percentage illuminated, and not whether it is waxing, waning etc-

import math
import datetime
#from DateToDayNumber import DateToDayNumber
from LeapYear import LeapYear
# calculate the phase of the moon
# method from "Practical Astronomy with your Calculator or Spreadsheet"
# Epoch is Jan 0 2010
dt = datetime

class Moon:
	def __init__(self, oDate):
		self.oDate = oDate


	def MoonPercent (self):
		#print("--------- start ---------")
		# define some constants
		epsilon_g = 279.557208 # ecliptic longitude at epoch 2010.0
		omegabar_g = 283.112438 # ecliptic longitude of pirigee at epoch 2010.0
		e = 0.016705 # eccentricity of orbit at epoch 2010.0 (sun)
		r_0 = 1.495985 *(10**8) # semi-major axis
		
		l_0 = 91.929336 # moons mean longitude at the epoch
		P_0 = 130.14076 # mean longitude of the perigee at the epoch
		N_0 = 291.682547 # mean longitude of the node at the epoch
		i_ = 5.145396 # inclination of Moons orbit at the epoch
		e_ = 0.0549 # eccentricity of the Moon's orbit at epoch
		a = 384401 # semi-major axis of Moons orbit at epoch
		theta_0 = 0.5181 # Moon's angular diameter at distance a from the earth at epoch
		pi_0 = 0.9507 # Moon's paralax at distance a from the earth
	
		epochyr = 2010
		#print(r_0)
		
		# find number of days elapsed since our epoch
		DayNumber = Moon.DateToDayNumber(self.oDate)
		#print("DayNumber:",DayNumber)
		#print("self.oDate.year: ",self.oDate.year)
		# add 365 days for every year since 2010 plus 1 extra day for each leapyear
		ExtraDays = (365 * abs(self.oDate.year - epochyr))
		#print(DayNumber)
		NoLeapYears=0
		# correct for the leapyears
		
		#sort out range from and to
		if self.oDate.year < epochyr +1:
			rangestart = self.oDate.year
			rangeend = epochyr +1
		else:
			rangestart = epochyr
			rangeend = self.oDate.year +1
		
		for i in range(rangestart,rangeend):
			#print("Year:",self.oDate.year)
			NoLeapYears += LeapYear(dt.date(i,1,1))
			#print(LeapYear(dt.date(i,1,1)))
		#print(ExtraDays,DayNumber,NoLeapYears)
		# calculate D
		if self.oDate.year < epochyr:
			D = DayNumber - (ExtraDays + NoLeapYears)
		else:
			D = DayNumber + (ExtraDays + NoLeapYears)
		#print("D: ",D)
	
		# calculate N = 360D/365.242191 and Make N(0,360]
		N = ((360 * D)/365.242191) % 360
		#print("N: ",N)
	
		# Find M_dot = N + epsilon-g − omegabar_g and add  360 if answer is negative
		M_dot = N + epsilon_g - omegabar_g
		if M_dot < 0:
			#print("M_dot befor adding 360: ",M_dot)
			M_dot +=360
		#print("M_dot: ",M_dot)
	
		# Find E_c = 360/pi e Sin(M_dot)
		E_c = 360/math.pi * e * math.sin(math.radians(M_dot))
		#print("E_c: ",E_c)
	
		# Find lamda_dot = N + E_c + epsilon_g and subtract 360 if result is greater than 360
		lamda_dot = N + E_c + epsilon_g
		if lamda_dot > 360:
			lamda_dot -=360
		#print("lamda_dot: ", lamda_dot)
	
		# we have now found the Sun's ecliptic longitude (lamda_dot) and mean anomaly (M_dot)
		# Find l = 13.176396D + l_0 adjust to (0,360]
		l = ((13.176396 * D) + l_0) % 360
		#print("l: ", l)
	
		# find M_m = l - 0.1114041D - P_0 adjust to (0,360]
		M_m = (l - (0.1114041 * D) - P_0) % 360
		#print("M_m: ", M_m)
		
		# Find N_ = N_0 -0.0529539 * D adjust to (0,360]
		N_ = (N_0 - (0.0529539 * D)) % 360
		#print("N_: ", N_)
	
		# Find E_v = 1.2739 Sin(2C-M_m); C = l - lamda_dot
		C = l - lamda_dot
		E_v = 1.2739 * math.sin(math.radians((2 * C) - M_m))
		#print("E_v: ",E_v)
	
		# Find A_e = 0.1858 Sin(M_dot)
		A_e = 0.1858 * math.sin(math.radians(M_dot))
		#print("A_e: ", A_e)	
	
		# Find A_3 = 0.37 Sin(M_dot)
		A_3 = 0.37 * math.sin(math.radians(M_dot))
		#print("A_3: ",A_3)
		
		# find the corrected anomaly Mdash_m = M_m + E_v - A_e - A_3
		Mdash_m = M_m + E_v - A_e - A_3
		#print("Mdash_m: ",Mdash_m)
	
		# Calculate E_c = 6.2886 Sin(Mdash_m)
		E_c = 6.2886 * math.sin(math.radians(Mdash_m))
		#print("E_c: ",E_c)
	
		# Calculate A_4 = 0.214 Sin(2 * Mdash_m)
		A_4 = 0.214 * math.sin(math.radians(2 * Mdash_m))
		#print("A_4: ",A_4)
	
		# calculate ldash = l + E_v + E_c - A_e + A_4
		ldash = l + E_v + E_c - A_e + A_4
		#print("ldash: ",ldash)	
	
		# calculate V = 0.6583 Sin [2(ldash - lamda_dot)]
		V = 0.6583 * math.sin(math.radians(2 * (ldash - lamda_dot)))
		#print("V: ",V)
	
		# find ldbldash = ldash + V
		ldbldash = ldash + V
		#print("ldbldash: ", ldbldash)
	
		# phase of moon is therefore F = 1/2(1-cos ((ldbldash - lamda_dot))
		F = 0.5 * (1 - math.cos(math.radians(ldbldash - lamda_dot)))
		#print("F: ", F)
	
	


		#print("---------- end -----------")
		return F
	
	##
	## get the phase of the moon such as waxing gibbous, new, full etc
	##
	## depends:
	## 	depends on MoonPercent
	## input:
	## 	oDate - Datetime object containing today's date
	## output:
	##	tuple containing Phase Name, Percent, Age
	##

	def MoonPhase(self):
		
		# see whether moon percent is increasing, decreasing or point of inflection. Get moon illumination for previous day [0], date [1] and following day [1]
		
		t = [0,0,0]
		#modify oDate property so that it looks one day before (t[0]), today (t[1]), and one day afterwards t[2].
		t[1] = Moon.MoonPercent(self)
		self.oDate=self.oDate+dt.timedelta(days=-1)	
		t[0] = Moon.MoonPercent(self)
		self.oDate=self.oDate+dt.timedelta(days=+2)	
		t[2] = Moon.MoonPercent(self)
		# reset self.oDate back to normal
		self.oDate=self.oDate+dt.timedelta(days=-2)	
		#print(t)
	
		if (t[0] <= t[1] <= t[2]):
			#print("increasing")
			# moon is waxing / getting brighter
			# check whether it is quarter, crescent or gibbous
			if 0 < t[1] <= 0.44:
				# up to a quarter then waxing crescent
				Phase = "Waxing Crescent"
			elif 0.44 < t[1] <= 0.56:
				# "exactly" a quarter - 1st quarter
				Phase = "First Quarter"
			elif 0.56 < t[1] < 1:
				# over half bright, but less than full - waxing gibbous
				Phase = "Waxing Gibbous"
	
		elif t[0] >= t[1] >= t[2]:
			#print("decreasing")
			# moon is waning / getting darker
			# check whether it is quarter, crescent or gibbous
			if 0 < t[1] <= 0.44:
				# up to a quarter then waxing crescent
				Phase = "Waning Crescent"
			elif 0.44 < t[1] <= 0.56:
				# "exactly" a quarter - 1st quarter
				Phase = "Last Quarter"
			elif 0.56 < t[1] < 1:
				# over half bright, but less than full - waxing gibbous
				Phase = "Waning Gibbous"
		
		elif t[0] >= t[1] <= t[2]:
			#print("lowest")
			# moon illumination is at the darkest therefore must be a new moon
			Phase = "New"
		elif t[0] <= t[1] >= t[2]:
			#print("highest")
			# moon illumination is at the brightest, therefore must be a full moon
			Phase = "Full"
	
		#print("Phase",Phase)
		#return {'Phase':Phase,'MoonPercent':MoonPercent(self.oDate)}
		return Phase

		####################################
		# DateToDayNumber
		# expects dt object being passed to it

	def DateToDayNumber (ooDate):
		#print("--------- start DateToDayNumber ---------")
		
		if ooDate.month <= 2:
			M = ooDate.month - 1
			M = M*(62 + LeapYear(ooDate))
			M = int(M/2)
	
		else:
			M = ooDate.month + 1
			M = int(M*30.6)
			M = M - (63 - LeapYear(ooDate))
	
		D = M + ooDate.day


		#print("---------- end DateToDayNumber -----------")

		return D

		
