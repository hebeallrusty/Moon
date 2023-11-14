# Moon Calculations
#
# Calculate rise and set times of the moon based on location.

# The heavy lifting is done by the EstimateMoon function, however MoonRiseSet should be the one invoked as this iterates the result to improve accurancy.
# Moon Rise and Set cannot be calculated without iterating due to the speed at which it moves in the sky.
#
# The method of calculation is based on those by Jean Meeus in Astronomical Algorithms and relevant chapters / equation locations are noted
#
# Moon Rise, Set or Transit can be obtained by calling MoonTime(YEAR,MONTH,DAY,Latitude,Longitude,Event) where event is either "Rise","Set" or "Transit" and Lat / Lon is decimal degrees


import math

def DegtoHMS(angle):
    HR = math.floor(angle / 15)
    MIN = math.floor(((angle / 15) - HR) * 60)
    SEC = ((((angle / 15) - HR) * 60) - MIN) * 60 
    
    return((HR,MIN,SEC))

def DegtoDMS(angle):
    
    if angle < 0:
        sign = (-1)
        angle = abs(angle)
       
    else:
        sign = 1
    
    DEG = math.floor(angle)
    MIN = math.floor(((angle - DEG) * 60))
    SEC = (((angle - DEG) * 60) - MIN) * 60
    
    return((DEG * sign,MIN,SEC))

def Hrs(DecimalHrs):
    HRS = math.floor(DecimalHrs)
    MIN = math.floor((DecimalHrs - HRS) * 60)
    SEC = math.floor(((((DecimalHrs - HRS) * 60)) - MIN) * 60)
    
    return((HRS, MIN,SEC))

def JulianDay(YEAR,MONTH,DAY):
    # calculates the Julian Day in respect of JD2000
    
    if MONTH <= 2:
        YEAR -= 1
        MONTH += 12
    
    A = math.floor(YEAR / 100)
    
    if YEAR >= 1582: # the year 
        B = 2 - A + math.floor(A / 4)
    else:
        B = 0
        
    return(math.floor(365.25 * (YEAR + 4716)) + math.floor(30.6001 * (MONTH + 1)) + DAY + B - 1524.5)

def CalendarDate(JDE):
    # from Meeus ch 7
    
    JD = JDE + 0.5
    
    Z = math.floor(JD)
    F = JD - Z
    
    #print("Z",Z)
    
    if Z < 2299161:
        A = Z
    else:
        alpha = math.floor((Z - 1867216.25)/36524.25)
        A = Z + 1 + alpha - math.floor(alpha / 4)
        #print("alpha",alpha)
        #print("A",A)
        
    B = A + 1524
    #print("B",B)
    C = math.floor((B - 122.1)/365.25)
    #print("C",C)
    D = math.floor(365.25 * C)
    #print("D",D)
    E = math.floor((B - D) / 30.6001)
    #print("E",E)
    
    DAY = B - D - math.floor(30.6001 * E)
    
    if E < 14:
        MONTH = E - 1
    else:
        MONTH = E - 13
        
    if MONTH > 2:
        YEAR = C - 4716
    else:
        YEAR = C - 4715
        
    TIME = Hrs(F * 24)
    
    return((YEAR,MONTH,DAY,TIME[0],TIME[1],TIME[2]))
    


def CalculateT(YEAR,MONTH,DAY):
    # calculates time T from Epoch JD2000 (JDE 2451545.0)
    
    return((JulianDay(YEAR,MONTH,DAY) - 2451545)/36525)

def Theta0(YEAR,MONTH,DAY):
    # aparent sidereal time at 0h UT on day D (from Meeus CH 12)
    # only valid for 0h UT of a date
    UT = (JulianDay(YEAR,MONTH,math.floor(DAY)) - 2451545) / 36525
    #print(UT)
    THETA0 = 100.46061837 + (36000.770053608 * UT) + (0.000387933 * UT * UT) - (UT * UT * UT / 38710000)
    
    return(THETA0)
    
def lrbMatrixSum(matrix,ecc,coeff,trig,D,M,Mdash,F):
    # used for summing up the terms involving l, r or b;
    # matrix is the list of coefficients and arguments
    # ecc is E - the eccentricity of the earths orbit (47.6 in Meeus)
    # coeff is the position of the coefficient for l, r or b, in the matrix (it'll be either 4 or 5 depending on which matrix is used)
    # trig is the function to perform - sin or cos
    # D, M, Mdash and F should be self explanitory and be previously calculated
    
    # positions in the matrix where D, M, Mdash and F can all be found
    Dpos = 0
    Mpos = 1
    Mdashpos = 2
    Fpos = 3
    
        # set up variable for summation
    s = 0
    
    # set up a generic function reference so that all sums can be calculated efficiently without mostly duplicating a def function
    if trig == "cos":
        trig = math.cos
    elif trig == "sin":
        trig = math.sin
    else:
        raise SystemExit('lrbMatrixSum requires trig argument to be either "cos" or "sin"') 
    
    # iterate through whole matrix adding up all the arguments and coefficients
    for row in matrix:
        #print(matrix[i])
        #print(row[Mpos])
        
        # as note above 47.6 in Meeus - if coeff of M is 1 or -1 then multiply coeff by E; if coeff of M is 2 or -2 then multiply coeff by (E*E)
        if abs(row[Mpos]) == 1:
            E = ecc
        elif abs(row[Mpos]) == 2:
            E = ecc * ecc
        else:
            E = 1
        #print(E)
        
        s += row[coeff] * E * trig(math.radians((D * row[Dpos]) + (M * row[Mpos]) + (Mdash * row[Mdashpos]) + (F * row[Fpos]))) 
    
    #print(s)
    return(s)

def psiepsilonMatrixSum(matrix,coeff,trig,D,M,Mdash,F,Ohmega):
    # used for summing up the terms involving delta-psi and delta-epsilon for;
    # matrix is the list of coefficients and arguments
    # coeff is the position of the coefficient for delta-psi or delta-epsilon, in the matrix (it'll be either 5 or 6 depending on which matrix is used)
    # trig is the function to perform - sin or cos
    # D, M, Mdash, F and Ohmega should be self explanitory and be previously calculated
    
    # positions in the matrix where D, M, Mdash and F can all be found
    Dpos = 0
    Mpos = 1
    Mdashpos = 2
    Fpos = 3
    Ohmegapos = 4
    
    # set up variable for sumation
    s = 0
    
    # set up a reference to the trig function to make function work for both delta-psi and delta-epsilon
    if trig == "cos":
        trig = math.cos
    elif trig == "sin":
        trig = math.sin
    else:
        raise SystemExit('lrbMatrixSum requires trig argument to be either "cos" or "sin"')
    
    # iterate through matrix adding all ters together
    for row in matrix:
        s += row[coeff] * trig(math.radians((D * row[Dpos]) + (M * row[Mpos]) + (Mdash * row[Mdashpos]) + (F * row[Fpos]) + (Ohmega * row[Ohmegapos])))
        
    return(s)

def PhaseCorrectionMatrixSum(matrix,M,Mdash,F,Ohmega):
    # calculate the sum of the corrections for Moon Phase
    Mpos = 0
    Mdashpos = 1
    Fpos = 2
    Ohmpos = 3
        
    # create list for the calculation in form [new moon, full moon, quarters]
    s = [0,0,0,0]
    
    #print("a:",a)
      
    for i in range(len(s)):
        for row in matrix:
            # sin(combination of M, Mdash, F and Ohmega) * moon phase coefficient 
            s[i] += math.sin(math.radians((M[i] * row[Mpos]) + (Mdash[i] * row[Mdashpos]) + (F[i] * row[Fpos]) + (Ohmega[i] * row[Ohmpos]))) * row[i+4]

    return(s)
        


def CalculateLdash(T):
    # Moon's mean longitude (Mean equinox of the date)
    return((218.3164477 + (481267.88123421 * T) - (0.0015786 * T * T) + (T * T * T / 538841) - (T * T * T * T / 65194000)) % 360)


def CalculateD(T):
    # Mean elongation of the Moon
    return((297.8501921 + (445267.1114034 * T) - (0.0018819 * T * T) + (T * T * T / 545868) - (T * T * T * T / 113065000)) % 360)

def CalculateM(T):
    # Sun's mean anomaly
    return((357.5291092 + (35999.0502909 * T) - (0.0001536 * T * T) + (T * T * T / 24490000)) % 360)

def CalculateMdash(T):
    # Moon's mean anomaly
    return((134.9633964 + (477198.8675055 * T) + (0.0087414 * T * T) + (T * T * T / 69699) - (T * T * T * T / 14712000)) % 360)

def CalculateF(T):
    # Moon's argument of latitude (mean distance of the Moon from its ascending node)
    return((93.2720950 + (483202.0175233 * T) - (0.0036539 * T * T) - (T * T * T / 3526000) + (T * T * T * T / 863310000)) % 360)

def CalculateOhmega(T):
    # Longitude of the ascending node of the moons mean orbit on the ecliptic
    return((125.04452 - (1934.136261 * T) + (0.0020708 * T * T) + (T * T * T / 450000) ) % 360)

def CalculateE(T):
    # calculate Earths Eccentricity of orbit around sum from Meeus 47.6
    return(1 - (0.002516 * T) - (0.0000074 * T * T)) 

def EpsilonPsi(T,D, M, Mdash, F):
    # calculate epsilon and psi based on ch 22 of Meeus
    # !!! units are scaled by 10000 and the results need to be divided by this to correct !!!

    Ohmega = CalculateOhmega(T)

    # matrix of arguments of D, M, Mdash, F and Ohmega and coefficients of delta-psi and delta-epsilon
    # [ [D, M, Mdash, F, Ohmega, Delta-psi, Delta-epsilon] ]
    psiepsilonMatrix = [ \
        [0,0,0,0,1,-171996 - (174.2 * T), 92025 + (8.9 * T)],
        [-2,0,0,2,2,-13187 - (1.6 * T),5736 - (3.1 * T)],
        [0,0,0,2,2,-2274 - (0.2 * T),977 - (0.5 * T)],
        [0,0,0,0,2,2062 + (0.2 * T), -895 + (0.5 * T)],
        [0,1,0,0,0,1426 - (3.4 * T), 54 - (0.1 * T)],
        [0,0,1,0,0,712 + (0.1 * T), -7],
        [-2,1,0,2,2,-517 + (1.2 * T),224 - (0.6 * T)],
        [0,0,0,2,1,-386 - (0.4 * T),200],
        [0,0,1,2,2,-301,129 - (0.1 * T)],
        [-2,-1,0,2,2,217 - (0.5 * T), -95 + (0.3 * T)],
        [-2,0,1,0,0,-158,0],
        [-2,0,0,2,1,129 + (0.1 * T),-70],
        [0,0,-1,2,2,123,-53],
        [2,0,0,0,0,63,0],
        [0,0,1,0,1,63 + (0.1 * T),-33],
        [2,0,-1,2,2,-59,26],
        [0,0,-1,0,1,-58 - (0.1 * T),32],
        [0,0,1,2,1,-51,27],
        [-2,0,2,0,0,48,0],
        [0,0,-2,2,1,46,-24],
        [2,0,0,2,2,-38,16],
        [0,0,2,2,2,-31,13],
        [0,0,2,0,0,29,0],
        [-2,0,1,2,2,29,-12],
        [0,0,0,2,0,26,0],
        [-2,0,0,2,0,-22,0],
        [0,0,-1,2,1,21,-10],
        [0,2,0,0,0,17 - (0.1 * T),0],
        [2,0,-1,0,1,16,-8],
        [-2,2,0,2,2,-16 + (0.1 * T),7],
        [0,1,0,0,1,-15,9],
        [-2,0,1,0,1,-13,7],
        [0,-1,0,0,1,-12,6],
        [0,0,2,-2,0,11,0],
        [2,0,-1,2,1,-10,5],
        [2,0,1,2,2,-8,3],
        [0,1,0,2,2,7,-3],
        [-2,1,1,0,0,-7,0],
        [0,-1,0,2,2,-7,3],
        [2,0,0,2,1,-7,3],
        [2,0,1,0,0,6,0],
        [-2,0,2,2,2,6,-3],
        [-2,0,1,2,1,6,-3],
        [2,0,-2,0,1,-6,3],
        [2,0,0,0,1,-6,3],
        [0,-1,1,0,0,5,0],
        [-2,-1,0,2,1,-5,3],
        [-2,0,0,0,1,-5,3],
        [0,0,2,2,1,-5,3],
        [-2,0,2,0,1,4,0],
        [-2,1,0,2,1,4,0],
        [0,0,1,-2,0,4,0],
        [-1,0,1,0,0,-4,0],
        [-2,1,0,0,0,-4,0],
        [1,0,0,0,0,-4,0],
        [0,0,1,2,0,3,0],
        [0,0,-2,2,2,-3,0],
        [-1,-1,1,0,0,-3,0],
        [0,1,1,0,0,-3,0],
        [0,-1,1,2,2,-3,0],
        [2,-1,-1,2,2,-3,0],
        [0,0,3,2,2,-3,0],
        [2,-1,0,2,2,-3,0]
        ]
    
    # in seconds of arc - needs converting from dms to decimal.
    delta_psi = psiepsilonMatrixSum(psiepsilonMatrix,5,"sin",D,M,Mdash,F,Ohmega) / 10000 
    delta_epsilon = psiepsilonMatrixSum(psiepsilonMatrix,6,"cos",D,M,Mdash,F,Ohmega) / 10000
    
    #print(delta_psi)
    #print(delta_epsilon)
    
    U = T / 100
    
    # from Meeus 22.3
    epsilon_0 = 23 + (26/60) + (21.448 / 3600) - ((4680.93 / 3600) * U ) - \
                (1.55 * U * U) + \
                (1999.25 * U * U * U) - \
                (51.38 * U * U * U * U) - \
                (249.67 * U * U * U * U * U) - \
                (39.05 * U * U * U * U * U * U) + \
                (7.12 * U * U * U * U * U * U * U) + \
                (27.87 * U * U * U * U * U * U * U * U) + \
                (5.79 * U * U * U * U * U * U * U * U * U) + \
                (2.45 * U * U * U * U * U * U * U * U * U * U)
    
    #print(epsilon_0)
    
    # from Meeus ch. 22
    epsilon = epsilon_0 + (delta_epsilon / 3600)
    
    #print(epsilon)
    #print("T",T,"D",D,"M",M,"Mdash",Mdash,"F",F,"Ohm",Ohmega)
    
    return((delta_psi,delta_epsilon,epsilon_0,epsilon))

def RAandDec(T):
    # get the Rise Ascension (alpha) and declination (delta) of the Moon
    
    # calculate angles L', D, M, M'
    
    Ldash = CalculateLdash(T)
    D = CalculateD(T)
    M = CalculateM(T)
    Mdash = CalculateMdash(T)
    F = CalculateF(T)

    # Calculate 3 further arguments
    A1 = (119.75 + (131.849 * T)) % 360
    A2 = (53.09 + (479264.290 * T)) % 360
    A3 = (313.45 + (481266.484 * T)) % 360
    
    #print(A1)
    #print(A2)
    #print(A3)
    
    # calculate Earths Eccentricity of orbit around sum from Meeus 47.6
    E = CalculateE(T) 
    
    #print(E)
    
    # matrix of coefficients and arguments for Sum l and Sum r (Table 47.A)
    # l and r matrix is as follows:
    # form is [D,M,Mdash,F, Sum l coefficient, Sum r coefficient]
    landrMatrix = [ \
        [0,0,1,0,6288774,-20905355],
        [2,0,-1,0,1274027,-3699111],
        [2,0,0,0,658314,-2955968],
        [0,0,2,0,213618,-569925],
        [0,1,0,0,-185116,48888],
        [0,0,0,2,-114332,-3149],
        [2,0,-2,0,58793,246158],
        [2,-1,-1,0,57066,-152138],
        [2,0,1,0,53322,-170733],
        [2,-1,0,0,45758,-204586],
        [0,1,-1,0,-40923,-129620],
        [1,0,0,0,-34720,108743],
        [0,1,1,0,-30383,104755],
        [2,0,0,-2,15327,10321],
        [0,0,1,2,-12528,0],
        [0,0,1,-2,10980,79661],
        [4,0,-1,0,10675,-34782],
        [0,0,3,0,10034,-23210],
        [4,0,-2,0,8548,-21636],
        [2,1,-1,0,-7888,24208],
        [2,1,0,0,-6766,30824],
        [1,0,-1,0,-5163,-8379],
        [1,1,0,0,4987,-16675],
        [2,-1,1,0,4036,-12831],
        [2,0,2,0,3994,-10445],
        [4,0,0,0,3861,-11650],
        [2,0,-3,0,3665,14403],
        [0,1,-2,0,-2689,-7003],
        [2,0,-1,2,-2602,0],
        [2,-1,-2,0,2390,10056],
        [1,0,1,0,-2348,6322],
        [2,-2,0,0,2236,-9884],
        [0,1,2,0,-2120,5751],
        [0,2,0,0,-2069,0],
        [2,-2,-1,0,2048,-4950],
        [2,0,1,-2,-1773,4130],
        [2,0,0,2,-1595,0],
        [4,-1,-1,0,1215,-3958],
        [0,0,2,2,-1110,0],
        [3,0,-1,0,-892,3258],
        [2,1,1,0,-810,2616],
        [4,-1,-2,0,759,-1897],
        [0,2,-1,0,-713,-2117],
        [2,2,-1,0,-700,2354],
        [2,1,-2,0,691,0],
        [2,-1,0,-2,596,0],
        [4,0,1,0,549,-1423],
        [0,0,4,0,537,-1117],
        [4,-1,0,0,520,-1571],
        [1,0,-2,0,-487,-1739],
        [2,1,0,-2,-399,0],
        [0,0,2,-2,-381,-4421],
        [1,1,1,0,351,0],
        [3,0,-2,0,-340,0],
        [4,0,-3,0,330,0],
        [2,-1,2,0,327,0],
        [0,2,1,0,-323,1165],
        [1,1,-1,0,299,0],
        [2,0,3,0,294,0],
        [2,0,-1,-2,0,8752]
        ]
    
    #print(landrMatrix)
    
    # matrix of coefficients and arguments for Sum b (Table 47.B)
    # b matrix is as follows:
    # form is [D,M,Mdash,F, Sum b coefficient]
    bMatrix = [ \
        [0,0,0,1,5128122],
        [0,0,1,1,280602],
        [0,0,1,-1,277693],
        [2,0,0,-1,173237],
        [2,0,-1,1,55413],
        [2,0,-1,-1,46271],
        [2,0,0,1,32573],
        [0,0,2,1,17198],
        [2,0,1,-1,9266],
        [0,0,2,-1,8822],
        [2,-1,0,-1,8216],
        [2,0,-2,-1,4324],
        [2,0,1,1,4200],
        [2,1,0,-1,-3359],
        [2,-1,-1,1,2463],
        [2,-1,0,1,2211],
        [2,-1,-1,-1,2065],
        [0,1,-1,-1,-1870],
        [4,0,-1,-1,1828],
        [0,1,0,1,-1794],
        [0,0,0,3,-1749],
        [0,1,-1,1,-1565],
        [1,0,0,1,-1491],
        [0,1,1,1,-1475],
        [0,1,1,-1,-1410],
        [0,1,0,-1,-1344],
        [1,0,0,-1,-1335],
        [0,0,3,1,1107],
        [4,0,0,-1,1021],
        [4,0,-1,1,833],
        [0,0,1,-3,777],
        [4,0,-2,1,671],
        [2,0,0,-3,607],
        [2,0,2,-1,596],
        [2,-1,1,-1,491],
        [2,0,-2,1,-451],
        [0,0,3,-1,439],
        [2,0,2,1,422],
        [2,0,-3,-1,421],
        [2,1,-1,1,-366],
        [2,1,0,1,-351],
        [4,0,0,1,331],
        [2,-1,1,1,315],
        [2,-2,0,-1,302],
        [0,0,1,3,-283],
        [2,1,1,-1,-229],
        [1,1,0,-1,223],
        [1,1,0,1,223],
        [0,1,-2,-1,-220],
        [2,1,-1,-1,-220],
        [1,0,1,1,-185],
        [2,-1,-2,-1,181],
        [0,1,2,1,-177],
        [4,0,-2,-1,176],
        [4,-1,-1,-1,166],
        [1,0,1,-1,-164],
        [4,0,1,-1,132],
        [1,0,-1,-1,-119],
        [4,-1,0,-1,115],
        [2,-2,0,1,107],
        ]
    
    # calculate the sum of r, l and b
    Sumr = lrbMatrixSum(landrMatrix,E,5,"cos",D,M,Mdash,F)
    Suml = lrbMatrixSum(landrMatrix,E,4,"sin",D,M,Mdash,F)
    Sumb = lrbMatrixSum(bMatrix,E,4,"sin",D,M,Mdash,F)
    
    #print(Sumr)
    
    # calculate additive terms to above sums   
    # additive to Sum l:    
    Suml += 3958 * math.sin(math.radians(A1)) + \
            1962 * math.sin(math.radians(Ldash - F)) + \
            318 * math.sin(math.radians(A2))
    #print(Suml)
    
    # additive to Sum b:    
    Sumb += -2235 * math.sin(math.radians(Ldash)) + \
            382 * math.sin(math.radians(A3)) + \
            175 * math.sin(math.radians(A1 - F)) + \
            175 * math.sin(math.radians(A1 + F)) + \
            127 * math.sin(math.radians(Ldash - Mdash)) - \
            115 * math.sin(math.radians(Ldash + Mdash))
    
    #print(Sumb)
    
    # calculate lambda, beta and delta
    MoonLambda = Ldash + (Suml / 1000000)
    MoonBeta = Sumb / 1000000
    MoonDelta = 385000.56 + (Sumr / 1000) # distance
    MoonPi = math.degrees(math.asin(6378.14 / MoonDelta)) % 360
    
    
    # get delta_psi and epsilon from calculations in Meeus Ch 22
    delta_psi,delta_epsilon,epsilon_0,epsilon = EpsilonPsi(T,D, M, Mdash, F)
    #print(epsilon)
    ApparentLambda = MoonLambda + (delta_psi / 3600) # in deg.
    
    #print(ApparentLambda)
    #print(MoonLambda)
    #print(MoonBeta)
    #print(MoonDelta)
    #print(MoonPi)
    
    # Meeus 13.3
    alpha = (math.degrees(math.atan2(math.sin(math.radians(ApparentLambda)) * math.cos(math.radians(epsilon)) - (math.tan(math.radians(MoonBeta)) * math.sin(math.radians(epsilon))), math.cos(math.radians(ApparentLambda)))))
    
    # Meeus 13.4
    sin_delta = ((math.sin(math.radians(MoonBeta)) * math.cos(math.radians(epsilon))) + (math.cos(math.radians(MoonBeta)) * math.sin(math.radians(epsilon)) * math.sin(math.radians(ApparentLambda))))
    #print("sin_delta:",sin_delta)
    if sin_delta < 0:
        delta = (180 - (math.degrees(math.asin(sin_delta))) % 180) * (-1)
    else:   
        delta = (math.degrees(math.asin(sin_delta))) % 180
    
    #print(alpha)
    #print("delta:",delta)
    
    return((alpha,delta,epsilon,delta_psi/3600,MoonPi))

def EstimateMoon(YEAR,MONTH,DAY,Latitude, Longitude,Event):
    # Meeus expects Longitude to be positive in the West and Negative in East
    # calculate Moon Rise and Set for iterative purposes - subsequent estimates should be obtained by increasing the DAY variable by fractions of a day i.e. DAY + (hours / 24)
    # for example DAY = 21; first estimate is 16:15, the DAY should then be entered as 21 + (16.25 / 24) = 21.67708
    
    # alpha is Right Ascention (RA)
    # delta is Declination (Dec)
    # T is Julian Day in respect of JD2000
    # Event is either "Rise","Set", or "Transit"
    
    # check if Event matches one of the following, and if it doesn't, halt execution further
    if ((Event != "Rise") + (Event != "Set") + (Event != "Transit")) == 3:
        raise SystemExit('Event must be one of "Rise","Set" or "Transit"')
    

    T = CalculateT(YEAR,MONTH,DAY)
    JD = JulianDay(YEAR,MONTH,DAY)
    #print(JD)
   
    # get Right Assension, declination, apparent parallax
    alpha,delta,epsilon,delta_psi,MoonPi = RAandDec(T)
    
    #print("alpha:", DegtoHMS(alpha), "delta:", DegtoDMS(delta), MoonPi)
    #print("alpha:", alpha, "delta:",delta, MoonPi)
    
    # standard altitude including apparent parallax (MoonPi)  
    h0 = (0.7275 * MoonPi) - (34/60)
    
    #DEBUG
    #h0 = 0.125
    
    #DEBUG
    #BigTheta_0 = 177.74208
    #h0 = -0.5667
    #alpha = 41.73129
    #delta = 18.44092
    
    #print(alpha)
    # from Meeus 15.1
    CosH0 = (math.sin(math.radians(h0)) - (math.sin(math.radians(Latitude)) * math.sin(math.radians(delta)))) / (math.cos(math.radians(Latitude)) * math.cos(math.radians(delta)))
    
    #print("CosH0",CosH0)
    
    # check if (-1 < CosH0 < 1) is - if it isn't, then the event doesn't happen (i.e. it won't rise or set)
    if abs(CosH0) > 1:
        return(False)
    
    
    H0 = (math.degrees(math.acos(CosH0))) % 180
    #print("H0", H0)
 
    nutation = delta_psi * math.cos(math.radians(epsilon))
    #print("nutation", nutation)
    
    #DEBUG
    BigTheta_0 = Theta0(YEAR,MONTH,DAY) - (nutation / 15)
    #print("BigTheta_0:",BigTheta_0)

    
    # from Meeus 15.2
    m = [0] * 3      
    m[0] = (alpha + Longitude - BigTheta_0) / 360 # transit
    m[1] = m[0] - (H0 / 360) # rising
    m[2] = m[0] + (H0 / 360) # setting
    
    #print("before correction m",m)
    
    # adjust to be between 0 and 1
    m = [i % 1 for i in m]
    
    #print("after correction m",m)
    
    # find sidereal time at Grenwhich in deg from
    theta_0 = [(BigTheta_0 + (360.985647 * i)) % 360 for i in m]
    #print("theta_0 (sidereal time)",theta_0)
    
   
    # DEBUG
    #alpha = 42.59324
    
    # find local hour angle of moon
    H = [((i) - (Longitude) - (alpha)) for i in theta_0]
    #print("H:",H)
    
    #Debug
    #H[0] = -0.05257
    #H[1] = -108.56577
    #H[2] = 108.52570
    #delta = 18.64229
    
    #print("H",H)
    
    # Moons altitude as Meeus 13.6 (ignore altitude for transit i.e. H[0]
    sinh = [(math.sin(math.radians(Latitude)) * math.sin(math.radians(delta))) + (math.cos(math.radians(Latitude)) * math.cos(math.radians(delta)) * math.cos(math.radians(i))) for i in H]
    
    # obtain h from sin h
    h = [math.degrees(math.asin(i)) for i in sinh]
    #print("h",h)

    
    # set up delta_m list for transit / rise and set
    delta_m = [0] * 3
    
    # H[0] for transit needs to be in range -180 to +180
    if H[0] < 0:
        sign = (-1)
        
    else:
        sign = 1
    
    delta_m[0] = (- (abs(H[0]) % 180) * sign) / 360
    
    for i in range(1,3):
        delta_m[i] = (h[i] - h0) / (360 * math.cos(math.radians(delta)) * math.cos(math.radians(Latitude)) * math.sin(math.radians(H[i])))
    
    #print("delta_m", delta_m)
    
    m = [(m[i] + delta_m[i]) * 24 for i in range(3) ]
    
    
    # return correct calculation:
    if Event == "Rise":
        return(m[1])
    elif Event == "Set":
        return(m[2])
    else:
        return(M[0])

    return((m[1],m[2],m[0]))



def MoonTime(YEAR,MONTH,DAY,Latitude,Longitude,Event):
    #
    # Longitude is positive west, negative east!!
    # correction to add to day is fractions of a day
    Times = 0
    TimesPrevDay = 0
    
    for i in range(5):
        # get the required time of the event
        Times = EstimateMoon(YEAR,MONTH,DAY + (Times / 24),Latitude,Longitude,Event)
        # also check yesterday's time of event
        TimesPrevDay = EstimateMoon(YEAR,MONTH, DAY - 1 + (TimesPrevDay / 24),Latitude,Longitude,Event)
        
        # if result is Circumpolar, then no event
        if Times == False:
            return(False)
    
    #print("Time",Times, "Time Previous Day",TimesPrevDay)
    
    # check the difference between yesterday's and today's event, if there is less that 5 mins of difference, then the event didn't happen
    if abs(Times - TimesPrevDay) < (5/60):
        return(False)
    else:
        return(Hrs(Times))
        
    #print("Rise:", Hrs(trise[0]), "Set:",Hrs(tset[1]))
    return(False)

def Phase(YEAR):
    # calculate all 4 phases at once by using lists
    # times off for Full and Last Quarter - phase correction matrix seems to nearly calculate approx the same as Meeus which could be rounding / implementation of language
    
    # YEAR to include fractional part to define month - beware that Month does not necessarily mean the event will occur in that month, and also event could happen twice in the month (blue moon for example)
    # also calculation is centric on New Month.
    
    # calculate k Meeus 49.2
    k_orig = ((YEAR) - 2000) * 12.3685
    
    # k is integer at New Moon, First quarter is +0.25; Full Moon is +0.5; Last Quarter is 0.75
    k = [0] * 4
    kphasedecimal = [0,0.5,0.25,0.75] # this looks a bit odd, but the calculations below work in the order of: New Moon, Full Moon, First Quarter and Last Quarter
    k = [math.floor(k_orig) + i for i in kphasedecimal ]
    #print("k",k)
    
    # calculate T Meeus 49.3
    T = [0] * 4
    T = [i / 1236.85 for i in k]

    #print("T",T)
    
    # Calculate JDE as Meeus 49.1
    JDE = [0] * 4
    JDE = [2451550.09766 + (29.530588861 * k[i]) + (0.00015437 * T[i] * T[i]) - (0.000000150 * T[i] * T[i] * T[i]) + (0.00000000073 * T[i] * T[i] * T[i] * T[i]) for i in range(4)]
    #print("JDE",JDE)

    # calculate E
    E = [0] * 4
    E = [CalculateE(i) for i in T]
    #print("E",E)
    
    # calculate M, Mdash, F, Ohmega - different to those used in Moon Rise / Set times
    
    # Suns mean anomaly at time JDE Meeus 49.4
    M = [0] * 4
    M = [(2.5534 + (29.10535670 * k[i]) - (0.0000014 * T[i] * T[i]) - (0.00000011 * T[i]* T[i] * T[i])) % 360 for i in range(4)] 
    #print("M",M)
    
    # Moons mean anomaly Meeus 49.5
    Mdash = [0] * 4
    Mdash = [(201.5643 + (385.81693528 * k[i]) + (0.0107582 * T[i] * T[i]) + (0.00001238 * T[i] * T[i] * T[i]) - (0.000000058 * T[i] * T[i] * T[i] * T[i])) % 360 for i in range(4)]
    #print("Mdash",Mdash)
    
    # Moon's argument of latitude Meeus 49.6
    F = [0] * 4
    F = [(160.7108 + (390.67050284 * k[i]) - (0.0016118 * T[i] * T[i]) - (0.00000227 * T[i] * T[i] * T[i]) + (0.000000011 * T[i] * T[i] * T[i] * T[i])) % 360 for i in range(4)]
    #print("F",F)
    
    # Longitude of the ascending node of the lunar orbit Meeus 49.7
    Ohmega = [0] * 4
    Ohmega = [(124.7746 - (1.56375588 * k[i]) + (0.0020672 * T[i] * T[i]) + (0.00000215 * T[i] * T[i] * T[i])) % 360 for i in range(4)]
    #print("Ohmega",Ohmega)
    
    
    #DEBUG
    #E = 1.0005753
    #M = 45.7375
    #Mdash = 95.3722
    #F = 120.9584
    #Ohmega = 207.3176
    
    # Planetary corrections - combining both into one matrix (14nr corrections, A1 = A[0] etc)
    # [coeff, argument]
    A = [0] * 4
    A = [[ \
        [0.000325,299.77 + (0.107408 * k[i]) - (0.009173 * T[i] * T[i]) ],
        [0.000165, 251.88 + (0.016321 * k[i])],
        [0.000164, 251.83 + (26.651886 * k[i])],
        [0.000126, 349.42 + (36.412478 * k[i])],
        [0.000110, 84.66 + (18.206239 * k[i])],
        [0.000062, 141.74 + (53.303771 * k[i])],
        [0.000060, 207.14 + (2.453732 * k[i])],
        [0.000056, 154.84 + (7.306860 * k[i])],
        [0.000047, 34.52 + (27.261239 * k[i])],
        [0.000042, 207.19 + (0.121824 * k[i])],
        [0.000040, 291.34 + (1.844379 * k[i])],
        [0.000037, 161.72 + (24.198154 * k[i])],
        [0.000035, 239.56 + (25.513099 * k[i])],
        [0.000023, 331.55 + (3.592518 * k[i])],
        ] for i in range(4)]
    #print("A",A)
    
    #DEBUG    
    #E = [1,1,1,1]
    
    # [ M, Mdash, F, Ohmega, New Moon, Full Moon, First & Last Quarter]    
    CorrectionsMatrix = [ \
        [0,1,0,0,-0.40720,-0.40614,-0.62801,-0.62801],
        [1,0,0,0, 0.17241 * E[0], 0.17302 * E[1],0.17172 * E[2],0.17172 * E[3] ],
        [0,2,0,0,0.01608,0.01614,0.00862,0.00862],
        [0,0,2,0,0.01039,0.01043,0.00804,0.00804],
        [-1,1,0,0,0.00739 * E[0], 0.00734 * E[1],0.00454 * E[2],0.00454 * E[3]],
        [1,1,0,0,-0.00514 * E[0], -0.00515 * E[1],-0.01183 * E[2],-0.01183 * E[3]],
        [2,0,0,0,0.00208 * E[0] * E[0], 0.00209 * E[1] * E[1],0.00204 * E[2] * E[2],0.00204 * E[3] * E[3]],
        [0,1,-2,0,-0.00111,-0.00111,-0.00180,-0.00180],
        [0,1,2,0,-0.00057,-0.00057,-0.00070,-0.00070],
        [1,2,0,0,0.00056 * E[0],0.00056 * E[1],0.00027 * E[2],0.00027 * E[3]],
        [0,3,0,0,-0.00042,-0.00042,-0.00040,-0.00040],
        [1,0,2,0,0.00042 * E[0], 0.00042 * E[1],0.00032 * E[2],0.00032 * E[3]],
        [1,0,-2,0,0.0038 * E[0], 0.00038 * E[1],0.00032 * E[2],0.00032 * E[3]],
        [-1,2,0,0,-0.00024 * E[0], -0.00024 * E[1],-0.00034 * E[2],-0.00034 * E[3]],
        [0,0,0,1,-0.00017,-0.00017,-0.00017,-0.00017],
        [2,1,0,0,-0.00007,-0.00007,-0.0028 * E[2] * E[2],-0.0028 * E[3] * E[3]],
        [0,2,-2,0,0.00004,0.00004,0.00002,0.00002],
        [3,0,0,0,0.00004,0.00004,0.00003,0.00003],
        [1,1,-2,0,0.00003,0.00003,0.00003,0.00003],
        [0,2,2,0,0.00003,0.00003,0.00004,0.00004],
        [1,1,2,0,-0.00003,-0.00003,-0.00004,-0.00004],
        [-1,1,2,0,0.00003,0.00003,0.00002,0.00002],
        [-1,1,-2,0,-0.00002,-0.00002,-0.00005,-0.00005],
        [1,3,0,0,-0.00002,-0.00002,-0.00002,-0.00002],
        [0,4,0,0,0.00002,0.00002,0,0],
        [-2,1,0,0,0,0,0.00004,0.00004]
            ]
    # DEBUG
    #N = []
    #for row in CorrectionsMatrix:
    #    N.append([row[0],row[1],row[2],row[3],row[5]])
    #print(N)
    
    
    # W for quarter phases only
    W = [0] * 2
    W = [0.00306 - (0.00038 * E[i] * math.cos(math.radians(M[i]))) + (0.00026 * math.cos(math.radians(Mdash[i]))) - (0.00002 * math.cos(math.radians(Mdash[i] - M[i]))) + (0.00002 * math.cos(math.radians(Mdash[i] + M[i]))) + (0.00002 * math.cos(math.radians(2 * F[i]))) for i in range(2,4)]
    #print("W",W)

    # calculate phase correction from the matrix above.
    ApparentPhase = PhaseCorrectionMatrixSum(CorrectionsMatrix,M,Mdash,F,Ohmega)
    
    #print("ApparentPhase",ApparentPhase)
       
    AdditionalCorrection = [0] * 4
    
    for i in range(4):
        for row in A[i]:
            # calculate the additional correction required
            AdditionalCorrection[i] += row[0] * math.sin(math.radians(row[1]))
    #print("AdditionalCorrection", AdditionalCorrection)
    
    # add the additional correction and W to the phases to get correct phase
    for i in range(4):
        ApparentPhase[i] = ApparentPhase[i] + AdditionalCorrection[i]
    
    # correction for first and last quarters
    ApparentPhase[2] += W[0]
    ApparentPhase[3] -= W[1]
    
    #print("Final ApparentPhase",ApparentPhase, "JDE",JDE[3] + ApparentPhase[3])
    
    # create dictionary - labels for dict
    PhaseLabel = ["New Moon","Full Moon","First Quarter", "Last Quarter"]
    
    return(dict(zip(PhaseLabel,[CalendarDate(JDE[i] + ApparentPhase[i]) for i in range(4)])))