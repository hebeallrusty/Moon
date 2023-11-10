# Moon Calculations
#
# Calculate rise and set times of the moon based on location.

# The heavy lifting is done by the EstimateMoon function, however MoonRiseSet should be the one invoked as this iterates the result to improve accurancy.
# Moon Rise and Set cannot be calculated without iterating due to the speed at which it moves in the sky.
#
# The method of calculation is based on those by Jean Meeus in Astronomical Algorithms and relevant chapters / equation locations are noted
#
# TODO - deal with circumpolar situations / days when the moon doesn't rise or set


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



def EpsilonPsi(T): #,D, M, Mdash, F):
    # calculate epsilon and psi based on ch 22 of Meeus
    # !!! units are scaled by 10000 and the results need to be divided by this to correct !!!
    
    # Mean elongation of the moon from the sun:
    D = (297.85036 + (445267.11480 * T) - (0.0019142 * T * T) + (T * T * T / 189474)) % 360
    
    # Mean anomaly of the Sun (Earth)
    M = (357.52772 + (35999.050340 * T) - (0.0001603 * T * T) - (T * T * T / 300000)) % 360
    
    # Mean anomaly of the Moon:
    Mdash = (134.96298 + (477198.867398 * T) + (0.0086972 * T * T) + (T * T * T / 56250) ) % 360
    
    # Moons argument of latitude:
    F = (93.27191 + (483202.017538 * T) - (0.0036825 * T * T) + (T * T * T / 327270) ) % 360
    
    # Longitude of the ascending node of the moons mean orbit on the ecliptic
    Ohmega = (125.04452 - (1934.136261 * T) + (0.0020708 * T * T) + (T * T * T / 450000) ) % 360

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
    
    #T = CalculateT(YEAR,MONTH,DAY)
    
    # calculate angles L', D, M, M'
    
    # Moon's mean longitude (Mean equinox of the date)
    Ldash = (218.3164477 + (481267.88123421 * T) - (0.0015786 * T * T) + (T * T * T / 538841) - (T * T * T * T / 65194000)) % 360
    
    # Mean elongation of the Moon
    D = (297.8501921 + (445267.1114034 * T) - (0.0018819 * T * T) + (T * T * T / 545868) - (T * T * T * T / 113065000)) % 360
    
    # Sun's mean anomaly
    M = (357.5291092 + (35999.0502909 * T) - (0.0001536 * T * T) + (T * T * T / 24490000)) % 360
    
    # Moon's mean anomaly
    Mdash = (134.9633964 + (477198.8675055 * T) + (0.0087414 * T * T) + (T * T * T / 69699) - (T * T * T * T / 14712000)) % 360
    
    # Moon's argument of latitude (mean distance of the Moon from its ascending node)
    F = (93.2720950 + (483202.0175233 * T) - (0.0036539 * T * T) - (T * T * T / 3526000) + (T * T * T * T / 863310000)) % 360
    
    #print(Ldash)
    #print(D)
    #print(M)
    #print(Mdash)
    #print(F)
    
    # Calculate 3 further arguments
    A1 = (119.75 + (131.849 * T)) % 360
    A2 = (53.09 + (479264.290 * T)) % 360
    A3 = (313.45 + (481266.484 * T)) % 360
    
    #print(A1)
    #print(A2)
    #print(A3)
    
    # calculate Earths Eccentricity of orbit around sum
    E = 1 - (0.002516 * T) - (0.0000074 * T * T)
    
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
    delta_psi,delta_epsilon,epsilon_0,epsilon = EpsilonPsi(T)#,D, M, Mdash, F)
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

def EstimateMoon(YEAR,MONTH,DAY,Latitude, Longitude):
    # Meeus expects Longitude to be positive in the West and Negative in East
    # calculate Moon Rise and Set
    
    # alpha is Right Ascention (RA)
    # delta is Declination (Dec)
    # T is Julian Day in respect of JD2000

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
    
    H0 = (math.degrees(math.acos(CosH0))) % 180
    #print("H0", H0)
    
    
    #DEBUG
    BigTheta_0 = Theta0(YEAR,MONTH,DAY) 
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
    
    nutation = delta_psi * math.cos(math.radians(epsilon))
    #print("nutation", nutation)
    
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
    
    

    return((m[1],m[2]))



def MoonRiseSet(YEAR,MONTH,DAY,Latitude,Longitude):
    #
    # Longitude is positive west, negative east!!
    # correction to add to day is fractions of a day
    trise = [0,0]
    tset = [0,0]
    
    for i in range(3):
        trise = EstimateMoon(YEAR,MONTH,DAY + (trise[0] / 24),Latitude,Longitude)
        tset = EstimateMoon(YEAR,MONTH,DAY + (tset[1] / 24),Latitude,Longitude)
        
    #print("Rise:", Hrs(trise[0]), "Set:",Hrs(tset[1]))
    
    return([Hrs(trise[0]),Hrs(tset[1])])

    