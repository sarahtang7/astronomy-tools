# LSPR Code - astrometry (returns RA and DEC of object)
## input file takes in RA, DEC, and pixel coordinates of comparison stars
## function takes in input file and pixel coordinates of asteroid (and boolean option of flattening# LSPR Code
# Sarah Tang
# 7/1/2019

import astropy as astropy
from astropy.io import ascii
import numpy as np
import math
from math import *

def LSPR(file, xUnknown, yUnknown, flatten):
    data = astropy.io.ascii.read(file) #read in the text file
    xStar = data[0][:] #column of x coordinates of the given reference star - list
    yStar = data[1][:] #column of y coordinates of the given reference star - list
    RAStar = data[2][:] #column of RA coordinates for the given reference star - list
    DECStar = data[3][:] #column of DEC coordinates for the given reference star - list

    #build sums for matrix used for RA and DEC (6 different)
    numStars = len(xStar)
    sumXStar = np.sum(xStar)
    sumYStar = np.sum(yStar)
    sumXStarSqr = np.sum(xStar**2)
    sumXYStar = np.sum(xStar*yStar)
    sumYStarSqr = np.sum(yStar**2)

    #convert RA sexagesimal J2000 to degrees
    i = 0
    RAStarDeg = []
    while i < len(RAStar):
        line = RAStar[i]
        RAhours = (float(line[:2]))+(float(line[3:5])/60)+(float(line[6:])/3600)
        RAStarDeg.append((RAhours/24)*360)
        i += 1

    #convert DEC sexagesimal J2000 to degrees
    j = 0
    DECStarDeg = []
    while j < len(DECStar):
        line = DECStar[j]
        DECdegrees = (float(line[1:3]))+(float(line[4:6])/60)+(float(line[7:])/3600)
        DECStarDeg.append(DECdegrees)
        j += 1

    #build matrix used for RA and DEC - same for flat and not flat
    midMatrix = np.array([[numStars, sumXStar, sumYStar],
                         [sumXStar, sumXStarSqr, sumXYStar],
                          [sumYStar,sumXYStar, sumYStarSqr]])

    #build sums for RA matrix
    sumRA = np.sum(RAStarDeg)
    sumRAXStar = np.sum(RAStarDeg * xStar)
    sumRAYStar = np.sum(RAStarDeg * yStar)

    #build RA matrix
    RAmatrix = np.array([[sumRA],
                         [sumRAXStar],
                         [sumRAYStar]])

    #invert matrix and multiply matrices
    inverseMidMatrix = np.linalg.inv(midMatrix)
    multMatrixRA = np.matmul(inverseMidMatrix, RAmatrix)
    b1 = float(multMatrixRA[0])
    a11 = float(multMatrixRA[1])
    a12 = float(multMatrixRA[2])

    #build sums for DEC matrix
    sumDEC = np.sum(DECStarDeg)
    sumDECXStar = np.sum(DECStarDeg * xStar)
    sumDECYStar = np.sum(DECStarDeg * yStar)

    #build DEC matrix
    DECmatrix = np.array([[sumDEC],
                          [sumDECXStar],
                          [sumDECYStar]])

    #multiply matrices for DEC
    multMatrixDEC = np.matmul(inverseMidMatrix, DECmatrix)
    b2 = float(multMatrixDEC[0])
    a21 = float(multMatrixDEC[1])
    a22 = float(multMatrixDEC[2])

    #print the 6 plate constants
    print("UNFLATTENED PLATE CONSTANTS")
    print("b1: ", b1, "deg") #plate constants for RA
    print("b2: ", b2, "deg") #plate constants for DEC
    print("a11: ", a11, "deg/pix") #plate constants for RA
    print("a12: ", a12, "deg/pix") #plate constants for RA
    print("a21: ", a21, "deg/pix") #plate constants for DEC
    print("a22: ", a22, "deg/pix") #plate constants for DEC

    #convert all to radians for FLATTEN
    m = 0
    RAStarRad = []
    while m < len(RAStarDeg):
        RAStarRad.append(RAStarDeg[m] * math.pi / 180)
        m += 1

    h = 0
    DECStarRad = []
    while h < len(DECStarDeg):
        DECStarRad.append(DECStarDeg[h] * math.pi / 180)
        h += 1

    #Uncertainty for reference stars
    k = 0
    RAerrorSum = 0
    DECerrorSum = 0
    while k < len(RAStar):
        RAfit = (b1 + a11*xStar[k] + a12*yStar[k])*3600 #calculate RA fit
        RAerrorSum +=((((RAStarDeg[k]*3600))-RAfit)**2)
        DECfit = (b2 + a21*xStar[k] + a22*yStar[k])*3600 #calculate DEC fit
        DECerrorSum +=((((DECStarDeg[k]*3600))-DECfit)**2)
        k += 1

    RAuncertainty = math.sqrt(abs(1/((len(RAStar)-3)))*RAerrorSum) #uncertainty for RA
    DECuncertainty = math.sqrt(abs(1/((len(DECStar)-3)))*DECerrorSum) #uncertainty for DEC

    print()
    print("UNFLATTENED UNCERTAINTY: ")
    print("RA uncertainty: ", RAuncertainty, "arcsec")
    print("DEC uncertainty: ", DECuncertainty, "arcsec")
    print()

    #to find RA and DEC of unknown
    bMatrix = np.array([[b1],
                        [b2]])
    aMatrix = np.array([[a11, a12],
                        [a21, a22]])
    xyMatrix = np.array([[xUnknown],
                         [yUnknown]])
    multipliedMat = np.matmul(aMatrix, xyMatrix)
    addMat = np.add(bMatrix, multipliedMat)
    RAobject = float(addMat[0])
    DECobject = float(addMat[1])

    #convert RAobject from decimal degrees to hh:mm:ss.ss
    RAobjecthr = RAobject/360*24
    RAobjectmin = (RAobjecthr % int(RAobjecthr))*60
    RAobjectsec = (RAobjectmin % int(RAobjectmin))*60

    #convert DECobject from decimal degrees to dd:mm:ss.s
    DECobjectDeg = DECobject
    DECobjectArcMin = (DECobjectDeg % int(DECobjectDeg))*60
    DECobjectArcSec = (DECobjectArcMin % int(DECobjectArcMin))*60

    print("UNFLATTENED RA AND DEC")
    print("RA: ", int(RAobjecthr), ":", int(RAobjectmin), ":", RAobjectsec)
    print("DEC: ", int(DECobjectDeg), ":", int(DECobjectArcMin), ":", DECobjectArcSec)
    print()

    #flatten RA and DEC of stars
    D = radians(np.sum(DECStarDeg)/len(DECStarDeg)) #average DEC of reference stars
    A = radians(np.sum(RAStarDeg)/len(RAStarDeg)) #average RA of reference stars
    flatRAvalues = []
    flatDECvalues = []
    val = 0
    L = 3.911/0.000009
    while val < len(RAStarDeg):
        H = (sin(DECStarRad[val]))*sin(D) + (cos(DECStarRad[val]))*cos(D)*cos(RAStarRad[val]-A)
        flatRA = (((cos(DECStarRad[val])*sin(RAStarRad[val]-A)) / H)) - (xStar[val]/L)
        flatRAvalues.append(flatRA)
        flatDEC = (((sin(DECStarRad[val])*cos(D)) - (cos(DECStarRad[val])*sin(D)*cos(RAStarRad[val]-A))) / H) - (yStar[val]/L)
        flatDECvalues.append(flatDEC)
        val += 1

    #build sums for RA matrix - FLATTEN
    if (flatten==True):
        sumRAflat = np.sum(flatRAvalues)
        sumRAXStarflat = np.sum(flatRAvalues * xStar)
        sumRAYStarflat = np.sum(flatRAvalues * yStar)

        #build RA matrix - FLATTEN
        RAmatrixflat = np.array([[sumRAflat],
                             [sumRAXStarflat],
                             [sumRAYStarflat]])
        
        #invert matrix and multiply matrices - FLATTEN
        inverseMidMatrixflat = np.linalg.inv(midMatrix)
        multMatrixRAflat = np.matmul(inverseMidMatrixflat, RAmatrixflat)
        b1flat = float(multMatrixRAflat[0])
        a11flat = float(multMatrixRAflat[1])
        a12flat= float(multMatrixRAflat[2])

        #build sums for DEC matrix - FLATTEN
        sumDECflat = np.sum(flatDECvalues)
        sumDECXStarflat = np.sum(flatDECvalues * xStar)
        sumDECYStarflat = np.sum(flatDECvalues * yStar)

        #build DEC matrix - FLATTEN
        DECmatrixflat = np.array([[sumDECflat],
                              [sumDECXStarflat],
                              [sumDECYStarflat]])

        #multiply matrices for DEC - FLATTEN
        multMatrixDECflat = np.matmul(inverseMidMatrixflat, DECmatrixflat)
        b2flat = float(multMatrixDECflat[0])
        a21flat = float(multMatrixDECflat[1])
        a22flat = float(multMatrixDECflat[2])

        #print the 6 plate constants - FLATTEN
        print("*************************")
        print()
        print("FLATTENED PLATE CONSTANTS")
        print("b1: ", b1flat, "deg") #plate constants for RA
        print("b2: ", b2flat, "deg") #plate constants for DEC
        print("a11: ", a11flat, "deg/pix") #plate constants for RA
        print("a12: ", a12flat, "deg/pix") #plate constants for RA
        print("a21: ", a21flat, "deg/pix") #plate constants for DEC
        print("a22: ", a22flat, "deg/pix") #plate constants for DEC
        print()

        #Uncertainty for reference stars - FLATTENED
        t = 0
        RAerrorSumflat = 0
        DECerrorSumflat = 0
        while t < len(RAStar):
            RAfitflat = ((b1flat + a11flat*xStar[t] + a12flat*yStar[t]) + (xStar[t]/L)) #calculate RA fit - flat, in degrees
            DECfitflat = ((b2flat + a21flat*xStar[t] + a22flat*yStar[t]) + (yStar[t]/L)) #calculate DEC fit - flat, in degrees
            delta = cos(D) - DECfitflat*sin(D)
            denom = math.sqrt(RAfitflat**2 + delta**2)
            
            RAfitunflattened = A + atan(RAfitflat/delta) #unflatten
            DECfitunflattened = atan((sin(D)+DECfitflat*cos(D)) / denom) #unflatten

            #convert to degrees
            RAfitunflattened = RAfitunflattened / math.pi * 180
            DECfitunflattened = DECfitunflattened / math.pi * 180
            
            RAerrorSumflat +=(((RAStarDeg[t])-RAfitunflattened)**2)
            DECerrorSumflat +=(((DECStarDeg[t])-DECfitunflattened)**2)
            t += 1

        RAuncertaintyflat = (math.sqrt(abs(1/((len(RAStarDeg)-3)))*RAerrorSumflat))*3600 #uncertainty for RA, convert to arcseconds
        DECuncertaintyflat = (math.sqrt(abs(1/((len(DECStarDeg)-3)))*DECerrorSumflat))*3600 #uncertainty for DEC, convert to arcseconds

        print("FLATTENED UNCERTAINTY:")
        print("RA uncertainty: ", RAuncertaintyflat, "arcsec")
        print("DEC uncertainty: ", DECuncertaintyflat, "arcsec")
        print()

        #to find RA and DEC of unknown
        bMatrixflat = np.array([[b1flat],
                            [b2flat]])
        aMatrixflat = np.array([[a11flat, a12flat],
                            [a21flat, a22flat]])
        xyMatrixflat = np.array([[xUnknown],
                             [yUnknown]])
        multipliedMatflat = np.matmul(aMatrixflat, xyMatrixflat)
        addMatflat = np.add(bMatrixflat, multipliedMatflat)
        
        RAobjectflat = float(addMatflat[0]) + (xUnknown/L)
        DECobjectflat = float(addMatflat[1]) + (yUnknown/L)

        delta = cos(D) - DECobjectflat*sin(D)
        denom = math.sqrt(RAobjectflat**2 + delta**2)
        
        RAobjunflat = A + atan(RAobjectflat/delta) 
        DECobjunflat = atan((sin(D)+DECobjectflat*cos(D)) / denom) 

        #convert RAobject from decimal degrees to hh:mm:ss.ss
        RAobjecthrflat = RAobjunflat/360*24 * (180/math.pi)
        RAobjectminflat = (RAobjecthrflat % int(RAobjecthrflat))*60
        RAobjectsecflat = (RAobjectminflat % int(RAobjectminflat))*60

        #convert DECobject from decimal degrees to dd:mm:ss.s
        DECobjectDegflat = DECobjunflat/math.pi*180
        DECobjectArcMinflat = (DECobjectDegflat % int(DECobjectDegflat))*60
        DECobjectArcSecflat = (DECobjectArcMinflat % int(DECobjectArcMinflat))*60

        print("FLATTENED RA AND DEC")
        print("RA: ", int(RAobjecthrflat), ":", int(RAobjectminflat), ":", RAobjectsecflat)
        print("DEC: ", int(DECobjectDegflat), ":", int(DECobjectArcMinflat), ":", DECobjectArcSecflat)

    return("")

##print(LSPR('LSPRtestinput1.txt', 484.35, 382.62, True))
print(LSPR('LSPRtestinput1.txt', 484.35, 382.62, False))

