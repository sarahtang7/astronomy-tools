# Aperture Photometry w/ Fractional Pixels
# Sarah Tang
# Created 7/5/2019

import math
from math import *
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

file = input("Enter the file name. ")
xPixelValue = int(input("Enter the object's x pixel value. "))
yPixelValue = int(input("Enter the object's y pixel value. "))
radius = int(input("Enter the radius of the circular aperture. "))
innerAnnulus = int(input("Enter the inner radius of the background annulus. "))
outerAnnulus = int(input("Enter the outer radius of the background annulus. "))

image = fits.getdata(file)

def getSquareMatrix():
    #splice a square matrix around the aperture
    squareMatrix = image[yPixelValue-radius-1:yPixelValue+radius, xPixelValue-radius-1:xPixelValue+radius]
    return(squareMatrix)

def getAnnSquare():
    annulusSquareMatrix = image[yPixelValue-outerAnnulus-1:yPixelValue+outerAnnulus, xPixelValue-outerAnnulus-1:xPixelValue+outerAnnulus]
    return(annulusSquareMatrix)

def getAnnulus():
    #ANNULUS WITH BORDER PIXELS ACCEPTED
    annulus = np.zeros((2*outerAnnulus+1, 2*outerAnnulus+1))
    numAnnulus = 0
    for rowNum in range(len(getAnnSquare())):
        for columnNum in range(len(getAnnSquare()[0])):
            distCenter = sqrt((columnNum-outerAnnulus)**2 + (rowNum-outerAnnulus)**2) #start four corner test
            if(distCenter<outerAnnulus and distCenter>innerAnnulus):
                annulus[rowNum, columnNum] = getAnnSquare()[rowNum, columnNum]
                numAnnulus += 1
            else:
                annulus[rowNum, columnNum] = 0
    return(annulus, numAnnulus)

def borderAccept(file, xPixelValue, yPixelValue, radius, fracPixels):
    image = fits.getdata(file)

    if fracPixels == True:
        ###FRACTIONAL PIXELS - FOR APERTURE :)
        quadrant = image[yPixelValue-radius-1:yPixelValue+radius, xPixelValue-radius-1:xPixelValue+radius] #square splice
        step = 0.001
        fracPixels = np.zeros((2*radius+1, 2*radius+1)) #same size as quadrant
        numAper = pi*(radius**2)
        for row in range(len(quadrant)):
            for col in range(len(quadrant[0])):
                tempSum = 0 #for each pixel
                currentCol = col - 0.5 + step - radius
                while currentCol <= col + 0.5 - radius:
                    currentRow = row - 0.5 + step - radius
                    while currentRow <= row + 0.5 - radius:
                        dist = sqrt(currentCol**2 + currentRow**2)
                        if dist <= radius:
                            tempSum += step**2
                        currentRow += step
                    currentCol += step
                fracPixels[row, col] = tempSum * quadrant[row, col]
        #print(fracPixels) #fracPixels replaces inApertureBorder

        sumAperBorder = 0
        for row in range(len(fracPixels)):
            for col in range(len(fracPixels[0])):
                    sumAperBorder += fracPixels[row, col]

        #calculate sum of counts in the annulus
        sumAnn = 0
        for row in range(len(getAnnulus()[0])):
            for col in range(len(getAnnulus()[0][0])):
                 sumAnn += getAnnulus()[0][row, col]

        #calculate signal (ADU) - WITH BORDER PIXELS
        sumSkyCountOnly = sumAnn / getAnnulus()[1] #ADU/pix
        signal = sumAperBorder - (numAper*sumSkyCountOnly) #ADU
        print()
        print("BORDER PIXELS ACCEPTED:")
        print("(with fractional pixels)")
        print("Signal in ADU (border pixels accepted): ", signal)

        #calculate SNR
        gain = 0.8 #e-/ADU
        signalElec = signal*gain #in e-
        skyElec = sumSkyCountOnly*gain #in e-/pix
        darkCurrent = 10 #e-/pix
        readNoise = 11 #e-/pix
        electronicNoise = readNoise**2 + (gain/sqrt(12))**2 #squared
        SNRnumer = sqrt(signalElec)
        SNRdenom = sqrt(1 + numAper*(1+(numAper/getAnnulus()[1]))*((skyElec+darkCurrent+electronicNoise)/signalElec))
        SNR = SNRnumer/SNRdenom
        print("SNR (border pixels accepted): ", SNR)

        #instrumental magnitude
        instrumMag = -2.5*log10(signal)
        print("Instrumental Magnitude: ", instrumMag)
        #instrumental magnitude error
        IMerror = 1.0857/SNR
        print("Instrumental Magnitude Error: ", IMerror)
        
    else: #no fractional pixels
        #BORDER PIXELS INCLUDED
        inApertureBorder = np.zeros((2*radius+1, 2*radius+1))
        numAper = 0
        for rowNum in range(len(getSquareMatrix())):
            
            for columnNum in range(len(getSquareMatrix()[0])):
                distLowerLeft = sqrt((columnNum-radius-0.5)**2 + (rowNum-radius+0.5)**2) #start four corner test
                distLowerRight = sqrt((columnNum-radius+0.5)**2 + (rowNum-radius+0.5)**2)
                distUpperLeft = sqrt((columnNum-radius-0.5)**2 + (rowNum-radius-0.5)**2)
                distUpperRight = sqrt((columnNum-radius+0.5)**2 + (rowNum-radius-0.5)**2)
                
                if(distLowerLeft<=radius or distLowerRight <= radius or distUpperLeft <= radius or distUpperRight <= radius):
                    inApertureBorder[rowNum, columnNum] = getSquareMatrix()[rowNum, columnNum]
                    numAper += 1
                    
                else:
                    inApertureBorder[rowNum, columnNum] = 0

        sumAperBorder = 0
        for row in range(len(inApertureBorder)):
            for col in range(len(inApertureBorder[0])):
                    sumAperBorder += inApertureBorder[row, col]

        #calculate sum of counts in the annulus
        sumAnn = 0
        for row in range(len(getAnnulus()[0])):
            for col in range(len(getAnnulus()[0][0])):
                 sumAnn += getAnnulus()[0][row, col]

        #calculate signal (ADU) - WITH BORDER PIXELS
        sumSkyCountOnly = sumAnn / getAnnulus()[1] #ADU/pix
        signal = sumAperBorder - (numAper*sumSkyCountOnly) #ADU
        print()
        print("BORDER PIXELS ACCEPTED:")
        print("Signal in ADU (border pixels accepted): ", signal)

        #calculate SNR
        gain = 0.8 #e-/ADU
        signalElec = signal*gain #in e-
        skyElec = sumSkyCountOnly*gain #in e-/pix
        darkCurrent = 10 #e-/pix
        readNoise = 11 #e-/pix
        electronicNoise = readNoise**2 + (gain/sqrt(12))**2 #squared
        SNRnumer = sqrt(signalElec)
        SNRdenom = sqrt(1 + numAper*(1+(numAper/getAnnulus()[1]))*((skyElec+darkCurrent+electronicNoise)/signalElec))
        SNR = SNRnumer/SNRdenom
        print("SNR (border pixels accepted): ", SNR)

        #instrumental magnitude
        instrumMag = -2.5*log10(signal)
        print("Instrumental Magnitude: ", instrumMag)
        #instrumental magnitude error
        IMerror = 1.0857/SNR
        print("Instrumental Magnitude Error: ", IMerror)
        
    return(" ")

def borderReject(file, xPixelValue, yPixelValue, radius):
    image = fits.getdata(file)
    #BORDER PIXELS NOT INCLUDED
    inApertureNoBorder = getSquareMatrix()
    numAperNoBor = 0
    for rowNum in range(len(getSquareMatrix())):
        
        for columnNum in range(len(getSquareMatrix()[0])):
            distLowerLeft = sqrt((columnNum-radius-0.5)**2 + (rowNum-radius+0.5)**2) #start four corner test
            distLowerRight = sqrt((columnNum-radius+0.5)**2 + (rowNum-radius+0.5)**2)
            distUpperLeft = sqrt((columnNum-radius-0.5)**2 + (rowNum-radius-0.5)**2)
            distUpperRight = sqrt((columnNum-radius+0.5)**2 + (rowNum-radius-0.5)**2)
            
            if(distLowerLeft<=radius and distLowerRight <= radius and distUpperLeft <= radius and distUpperRight <= radius):
                inApertureNoBorder[rowNum, columnNum] = getSquareMatrix()[rowNum, columnNum]
                numAperNoBor+=1
                
            else:
                inApertureNoBorder[rowNum, columnNum] = 0
 
    sumAperNB = 0
    for row in range(len(inApertureNoBorder)):
        for col in range(len(inApertureNoBorder[0])):
                sumAperNB += inApertureNoBorder[row, col]

    #calculate sum of counts in the annulus
    sumAnn = 0
    for row in range(len(getAnnulus()[0])):
        for col in range(len(getAnnulus()[0][0])):
             sumAnn += getAnnulus()[0][row, col]
             
    ###calculate signal (ADU) - W/O BORDER PIXELS
    sumSkyCountNB = sumAnn / getAnnulus()[1] #ADU/pix
    signalNB = sumAperNB - (numAperNoBor*sumSkyCountNB) #ADU
    print("BORDER PIXELS REJECTED:")
    print("Signal in ADU (border pixels rejected): ", signalNB)

    #calculate SNR
    gain = 0.8 #e-/ADU
    signalElecNB = signalNB*gain #in e-
    skyElecNB = sumSkyCountNB*gain #in e-/pix
    darkCurrent = 10 #e-/pix
    readNoise = 11 #e-/pix
    electronicNoise = readNoise**2 + (gain/sqrt(12))**2 #squared
    SNRnumerNB = sqrt(signalElecNB)
    SNRdenomNB = sqrt(1 + numAperNoBor*(1+(numAperNoBor/getAnnulus()[1]))*((skyElecNB+darkCurrent+electronicNoise)/signalElecNB))
    SNRnobor = SNRnumerNB/SNRdenomNB
    print("SNR (border pixels accepted): ", SNRnobor)

    #instrumental magnitude
    instrumMagNB = -2.5*log10(signalNB)
    print("Instrumental Magnitude: ", instrumMagNB)
    #instrumental magnitude error
    IMerrorNB = 1.0857/SNRnobor
    print("Instrumental Magnitude Error: ", IMerrorNB)
    return("")

print(borderAccept(file, xPixelValue, yPixelValue, radius, True)) #if you want fractional pixels --> True

print(borderReject(file, xPixelValue, yPixelValue, radius))

