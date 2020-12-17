# FINDS CENTROID OF OBJECT
# Sarah Tang
# Created 7/2019

import math
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

def findCentroid(fileName, xCoor, yCoor, radius):
    image = fits.getdata(fileName)
    #print(image[153:156, 350:353])

    #GET X POSITION
    i = xCoor-radius
    sizeOfCol = xCoor+radius+1
    sumForX = 0
    numerator = 0
    while i < sizeOfCol:
        eachCol = 0
        for column in image:
            eachCol = eachCol + column[i]
        numerator += eachCol*i
    ##    print(eachCol)
        sumForX += eachCol
        i += 1
    numeratorForX = numerator

    #GET Y POSITION
    sizeOfRow = yCoor+radius+1
    j = yCoor-radius
    sumForY = 0
    numeratorY = 0
    while j < sizeOfRow:
        eachRow = 0
        for row in image:
            eachRow += np.sum(j,)
        numeratorY += eachRow*j
    ##    print(eachRow)
        sumForY += eachRow
        j += 1
    numeratorForY = numeratorY

    #TO GET X
    xValueCentroid = numeratorForX / sumForX

    #TO GET Y
    yValueCentroid = numeratorForY / sumForY

    #FOR X ERROR
    sizeOfCol = xCoor+radius+1 
    i = xCoor-radius
    sumForX = 0
    numerator = 0
    k = 0
    numeratorErrX = 0
    while i < sizeOfCol:
        eachCol = 0
        for column in image:
            eachCol = eachCol + column[i]
        numeratorErrX += (np.sum(image[yCoor-radius:yCoor+radius+1, i:i+1]))*((i-xValueCentroid)**2)
        numerator += eachCol*i
    ##    print(eachCol)
        sumForX += np.sum(image[yCoor-radius:yCoor+radius+1, i:i+1])
        i += 1
    sumX = sumForX
    numeratorForX = numerator
    numeratorErrorX = numeratorErrX
    errX = math.sqrt(numeratorErrorX/(sumX-1))
    errorX = errX / math.sqrt(sumX)
    ##    print(errorX)

    #FOR Y ERROR
    sizeOfRow = yCoor+radius+1
    j = yCoor-radius
    sumForY = 0
    numeratorY = 0
    h = 0
    numeratorErrY = 0
    while j < sizeOfRow:
        eachRow = 0
        for row in image:
            eachRow += np.sum(j,)
        numeratorY += eachRow*j
        numeratorErrY += (np.sum(image[j:j+1, xCoor-radius:xCoor+radius+1]))*((j-yValueCentroid)**2)
    ##    print(eachRow)
        sumForY += np.sum(image[j:j+1, xCoor-radius:xCoor+radius+1])
        j += 1
    sumY = sumForY
    numeratorForY = numeratorY
    numeratorErrorY = numeratorErrY
    errY = math.sqrt(numeratorErrorY/(sumY-1))
    errorY = errY / math.sqrt(sumY)
    ##print(errorY)

    return(xValueCentroid, yValueCentroid, errorX, errorY)

print("Centroid (x, y), Uncertainty (x, y): ", findCentroid("sampleimage.fits", 351, 154, 1))
