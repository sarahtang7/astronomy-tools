# FINDS MAGNITUDE OF OBJECT USING COMPARISON STARS
# Sarah Tang
# Created 7/2019

import numpy as np
import math
from astropy.io import fits
import matplotlib.pyplot as plt

def magObj(filename, xCoorStar, yCoorStar, widthStar, magStar, xCoorAst, yCoorAst, widthAst, xCoorBlank, yCoorBlank):
    image = fits.getdata(filename)

    # FOR STAR
    starBox = image[yCoorStar-(int(widthStar/2)): yCoorStar+(int(widthStar/2))+1, xCoorStar-(int(widthStar/2)): xCoorStar+(int(widthStar/2))+1]
    starAndSky = 0
    for row in starBox:
        for column in row:
            starAndSky += column

    blankSkyBox = image[yCoorBlank-1:yCoorBlank+2, xCoorBlank-1:xCoorBlank+2]
    avgSkySum = 0
    for row in blankSkyBox:
        for column in row:
            avgSkySum += column
    avgSky = avgSkySum / 9

    signalStar = starAndSky - (avgSky*(widthStar**2))
    constant = magStar + 2.5*math.log10(signalStar)

    # FOR TARGET OBJECT
    asteroidBox = image[yCoorAst-(int(widthAst/2)): yCoorAst+(int(widthAst/2))+1, xCoorAst-(int(widthAst/2)): xCoorAst+(int(widthAst/2))+1]
    astAndSky = 0
    for row in asteroidBox:
        for column in row:
            astAndSky += column

    signalAst = astAndSky - (avgSky*(widthAst**2))
    magObj = -2.5*math.log10(signalAst) + constant

    # Display original .fits file
    plt.imshow(image)
    plt.gray()
    plt.show()

    return magObj

print("Magnitude of the Object (using first star): ", magObj("sampleimage.fits", 173, 342, 5, 15.26, 351, 154, 3, 200, 200))
print("Magnitude of the Object (using second star): ", magObj("sampleimage.fits", 355, 285, 5, 16.11, 351, 154, 3, 200, 200))
