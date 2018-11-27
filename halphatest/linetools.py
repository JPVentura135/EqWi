# Jean-paul Ventura
# January 12, 2017
# Modified: February 20th, 2016

from astropy.io import fits
import numpy as np


# Calculate equivalent width
def eqwidth(filename,wavelength1,wavelength2,wavelength3,wavelength4):
    '''

    This function measures the equivalent width of spectral emission/absorption
    features.

    ====================
    Function parameters:
    ====================

    filename 	- SDSS spectra .fits filename entered as a string. You may
                  specify entire path as a string or if file is located in the
                  current working directory, simply the filename
                  (e.g. 'plate-fiber-mjd.fits')

    wavelength1 - wavelength value where continuum region to the left of the
                  spectral feature begins.

    wavelength2 - wavelength value where continuum region to the left of the
                  spectral feature ends.

    wavelength3 - wavelength value where continuum region to the right of the
                  spectral feature begins.

    wavelength4 - wavelength value where continuum region to the right of the
                  spectral feature ends.

    =================
    Function returns:
    =================

    The measured equivalent width of the spectral emission/absorption feature

    '''

    # Open .fits file
    hdu = fits.open(filename)

    
    # Define wavelength and flux arrays by accessing .fits sloan data
    lmbda = 10**np.array(hdu[1].data['loglam'])
    flux = hdu[1].data['flux']


    # Normalize flux by dividing flux array by mean flux value over relatively quiet
    # continuum region
    indx = np.where((lmbda >= 8250) & (lmbda <= 8350))
    avgval = np.mean(flux[indx])
    nflux = np.array(flux)/avgval  # changed to np array, test to see if it
                                   # makes a difference...

    # Index for the entire line feature (line + wings)
    lineindex = np.where((lmbda >= wavelength1) & (lmbda <= wavelength4))

    # Index for the line feature 'wings'
    continuumindex = np.where((lmbda >= wavelength1) & (lmbda <= wavelength2)) | ((lmbda >= wavelength3) & (lmbda <= wavelength4))

    # Isolate the wavelength region of the line feature 'wings'.
    x1 = lmbda[continuumindex]

    #Isolate the normalized flux over wavelength region of the feature wings.
    y1 = nflux[continuumindex]

    # get 1st degree polynomial coefficients to wing region flux to create continuum baseline
    m,c = np.polyfit(x1,y1,1)

    #print m,c

    #find error the continuum baseline.
    rmse = np.sqrt( np.mean( (m*x1 + c - y1)**2 ) )

    # reassign wavelength and flux values to correspond to entire line feature  region
    lmbda = lmbda[lineindex]
    flux  = nflux[lineindex]

    # assign infinitesimal/wavelength (wavelength sampling rate) to variable.
    dlmbda = lmbda[2]-lmbda[1]

    # Create flux array corresponding to continuum baseline
    fitline = m * lmbda + c

    ewindex = np.where((lmbda >= wavelength2) & (lmbda <= wavelength3))
    eqwi = np.sum(1-np.array(nflux[ewindex]) ) * dlmbda

    # For error
    ewidth, error = mcerror(lmbda,nflux,wavelength2,wavelength3,m,c,rmse)

    return eqwi , error


def mcerror(lmbda,nflux,wavelength1,wavelength2,m,c,rmse):
    trials  = 3000
    guess   = []
    ewindex = np.where((lmbda >= wavelength1) & (lmbda <= wavelength2))
    dlmbda  = lmbda[2]-lmbda[1]

    for i in range(trials):
        fitline = m * lmbda + c + np.uniform(-1,1) * rmse
        fx = np.sum(1-nflux[ewindex]) * dlmbda
        guess.append(fx)

    return np.mean(guess), np.std(guess)

#Calculate line height of spectral emission feature
def line_height(filename,wavelength1,wavelength2,wavelength3,wavelength4):

    # Open sloan .fits file
    hdu = fits.open(filename)

    #Construct wavelength (lmbda) and flux arrays from data
    lmbda = 10**np.array(hdu[1].data['loglam'])
    flux = hdu[1].data['flux']

    # Normalize flux by dividing flux array by mean flux value over relatively quiet
    # continuum region
    indx = np.where((lmbda > 8250) & (lmbda < 8350))
    avgval = flux[indx].mean()
    nflux = flux/avgval

    # designate fundamental and continuum indices
    fuind = (lmbda >= wavelength1) & (lmbda <= wavelength4)
    coind = ((lmbda >= wavelength1) & (lmbda <= wavelength2)) | ((lmbda >= wavelength3) & (lmbda <= wavelength4))

    lmbda = lmbda[fuind]
    flux = nflux[fuind]
    x1 = lmbda[coind]
    y1 = nflux[coind]
    m,c = np.polyfit(x1,y1,1)
    #print m,c
    rms = sqrt( mean( (m*x1 + c - y1)**2 ) )


    lnheight = max(y1)-min(y1)

    return lnheight
