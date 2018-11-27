# Jean-paul Ventura
# January 12, 2017
# Modified: February 20th, 2016

from astropy.io import fits
import numpy as np


# Calculate equivalent width
def eqwidth(filename, wavelength1, wavelength2, wavelength3, wavelength4,
            RV_correct):

    """
    Measures the equivalent width of spectral emission/absorption features.

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

    """
    # Open .fits file

    hdu    = fits.open(filename)
    prihdr = hdu[0]
    data   = prihdr.header

    # Define wavelength and flux arrays by accessing .fits sloan data

    lmbda  = []
    flux   = prihdr.data[0]
    for i in range(0, prihdr.header['NAXIS1']):
        wavelengths = 10**np.array(prihdr.header['COEFF0'] + (prihdr.header['COEFF1'] * i))
        lmbda.append(wavelengths)

    lmbda   = np.array(lmbda)

    # RV correct the wavelength array.

    ha_line = 6562.44   # in Angstroms, SDSS DR7 M dwarf designation for Hα

    c       = 299792.0 # km/s

    obs_lmbda      = ha_line * ((RV_correct/c) + 1 )

    lmbdashift     = ha_line - obs_lmbda

    lmbda          = lmbda + lmbdashift

    # Normalize flux by dividing flux array by mean flux value over relatively quiet continuum region

    normindex = np.where((lmbda >= 8250) & (lmbda <= 8350))
    avgval    = np.mean(flux[normindex])
    nflux     = flux/avgval  # changed to np array, test to see if it
                             # makes a difference...

    # Index for the entire line feature (line + wings)
    # lineindex = np.where((lmbda >= wavelength1) & (lmbda <= wavelength4))
    # Index for the line feature 'wings'
    continuumindex = np.where(((lmbda >= wavelength1) & (lmbda <= wavelength2)) | ((lmbda >= wavelength3) & (lmbda <= wavelength4)))


    # Isolate the wavelength region of the line feature 'wings'.
    x1  = lmbda[continuumindex]

    #Isolate the normalized flux over wavelength region of the feature wings.
    y1  = nflux[continuumindex]

    # get 1st degree polynomial coefficients to wing region flux to create continuum baseline
    m,c = np.polyfit(x1,y1,1)

    #print m,c

    # find error the continuum baseline.
    rmse = np.sqrt( np.mean( np.array((m*x1 + c) - y1)**2 ) )

    # reassign wavelength and flux values to correspond to entire line feature  region
    # lmbda = lmbda[lineindex]
    # lineflux  = nflux[lineindex]

    # assign infinitesimal/wavelength (wavelength sampling rate) to variable.
    # dlmbda = lmbda[2]-lmbda[1]

    # Create flux array corresponding to continuum baseline
    fitline = m * lmbda + c

    lineindex = np.where((lmbda >= wavelength2) & (lmbda <= wavelength3))
    #eqwi = -1.0 * np.sum((1 - nflux[lineindex]/np.mean(y1))*dlmbda)
    eqwi = -1.0 * np.trapz(1 - nflux[lineindex]/np.mean(y1), lmbda[lineindex])

    # For error
    ewidth, error = mcerror(lmbda,nflux,wavelength2,wavelength3,m,c,rmse,y1)

    hdu.close()

    return eqwi , error


def mcerror(lmbda,nflux,wavelengtha,wavelengthb,m,c,rmse,contAB):
    trials  = 1000
    guess   = []
    ewindex = np.where((lmbda >= wavelengtha) & (lmbda <= wavelengthb))
    dlmbda  = lmbda[2]-lmbda[1]

    for i in range(trials):
        fitline = m * lmbda + c + np.random.normal(-1,1) * rmse
        fx = -1.0 * np.trapz(1 - (nflux[ewindex]/np.mean(contAB)),lmbda[ewindex]) - fitline
        #fx = -1.0 * np.sum(1 - nflux[ewindex]/np.mean(fitline))
        guess.append(fx)

    return np.mean(guess), np.std(guess)

#Calculate line height of spectral emission feature
def line_height(filename,wavelength1,wavelength2,wavelength3,wavelength4):

    # Open sloan .fits file
    hdu = fits.open(filename)
    prihdr = hdu[0]        # **check indices for error vector** #
    data = prihdr.header

    # Define wavelength and flux arrays by accessing .fits sloan data
    lmbda = []
    flux  = prihdr.data[0]
    for i in range(0,prihdr.header['NAXIS1']):
        wlenvalue = 10**np.array(prihdr.header['COEFF0'] + (prihdr.header['COEFF1'] * i))
        lmbda.append(wlenvalue)

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
