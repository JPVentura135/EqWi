#!/usr/bin/python

"""
David Rodriguez
Sep 6, 2012
Load up a spectrum and plot it
Optimized for WiFeS
Modified: June 10, 2013
Use new method to estimate EW errors
"""

from astropy.io import fits
from pylab import *
import asciitable
from random import random,uniform

# Calculate equivalent width
def ewidth(wav,dat,w1,w2,w3,w4):
	wav = array(wav)
	dat = array(dat)

	fuind = (wav >= w1) & (wav <= w4)
	coind = ((wav >= w1) & (wav <= w2)) | ((wav >= w3) & (wav <= w4))

	x1 = wav[coind]
	y1 = dat[coind]
	m,c = polyfit(x1,y1,1)
	#print m,c
	rms = sqrt( mean( (m*x1 + c - y1)**2 ) )

	wav = wav[fuind]
	dat = dat[fuind]
	fitl = m*wav + c
	ndat = dat / fitl
	dlam = wav[2]-wav[1]

	ewind = (wav >= w2) & (wav <= w3)
	fx = sum(1-ndat[ewind]) *dlam

	# For error
	fx1, err = mcerr(wav,dat,w2,w3,m,c,rms)

	return fx, err

# =========================

def mcerr(wav,dat,w1,w2,m,c,rms):
	siz = 3000
	guess = []
	dlam = wav[2]-wav[1]

	for i in range(siz):
	  fitl = m*wav + c + uniform(-1,1)*rms
	  ndat = dat / fitl
	  ewind = (wav >= w1) & (wav <= w2)
	  fx = sum(1-ndat[ewind]) *dlam
	  guess.append(fx)

	return mean(guess),std(guess)

# =========================

#num = 174

#d = asciitable.read('2012Sepwifes_Rodriguez.log')
d = asciitable.read('2012Sepwifes_Rodriguez_3.log', names=['file','name','exp','ra','dec','am','ut','type'])

filenums = d['file']
names = d['name']

index = where(filenums == num)
name = names[index]

file = 'sp'+str(num)+'rx.fits'

hdu = fits.open(file)
data = hdu[0].data

crval1 = hdu[0].header['CRVAL1']
cdelt1 = hdu[0].header['CDELT1']
naxis1 = hdu[0].header['NAXIS1']

arr1 = array([float(i) for i in range(naxis1)])

wave = crval1 + arr1*cdelt1

hdu.close()

# Halpha
w1 = 6548
w2 = 6556.
w3 = 6570
w4 = 6580
st = 2 # spread to consider
cc = 0.002

#value = eqwidth(wave,data,w1,w2,w3,w4)
#print value
#val1, err1 = mcew(wave,data,w1,w2,w3,w4,st,cc)
val1, err1 = ewidth(wave,data,w1,w2,w3,w4)
print val1, err1


fuind = (wave >= w1) & (wave <= w4)
coind = ((wave >= w1) & (wave <= w2)) | ((wave >= w3) & (wave <= w4))

x1 = wave[coind]
y1 = data[coind]
m,c = polyfit(x1,y1,1)

wav1 = wave[fuind]
dat1 = data[fuind]
fitl1 = m*wav1 + c
ndat1 = dat1 / fitl1


fuind = (wave >= 6540) & (wave <= 6750)
data = data[fuind]
wave = wave[fuind]

plot(wave,data,"k")
plot(wav1,dat1,"b")
plot(wav1,fitl1,"r")

axvspan(w2,w3,facecolor="0.8")



# Li
w1 = 6700.5
w2 = 6706.
w3 = 6711.5
w4 = 6713.
st = 1.
cc = 0.000

#val2, err2 = mcew(wave,data,w1,w2,w3,w4,st,cc)
val2, err2 = ewidth(wave,data,w1,w2,w3,w4)
print val2, err2


fuind = (wave >= w1) & (wave <= w4)
coind = ((wave >= w1) & (wave <= w2)) | ((wave >= w3) & (wave <= w4))

x1 = wave[coind]
y1 = data[coind]
m,c = polyfit(x1,y1,1)

wav2 = wave[fuind]
dat2 = data[fuind]
fitl2 = m*wav2 + c
ndat2 = dat2 / fitl2


plot(wav2,dat2,"b")
plot(wav2,fitl2,"r")
axvspan(w2,w3,facecolor="0.8")

tout = '%s: %2.3f+/-%0.3f, %1.3f+/-%0.3f' % (name[0], val1, err1, val2, err2)

xlabel('Wavelength (Ang)')
ylabel('Flux')
title(tout)
ranges = axis()
x = [ranges[0],ranges[1]]
x = [6540,6750]
y = [ranges[2],ranges[3]]
#y = [-70,-20]
#x.reverse()
ranges = concatenate((x,y))
axis(ranges)
grid()

outname = name[0]+'_'+str(num)+'.png'
#savefig(outname)

show()
