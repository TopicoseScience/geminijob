import numpy as np
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pyfits
from astropy.io import fits
## lendo os fits

def img():
    h1=pyfits.open('1.fits')
    h2=pyfits.open('2.fits')
    h1data=h1[1].data
    h2data=h2[1].data
    xccd=h2data.field('x_ccd')
    yccd=h2data.field('y_ccd')
## plotando a imagem com catalogo
    plt.imshow(h1data, cmap='gray', norm=LogNorm())
    plt.colorbar()
    plt.plot(xccd, yccd, 'rx')
    plt.show()


## plotando as imagens de spectro
def get_wstart(ref, wave_ref, wave_per_pixel):
    """
    Obtain the starting wavelength of a spectrum.

    Parameters
    ----------

    ref: int,
        Reference pixel.
    
    wave_ref: float,
        Coordinate at reference pixel.

    wave_per_pixel: float,
        Coordinate increase per pixel.

    Returns
    -------

    wstart: float,
        Starting wavelength.
    """
    
    return wave_ref - ((ref-1) * wave_per_pixel)

def get_wavelength(start_wave, wave_per_pixel, size):
    """
    Obtain an array of wavelengths according to input values.

    Parameters
    ----------

    start_wave: float,
        Starting wavelength.

    wave_per_pixel: float,
        Wavelength per pixel.

    size: int,
        Size of array.

    Returns
    -------

    wave_array: numpy.ndarray,
        Wavelength array
    """
    
    return np.array([start_wave + i*wave_per_pixel for i in range(size)])


def load_spectrum(fname,n):
    """
    Loads the spectrum in FITS format to a numpy.darray.

    Parameters
    ----------

    fname: str,
        File name of the FITS spectrum.

    Returns
    -------

    spectrum: ndarray,
        Spectrum array with wavelength and flux.
    """
    
    # Load spectrum
    spec_FITS = pyfits.open(fname)
    #Load flux
    flux = spec_FITS[n].data
    
    #Obtain parameters for wavelength determination from header
    ref_pixel = spec_FITS[n].header['CRPIX1']       #Reference pixel
    coord_ref_pixel = spec_FITS[n].header['CRVAL1'] #Wavelength at reference pixel
    wave_pixel = spec_FITS[n].header['CDELT1']      #Wavelength per pixel
    
    #Get starting wavelength
    wstart = get_wstart(ref_pixel, coord_ref_pixel, wave_pixel)
    
    #Obtain array of wavelength
    wave = get_wavelength(wstart, wave_pixel, len(flux))
    
    return np.dstack((wave, flux))[0]

def dospec():
    full_spec = load_spectrum('2.fits',2)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'r-')
    full_spec = load_spectrum('2.fits',3)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'b-')
    full_spec = load_spectrum('2.fits',4)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'g-')
    full_spec = load_spectrum('2.fits',5)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'y-')
    full_spec = load_spectrum('2.fits',6)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'c-')
    full_spec = load_spectrum('2.fits',7)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'm-')
    full_spec = load_spectrum('2.fits',8)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'k-')
    full_spec = load_spectrum('2.fits',9)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'w-')
    full_spec = load_spectrum('2.fits',10)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'r-')
    full_spec = load_spectrum('2.fits',11)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'b-')
    full_spec = load_spectrum('2.fits',12)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'g-')
    full_spec = load_spectrum('2.fits',13)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'y-')
    full_spec = load_spectrum('2.fits',14)
    plt.plot(full_spec[:, 0], full_spec[:, 1], 'c-')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.show()
