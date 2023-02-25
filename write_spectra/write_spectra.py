import os
from astropy.table import Table
from specutils import Spectrum1D
from header_function import *
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.nddata import StdDevUncertainty
from astropy.io import fits
import pprint


absolute_path = os.path.dirname(__file__)
NIRSpec_filename = 'VHS1256b_futz.fits'
fits_filename_NIR = os.path.join(absolute_path, NIRSpec_filename)


#  Header information about both NIRSpec and MIRI Spectra
JWST_header_dict = {
    'RA': '194.007637',
    'dec': '-12.957692',
    'generated_history':  'Data Reduced by JWST pipeline version 1.7.2, CRDS version 11.16.12, CRDS context jwst_0977.pmap',
    'object': 'VHS1256b',
    'telescope': 'JWST'
     }


#################Generate NIRSpec Files ##################

# Input your NIRSpec wavelength, flux, and error data here
# Be sure to include units
NIRspec_wavelength = np.arange(0, 5, 0.5) * u.um
NIRspec_flux = [0.15, .15, .07, .05, .27, .022, .06, .037, .02, .2] * u.Jy
NIRspec_error = [0.0014, .006, .0088, .009, .0097, .0098, .007, .0013, .0012, .009] * u.Jy

NIR_table = Table({'wavelength': NIRspec_wavelength,
                   'flux': NIRspec_flux,
                   'flux_uncertainty': NIRspec_error})

spec1d_NIR = Spectrum1D(spectral_axis=NIRspec_wavelength,
                        flux=NIRspec_flux,
                        uncertainty=StdDevUncertainty(NIRspec_error))

# Header Info specific to NIRSpec spectra
NIRSpec_header_dict = {
    'instrument': 'NIRSpec',
    'bandpass': 'NIR'
     }


NIRSpec_dict = {**JWST_header_dict, **NIRSpec_header_dict}  # combine dictionaries to create headers

# Turns dictionary ino properly formatted FITS header
header = compile_header(NIRspec_wavelength, **NIRSpec_dict)

print('BEFORE:')
print(repr(header))
spec1d_NIR.meta['header'] = header

# for key, val in spec1d_NIR.meta.items():
#    print(key, val, type(val))

# Make a dedicated HDU for the header and spectrum
# hdu0 = fits.PrimaryHDU(header=header)
# hdu1 = fits.BinTableHDU(data=NIR_table)
# hdu1.header['EXTNAME'] = 'SPECTRUM'
# hdu1.header.set('OBJECT', NIRSpec_dict['object_name'], 'Object Name')

# Combine the two HDUs into a multi-extension FITS (MEF)
# spectrum_mef = fits.HDUList([hdu0, hdu1])  # hdu0 is header and hdu1 is data

# Write out the MEF to a file
#spectrum_mef.writeto(fits_filename_NIR, output_verify="exception")
try:
    spec1d_NIR.write(fits_filename_NIR, format='tabular-fits', overwrite=True)
    print(f'Wrote out {fits_filename_NIR}')
except:
    raise IOError

# This is extra code which demonstrates how to read the file back in using specutils
spec1d_NIR_rt = Spectrum1D.read(fits_filename_NIR, format='tabular-fits')

print('META:')
header_NIR_meta = spec1d_NIR_rt.meta['header']
pprint.pprint(header_NIR_meta)

print('FITS:')
header_NIR_fits = fits.getheader(fits_filename_NIR, 1)
print(repr(header_NIR_fits))


# name_NIR = header_NIR_meta['OBJECT']  # get object name from the header
# print(name_NIR)

# Plot spectra
ax = plt.subplots()[1]
ax.plot(spec1d_NIR.spectral_axis, spec1d_NIR.flux)
ax.errorbar(spec1d_NIR.spectral_axis.value, spec1d_NIR.flux.value, yerr=spec1d_NIR.uncertainty.array, fmt='-')
ax.set_xlabel(f"Dispersion ({spec1d_NIR.spectral_axis.unit})")
ax.set_ylabel(f"Flux ({spec1d_NIR.flux.unit})")
# plt.title(f"{name_NIR} {header_NIR_meta['TELESCOP']} {header_NIR_meta['INSTRUME']}")
# plt.savefig(f'{absolute_path}/{name_NIR}.png')