'''
Astropy's FITS file format analyzer is actually extremely useful.
This file just shows some examples of its use. In the future, it
might be required to work with FITS files, so I will just leave
this here.
'''

import numpy as np
import math
from matplotlib import pyplot as plt
import astropy.io.fits as FITS
import astropy
import astropy.cosmology as cosmo
from astropy.cosmology import LambdaCDM
dark_cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
import csv
import pandas as pd
#import vtk
#import vtk.util.numpy_support as nps

def readingVotable(fileName):
    from astropy.io.votable import parse_single_table
    dataTable = parse_single_table(fileName).to_table(use_names_over_ids=True)
    dataTableFrame = dataTable.to_pandas()
    print dataTableFrame.columns
    return dataTableFrame

def readingFits(fileName):
    fits_file = FITS.open('DR10Q_v2.fits', memmap=True)
    print(fits_file.info())
    print "\n"
    # fits_file[0] is the primary table
    fits_file_headr_0 = fits_file[0].header
    fits_file_data_0 = fits_file[0].data
    fits_file_headr_1 = fits_file[1].header
    fits_file_data_1 = fits_file[1].data
    fits_file_columns_1 = fits_file[1].columns
    #fits_file_columns_2 = fits_file[2].columns
    cols = fits_file_columns_1
    #print(cols.names)
    fits_file.close()
    return fits_file_data_1
fits_file_data_1 = readingVotable("II_326_zcatrev-161102_a.vot")

"""zipped = zip(fits_file_data_1['RA'],  fits_file_data_1['DEC'],
        fits_file_data_1['Z_VI'],  fits_file_data_1['MI'],
        fits_file_data_1['DGMI'],
        fits_file_data_1['ALPHA_NU'],fits_file_data_1['FWHM_CIV'],
        fits_file_data_1['REWE_CIV'], fits_file_data_1['FWHM_CIII'],
        fits_file_data_1['REWE_CIII'], fits_file_data_1['FWHM_MGII'],
        fits_file_data_1['REWE_MGII'])
np.savetxt('original.csv', zipped,
        fmt='%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f')"""
RAJ = fits_file_data_1['RAJ2000'] # RA in the fits file
DEJ = fits_file_data_1['DEJ2000']#DC in the fits file
zb = fits_file_data_1['zsp']# Z_VI in the fits file
#FWCIV = fits_file_data_1['FWHM_CIV']
#iband = fits_file_data_1['imag'] #'MI' in the fits file

max_zb = 0
maxCo  = 0
maxLum = 0

x_no_change = []
x_co = []
x_co_scaled = []
x_lum = []
x_lum_scaled = []
y_no_change = []
y_co = []
y_co_scaled = []
y_lum = []
y_lum_scaled = []
z_no_change = [] 
z_co = [] 
z_co_scaled = []
z_lum = []
z_lum_scaled = []
co_zb = []
co_zb_scaled = []
lum_zb = []
lum_zb_scaled = []
result = []
for i in range(0,len(zb)):
	co_zb.append( dark_cosmo.comoving_distance(zb[i]).value) 
	lum_zb.append(dark_cosmo.luminosity_distance(zb[i]).value)

	if zb[i] > max_zb:
	   max_zb = zb[i]
	if co_zb[i] > maxCo:
	   maxCo = co_zb[i]
	if lum_zb[i] > maxLum:
	   maxLum = lum_zb[i] 

	'''
	This block of code calculates the X, Y, and Z coordinates of the
	galaxies described in the CSV file we are reading. We calculate
	distance based purely on using redshift as a distance measure,
	in addition to comoving and luminosity distance measures.

	These are spherical coordinate calculations. Explanations for the
	measures RAJ and DEJ are available online, but they are essentially
	the angles used in such a coordinate system. The comoving distance,
	the luminosity distance, and zb are used as the radial measures. 

	'''
	 ########################################################################
	x_no_change.append(zb[i] * math.sin(math.radians(90.0 - DEJ[i])) * math.cos(math.radians(RAJ[i])))

	x_co.append(co_zb[i] * math.sin(math.radians(90.0 - DEJ[i])) * math.cos(math.radians(RAJ[i])))

	x_lum.append(lum_zb[i] * math.sin(math.radians(90.0 - DEJ[i])) * math.cos(math.radians(RAJ[i])))

	########################################################################

	y_no_change.append(zb[i] * math.sin(math.radians(90.0 - DEJ[i])) * math.sin(math.radians(RAJ[i])))

	y_co.append(co_zb[i] * math.sin(math.radians(90.0 - DEJ[i])) * math.sin(math.radians(RAJ[i])))

	y_lum.append(lum_zb[i] * math.sin(math.radians(90.0 - DEJ[i])) * math.sin(math.radians(RAJ[i])))
	########################################################################

	z_no_change.append(zb[i]  * math.cos(math.radians(90.0 - DEJ[i])))

	z_co.append(co_zb[i]  * math.cos(math.radians(90.0 - DEJ[i])))

	z_lum.append(lum_zb[i] * math.cos(math.radians(90.0 - DEJ[i])))

	########################################################################

	''' 
	This is the black hole mass calculation. Lots
	of magic numbers. This calculation should be
	documented elsewhere in this repository.
	
	
	if FWCIV[i] > 0:
	    	FWCIV[i] = np.float64(FWCIV[i] / 1000);
	    	iband[i] = np.float64(iband[i] + 1.486)        

	    	L1450 = 10**(-.4 * (iband[i] - 4.85))
	    	L1450 = L1450 * 3.826 * 10**(33)
	    	L1450 = L1450 / 10**(44)
		#print(L1450,'fwciv:',FWCIV[i], 'i:',i)
	
	    	result.append(math.log(FWCIV[i]**2 * (L1450**(.53)), 10) + 6.66)
	else:
		result.append(float('nan'))
    		#calc.append(result)
'''
	       
	  ########################################################################
for i in range(0,len(co_zb)):
	co_zb_scaled.append(co_zb[i] / maxCo)
	lum_zb_scaled.append(lum_zb[i] / maxLum)
	x_co_scaled.append(co_zb_scaled[i] * math.sin(math.radians(90 - DEJ[i])) * math.cos(math.radians(RAJ[i])))
	x_lum_scaled.append(lum_zb_scaled[i] * math.sin(math.radians(90 - DEJ[i])) * math.cos(math.radians(RAJ[i])))
	y_co_scaled.append(co_zb_scaled[i] * math.sin(math.radians(90 - DEJ[i])) * math.sin(math.radians(RAJ[i])))
   	y_lum_scaled.append(lum_zb_scaled[i] * math.sin(math.radians(90 - DEJ[i])) * math.sin(math.radians(RAJ[i]))) 
	z_co_scaled.append(co_zb_scaled[i] * math.cos(math.radians(90 - DEJ[i])))
   	z_lum_scaled.append(lum_zb_scaled[i] * math.cos(math.radians(90 - DEJ[i])))

"""col1 = FITS.Column(name='co_zb', format='E', array=co_zb)
col2 = FITS.Column(name='lum_zb', format='E', array=lum_zb)
col3 = FITS.Column(name='co_zb_scaled', format='E', array=co_zb_scaled)
col4 = FITS.Column(name='lum_zb_scaled', format='E', array=lum_zb_scaled)
col5 = FITS.Column(name='x_no_change', format='E', array=x_no_change)
col6 = FITS.Column(name='x_co', format='E', array=x_co)
col7 = FITS.Column(name='x_lum', format='E', array=x_lum)
col8 = FITS.Column(name='y_no_change', format='E', array=y_no_change)
col9 = FITS.Column(name='y_co', format='E', array=y_co)
col10 = FITS.Column(name='y_lum', format='E', array=y_lum)
col11 = FITS.Column(name='z_no_change', format='E', array=z_no_change)
col12 = FITS.Column(name='z_co', format='E', array=z_co)
col13 = FITS.Column(name='z_lum', format='E', array=z_lum)
col14 = FITS.Column(name='x_co_scaled', format='E', array=x_co_scaled)
col15 = FITS.Column(name='x_lum_scaled', format='E', array=x_lum_scaled)
col16 = FITS.Column(name='y_co_scaled', format='E', array=y_co_scaled)
col17 = FITS.Column(name='y_lum_scaled', format='E', array=y_lum_scaled)
col18 = FITS.Column(name='z_co_scaled', format='E', array=z_co_scaled)
col19 = FITS.Column(name='z_lum_scaled', format='E', array=z_lum_scaled)
col20 = FITS.Column(name='BHM', format='E', array=result)


cols = FITS.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20])
tbhdu = FITS.BinTableHDU.from_columns(cols)
prihdr = FITS.Header()
prihdr['OBSERVER'] = 'Ayat Mohammed'
prihdr['COMMENT'] = "new FITS file output"
prihdu = FITS.PrimaryHDU(header=prihdr)
thdulist = FITS.HDUList([prihdu, tbhdu])
thdulist.writeto('DR10Q_v2_results.fits')"""


zipped = zip(x_no_change, y_no_change, z_no_change, x_co, y_co, z_co, x_lum,  y_lum, z_lum, x_co_scaled, y_co_scaled, z_co_scaled, x_lum_scaled, y_lum_scaled, z_lum_scaled, co_zb, lum_zb, co_zb_scaled, lum_zb_scaled, fits_file_data_1['RAJ2000'], fits_file_data_1['DEJ2000'],
        fits_file_data_1['zsp'],fits_file_data_1['imag'])
""",
        fits_file_data_1['DGMI'],
        fits_file_data_1['ALPHA_NU'],fits_file_data_1['FWHM_CIV'],
        fits_file_data_1['REWE_CIV'], fits_file_data_1['FWHM_CIII'],
        fits_file_data_1['REWE_CIII'], fits_file_data_1['FWHM_MGII'],
        fits_file_data_1['REWE_MGII'])""" 
np.savetxt('votable_all_result_votable.csv', zipped, fmt='%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f')
"""
zipped = zip(x_no_change, y_no_change, z_no_change)
np.savetxt('no_change.csv', zipped, fmt='%f,%f,%f')

zipped2 = zip(x_co_scaled, y_co_scaled, z_co_scaled)
np.savetxt('comoving_scaled.csv', zipped2, fmt='%f,%f,%f')

zipped3 = zip(x_lum_scaled, y_lum_scaled, z_lum_scaled)
np.savetxt('lum_scaled.csv', zipped3, fmt='%f,%f,%f')

zipped4 = zip(co_zb, lum_zb, co_zb_scaled, lum_zb_scaled)
np.savetxt('zb.csv', zipped4, fmt='%f,%f,%f,%f')


print(FWCIV[2])
print(iband[0])
"""
#np.savetxt('votableBHM.csv', result, fmt='%f')
"""with open('headr_1.csv', 'wb') as csvfile:
     header_0_w = csv.writer(csvfile, delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
     for index in range(0,len(fits_file_headr_1)):
        header_0_w.writerow([fits_file_headr_1[index]])"""

"""with open('headr_1_keys.csv', 'wb') as csvfile:
     header_0_w = csv.writer(csvfile, delimiter=' ',quotechar='|',quoting=csv.QUOTE_MINIMAL)
     #for index in range(0,len(fits_file_headr_0)):
     header_0_w.writerow([fits_file_headr_1.keys()])"""       

"""with open('data_test.csv', 'wb') as csvfile:
     data_1_w = csv.writer(csvfile, delimiter=' ',quotechar='|',quoting=csv.QUOTE_MINIMAL)
     for index in range(0,len(fits_file_data_1)):
        data_1_w.writerow([fits_file_data_1[index][1],fits_file_data_1[index][2],fits_file_data_1[index][8],fits_file_data_1[index][26],
            fits_file_data_1[index][31], fits_file_data_1[index][32],
            fits_file_data_1[index][33], fits_file_data_1[index][34],
            fits_file_data_1[index][35], fits_file_data_1[index][36],
            fits_file_data_1[index][43], fits_file_data_1[index][44],
            fits_file_data_1[index][45], fits_file_data_1[index][46],
            fits_file_data_1[index][47], fits_file_data_1[index][48],
            fits_file_data_1[index][49], fits_file_data_1[index][72],
            fits_file_data_1[index][73], fits_file_data_1[index][75],
            fits_file_data_1[index][76], fits_file_data_1[index][95],
            fits_file_data_1[index][96], fits_file_data_1[index][108],
            fits_file_data_1[index][109], fits_file_data_1[index][112],
            fits_file_data_1[index][113], fits_file_data_1[index][116],
            fits_file_data_1[index][117], fits_file_data_1[index][120],
            fits_file_data_1[index][121]])

"""

#numpy.savetxt("foo.csv", fits_file_headr_0, delimiter=",")
#my_df = pd.DataFrame(fits_file_headr_0)
#my_df.to_csv('my_csv.csv', index=False, header=False)
#scidata1 = fits_file[1].data


#header1  = fits_file[1].header
#print header1
#scidata2 = fits_file[2].data

#print len(scidata1)
#print len(scidata1[0])
print "\n"

#print scidata1[0][0]
#print scidata1[87821][0]

#print len(scidata2)
#print len(scidata2[0])

#print scidata2[0][0]
#fits_file.close()

