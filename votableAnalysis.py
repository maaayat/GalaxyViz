'''
Astropy's FITS file format analyzer is actually extremely useful.
This file just shows some examples of its use. In the future, it
might be required to work with FITS files, so I will just leave
this here.
'''

import numpy as np
import math
#from matplotlib import pyplot as plt
from astropy.io.votable import parse_single_table
#import astropy
import astropy.cosmology as cosmo
from astropy.cosmology import LambdaCDM
dark_cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
import csv
table_read = parse_single_table("II_326_zcatrev-161102_a.vot").to_table(use_names_over_ids=True)
print dir(table)
print table.pprint #prints out a preview of the table
print table.info # prints out information about the colums in the table
a = table.to_pandas() # convert table to a pandas dataframe and now you can access the dataframe elements



