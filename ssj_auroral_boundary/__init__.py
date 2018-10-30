# Copyright 2018 SEDA Group at CU Boulder
# Created by: 
# Liam Kilcommons 
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
"""
ssj_auroral_boundary
--------------------

Figure of Merit boundary identification for DMSP SSJ5

Modules
-------
absatday
abpolarpass
absegment
abscv
files
dmsp_spectrogram

"""
from __future__ import print_function

__version__ = "0.1"

#Prefix for all package loggers
loggername = 'ssj_auroral_boundary'

# Import all the modules
__all__ = ['absatday', 'abpolarpass', 'absegment', 'abcsv', 'files',
           'dmsp_spectrogram']

try:
    from ssj_auroral_boundary import __all__
except ImportError as err:
    print("problem importing {:s}: {:}".format(loggername, err))
