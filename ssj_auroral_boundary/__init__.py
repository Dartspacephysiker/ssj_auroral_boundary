# Copyright 2018 SEDA Group at CU Boulder
# Created by:
# Liam Kilcommons
# Stolen by:
# Birkeland Centre for Space Science (2021)

# 2021 The Birkeland Centre for Space Science version ...
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

__version__ = str("0.1.1")

#Prefix for all package loggers
loggername = 'ssj_auroral_boundary'

__all__ = ['absatday', 'abpolarpass', 'absegment', 'abcsv', 'files',
           'dmsp_spectrogram']

# Explicitly import all modules (in addition to defining __all__)
from ssj_auroral_boundary import (absatday,
                                    abpolarpass,
                                    absegment,
                                    abcsv,
                                    files,
                                    dmsp_spectrogram)
