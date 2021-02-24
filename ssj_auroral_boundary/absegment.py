# Copyright 2018 SEDA Group at CU Boulder
# Created by:
# Liam Kilcommons
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
import numpy as np
import logging
from ssj_auroral_boundary import loggername,loggername_poes

class absegment(object):
    """ A class for encapulating one 'segment', which is a region of
    over-threshold flux

    """
    def __init__(self, polarpass, ind_seg_start, ind_seg_end):
        """ Initialize segment class object

        Parameters
        ----------
        polarpass
        ind_seg_start
        ind_seg_end

        """

        self.log = logging.getLogger(loggername + '.' + self.__class__.__name__)
        self.polarpass = polarpass
        self.si = ind_seg_start
        self.ei = ind_seg_end
        self.area = np.nansum(self['intflux']) # Total flux in this segment
        self.area_uncert = np.nanmean(self['total_flux_std']) # Avg Relative Unc
        self.twidth = abs(self['uts'][-1] - self['uts'][0]) #Width in seconds
        self.log.debug(('Segment from UT second '
                        +'%.1f - %.1f ' % (self['uts'][0],self['uts'][-1])
                        +'(%.1f seconds) ' % (self.twidth)
                        +'had area (total flux) %.2g' % (self.area)))

    def __str__(self):
        """Construct the output string
        """
        outstring = "UTS: %.1f - %.1f (DT=" % (self['uts'][0], self['uts'][-1])
        outstring += "%.1f s), area (totflux): %.2g" % (self.twidth, self.area)
        return outstring

    def __getitem__(self,var):
        """Call the parent's getitem for the variable, and then subscript it
        """
        return self.polarpass[var][self.si:self.ei]


class absegment_poes(object):
    """ A class for encapulating one 'segment', which is a region of
    over-threshold flux

    """
    def __init__(self, polarpass, ind_seg_start, ind_seg_end):
        """ Initialize segment class object

        Parameters
        ----------
        polarpass
        ind_seg_start
        ind_seg_end

        """

        self.log = logging.getLogger(loggername_poes + '.' + self.__class__.__name__)
        self.polarpass = polarpass
        self.si = ind_seg_start
        self.ei = ind_seg_end
        self.area = np.nansum(self['intflux']) # Total flux in this segment
        self.area_uncert = np.nanmean(self['total_flux_std']) # Avg Relative Unc
        self.twidth = abs(self['uts'][-1] - self['uts'][0]) #Width in seconds
        self.log.debug(('Segment from UT second '
                        +'%.1f - %.1f ' % (self['uts'][0],self['uts'][-1])
                        +'(%.1f seconds) ' % (self.twidth)
                        +'had area (total flux) %.2g' % (self.area)))

    def __str__(self):
        """Construct the output string
        """
        outstring = "UTS: %.1f - %.1f (DT=" % (self['uts'][0], self['uts'][-1])
        outstring += "%.1f s), area (totflux): %.2g" % (self.twidth, self.area)
        return outstring

    def __getitem__(self,var):
        """Call the parent's getitem for the variable, and then subscript it
        """
        return self.polarpass[var][self.si:self.ei]

