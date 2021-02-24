# Copyright 2018 SEDA Group at CU Boulder
# Created by:
# Liam Kilcommons
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
import os
import logging

from geospacepy import special_datetime

from ssj_auroral_boundary import loggername,loggername_poes
from ssj_auroral_boundary.abpolarpass import abpolarpass,abpolarpass_poes
from ssj_auroral_boundary.abcsv import abcsv,abcsv_poes

class absatday(object):
    """Class for one satellite-day of SSJ data (one CDF file)

    Implements __getitem__ interface to access any data marked with
        (from CDF) below.

    Attributes
    ----------
    cdffn : str
        Full file path to CDF file
    satnum : int
        DMSP number from CDF
    log : logging.logger
        Logger for this instance
        Expects a root logger with name 'ssj_auroral_boundary'
        in the calling function, otherwise won't display anything
    writecsv : bool
        Write out the CSV
    csv : abcsv.abcsv
        Instance for adding lines to CSV
    imgdir : str
        Path to write image files to if make_plot is True or plot_failed is True
        (see constructor for defualt behavior and fallbacks)
    make_plot : bool
           Plot successful passes
    plot_failed : bool
        Also plot unsuccessful passes
    cdf : pycdf.CDF
        Open CDF file
    time : np.ndarray, shape=(86400,1)
        Timestamp of SSJ data as datetimes (Universal Time)
    uts : np.ndarray, shape=(86400,1)
        Timestamp of SSJ data as Second of Day (Universal Time)
    hod : np.ndarray, shape=(86400,1)
        Timestamp of SSJ data as Hour of Day (Universal Time)
    mlat : np.ndarray, shape=(86400,1)
        Magnetic latitude AACGM (NOAA CDF v1.1.2) / Apex (CU internal SSJ CDF)
    mlt : np.ndarray, shape=(86400,1)
        Magnetic local time AACGM (NOAA CDF v1.1.2) / Apex (CU internal SSJ CDF)
    counts : np.ndarray, shape=(86400,19)
        Detector counts for all SSJ channels (from CDF)
    diff_flux : np.ndarray, shape=(86400,19)
        Differential electron energy energy flux (from CDF)
    diff_flux_std : np.ndarray, shape=(86400,19)
        Relative error in differential electron energy flux (from CDF)
    total_flux : np.ndarray, shape=(86400,1)
        Total energy flux intergrated across all 19 channels (from CDF)
    total_flux_std : np.ndarray, shape=(86400,1)
        Relative error in total electron energy flux (from CDF)
    channel_energies : np.ndarray, shape=(19,)
        Energy in electron volts (from CDF)
    xings : list
        Indices into data arrays of equator crossings
    polarpasses : list
        List of abpolarpass.abpolarpass obj for each polar crossing (half orbit)
    """
    def __init__(self, cdffile, imgdir=None, make_plot=True, plot_failed=False,
                 csvdir=None, writecsv=True, csvvars=['mlat', 'mlt']):
        """Constructor for absatday

        Parameters
        ----------
        cdffile : str
            DMSP SSJ CDF file (probably from NASA CDAWeb)
        imgdir : str, optional
            Path to dump boundary identification images to (must exist)
            If None looks for environment variable DMSP_DIR_ABIMG
        make_plot : bool, optional
            Plot of each successful identification (the default is True)
        plot_failed : bool, optional
            Also plot unsuccesful passes (the default is False)
        csvdir : str, optional
            Directory to dump CSV files to (must exist)
            If None looks for environment variable DMSP_DIR_ABCSV
            If still fails raises RuntimeError
        writecsv : bool, optional
            Write a CSV file of boundary identifications (the default is True)
        csvvars : list, optional
            List of optional variables to include in each line of the CSV file.
            See abcsv for more details.  (default=['mlat', 'mlt'])

        """

        try:
            from spacepy import pycdf
        except Exception as e:
            print('Failed to import spacepy.pycdf; likely CDF C library was not found')
            print('this is a fatal error unless you are building documentation')
            print(e)

        self.log = logging.getLogger(loggername+'.'+self.__class__.__name__)
        self.cdf = pycdf.CDF(cdffile)
        #Parse out spacecraft so we know how to handle J4/J5 differences
        if 'dmsp' in cdffile:
            # Get the spacecraft number from the filename
            self.satnum = int(os.path.split(cdffile)[-1].split('dmsp-f')[-1][:2])
            self.log.info("Satellite number determined to be "
                          + "{:d}".format(self.satnum))
        else:
            raise RuntimeError(('Unexpected CDF filename {:s}, '.format(cdffile)
                               +'could not parse out DMSP number' ))
        self.cdffn = cdffile
        self.make_plot = make_plot # Make plots of passes T/F
        self.plot_failed = plot_failed #Plot failed identifications also T/F
        self.writecsv = writecsv # Write pass identifications to a file
        self.time = self.cdf['Epoch'][:]
        self.uts = special_datetime.datetimearr2sod(self.time)
        self.hod = self.uts/3600.
        self.diff_flux = self.cdf['ELE_DIFF_ENERGY_FLUX'][:]
        self.diff_flux_std = self.cdf['ELE_DIFF_ENERGY_FLUX_STD'][:]
        self.total_flux = self.cdf['ELE_TOTAL_ENERGY_FLUX'][:]
        #The uncertainty in the CDF is relative
        self.total_flux_std = self.cdf['ELE_TOTAL_ENERGY_FLUX_STD'][:]

        #Handle filtering out any data without enough counts
        countthresh = 2.
        self.counts = (self.cdf['ELE_COUNTS_OBS'][:] \
                        -self.cdf['ELE_COUNTS_BKG'][:])

        #Zero out any dubious fluxes
        self.diff_flux[self.counts <= countthresh] = 0.0

        latvar,ltvar = 'SC_APEX_LAT','SC_APEX_MLT'
        if latvar not in self.cdf or ltvar not in self.cdf:
            #v1.1.3
            self.log.warning(('Unable to find APEX latitude or local '
                       + 'time variables in CDF file. Falling '
                       +'back to AACGM magnetic coordinates'))
            latvar,ltvar = 'SC_AACGM_LAT','SC_AACGM_LTIME'

        self.mlat = self.cdf[latvar][:]
        self.mlt = self.cdf[ltvar][:]
        self.channel_energies = self.cdf['CHANNEL_ENERGIES'][:]
        self.xings = self.simple_passes(self.mlat)
        self.polarpasses = []

        #Look for environemnt variables to define paths if no paths provided
        imgdir = self.if_none_use_envvar(imgdir,'DMSP_DIR_ABIMG')
        if imgdir is None:
            raise RuntimeError('No image dir passed & no '
                                           + 'DMSP_DIR_ABIMG envvar')
        self.imgdir = imgdir

        csvdir = self.if_none_use_envvar(csvdir, 'DMSP_DIR_ABCSV')
        if csvdir is None:
            raise RuntimeError('No csv dir passed & no '
                               + 'DMSP_DIR_ABCSV envvar')

        cdffn_noext = os.path.splitext(os.path.split(cdffile)[-1])[0]
        csvfile = cdffn_noext + '_boundaries.csv'
        self.csv = abcsv(csvdir, csvfile, cdffile, csvvars=csvvars,
                         writecsv=self.writecsv)

        #Start processing the polar passes one by one
        for i in range(len(self.xings)-1):
            newpass = abpolarpass(self, self.xings[i],
                                              self.xings[i+1]-1)
            self.polarpasses.append(newpass)

    def if_none_use_envvar(self,checkvar,envvar):
        """Check for environment variable envvar if checkvar is None
        """
        if checkvar is None and envvar in os.environ:
                return os.environ[envvar]
        elif checkvar is None:
            return None
        else:
            return checkvar

    def simple_passes(self,latitude):
        """Finds all the equator crossings"""
        npts = len(latitude.flatten())
        entered_north = []
        entered_south = []

        for k in range(1,npts):
            #poleward crossing
            if latitude[k-1] < 0. and latitude[k] >= 0.:
                entered_north.append(k)
                self.log.info("Entered Northern Hemisphere: ind"
                              + ":%d,lat:%.3f" % (k, latitude[k]))
            elif latitude[k-1] > 0. and latitude[k] <= 0.:
                entered_south.append(k)
                self.log.info("Entered Southern Hemisphere: ind"
                              + ":%d,lat:%.3f" % (k, latitude[k]))

        xings = entered_north+entered_south
        xings.sort()
        return xings

    def __getitem__(self,var):
        if hasattr(self,var):
            return getattr(self,var)
        elif var in self.cdf:
            return self.cdf[var][:]
        else:
            self.log.error(("Non-existent variable %s " % (str(var))
                            + "requested through getattr. Returning None"))
            return None


class absatday_poes(object):
    """Class for one satellite-day of TED data (one NC file)

    Implements __getitem__ interface to access any data marked with
        (from NC) below.

    Attributes
    ----------
    ncfn : str
        Full file path to NC file
    satnum : int
        DMSP number from NC
    log : logging.logger
        Logger for this instance
        Expects a root logger with name 'ted_auroral_boundary'
        in the calling function, otherwise won't display anything
    writecsv : bool
        Write out the CSV
    csv : abcsv.abcsv_poes
        Instance for adding lines to CSV
    imgdir : str
        Path to write image files to if make_plot is True or plot_failed is True
        (see constructor for defualt behavior and fallbacks)
    make_plot : bool
           Plot successful passes
    plot_failed : bool
        Also plot unsuccessful passes
    cdf : pycdf.NC
        Open NC file
    time : np.ndarray, shape=(86400,1)
        Timestamp of TED data as datetimes (Universal Time)
    uts : np.ndarray, shape=(86400,1)
        Timestamp of TED data as Second of Day (Universal Time)
    hod : np.ndarray, shape=(86400,1)
        Timestamp of TED data as Hour of Day (Universal Time)
    mlat : np.ndarray, shape=(86400,1)
        Magnetic latitude AACGM (POES/MetOp NC v1.1.2) / Apex (CU internal TED NC)
    mlt : np.ndarray, shape=(86400,1)
        Magnetic local time AACGM (POES/MetOp NC v1.1.2) / Apex (CU internal TED NC)
    counts : np.ndarray, shape=(86400,19)
        Detector counts for all TED channels (from NC)
    diff_flux : np.ndarray, shape=(86400,19)
        Differential electron energy energy flux (from NC)
    diff_flux_std : np.ndarray, shape=(86400,19)
        Relative error in differential electron energy flux (from NC)
    total_flux : np.ndarray, shape=(86400,1)
        Total energy flux intergrated across all 19 channels (from NC)
    total_flux_std : np.ndarray, shape=(86400,1)
        Relative error in total electron energy flux (from NC)
    channel_energies : np.ndarray, shape=(19,)
        Energy in electron volts (from NC)
    xings : list
        Indices into data arrays of equator crossings
    polarpasses : list
        List of abpolarpass.abpolarpass obj for each polar crossing (half orbit)
    """
    def __init__(self, ncfile, imgdir=None, make_plot=True, plot_failed=False,
                 csvdir=None, writecsv=True, csvvars=['mlat', 'mlt']):
        """Constructor for absatday

        Parameters
        ----------
        ncfile : str
            POES/MetOp TED NC file (probably from NASA CDAWeb)
        imgdir : str, optional
            Path to dump boundary identification images to (must exist)
            If None looks for environment variable POES_DIR_ABIMG
        make_plot : bool, optional
            Plot of each successful identification (the default is True)
        plot_failed : bool, optional
            Also plot unsuccesful passes (the default is False)
        csvdir : str, optional
            Directory to dump CSV files to (must exist)
            If None looks for environment variable POES_DIR_ABCSV
            If still fails raises RuntimeError
        writecsv : bool, optional
            Write a CSV file of boundary identifications (the default is True)
        csvvars : list, optional
            List of optional variables to include in each line of the CSV file.
            See abcsv_poes for more details.  (default=['mlat', 'mlt'])

        """

        try:
            import xarray as xr
        except Exception as e:
            print('Failed to import xarray')
            print('this is a fatal error unless you are building documentation')
            print(e)

        self.log = logging.getLogger(loggername_poes+'.'+self.__class__.__name__)
        self.ds = xr.open_dataset(ncfile)

        # self.cdf = pycdf.NC(ncfile)

        #Parse out spacecraft so we know how to handle J4/J5 differences
        if 'poes' in ncfile:
            # Get the spacecraft number from the filename
            self.satnum = int(os.path.split(ncfile)[-1].split('poes_')[-1][1:3])
            self.log.info("Satellite number determined to be "
                          + "{:d}".format(self.satnum))
        else:
            raise RuntimeError(('Unexpected NC filename {:s}, '.format(ncfile)
                               +'could not parse out POES number' ))

        import pandas as pd
        import numpy as np
        from datetime import datetime,timedelta
        import apexpy
        from pyamps.mlt_utils import mlon_to_mlt
        apexrefh = 110


        self.ncfn = ncfile
        self.make_plot = make_plot # Make plots of passes T/F
        self.plot_failed = plot_failed #Plot failed identifications also T/F
        self.writecsv = writecsv # Write pass identifications to a file

        # self.time = self.cdf['Epoch'][:]
        # self.time = pd.DatetimeIndex([datetime(yr,1,1)+timedelta(days=day-1,milliseconds=msec) for yr,day,msec in zip(self.ds['year'],
        #                                                                                                               self.ds['day'],
        #                                                                                                               self.ds['msec'])])
        self.time = pd.DatetimeIndex([datetime(yr,1,1)+timedelta(days=day-1,milliseconds=msec) for yr,day,msec in \
                                      zip(self.ds['year'].values,
                                          np.float64(self.ds['day'].values),
                                          np.float64(self.ds['msec'].values))])

        self.uts = special_datetime.datetimearr2sod(self.time)
        self.hod = self.uts/3600.

        # self.diff_flux = self.cdf['ELE_DIFF_ENERGY_FLUX'][:]
        # self.diff_flux_std = self.cdf['ELE_DIFF_ENERGY_FLUX_STD'][:]
        # self.total_flux = self.cdf['ELE_TOTAL_ENERGY_FLUX'][:]

        self.total_flux = self.ds['ted_ele_tel0_hi_eflux'].values
        # self.diff_flux_std = self.ds['ted_ele_tel0_hi_eflux_error'].values 

        #The uncertainty in the NC is relative
        # ted_ele_tel0_hi_eflux_error : 'TED electron (1-20 keV) 0 deg telescope energy flux percent error'
        # self.total_flux_std = self.cdf['ELE_TOTAL_ENERGY_FLUX_STD'][:]
        self.total_flux_std = self.ds['ted_ele_tel0_hi_eflux_error'].values

        # Apex coordinates
        apexreftime = self.time[0]
        a = apexpy.Apex(date=apexreftime,refh=apexrefh)

        #Assume that poes netCDF variables 'lat' and 'lon' are geodetic, not geocentric
        self.mlat, self.mlon = a.geo2apex(self.ds['lat'].values,
                                          self.ds['lon'].values,
                                          self.ds['alt'].values)
        self.mlt = mlon_to_mlt(self.mlon, self.time, apexreftime.year)

        
        #Handle filtering out any data without enough counts
        print("WARNING! KC algorithm zeroes out low counts. What do we do?")
        print("countthresh = 2.")
        print("self.counts = (self.cdf['ELE_COUNTS_OBS'][:] ")
        print("               -self.cdf['ELE_COUNTS_BKG'][:])")


        #Zero out any dubious fluxes
        print("")
        print("self.diff_flux[self.counts <= countthresh] = 0.0")

        # latvar,ltvar = 'SC_APEX_LAT','SC_APEX_MLT'
        # if latvar not in self.cdf or ltvar not in self.cdf:
        #     #v1.1.3
        #     self.log.warning(('Unable to find APEX latitude or local '
        #                + 'time variables in NC file. Falling '
        #                +'back to AACGM magnetic coordinates'))
        #     latvar,ltvar = 'SC_AACGM_LAT','SC_AACGM_LTIME'

        # self.mlat = self.cdf[latvar][:]
        # self.mlt = self.cdf[ltvar][:]

        # self.channel_energies = self.cdf['CHANNEL_ENERGIES'][:]
        print("WARNING! What did KC algorithm want with 'CHANNEL_ENERGIES'?")

        self.xings = self.simple_passes(self.mlat)
        self.polarpasses = []

        #Look for environemnt variables to define paths if no paths provided
        imgdir = self.if_none_use_envvar(imgdir,'POES_DIR_ABIMG')
        if imgdir is None:
            raise RuntimeError('No image dir passed & no '
                                           + 'POES_DIR_ABIMG envvar')
        self.imgdir = imgdir

        csvdir = self.if_none_use_envvar(csvdir, 'POES_DIR_ABCSV')
        if csvdir is None:
            raise RuntimeError('No csv dir passed & no '
                               + 'POES_DIR_ABCSV envvar')

        ncfn_noext = os.path.splitext(os.path.split(ncfile)[-1])[0]
        csvfile = ncfn_noext + '_boundaries.csv'
        self.csv = abcsv_poes(csvdir, csvfile, ncfile, csvvars=csvvars,
                              writecsv=self.writecsv)

        #Start processing the polar passes one by one
        for i in range(len(self.xings)-1):
            newpass = abpolarpass_poes(self, self.xings[i],
                                              self.xings[i+1]-1)
            self.polarpasses.append(newpass)

    def if_none_use_envvar(self,checkvar,envvar):
        """Check for environment variable envvar if checkvar is None
        """
        if checkvar is None and envvar in os.environ:
                return os.environ[envvar]
        elif checkvar is None:
            return None
        else:
            return checkvar

    def simple_passes(self,latitude):
        """Finds all the equator crossings"""
        npts = len(latitude.flatten())
        entered_north = []
        entered_south = []

        for k in range(1,npts):
            #poleward crossing
            if latitude[k-1] < 0. and latitude[k] >= 0.:
                entered_north.append(k)
                self.log.info("Entered Northern Hemisphere: ind"
                              + ":%d,lat:%.3f" % (k, latitude[k]))
            elif latitude[k-1] > 0. and latitude[k] <= 0.:
                entered_south.append(k)
                self.log.info("Entered Southern Hemisphere: ind"
                              + ":%d,lat:%.3f" % (k, latitude[k]))

        xings = entered_north+entered_south
        xings.sort()
        return xings

    def __getitem__(self,var):
        if hasattr(self,var):
            return getattr(self,var)
        elif var in self.cdf:
            return self.cdf[var][:]
        else:
            self.log.error(("Non-existent variable %s " % (str(var))
                            + "requested through getattr. Returning None"))
            return None
