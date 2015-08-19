# -*- coding: utf-8 -*-
"""
Turning Points, Rainflow Counting and Equivalent Load
=====================================================

sig2ext and rfc_astm are from Nieslony's rainflow counting algorithm. Originally
written for Matlab and freely available under a BSD license at the
Matlab Central: http://www.mathworks.com/matlabcentral/fileexchange/3026

findtp, _findcross are extracted from the wafo toolbox, which is available
at: http://code.google.com/p/pywafo/

The fatigue damage is calculated via the Palmgren-Miner rule and is used to
express an equivalent load.

Created on Wed Feb 16 16:53:18 2011
@author: David Verelst
"""

from time import time
import os
import array
import numpy as np
import scipy.io
import rainflowlib
import warnings
from matplotlib import pyplot as plt
#import wafo

class LoadResults:
    """
    Read a HAWC2 result data file
    =============================

    sig = LoadResults(file_path, file_name)

    This class is called like a function:
    LoadResults() will read the specified file upon object initialization.

    Parameters
    ----------
    file_path : string
        Specify the path to the results file folder.

    file_name : string
        HAWC2 result file name, extension .htc, .sel, .dat, .log are removed.

    debug : boolean, optional
        False by default. If set to True, some addiation information will be
        printed while running.

    Member variables
    ----------------
    sig.sig[timeStep,channel] : array-like
        complete result file in a numpy array

    sig.ch_details[channel,(0=ID; 1=units; 2=description)] : array-like

    sig.error_msg : array-like
        is 'none' if everything went OK, otherwise it holds the error(s)

    Created by David Verelst
    """

    # start with reading the .sel file, containing the info regarding
    # how to read the binary file and the channel information
    def __init__(self, file_path, file_name, debug=False):

        self.debug = debug

        # timer in debug mode
        if self.debug:
            start = time()

        self.file_path = file_path
        # remove .log, .dat, .sel extensions who might be accedental left
        file_name = file_name.replace('.htc', '')
        file_name = file_name.replace('.sel', '')
        file_name = file_name.replace('.dat', '')
        file_name = file_name.replace('.log', '')
        self.file_name = file_name
        self.read_sel()
        # continue if the file has been succesfully read
        if self.error_msg == 'none':
            # load the channel id's and scale factors
            scale_factors = self.data_sel()
            # read the binary file
            if self.FileType == 'BINARY':
                self.read_bin(scale_factors)
            # read the ASCII file
            elif self.FileType == 'ASCII':
                self.read_ascii()
            else:
                print '========================================================'
                print 'unknown file type: ' + self.FileType
                print '========================================================'
                self.error_msg = 'error: unknown file type'
                self.sig = []

        if self.debug:
            stop = time() - start
            print 'time to load HAWC2 file:', stop, 's'

    def read_sel(self):
        # anticipate error on file reading
        try:
            # open file, read and close
            go_sel = self.file_path + self.file_name + '.sel'
            FILE = open(go_sel, "r")
            self.lines = FILE.readlines()
            FILE.close()
            self.error_msg = 'none'

        # error message if the file does not exists
        except:
            # print 26*' ' + 'ERROR'
            print 50*'='
            print self.file_path
            print self.file_name + '.sel could not be found'
            print 50*'='
            self.error_msg = 'error: file not found'

    def data_sel(self):
        # increase precision
        # D.getcontext().prec = 50

        # scan through all the lines in the file
        line_nr = 1
        # channel counter for ch_details
        ch = 0
        for line in self.lines:
            # on line 9 we can read following paramaters:
            if line_nr == 9:
                # remove the end of line character
                line = line.replace('\n','')

                settings = line.split(' ')
                # delete all empty string values
                for k in range(settings.count('')):
                    settings.remove('')

                # and assign proper values with correct data type
                self.N = int(settings[0])
                self.Nch = int(settings[1])
                self.Time = float(settings[2])
                # there are HEX values at the end of this line...
                # On Linux they will show up in the last variable, so don't inc
                if os.name == 'posix':
                    nrchars = len(settings[3])-1
                elif os.name == 'nt':
                    nrchars = len(settings[3])
                else:
                    raise UserWarning, \
                    'Untested platform:', os.name
                settings[3] = settings[3][0:nrchars]
                self.FileType = settings[3]
                self.Freq = self.N/self.Time

                # prepare list variables
                self.ch_details = np.ndarray(shape=(self.Nch,3),dtype='<U100')
                # it seems that float64 reeds the data correctly from the file
                scale_factors = scipy.zeros(self.Nch, dtype='Float64')
                #self.scale_factors_dec = scipy.zeros(self.Nch, dtype='f8')
                i = 0

            # starting from line 13, we have the channels info
            if line_nr > 12:
                # read the signal details
                if line_nr < 13 + self.Nch:
                    # remove leading and trailing whitespaces from line parts
                    self.ch_details[ch,0] = line[12:43].strip() # chID
                    self.ch_details[ch,1] = line[43:54].strip() # chUnits
                    self.ch_details[ch,2] = line[54:-1].strip() # chDescr
                    ch += 1
                # read the signal scale parameters for binary format
                elif line_nr > 14 + self.Nch:
                    scale_factors[i] = line
                    # print scale_factors[i]
                    #self.scale_factors_dec[i] = D.Decimal(line)
                    i = i + 1
                # stop going through the lines if at the end of the file
                if line_nr == 2*self.Nch + 14:
                    self.scale_factors = scale_factors

                    if self.debug:
                        print 'N       ', self.N
                        print 'Nch     ', self.Nch
                        print 'Time    ', self.Time
                        print 'FileType', self.FileType
                        print 'Freq    ', self.Freq
                        print 'scale_factors', scale_factors.shape

                    return scale_factors
                    break

            # counting the line numbers
            line_nr = line_nr + 1

    def read_bin(self, scale_factors):
        # if there is an error reading the binary file (for instance if empty)
        try:
            # read the binary file
            go_binary = self.file_path + self.file_name + '.dat'
            FILE = open(go_binary, mode='rb')

            # create array, put all the binary elements as one long chain in it
            binvalues = array.array('h')
            binvalues.fromfile(FILE, self.N * self.Nch)
            FILE.close()
            # convert now to a structured numpy array
            # sig = np.array(binvalues, np.float)
            sig = np.array(binvalues)

            if self.debug: print self.N, self.Nch, sig.shape

            # reshape the array to 2D and transpose (Fortran to C array)
            sig = sig.reshape((self.Nch, self.N)).T

            # create diagonal vector of size (Nch,Nch)
            dig = np.diag(scale_factors)
            # now all rows of column 1 are multiplied with dig(1,1)
            sig = np.dot(sig,dig)
            self.sig = sig
            # 'file name;' + 'lnr;msg;'*(len(MsgList)) + '\n'
        except:
            self.sig = []
            self.error_msg = 'error: reading binary file failed'
            print '========================================================'
            print self.error_msg
            print '========================================================'

    def read_ascii(self):

        try:
            go_ascii = self.file_path + self.file_name + '.dat'
            self.sig = np.fromfile(go_ascii, dtype=np.float32, sep='  ')
            self.sig = self.sig.reshape((self.N, self.Nch))
        except:
            self.sig = []
            self.error_msg = 'error: reading ascii file failed'
            print '========================================================'
            print self.error_msg
            print '========================================================'

class Fatigue:
    """
    Fatigue Damage
    ==============

    Wrapper class for calculating turning points (sig2ext or findtp), counting
    rainflow cycles (rfc_astm), dividing them into bins (rfc_hist) and calculate
    the accumulated damage with the Palmgren-Miner rule.

    The output of the fatigue analysis is the amplitude that, given the number
    of cycles as defined in neq or EQ number, results in the same amount of
    fatigue damage according to the Palmgren-Miner rule.

    Paramters
    ---------

    clsn : int, optional
        Devide sig into classes before searching for the extremes

    nrbins : int, optional
        Divide the rainflow counted amplitudes in a number of equally spaced
        bins.

    tp_method : {'wafo', 'nieslony'}, optional
        The turning points can be determined by Nieslony's or wafo's approach.
        Both give the same result although the wafo method is faster.

    material: float or array-like, optional
        Material parameter for the Palmgren-Miner damage accumulation.

    neq : int, optional
        Reference number of cycles for the Palmgren-Miner damage accumulation.

    Returns
    -------

    S1 : array-like or float
        Amplitude of the load that will result in a similar amount of fatigue
        damage (according to Palmgren-Miner) considering neq number of cycles.

    Example
    -------

    >>> f = Fatigue()
    >>> S1 = f(sig.sig[:,12])

    """
    def __init__(self, tp_method='wafo', clsn=None, material=12., nrbins=46,
                 neq = 1.):
        """
        """
        # default parameters
        self.tp_method = tp_method
        self.clsn = clsn
        self.material = material
        self.neq = neq
        self.nrbins = nrbins

    def __call__(self, sig):
        """
        """
        self.sig_rf = rainflow(sig, tp_method=self.tp_method, clsn=self.clsn)

        self.hist, self.bin_edges, self.bin_avg = \
                    rfc_hist(self.sig_rf, nrbins=self.nrbins)

        S1 = palmgren_mimer(self.hist, self.bin_avg,
                                 material=self.material, neq=self.neq)

        return S1

def rfc_hist(sig_rf, nrbins=46):
    """
    Histogram of rainflow counted cycles
    ====================================

    hist, bin_edges, bin_avg = rfc_hist(sig, nrbins=46)

    Divide the rainflow counted cycles of a signal into equally spaced bins.
    Combine

    Parameters
    ----------

    sig_rf : array-like
        As outputted by rfc_astm or rainflow

    nrbins : int, optional
        Divide the rainflow counted amplitudes in a number of equally spaced
        bins.

    Output
    ------

    hist : array-like
        Counted rainflow cycles per bin, has nrbins elements

    bin_edges : array-like
        Edges of the bins, has nrbins+1 elements.

    bin_avg : array-like
        Average rainflow cycle amplitude per bin, has nrbins elements.
    """

    # convert to half cycles
    # select only full cycles
    i_full = (sig_rf[:,2] == 1.0).nonzero()[0]
    # and duplicate the full cycles (2*half = 1*full)
    rf_half = np.r_[sig_rf[:,0], sig_rf[i_full,0]]

#    rf_half = sig_rf

    # the Matlab approach is to divide into 46 bins
    bin_edges = np.linspace(0,1, num=nrbins+1)*rf_half.max()
    hist = np.histogram(rf_half, bins=bin_edges)[0]
    # calculate the average per bin
    hist_sum = np.histogram(rf_half, weights=rf_half, bins=bin_edges)[0]
    # replace zeros with one, to avoid 0/0
    hist_ = hist.copy()
    hist_[(hist == 0).nonzero()] = 1.0
    # since the sum is also 0, the avg remains zero for those whos hist is zero
    bin_avg = hist_sum / hist_

    return hist, bin_edges, bin_avg


def palmgren_mimer(hist, bin_avg, material=12., neq=1.):
    """
    Fatigue Damage With Palmgren-Miner Rule
    =======================================

    S1 = palmgren_mimer(hist, bin_avg, material=12., neq=1.)

    Paramters
    ---------

    hist : array-like
        Number of half cycles per bin, as given by rfc_hist().

    bin_avg : array-like
        Average value for each bin of the rainflow cycle amplitudes.

    material: float or array-like, optional
        Material parameter for the Palmgren-Miner damage accumulation.

    neq : int, optional
        Reference number of cycles for the Palmgren-Miner damage accumulation.

    Output
    ------

    S1 : array-like or float
        Amplitude of the load that will result in a similar amount of fatigue
        damage (according to Palmgren-Miner) considering neq number of cycles.

    """
    # calculate the equivalent load with palmgren_mimer
    # for the range of the amplitudes, take the average of each bin
    # since hist counts the number of half cycles, divide by 2 to convert back
    # to full cycles
    S1 = np.power(np.sum(0.5*hist*np.power(bin_avg, material))/neq, 1./material)
    return S1


def rainflow(sig, tp_method='wafo', clsn=None):
    """
    sig_rfc = rainflow(sig, tp_method='wafo', clsn=None)

    Convinience function for turning points sig2ext or findtp and rfc_astm

    Paramters
    ---------

    sig : array-like
        Time history of loading. Can be one-dimensional (no given time history)
        or two-dimensional (first column is time, second column is the signal).
        All other columns will be ignored. Time is also ignore in combination
        with tp_method=wafo.

    tp_method : {'wafo', 'nieslony'}, optional
        The turning points can be determined by Nieslony's or wafo's approach.
        Both give the same result although the wafo method is faster.

    clsn : int
        Input for the Nieslony's determination of the turning points in sig2ext.
        Devide sig into classes before searching for the extremes.

    Output
    ------

    sig_rfc : array-like
        rainflow cycles: array nx3
        sig_rfc[:,0] Cycles amplitude
        sig_rfc[:,1] Cycles mean value
        sig_rfc[:,2] Number of cycles (0.5 or 1.0)
        If time is given in sig, sig_rfc has two additional columns:
        sig_rfc[:,3] Begining time (when input includes dt or extt data)
        sig_rfc[:,4] Cycle period (when input includes dt or extt data)

    """
    # check input data validity
    if not type(sig).__name__ == 'ndarray':
        raise TypeError, 'sig should be ndarray, not: ' + type(sig).__name__

    elif len(sig.shape) not in (1,2):
        raise TypeError, 'sig should be 1D or 2D, not: ' + str(len(sig.shape))

    dt = None
    if len(sig.shape) == 2:
        if sig.shape[1] > 1:
            # time on first column, signal on second, ignore all others
            dt = sig[:,0]
            sig = sig[:,1]

    # find turning points
    if tp_method == 'wafo':
        # give warning if combined with clsn and/or dt
        if dt is not None:
            raise Warning, 'tp_method=wafo does not calculate timings of the '\
                + 'turning points.'
        if clsn is not None:
            raise Warning, 'tp_method=wafo has no support for the clsn.'

        ext = sig[findtp(sig)]
        # the wafo method does calculate the timings of the turning points
        exttime = None
    else:
        ext, exttime = sig2ext(sig, dt=dt, clsn=clsn)

    sig_rfc = rfc_astm(ext, exttime=exttime)

    return sig_rfc


def rfc_astm(sig_tp, exttime=None):
    """
    Rainflow Counting by Nieslony
    =============================

    sig_rfc = rfc_astm(sig_tp, exttime=None)

    Paramters
    ---------

    sig_tp : array-like, 1D array holding the signals turning points

    exttime : array-like
        signal time, vector nx1, exact time of occurrence of turning points.

    Returns
    -------

    sig_rfc : array-like
        rainflow cycles: array nx3
        sig_rfc[:,0] Cycles amplitude
        sig_rfc[:,1] Cycles mean value
        sig_rfc[:,2] Number of cycles (0.5 or 1.0)
    If exttime is given, sig_rfc has two additional columns :
        sig_rfc[:,3] Begining time (when input includes dt or extt data)
        sig_rfc[:,4] Cycle period (when input includes dt or extt data)

    References
    ----------

    Adam Niesłony, “Determination of fragments of multiaxial service loading
    strongly influencing the fatigue of machine components,”
    Mechanical Systems and Signal Processing 23, no. 8 (2009): 2712-2721.

    and is based on the following standard:
    ASTM E 1049-85 (Reapproved 1997), Standard practices for cycle counting in
    fatigue analysis, in: Annual Book of ASTM Standards, vol. 03.01, ASTM,
    Philadelphia, 1999, pp. 710–718.

    Copyright (c) 1999-2002 by Adam Nieslony
    Ported to Python by David Verelst
    """
    # in case the rainflow library failes, return None
    if rainflowlib is None:
        return None

    # check the input data

    # if exttime is not defined, use rf3
    if exttime is None:
        sig_rfc, cnr = rainflowlib.rf3(sig_tp)

    # exttime is an array with all the time occurences
    elif type(exttime).__name__ in ('ndarray'):
        sig_rfc, cnr = rainflowlib.rf5(sig_tp, exttime)
    else:
        raise TypeError, 'exttime should be ndarray, not: ' \
                + type(exttime).__name__

    # the sig_rfc was constructed too big in rainflow.rf3, so
    # reduce the sig_rfc array as done before by the mx and mex c function
    n = len(sig_rfc)
    # sig_rfc2 = sig_rfc[0:n-cnr[0],:]
    sig_rfc = sig_rfc.__getslice__(0,n-cnr[0])
    return sig_rfc


def sig2ext(sig, dt=None, clsn=None, debug=False):
    """
    SIG2EXT
    =======

    sig2ext(sig, dt=None, clsn=None)

    Search for local extrema in a time history signal

    Parameters
    ----------

    sig : array-like, time history of the loading

    dt : int, float or array-like, optional

        sampling or if array-like timings corresponding to sig (equal lenghts)

    clsn : int, optional
        Devide sig into classes before searching for the extremes

    Output
    ------

    ext : array-like
        Turning points of the input signal sig

    exttime : array-like or None
        Occurence for each turning point in ext. None if dt isn't defined

    By Adam Nieslony
    Revised, 10-Nov-2009
    Visit the MATLAB Central File Exchange for latest version

    Copyright (c) 1999-2002 by Adam Nieslony
    Ported to Python by David Verelst
    """

    # NOTE: array.__getslice__(a,b) is more efficient compared to array[a:b]:
    # it is twice as fast. See debug.py for the results and code.

    # Should there be a time analysis as well? Do not assume default value
    # for time step, but force to specify dt
    if dt is None:
        TimeAnalyse = False
    elif type(dt).__name__ in ('int', 'float', 'float64', 'float32', 'ndarray'):
        TimeAnalyse = True
    else:
        msg = 'dt should be either int, float or ndarray and not' \
            + type(dt).__name__
        raise TypeError, msg
    # make sure the signal is 1D or has one column
    if type(sig).__name__ in ('ndarray'):
        if len(sig.shape) == 1:
            pass
        elif len(sig.shape) == 2:
            if sig.shape[1] > 1:
                msg = 'sig should be 1D or have 1 column! Wrong shape:' \
                        + str(sig.shape)
                raise TypeError, msg
        else:
            msg = 'sig should be 1D or have 1 column! Wrong shape:' \
                        + str(sig.shape)
            raise TypeError, msg

    # Force to 1D matrix
    sig = sig.ravel()

    # If clsn is given, divide into classes
    if clsn is not None:
        clsn = round(clsn)

        smax = np.max(sig)
        smin = np.min(sig)

        sig = clsn*((sig-smin)/(smax-smin))
        sig = np.fix(sig)
        # find all values in sig equal to clsn
        # return array with 0 where sig is not clsn and clsn on the rest
        temp = sig.__eq__(clsn)
        # take all the non zero elements only and do -1
        sig[temp.nonzero()] = clsn - 1.0
        # some kind of scaling...
        sig = (smax-smin)/(clsn-1)*sig+smin

    # Create a binary vector where 1 means equality of the extreme or,
    # It is recognized that the first and last point is an ext
    w1 = np.diff(sig)
    # Multiply elements n and n+1, evaluate if they are <= 0
    w11 = (w1[0:len(w1)-1]*w1[1:len(w1)]).__le__(0)
    # and get the indices of the elements who comply the <=0 condition
    # nonzero() returns a tuple for each dimension
    w = np.r_[1, w11, 1]
    wi = w.nonzero()[0]
    ext = sig[wi]

    if TimeAnalyse:
        if type(dt).__name__ in ('int', 'float', 'float64', 'float32'):
            # evaluates to true if w equals to one
            # return indices of the true elements, if false, return zero
            exttime = w.__eq__(1).nonzero()[0]*dt
            # no need to do -1 as in Matlab, because indices are 0...n
        else:
            exttime = dt[wi]

    # Removing triple values
    w1 = np.diff(ext)
    p1 = w1[0:len(w1)-1].__eq__(0.0)
    p2 = w1[1:len(w1)].__eq__(0.0)

    # both p1 and p2 should be equal to each other
    w11 = p1*p2

    # using the slice method is 2x as fast, but it is very fast anyway
    w = scipy.zeros((len(w11)+2), dtype=bool)
    w.__setitem__(0,False)
    w.__setitem__(len(w)-1,False)
    w.__setslice__(1,len(w)-1,w11)
    w = w.__invert__()
    # or in a compact way this looks like
    # w = np.array(np.r_[0, w11, 0],dtype=bool).__invert__()

    # and back to indices instead of booleans
    w = w*1
    # select indices of the ones who comply the condtion
    wi = w.nonzero()[0]
    ext = ext[wi]

    if TimeAnalyse:
        exttime = exttime[wi]

    # Removes double value and move the time to the middle
#    p1 = ext.__getslice__(0,len(ext)-1)
#    p2 = ext.__getslice__(1,len(ext))
    w11 = np.equal(ext[0:len(ext)-1],ext[1:len(ext)])
    w = np.array(np.r_[0, w11], dtype=bool).__invert__()
    wi = w.nonzero()[0]
    ext = ext[wi]

    if TimeAnalyse:
        # w1=(exttime(2:end)-exttime(1:end-1))./2;
        end = len(exttime)

#        p1 = exttime.__getslice__(0,len(exttime)-1)
#        p2 = exttime.__getslice__(1,len(exttime))
        p1 = exttime[0:end-1]
        p2 = exttime[1:end]
        w1 = (p2 - p1) / 2.0
#        p3 = w.__getslice__(1,len(w))
        p3 = w[1:len(w)]
        exttime = np.r_[p1 + w1*(p3.__invert__()), exttime[end-1]]
        exttime = exttime[wi]

    # Again, check the extremes
    # at this point there might not be enough points left
    if len(ext) > 2:
        w1 = np.diff(ext)
        w11 = ((w1[0:-1]*w1[1:len(w1)]).__le__(0))
        w = np.r_[1, w11, 1]
        wi = w.nonzero()[0]
        ext = ext[wi]
        if TimeAnalyse:
            exttime = exttime[wi]

    try:
        type(exttime)
        return ext, exttime

    except:
        # set ext and exttime to None in case they are not asked
        return ext, None


def findtp(x):
    """
    Find Turning Points
    ===================

    ind = findtp(x)

    This is a condensed version withouth the wave filtering height h
    of the wafo object TimeSeries.turning_points

    Paramters
    ---------
    x : array-like
        1D time history of loading

    Returns
    -------
    ind : array-like
        1D array indices for the turning points

    Wafo function
    """
    n = len(x)

    # force it into a 1D array
    xn = np.atleast_1d(x).ravel()

    # element wise difference, convert to (-1, 0 or 1), and force to int8
    xn = np.int8(np.sign(np.diff(xn)))
    # indices to minima and maxima in the original sequence x
    ind = _findcross(xn) + 1

    # the Nieslony approach always put the first loading point as the first
    # turning point. The original wafo method does something slightly different.
    # Here the Nieslony's approach is followed
    if x[ind[0]] != x[0]:
        # make sure the first turning point is the first of the signal
        ind = np.r_[0, ind, n - 1]
    else:
        # only add the last point of the signal as well
        ind = np.r_[ind, n - 1]

    return ind


def _findcross(xn):
    """
    _findcross(xn)

    Return indices to zero up and downcrossings of a vector

    Wafo function
    """
    # try to load the findcross function in the C library first
    if rainflowlib is not None:
        ind, m = rainflowlib.findcross(xn, 0.0)
        return ind[:m]

    # if the C library failes, fall back to this python version
    n = len(xn)
    iz, = (xn == 0).nonzero()
    if len(iz) > 0:
        # Trick to avoid turning points on the crossinglevel.
        if iz[0] == 0:
            if len(iz) == n:
                warnings.warn('All values are equal to crossing level!')
                return np.zeros(0, dtype=np.int)

            diz = np.diff(iz)
            if len(diz) > 0 and (diz > 1).any():
                ix = iz[(diz > 1).argmax()]
            else:
                ix = iz[-1]

            #x(ix) is a up crossing if  x(1:ix) = v and x(ix+1) > v.
            #x(ix) is a downcrossing if x(1:ix) = v and x(ix+1) < v.
            xn[0:ix + 1] = -xn[ix + 1]
            iz = iz[ix + 1::]

        for ix in iz.tolist():
            xn[ix] = xn[ix - 1]

    #% indices to local level crossings ( without turningpoints)
    ind, = (xn[:n - 1] * xn[1:] < 0).nonzero()
    return ind

def test_loadresults():
    """
    Compare HAWC2 result files in binary and ascii formats
    """
    atol = 0.001
    rtol = 0.01
    print '\nCompare bin and ascii hawc2 file formats, tolerances are set to:'
    print 'abs tol = ', atol, 'rel tol =', rtol
    # the ascii file
    file_path = 'data/'
    file_name = 'hawc2test_ascii'
    sig_ascii = LoadResults(file_path, file_name)

    # the binary file
    file_name = 'hawc2test_bin'
    sig_bin = LoadResults(file_path, file_name)

    # select which channels are equal to each other, given an absolute and
    # relative difference: absolute(a - b) <= (atol + rtol * absolute(b))
#    p=0
    wrong_chan = []
    for k in xrange(sig_bin.sig.shape[1]):
        g = np.allclose(sig_ascii.sig[:,k],sig_bin.sig[:,k],atol=atol,rtol=rtol)
        if not g:
#            print str(k).ljust(5),
            wrong_chan.append(k)
        else:
            pass
#            print str(g).ljust(5),
#        p+=1
#        if p > 10:
#            p=0
#            print ''

    if len(wrong_chan) == 0:
        print 'ASCII and BINARY files are within tolerances:'
    else:
        print '\nAbout the channels who are NOT wihtin the tolerances:'
        print 'show: channel index, number of innacurate points, channel range'
        print 'and compare the first 7 inacurate occurences of bin and ascii'
        np.set_printoptions(precision=4, suppress=True)
        rtol = 0.01
        for k in wrong_chan:
            b = sig_ascii.sig[:,k]
            err = np.abs(sig_bin.sig[:,k] - b).__ge__(np.abs(b)*rtol + atol)
            erri = err.nonzero()[0]
            print 'channeli:',k, len(erri), 'range:', np.abs(b.max()-b.min())
            if len(erri) < 7:
                print sig_bin.sig[erri,k]
                print b[erri]
            else:
                print sig_bin.sig[erri[0:7],k]
                print b[erri[0:7]]


def test_rainflow():
    """
    Compare the python results with original Matlab implementation from Nieslony
    """
    np.set_printoptions(precision=4, suppress=True)
    # --------------------------------------------------------------------------
    # test if Python and Matlab results are the same
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # load the reference data from the Matlab .mat file, no time signal used
    matfile = 'data/rainflow_matlab.mat'
    data = scipy.io.loadmat(matfile, matlab_compatible=True)

    sig = data['signal']

#    file_name='s1_steady_free_z_yawfreenoda_c0_y-10_ab00_12.5ms_sh_0.2tow_sh_0'
#    file_path = 'data/'
#    sig = LoadResults(file_path, file_name)
#    sig = sig.sig[:,[0,14,15,16,35,36,37]]

    # the one with the timings
    ext_mat = data['ext'][:,0]
    exttime_mat = data['exttime'][:,0]
    rf_mat = data['rf'].transpose()
    dt = data['dt'][:,0]
    clsn = data['clsn']
    # no timings
    tp2_mat = data['sig_tp2'][:,0]
    rfc2_mat = data['sig_rf2'].transpose()
    tp3_mat = data['sig_tp3'][:,0]
    rfc3_mat = data['sig_rf3'].transpose()
    tp4_mat = data['sig_tp4'][:,0]
    rfc4_mat = data['sig_rf4'].transpose()
    tp5_mat = data['sig_tp5'][:,0]
    rfc5_mat = data['sig_rf5'].transpose()
    tp6_mat = data['sig_tp6'][:,0]
    rfc6_mat = data['sig_rf6'].transpose()
    tp7_mat = data['sig_tp7'][:,0]
    rfc7_mat = data['sig_rf7'].transpose()
    # original s1_steady_free_z_yawfreenoda_c0_y-10_ab00_12.5ms_sh_0.2tow_sh_0
    # channels: 0,14,15,16,35,36,37 : time, blade root M, tower base M

    # --------------------------------------------------------------------------
    print 'all tests should evaluate to True!!'
    # Python, no time version
    print '----> rf3'
    tp2, tp2_time = sig2ext(sig[:,1])
    rfc2 = rfc_astm(tp2)
    rfc2_w = rainflow(sig[:,1], tp_method='nieslony')

    tp3, tp3_time = sig2ext(sig[:,2])
    rfc3 = rfc_astm(tp3)
    rfc3_w = rainflow(sig[:,2], tp_method='nieslony')

    tp4, tp4_time = sig2ext(sig[:,3])
    rfc4 = rfc_astm(tp4)
    rfc4_w = rainflow(sig[:,3], tp_method='nieslony')

    tp5, tp5_time = sig2ext(sig[:,4])
    rfc5 = rfc_astm(tp5)
    rfc5_w = rainflow(sig[:,4], tp_method='nieslony')

    tp6, tp6_time = sig2ext(sig[:,5])
    rfc6 = rfc_astm(tp6)
    rfc6_w = rainflow(sig[:,5], tp_method='nieslony')

    tp7, tp7_time = sig2ext(sig[:,6])
    rfc7 = rfc_astm(tp7)
    rfc7_w = rainflow(sig[:,6], tp_method='nieslony')

#    print 'compare all shapes of the results:'
#    print '-----> tp'
    print tp2.shape == tp2_mat.shape, tp3.shape == tp3_mat.shape,
    print tp4.shape == tp4_mat.shape, tp5.shape == tp5_mat.shape,
    print tp6.shape == tp6_mat.shape, tp7.shape == tp7_mat.shape
#    print '-----> rfc'
    print rfc2.shape == rfc2_mat.shape, rfc3.shape == rfc3_mat.shape,
    print rfc4.shape == rfc4_mat.shape, rfc5.shape == rfc5_mat.shape,
    print rfc6.shape == rfc6_mat.shape, rfc7.shape == rfc7_mat.shape
#    print '-----> rfc_w'
    print rfc2_w.shape == rfc2_mat.shape, rfc3_w.shape == rfc3_mat.shape,
    print rfc4_w.shape == rfc4_mat.shape, rfc5_w.shape == rfc5_mat.shape,
    print rfc6_w.shape == rfc6_mat.shape, rfc7_w.shape == rfc7_mat.shape

#    print 'compare all arrays with np.allclose on standard atol and rtol'
#    print '-----> tp'
    print np.allclose(tp2, tp2_mat), np.allclose(tp3, tp3_mat),
    print np.allclose(tp4, tp4_mat), np.allclose(tp5, tp5_mat),
    print np.allclose(tp6, tp6_mat), np.allclose(tp7, tp7_mat)
#    print '-----> rfc'
    print np.allclose(rfc2, rfc2_mat), np.allclose(rfc3, rfc3_mat),
    print np.allclose(rfc4, rfc4_mat), np.allclose(rfc5, rfc5_mat),
    print np.allclose(rfc6, rfc6_mat), np.allclose(rfc7, rfc7_mat)

    # Python, time version
    print '----> rf5'
    tp, tp_time = sig2ext(sig[:,1], dt=dt, clsn=clsn, debug=False)
    rfc = rfc_astm(tp, exttime=tp_time)
#    print 'shapes:',
    print tp.shape==ext_mat.shape, tp_time.shape==exttime_mat.shape,
    print rfc.shape==rf_mat.shape
#    print 'allclose:',
    print np.allclose(tp, ext_mat),
    print np.allclose(tp_time,exttime_mat), np.allclose(rfc,rf_mat)

    # compare turning points from wafo and sig2ext
    print '----> tp wafo and Nieslony'
    tpw2 = sig[findtp(sig[:,1]),1]
#    print tpw2.shape, tp2.shape, tpw2[0:3], tp2[0:3]

    tpw3 = sig[findtp(sig[:,2]),2]
#    print tpw3.shape, tp3.shape, tpw3[0:3], tp3[0:3]

    tpw4 = sig[findtp(sig[:,3]),3]
#    print tpw4.shape, tp4.shape, tpw4[0:3], tp4[0:3]

    tpw5 = sig[findtp(sig[:,4]),4]
#    print tpw5.shape, tp5.shape, tpw5[0:3], tp5[0:3]

    tpw6 = sig[findtp(sig[:,5]),5]
#    print tpw6.shape, tp6.shape, tpw6[0:3], tp6[0:3]

    tpw7 = sig[findtp(sig[:,6]),6]
    print np.allclose(tpw2, tp2),np.allclose(tpw3, tp3),np.allclose(tpw4, tp4),
    print np.allclose(tpw5, tp5),np.allclose(tpw6, tp6),np.allclose(tpw7, tp7)

    # testing the data check in rainflow
#    rfc00 = rainflow(sig[:,3], tp_method='wafo',clsn=5)

    # who is faster, tp wafo or Nieslony? Do in a IPython shell
#    %timeit sig[findtp(sig[:,6]),6]
#    %timeit sig2ext(sig[:,1])
    # and the winner is: WAFO!
    # wafo   : 10000 loops, best of 3: 179 us per loop
    # sig2ext:  1000 loops, best of 3: 437 us per loop

#    code_ini = "import numpy as np;from fatigue import *"
#    code_loop = "np.array(np.r_[0, gg, 0])"
#    t = timeit.Timer(code_loop,code_ini)
#    print t.timeit(number=50)/50.


def plot_hist():
    """
    """

    pass


if __name__ == '__main__':

    # used this test to compare the python results with Matlab
    test_rainflow()

#    # compare bin and ascii formats
#    test_loadresults()
#
    # Example
    # -------------------------------------------------------------------------
    # load the reference data from the HAWC2 result file
    # -------------------------------------------------------------------------
#    file_name='hawc2test_bin'
    file_name = '4ms_ti0.18_s11_tshad2_a1b1'
    file_path = 'data/'
    sig = LoadResults(file_path, file_name)
    # blade root moments: chi 12 13 14

    # -------------------------------------------------------------------------
    # rainflow counting
    tp_method = 'astm'
    sig_rf = rainflow(sig.sig[:,19], tp_method=tp_method, clsn=None)
    # the rainflow counted histogram
    hist, bin_edges, bin_avg = rfc_hist(sig_rf, nrbins=46)

#    plt.figure()
#    plt.plot(sig.sig[:,0], sig.sig[:,19])
#    plt.grid(True)
#    plt.show()

    neq = 600.0
    textbox = '%s\n' % tp_method
    textbox += 'EQ number=%1.1f' % neq
    for k in [3.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
        textbox += '\n'
        # equivalent load
        S = palmgren_mimer(hist, bin_avg, material=k, neq=neq)
        textbox += 'EQ(m=% 2.0f): %11.4f' % (k, 2.0*S)

    bbox = dict(boxstyle="round", edgecolor=(1., 0.5, 0.5),
                    facecolor=(1., 0.8, 0.8),)

    fig = plt.figure()
    plt.hist(hist, bins=46)
    plt.text(700, 25, textbox, fontsize=12, verticalalignment='bottom',
                 horizontalalignment='right', bbox=bbox)
    plt.grid(True)
    plt.show()

    # -------------------------------------------------------------------------
    # rainflow counting
    tp_method = 'wafo'
    sig_rf = rainflow(sig.sig[:,19], tp_method=tp_method, clsn=None)
    # the rainflow counted histogram
    hist, bin_edges, bin_avg = rfc_hist(sig_rf, nrbins=46)

#    plt.figure()
#    plt.plot(sig.sig[:,0], sig.sig[:,19])
#    plt.grid(True)
#    plt.show()

    neq = 600.0
    textbox = '%s\n' % tp_method
    textbox += 'EQ number=%1.1f' % neq
    for k in [3.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
        textbox += '\n'
        # equivalent load
        S = palmgren_mimer(hist, bin_avg, material=k, neq=neq)
        textbox += 'EQ(m=% 2.0f): %11.4f' % (k, 2.0*S)

    bbox = dict(boxstyle="round", edgecolor=(1., 0.5, 0.5),
                    facecolor=(1., 0.8, 0.8),)

    fig = plt.figure()
    plt.hist(hist, bins=46)
    plt.text(700, 25, textbox, fontsize=12, verticalalignment='bottom',
                 horizontalalignment='right', bbox=bbox)
    plt.grid(True)
    plt.show()

    # -------------------------------------------------------------------------
    # REQUIRES THE PYWAFO TOOLBOX TO BE INSTALLED
#    # wafo only approach
#    sig_ts = wafo.objects.mat2timeseries(sig.sig[:,[0,19]])
#    sig_tp = sig_ts.turning_points(h=0.0, wavetype=None)
#    # h=0 means no filtering of cycles with a certain height
#    # is this equivalent to whole cycles??
#    sig_cp = sig_tp.cycle_pairs(h=0, kind='max2min', method='clib')
#    # the amplitude is the (min - max) /2
#    sig_cp_Amp = sig_cp.amplitudes()
#    hist_cp = np.histogram(sig_cp_Amp, bins=46, range=(0,sig_cp_Amp.max()))
#    Nlife = np.append(hist_cp[0], np.array([0]))
#
#    neq = 600.0
#    textbox = '%s\n' % tp_method
#    textbox += 'EQ number=%1.1f' % neq
#    for Mat in [3.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
#        textbox += '\n'
#        # equivalent load
#        S = np.power((0.5*Nlife*np.power(hist_cp[1],Mat)).sum()/neq,1./Mat)
##        S = palmgren_mimer(hist, bin_avg, material=k, neq=neq)
#        textbox += 'EQ(m=% 2.0f): %11.4f' % (k, 2.0*S)
#
#    fig = plt.figure()
#    plt.hist(Nlife, bins=46, range=(0,sig_cp_Amp.max()))
#    plt.text(sig_cp_Amp.max()/3.0, 25, textbox, fontsize=12, verticalalignment='bottom',
#                 horizontalalignment='right', bbox=bbox)
#    plt.grid(True)
#    plt.show()

    # -------------------------------------------------------------------------
#    # or via the Fatigue class
#    f = Fatigue(material=14)
#    f.material = 1
#    S1 = f(sig.sig[:,12])

    from fatigue_pydap import rainflow_astm_wrapper
    import fatigue_pydap as fatigue


    rainflow_func = rainflow_astm_wrapper
    tp_method = 'pydap'
    sig_rf = rainflow_func(sig.sig[:,19])
    # place a text box in upper left in axes coords
    neq = 600.0
    textbox = '%s\n' % tp_method
    textbox += 'EQ number=%1.1f' % neq
    for m in [3.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
        textbox += '\n'
        # equivalent load
        S = fatigue.eq_load(sig.sig[:,19], m=m, neq=neq,
                                  rainflow_func=rainflow_func)
        textbox += 'EQ(m=% 2.0f): %11.4f' % (m, S[0])

    fig = plt.figure()
    plt.hist(sig_rf, bins=46)
    plt.text(1200, 600, textbox, fontsize=12, verticalalignment='bottom',
                 horizontalalignment='right', bbox=bbox)
    plt.grid(True)
    plt.show()
