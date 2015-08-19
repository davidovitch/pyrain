'''
Created on 04/03/2013

@author: mmpe
'''


from cython_compile import cython_import
from import_cython import import_cython
import numpy as np


def rfc_hist(sig_rf, nrbins=46):
    """
    Histogram of rainflow counted cycles
    ====================================

    hist, bin_edges, bin_avg = rfc_hist(sig, nrbins=46)

    Divide the rainflow counted cycles of a signal into equally spaced bins.

    Created on Wed Feb 16 16:53:18 2011
    @author: David Verelst
    Modified 10.10.2011 by Mads M Pedersen to elimintate __copy__ and __eq__

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

    rf_half = sig_rf

    # the Matlab approach is to divide into 46 bins
    bin_edges = np.linspace(0, 1, num=nrbins + 1) * rf_half.max()
    hist = np.histogram(rf_half, bins=bin_edges)[0]
    # calculate the average per bin
    hist_sum = np.histogram(rf_half, weights=rf_half, bins=bin_edges)[0]
    # replace zeros with one, to avoid 0/0
    hist_ = hist.copy()
    hist_[(hist == 0).nonzero()] = 1.0
    # since the sum is also 0, the avg remains zero for those whos hist is zero
    bin_avg = hist_sum / hist_

    return hist, bin_edges, bin_avg


def check_signal(signal):
    # check input data validity
    if not type(signal).__name__ == 'ndarray':
        raise TypeError('signal should be ndarray, not: ' + type(signal).__name__)

    elif len(signal.shape) not in (1, 2):
        raise TypeError('signal should be 1D or 2D, not: ' + str(len(signal.shape)))

    if len(signal.shape) == 2:
        if signal.shape[1] > 1:
            raise TypeError('signal should have one column only, not: ' + str(signal.shape[1]))


def rainflow_windap_wrapper(signal):
    levels = 255.
    thresshold = (255. / 50)
    check_signal(signal)
    #type <double> is reuqired by <find_extreme> and <rainflow>
    signal = signal - np.min(signal)
    if np.max(signal) > 0:
        gain = float(np.max(signal)) / levels
        signal = signal / gain
        signal = np.round(signal).astype(np.int)

        #Import peak_throug(find extremes) and pair_range(rainflowcount).
        #If possible the module is compiled using cython otherwise the python implementation is used

        cython_import('common.fatigue.peak_trough')
        from peak_trough import peak_trough


        #Convert to list of local minima/maxima where difference > thresshold
        sig_ext = peak_trough(signal, thresshold)

        import_cython("common.fatigue.pair_range")
        from pair_range import pair_range
        #rainflow count
        sig_rfc = pair_range(sig_ext)

        #This is not smart but corresponds to the implementation in Windap
        sig_rfc = np.array(sig_rfc)
        #sig_rfc = sig_rfc[np.where(sig_rfc>=thresshold)[0]]
        sig_rfc = np.round(sig_rfc / thresshold) * gain * thresshold
        #sig_rfc = np.round(np.array(sig_rfc)/thresshold)*gain*thresshold
        #return sig_rfc[np.where(sig_rfc>=gain*thresshold)[0]]
        return sig_rfc



def rainflow_astm_wrapper(signal):
    check_signal(signal)

    # type <double> is reuqired by <find_extreme> and <rainflow>
    signal = signal.astype(np.double)

    # Import find extremes and rainflow.
    # If possible the module is compiled using cython otherwise the python implementation is used
    cython_import('rainflowcount_astm')
    from rainflowcount_astm import find_extremes, rainflow_astm

    # Remove points which is not local minimum/maximum
    sig_ext = find_extremes(signal)

    # rainflow count
    sig_rfc = rainflow_astm(sig_ext)

    return np.array(sig_rfc)


def eq_load(signal, no_bins=46, m=[3, 4, 6, 8, 10, 12], neq=1, rainflow_func=rainflow_windap_wrapper):
    sig_rf = rainflow_func(signal)
    if sig_rf is None:
        return []
    hist_data, x, bin_avg = rfc_hist(sig_rf, no_bins)

    m = np.atleast_1d(m)

    eq = []
    for i in range(len(m)):
        eq.append(np.power(np.sum(0.5 * hist_data * np.power(bin_avg, m[i])) / neq, 1. / m[i]))
    return eq

# print eq_load(np.array([-2.0,0.0,1.0,0.0,-3.0,0.0,5.0,0.0,-1.0,0.0,3.0,0.0,-4.0,0.0,4.0,0.0,-2.0]),no_bins=50,neq=17)
