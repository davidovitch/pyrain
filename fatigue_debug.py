# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 14:00:40 2011

@author: dave
"""

import numpy as np
import scipy

def sig2ext_debug(sig, dt=None, clsn=None, debug=False):
    """
    same as sig2ext, but includes the debug process
    """
    # NOTE: array.__getslice__(a,b) is mentioned to be more efficient compared to
    # array[a:b]. However, using timeit it seems timings are comparable and very 
    # fast anyway
#    code_ini = "import numpy as np;import scipy;gg=np.random.rand(30000)"
#    code_loop = "np.array(np.r_[0, gg, 0])"
#    t = timeit.Timer(code_loop,code_ini)
#    print t.timeit(number=50)/50.
#    
#    code_loop = "w = scipy.zeros((len(gg)+2));" \
#                + "w.__setitem__(0,False);" \
#                + "w.__setitem__(len(w)-1,False);" \
#                + "w.__setslice__(1,len(w)-1,gg)"
#    t = timeit.Timer(code_loop,code_ini)
#    print t.timeit(number=50)/50.
    # the results
    # compact notion: 8.43191146851e-05
    # using __slice : 4.41217422485e-05
    # for IPython console
    # %timeit np.random.rand(30000,1)
    
    # Should there be a time analysis as well? Do not assume default value
    # for time step, but force to specify dt
    if dt is None:
        TimeAnalyse = False
    elif type(dt).__name__ in ('int', 'float', 'ndarray', 'float64', 'float32'):
        TimeAnalyse = True
    else:
        msg = 'dt should be either int, float or ndarray and not' \
            + type(dt).__name__
        raise TypeError, msg
    
    # Force to 1D matrix
    sig = sig.ravel()
    
    # If clsn is given, divide into classes
    if clsn is not None:
        clsn = round(clsn)
        
        smax = np.max(sig)
        smin = np.min(sig)
        
        if debug:
            matfile = 'data/sig_scaled.mat'
            data = scipy.io.loadmat(matfile, matlab_compatible=True)
#            sig1 = data['sig1'][:,0]
#            sig2 = data['sig2'][:,0]
#            sig3 = data['sig3'][:,0]
            sig4 = data['sig4'][:,0]
            smaxm = data['smax']
            sminm = data['smin']
            print smax, smaxm
            print smin, sminm
        
        sig = clsn*((sig-smin)/(smax-smin))
        sig = np.fix(sig)
        # find all values in sig equal to clsn
        # return array with 0 where sig is not clsn and clsn on the rest
        temp = sig.__eq__(clsn)*clsn
        # take all the non zero elements only and do -1
        sig[temp.nonzero()] = clsn - 1.0
        # some kind of scaling...
        sig = (smax-smin)/(clsn-1)*sig+smin
        
        if debug:
#            print np.allclose(sig, sig1),
#            print np.allclose(sig, sig2),
#            print np.allclose(sig, sig3)
            print 'sigs scales ==', np.allclose(sig, sig4),
            matfile = 'data/sig_scaled_final.mat'
            data = scipy.io.loadmat(matfile, matlab_compatible=True)
            sig_final = data['sig_scaled'][:,0]
            print np.allclose(sig, sig_final),
            print np.allclose(sig4, sig_final)
    
    # Create a binary vector where 1 means equality of the extreme or,
    # It is recognized that the first and last point is an ext
    w1 = np.diff(sig)
    # Multiply elements n and n+1, evaluate if they are <= 0
    # w = logical([1;(w1(1:end-1).*w1(2:end))<=0;1]);
    w11 = (w1[0:len(w1)-1]*w1[1:len(w1)]).__le__(0)
    # and get the indices of the elements who comply the <=0 condition
    # nonzero() returns a tuple for each dimension
    w = np.r_[1, w11, 1]
    # matlab returns w_i = 1, w_i=0 is not returned
    wi = w.nonzero()[0]
    ext = sig[wi]
    
    if TimeAnalyse:
        if type(dt).__name__ in ('int', 'float', 'float64', 'float32'):
            # find: returns indices if those who are w == 1
            # exttime=(find(w==1)-1).*dt;
            
            # evaluates to true if w equals to one
            # return indices of the true elements, if false, return zero
            exttime = w.__eq__(1).nonzero()[0]*dt
            # no need to do -1 as in Matlab, because indices are 0...n
        else:
            exttime = dt[wi]
    
    if debug:
        try: print 0, ext.shape, w.shape, exttime.shape
        except: print 0, w.shape
        matfile = 'data/sig2ext_checkpoint0.mat'
        data = scipy.io.loadmat(matfile, matlab_compatible=True)
        wm = data['w'][:,0]
        w1m = data['w1'][:,0]
        extm = data['ext'][:,0]
        print 'checkpoint0', np.allclose(w1,w1m), np.allclose(w,wm),
        print np.allclose(ext,extm)
        print w1.dtype, w1m.dtype, w.dtype, wm.dtype, ext.dtype, extm.dtype
    
    # w =~ logical([0; w1(1:end-1)==0 & w1(2:end)==0; 0]);
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
    
    if debug:
        try: print 1, ext.shape, w.shape, exttime.shape
        except: print 1, w.shape
        matfile = 'data/sig2ext_checkpoint1.mat'
        data = scipy.io.loadmat(matfile, matlab_compatible=True)
        wm = data['w'][:,0]
        w1m = data['w1'][:,0]
        extm = data['ext'][:,0]
        print 'checkpoint1', np.allclose(w1,w1m), np.allclose(w,wm),
        print np.allclose(ext,extm)
        print w1.dtype, w1m.dtype, w.dtype, wm.dtype, ext.dtype, extm.dtype
        matfile = 'data/sig2ext_checkpoint1_p1p2.mat'
        data = scipy.io.loadmat(matfile, matlab_compatible=True)
        p1m = data['p1'][:,0]
        p2m = data['p2'][:,0]
        p3m = data['p3'][:,0]
        print 'ps',np.allclose(p1,p1m), np.allclose(p2,p2m),np.allclose(w11,p3m)
        print p1.dtype, p1m.dtype, p2.dtype, p2m.dtype, w11.dtype, p3m.dtype
    
    # w=~logical([0; ext(1:end-1)==ext(2:end)]);
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
        # exttime=[exttime(1:end-1)+w1.*~w(2:end); exttime(end)];
        exttime = np.r_[p1 + w1*(p3.__invert__()), exttime[end-1]]
        exttime = exttime[wi]
    
    if debug:
        try:
            print 2, ext.shape, w.shape, exttime.shape
        except:
            print 2, w.shape
    
    # Again, check the extremes
    # at this point there might not be enough points left
    if len(ext) > 2:
        w1 = np.diff(ext)
        # w = logical([1; w1(1:end-1).*w1(2:end) < 0 ; 1]);
        w11 = ((w1[0:-1]*w1[1:len(w1)]).__le__(0))
        w = np.r_[1, w11, 1]
        wi = w.nonzero()[0]
        ext = ext[wi]
        if TimeAnalyse:
            exttime = exttime[wi]
    
    if debug:
        try:
            print 3, ext.shape, w.shape, exttime.shape
        except:
            print 3, w.shape
    
    try:
        type(exttime)
        return ext, exttime
    
    except:
        # set ext and exttime to None in case they are not asked
        return ext, None