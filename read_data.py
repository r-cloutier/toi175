from imports import *

def bin_data(t, f, ef, dt=.2, include_edges=False, tfull=np.zeros(0)):
    # check that the desired binning is coarser than the input sampling
    if np.diff(t).mean() > dt:
        if include_edges:
            assert tfull.size > 1
            tbin = np.append(np.append(tfull.min(), t), tfull.max())
            fbin = np.append(np.append(np.median(f[:10]), f), \
                             np.median(f[-10:]))
            efbin = np.append(np.append(np.median(ef[:10]), ef),
                              np.median(ef[-10:]))
            return tbin, fbin, efbin
        else:
            return t, f, ef
    
    Nbounds = int(np.floor((t.max()-t.min()) / dt))
    tbin, fbin, efbin = np.zeros(Nbounds-1), np.zeros(Nbounds-1), \
                        np.zeros(Nbounds-1)
    for i in range(Nbounds-1):
        inds = np.arange(t.size/Nbounds) + i*t.size/Nbounds
        tbin[i]  = t[inds].mean()
        fbin[i]  = np.median(f[inds])
        efbin[i] = np.mean(ef[inds]) / np.sqrt(inds.size)
    if include_edges:
        assert tfull.size > 1
        tbin = np.append(np.append(tfull.min(), tbin), tfull.max())
        fbin = np.append(np.append(np.median(f[:10]), fbin), np.median(f[-10:]))
        efbin = np.append(np.append(np.median(ef[:10]), efbin),
                          np.median(ef[-10:]))
    return tbin, fbin, efbin


def read_HARPS():
    # RVs
    bjd, rv, erv = np.loadtxt('input_data/TOI-175_050419_NAIRA.rdb',
                              skiprows=2).T
    bjd += 24e5
    print rv.min()*1e3
    rv -= rv.mean()
    rv *= 1e3
    erv *= 1e3

    # ccf params
    _,fwhm,bis,_,_ = np.loadtxt('input_data/TOI-175_050419_CCFparam.rdb',
                                skiprows=2).T

    # activity indices
    p = np.loadtxt('input_data/TOI-175_ActIndex.rdb',
                   skiprows=2)
    _,Halpha,eHalpha,Hbeta,eHbeta,Hgamma,eHgamma,NaD,eNaD,Sindex,eSindex = p.T
    
    return bjd, rv, erv, fwhm, bis, Halpha, eHalpha, Hbeta, eHbeta, Hgamma, eHgamma, NaD, eNaD, Sindex, eSindex


    
