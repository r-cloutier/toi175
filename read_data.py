from imports import *

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


    
