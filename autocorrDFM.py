# compute the autocorrelation time tau of parallel MCMC chains using the
# formalism present by Dan Foreman-Mackey here:
# https://dfm.io/posts/autocorr/
from imports import *

#def autocorr_new(y, c=5.0):
def autocorr_time_DFM(sampler_chain, c=5.):
    '''Use the DFM method of computing the autocorrelation time of an input 
    time-series. y has dimensions [nchains (nwalkers), nsteps, nparams]'''
    y = sampler_chain
    assert len(y.shape) == 3
    Nchains, Nsamples, Nparams = y.shape
    f = np.zeros_like(y)
    for i in range(Nparams):
        for j in range(Nchains):
            f[j,:,i] = autocorr_func_1d(y[j,:,i]) / Nchains / Nparams
    f = np.nansum(np.nansum(f,0),1)
    assert f.size == Nsamples
    taus = 2.0*np.cumsum(f)-1.0
    window = auto_window(taus, c)
    return taus[window], Nsamples


def auto_window(taus, c):
    '''from Sokal (1989)'''
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1


def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i


def autocorr_func_1d(x, norm=True):
    '''Compute the autocorrelation function of an input time-series.'''
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2*n)
    acf = np.fft.ifft(f * np.conjugate(f))[:len(x)].real
    acf /= 4*n
    
    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf
