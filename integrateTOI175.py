from DynamicalSim import *
import numpy as np
import pylab as plt
import rvs, sys, glob, os
from setupTOI175sim import setupTOI175simv2
import cPickle as pickle


global incb, Rs, Ms
Rs, Ms = .3137, .3124

class TOI175integration:

    def __init__(self, folder, index, index2, Nyrs=1e6, Nout=500):

	self.index, self.index2, self.Nyrs, self.Nout = int(index), int(index2), int(Nyrs), int(Nout)
	self.Ms = 0.
	while self.Ms <= 0:
	    self.Ms = Ms + np.random.randn() * .0032
	self.bjd0 = 2458354.101878819  # first TESS photometry epoch
	self.DONE = False	

	# Setup simulation
	self.folder = folder
        self.outname = '%s/SimArchive/archived%.4d_%.4d'%(self.folder, self.index, self.index2)
	sim = setupsim(self, self.outname)
	self.mps,self.sma0,self.ecc0,self.inc0,self.omega0,self.Omega0,self.theta0 = get_initial_parameters(sim)
	self.inc0 += 90
	#self.Rhill = get_Rhill_init(self.mps, self.Ms, self.sma0)
	self.pickleobject()

	# Integrate simulation
	print 'Integrating system...'
	self.bjds,self.RVs,self.smas,self.eccs,self.incs,self.dist,self.stable = \
		integrate_sim(np.linspace(0, rvs.yrs2days(self.Nyrs), self.Nout) + self.bjd0, sim)
	self.bs = np.zeros(self.incs.shape)
	for i in range(self.bjds.size):
	    self.bs[i] = rvs.impactparam_inc(rvs.AU2m(self.smas[i])/rvs.Rsun2m(Rs), self.incs[i]+90, ecc=self.eccs[i])
	self.DONE = True
	self.pickleobject()


    def pickleobject(self):
	f = open('%s/pickles/TOI175simulation%.4d_%.4d'%(self.folder, self.index, self.index2), 'wb')
	pickle.dump(self, f)
	f.close()


def setupsim_corr(self, outname, N=50):
    # Sample eccentircities logarithmically
    eccbs = np.logspace(-3, -1, N)
    ecccs = np.logspace(-3, -1, N)
    #eccbs = np.linspace(1e-3, .59, N)
    #ecccs = np.linspace(1e-3, .43, N)
    eccs = np.array([ecccs[self.index-1], eccbs[self.index-1]])
    eccs += np.random.randn(eccs.size) * eccs/10.
    while np.any(eccs<0) or np.any(eccs>=1):
	eccs = np.array([ecccs[self.index-1], eccbs[self.index-1]])
    	eccs += np.random.randn(eccs.size) * eccs/10.
    incs_deg = np.array([sampleK218cincdeg(0.06, eccs[0]), incb+np.random.randn()*.0084])
    # How often to save snapshots
    interval_yrs = self.Nyrs / self.Nout
    return setupK2_18simv2(self.bjd0, self.Ms, eccs, incs_deg, outname, int(interval_yrs))


def setupsim(self, outname):
    # Sample eccentricities logarithmically
    eccupperlim = 1.  # .1
    #eccs = np.array([10**(np.random.uniform(-3,np.log10(eccupperlim)))])
    #eccs = np.append(eccs, 10**(np.random.uniform(-3,np.log10(eccupperlim))))
    #eccs = np.random.uniform(1e-3,eccupperlim,2)
    eccs = np.random.uniform(0, eccupperlim, 3)
    #while np.any(eccs<0) or np.any(eccs>=1):
    #	eccs = np.array([10**(np.random.uniform(-3,-1))])
    #	eccs = np.append(eccs, 10**(np.random.uniform(-3,-1)))
    incs_deg = np.array([88.8, 89.4, 88.6])  # Kostov 2019 plus my stellar params
    incs_deg += np.array([np.random.randn()*.73, np.random.randn()*.34, np.random.randn()*.22]) # Kostov 2019
    # How often to save snapshots
    interval_yrs = self.Nyrs / self.Nout
    return setupTOI175simv2(self.bjd0, self.Ms, eccs, incs_deg, outname, int(interval_yrs))


def loadpickle(fname):
    f = open(fname, 'rb')
    self = pickle.load(f)
    f.close()
    return self

def sampleK218cincdeg(smac, eccc):
    '''Sample the orbital inclination from its eccentrcitity (<i^2>=<e^2>/4)
    ensuring that the planet does not transit. 
    smac in AU.'''
    b, ind = 0, 0
    while abs(b) < 1:
	incc = incb + np.rad2deg(np.sqrt(eccc**2/4.)) * np.random.choice([-1,1])
	# Add noise to inc if cannot find solution
	if ind > 100:
	    incc += np.random.randn() * ind/100.
    	b = rvs.impactparam_inc(smac/rvs.m2AU(rvs.Rsun2m(Rs)), incc)
	ind += 1
    return incc


def sampleK218incdeg_uncorr(incbb, sig=1.5):
    '''Sample incc from a Gaussian but reject inclinations that result in a transit.'''
    #sig = 2.2 * 3  # 3sigma on mode of mutual inclination distribution from Fabrycky+2014
    smac = .06
    incc = incbb + np.random.randn()*sig
    b = rvs.impactparam_inc(smac/rvs.m2AU(rvs.Rsun2m(Rs)), incc)
    while abs(b) < 1:
	incc = incbb + np.random.randn()*sig
    	b = rvs.impactparam_inc(smac/rvs.m2AU(rvs.Rsun2m(Rs)), incc)
    return incc


def RHill_from_sim(sim):
    '''Compute the Hill radius from a simulation with 1 star plus 2 planets.'''
    ps = sim.particles
    star, p1, p2 = ps['star'], ps[1], ps[2]
    return ((p1.m+p2.m)/(3.*star.m))**(1./3) * (p1.a+p2.a)/2.   # units of a (e.g. AU)


def get_Rhill_init(mps, Ms, smas):
    '''mps in Msun, Ms in Msun, smas in AU.'''
    return (mps.sum()/(3.*Ms))**(1./3) * smas.sum()/2.


def get_initial_parameters(sim):
    '''Return the initial planetary parameters.'''
    ps = sim.particles
    mps, smas, eccs, incs = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)
    omegas, Omegas, thetas = np.zeros(3), np.zeros(3), np.zeros(3)
    for i in range(3):
	mps[i] = ps[i+1].m
        smas[i] = ps[i+1].a
        eccs[i] = ps[i+1].e
        incs[i] = np.rad2deg(ps[i+1].inc)
        omegas[i] = np.rad2deg(ps[i+1].omega)
        Omegas[i] = np.rad2deg(ps[i+1].Omega)
        thetas[i] = np.rad2deg(ps[i+1].theta)
    return mps, smas, eccs, incs, omegas, Omegas, thetas


def do_i_run_this_simulation(fname):
    if not os.path.isfile(fname):
	return True
    else:
	d = loadpickle(fname)
	return False if d.DONE else True


if __name__ == '__main__':
    # Neccentricities == 10 in log space
    # so index in [1,1e2] with 10 simulations per ecc value
    Nyrs, Nsim_per_ecc, Nout = 1e6, 100, 1e3
    folder = 'pickles_uncorr_ALLecc_lin_toi175'
    index = int(sys.argv[1])
    for i in range(1,Nsim_per_ecc+1):
	if do_i_run_this_simulation('%s/SimArchive/archived%.4d_%.4d'%(folder, index, i)):
    	    self = TOI175integration(folder, index, i, Nyrs=Nyrs, Nout=Nout)
