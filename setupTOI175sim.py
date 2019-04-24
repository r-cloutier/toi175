import rebound
import numpy as np
import rvs
from PyAstronomy.pyasl import foldAt


def setupTOI175simv2(bjd0, Ms, eccs, incs_deg, outname, interval_yrs, 
		     random=True, DeltaOmegamax=180):
    '''Create a simulation of TOI175 with custom input for the 
    planetary parameters in the form of 1d arrays.'''
    # Initialize
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.units = ('AU','Msun','yr')
    sim.automateSimulationArchive("%s.bin"%outname, interval=interval_yrs)

    # Get keplerian parameters
    assert eccs.size == 3
    assert incs_deg.size == 3
    Ps = np.array([2.2531, 3.6904, 7.4512])
    if random:
        Ps += np.array([np.random.randn()*4e-4, np.random.randn()*3e-4, np.random.randn()*8e-4])
	while np.any(Ps <= 0):
	    Ps = np.array([2.2531, 3.6904, 7.4512]) + np.array([np.random.randn()*4e-4, np.random.randn()*3e-4, np.random.randn()*8e-4])

    T0s = np.array([1366.1708, 1367.2752, 1362.7376]) + 2457000
    if random:
	T0s += np.array([np.random.randn()*1e-4, np.random.randn()*6e-4, np.random.randn()*9e-4])

    # sample mass of the smallest planet directly from the PDF since it's not gaussian (i.e. a non-detection)
    mp_175d03 = np.random.choice(np.load('toi175d03_2d25_planetmass_PDF.npy')) + np.random.randn()*5e-4
    mps = np.array([mp_175d03, 2.6, 2.4])
    if random:
    	mps += np.array([0, np.random.randn()*0.4, np.random.randn()*0.6])
    	while np.any(mps <= 0):
	    mps = np.array([0.34, 2.6, 2.4]) + np.array([np.random.randn()*.42, np.random.randn()*0.4, np.random.randn()*0.6])

    smas = rvs.m2AU(rvs.semimajoraxis(Ps, Ms, mps))
    mps = rvs.kg2Msun(rvs.Mearth2kg(mps))
    nplanets = Ps.size
    omegas = np.random.uniform(-np.pi, np.pi, nplanets)
    Omegas = np.random.uniform(-np.pi, np.pi, nplanets)
    thetas = 2*np.pi * foldAt(bjd0, Ps, T0s) - np.pi

    # Add star
    sim.add(m=Ms, hash='star')

    # Add planets
    for i in range(nplanets):
        sim.add(m=mps[i], a=smas[i], inc=np.deg2rad(incs_deg[i]-90), e=eccs[i],
                omega=omegas[i], Omega=Omegas[i], theta=thetas[i])

    sim.move_to_com()

    RHill1 = (mps[:2].sum()/(3.*Ms))**(1./3) * smas[:2].mean()
    RHill2 = (mps[1:].sum()/(3.*Ms))**(1./3) * smas[1:].mean()
    sim.exit_min_distance = float(np.max([RHill1,RHill2]))

    return sim
