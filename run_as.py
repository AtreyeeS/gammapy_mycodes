import numpy as np
import os
import astropy
import regions
import sherpa
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.table import vstack as vstack_table
from regions import CircleSkyRegion
from gammapy.data import DataStore, ObservationList
from gammapy.data import ObservationStats, ObservationSummary
from gammapy.background.reflected import ReflectedRegionsBackgroundEstimator
from gammapy.utils.energy import EnergyBounds
from gammapy.spectrum import SpectrumExtraction, SpectrumObservation, SpectrumFit, SpectrumResult
#from gammapy.spectrum import SimulationRealBkg
from gammapy.spectrum.models import PowerLaw, ExponentialCutoffPowerLaw, LogParabola
from gammapy.spectrum import FluxPoints, SpectrumEnergyGroupMaker, FluxPointEstimator
from gammapy.image import SkyImage
from gammapy.spectrum import models

from simulation_realbkg import *

#from gammapy.irf import EnergyDispersion, EnergyDispersion2D, EffectiveAreaTable, EffectiveAreaTable2D

# Setup the logger
import logging
logging.basicConfig()
log = logging.getLogger('gammapy.spectrum')
log.setLevel(logging.WARNING)


#choose the obsevartion
name="Crab"
datastore = DataStore.from_dir("$HESS_DATA") 
src=SkyCoord.from_name(name)
sep=SkyCoord.separation(src,datastore.obs_table.pointing_radec)
srcruns=(datastore.obs_table[sep<2.0*u.deg]) 
obsid=srcruns['OBS_ID'].data
mylist=datastore.obs_list((obsid[0],))

# Define obs parameters
lo_threshold = 0.1 * u.TeV
hi_threshold = 60 * u.TeV

# Define spectral model
index = 3.0 * u.Unit('')
amplitude = 2.5 * 1e-11 * u.Unit('cm-2 s-1 TeV-1')
reference = 1 * u.TeV
model = PowerLaw(index=index, amplitude=amplitude, reference=reference)

# for one obs only

n_obs=1
seeds = np.arange(n_obs)
sim = SimulationRealBkg(source_model=model, obsrun=mylist[0], obspos=src)
sim.run(seeds)

obs=sim.result[0]
fit = SpectrumFit(obs, model=model, stat='wstat')
fit.model.parameters['index'].value = 2
fit.run()

fit.est_errors()

print fit.result[0]

#now run in a loop

n_obs=20
mylist=datastore.obs_list(obsid[0:n_obs])
seeds = np.arange(n_obs)
sims=[]
for i in range(n_obs):
    sim = SimulationRealBkg(source_model=model, obsrun=mylist[i], obspos=src)
    sim.run(seeds)
    sims.append(sim.result[0])

n_on = [obs.total_stats.n_on for obs in sims]
n_off = [obs.total_stats.n_off for obs in sims]
excess = [obs.total_stats.excess for obs in sims]
"""
fix, axes = plt.subplots(1,3, figsize=(12, 4))
axes[0].hist(n_on)
axes[0].set_xlabel('n_on')
axes[1].hist(n_off)
axes[1].set_xlabel('n_off')
axes[2].hist(excess)
axes[2].set_xlabel('excess')
"""
best_fit_index = []
best_fit_flux = []

for obs in sims:
    fit = SpectrumFit(obs, model=model, stat='wstat')
    fit.model.parameters['index'].value = 2.0
    fit.run()
    best_fit_index.append(fit.result[0].model.parameters['index'].value)
    best_fit_flux.append(fit.result[0].model.parameters['amplitude'].value)

fix, axes = plt.subplots(1,3)
axes[0].hist(best_fit_index)
axes[0].set_xlabel('index')
axes[1].hist(best_fit_flux)
axes[1].set_xlabel('amplitude')
axes[2].plot(best_fit_index,best_fit_flux,"ro")
#print('best_fit_index:', best_fit_index)
plt.show()
