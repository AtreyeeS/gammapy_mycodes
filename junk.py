# coding: utf-8
import gammapy
import gammapy
import numpy as np
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
from gammapy.spectrum.models import PowerLaw, ExponentialCutoffPowerLaw, LogParabola
from gammapy.spectrum import FluxPoints, SpectrumEnergyGroupMaker, FluxPointEstimator
from gammapy.image import SkyImage
from gammapy.spectrum import models
# Setup the logger
import logging
logging.basicConfig()
log = logging.getLogger('gammapy.spectrum')
log.setLevel(logging.WARNING)

datastore = DataStore.from_dir("$HESS_DATA") 
crab=SkyCoord.from_name("Crab")
sep=SkyCoord.separation(crab,datastore.obs_table.pointing_radec)
sep
sep1=sep<2.0*u.deg
sep1
sep=SkyCoord.separation(crab,datastore.obs_table.pointing_radec)
crabrun=(datastore.obs_table[sep<2.0*u.deg]) 
obsid=crabrun['OBS_ID'].data
obsid
obsid[0]
fakerun=datastore.obs_list(obsid[0])
fakerun=datastore.obs_list(obsid[0:2])
fakerun=ObservationList(obsid[0])
fakerun=datastore.obs(obsid[0])
# Define obs parameters
livetime = 1.0 * u.hr
offset = 0.01 * u.deg
lo_threshold = 0.1 * u.TeV
hi_threshold = 60 * u.TeV
# Define spectral model
index = 2.0 * u.Unit('')
amplitude = 2.5 * 1e-12 * u.Unit('cm-2 s-1 TeV-1')
reference = 1 * u.TeV
model = PowerLaw(index=index, amplitude=amplitude, reference=reference)
get_ipython().set_next_input(u'aeff=fakerun.aeff');get_ipython().magic(u'pinfo fakerun.aeff')
aeff=fakerun.aeff.to_effective_area_table(offset=offset)
edisp=fakerun.edisp.to_energy_dispersion(offset=offset)

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from gammapy.irf import EnergyDispersion, EnergyDispersion2D, EffectiveAreaTable, EffectiveAreaTable2D
from gammapy.irf import EnergyDependentMultiGaussPSF 
from gammapy.spectrum import SpectrumSimulation, SpectrumFit
from gammapy.spectrum.models import PowerLaw

import os
edisp=fakerun.edisp.to_energy_dispersion(offset=offset)
aeff=fakerun.aeff.to_effective_area_table(offset=offset)
aeff.lo_threshold = lo_threshold
aeff.hi_threshold = hi_threshold
get_ipython().magic(u'pinfo fakerun.bkg')
inubilutbkh=fakerun.bkg()
inubilutbkh=fakerun.bkg.evaluate_at_offset(offset=offset)
inubuiltbkg=fakerun.bkg.evaluate_at_offset(offset=offset)
inbuiltbkg=fakerun.bkg.evaluate_at_offset(offset=offset)
fig, axes = plt.subplots(1, 2, figsize=(12, 6)

)
edisp.plot_matrix(ax=axes[0])
aeff.plot(ax=axes[1])
plt.show()
aeff.lo_threshold = lo_threshold
aeff.hi_threshold = hi_threshold
sim = SpectrumSimulation(aeff=aeff, edisp=edisp, source_model=model, livetime=livetime)
sim.simulate_obs(seed=42, obs_id=0)
sim.obs.peek()
print sim.obs
plt.show()
amplitude = 2.5 * 1e-10 * u.Unit('cm-2 s-1 TeV-1')
sim = SpectrumSimulation(aeff=aeff, edisp=edisp, source_model=model, livetime=livetime)
sim.simulate_obs(seed=42, obs_id=0)
sim.obs.peek()
print sim.obs
plt.show()
amplitude
model = PowerLaw(index=index, amplitude=amplitude, reference=reference)
sim = SpectrumSimulation(aeff=aeff, edisp=edisp, source_model=model, livetime=livetime)
sim.simulate_obs(seed=42, obs_id=0)
sim.obs.peek()
print sim.obs
plt.show()
get_ipython().magic(u'save 1-40 sim-bkg,py')
get_ipython().magic(u'save 1-40 sim-bkg.py')
get_ipython().magic(u'save 1-40 sim-bkg.py')
get_ipython().magic(u'save sim-bkg.py 1-40')
crab
allsky_mask = SkyImage.read(EXCLUSION_FILE)
exclusion_mask = allsky_mask.cutout(
    position=crab,
    size=Angle('6 deg'),
)

EXCLUSION_FILE = '$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits'
allsky_mask = SkyImage.read(EXCLUSION_FILE)
exclusion_mask = allsky_mask.cutout(
    position=crab,
    size=Angle('6 deg'),
)
from gammapy.background.reflected import ReflectedRegionsBackgroundEstimator
get_ipython().magic(u'pinfo ReflectedRegionsBackgroundEstimator')
on_region=CircleSkyRegion(crab,0.15 * u.deg)
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=fakerun, on_region=on_region, exclusion_mask = exclusion_mask)
background_estimator.run()
get_ipython().magic(u'pinfo background_estimator.run')
get_ipython().magic(u'pinfo background_estimator.result')
fakerun.obs_id
background_estimator.process()
background_estimator.process(fakerun)
background_estimator.result
background_estimator.result[0]
background_estimator.result
print background_estimator.result
plt.figure(figsize=(8,8))
background_estimator.plot()
len(fakerun)
background_estimator.obs_list
print background_estimator.obs_list
print len(background_estimator.obs_list)
print len(background_estimator.obs_list)
fakerun1=datastore.obs(obsid[0:2])
fakerun1=datastore.obs_list(obsid[0:2])
edisp1=fakerun1.edisp.to_energy_dispersion(offset=offset)
fakerun1.make_mean_edisp()
get_ipython().magic(u'pinfo fakerun1.make_mean_edisp')
fakerun1.make_mean_edisp(offset=offset)
get_ipython().magic(u'pinfo fakerun1.make_mean_edisp')
fakerun1.make_mean_edisp(crab)
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=fakerun1, on_region=on_region, exclusion_mask = exclusion_mask)
background_estimator.process(fakerun)
background_estimator.run(fakerun1)
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=fakerun1, on_region=on_region, exclusion_mask = exclusion_mask)
background_estimator.run()
print(background_estimator.result[0])
plt.figure(figsize=(8,8))
background_estimator.plot()
plt.show()
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=fakerun1, on_region=on_region, exclusion_mask = exclusion_mask)
background_estimator.process()
background_estimator.process(fakerun)
background_estimator.result
background_estimator.run()
background_estimator.result[0]
background_estimator.result
print background_estimator.result
print background_estimator.result[0]
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=fakerun, on_region=on_region, exclusion_mask = exclusion_mask)

background_estimator.process(fakerun)
print background_estimator.result[0]
print background_estimator.result
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=fakerun, on_region=on_region, exclusion_mask = exclusion_mask)

background_estimator.process(fakerun)
background_estimator.plot()
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=fakerun, on_region=on_region, exclusion_mask = exclusion_mask)
background_estimator.run()
background_estimator.obs_list
print background_estimator.obs_list
for obs in  background_estimator.obs_list:
    print 1
    
print background_estimator.obs_list
print len(background_estimator.obs_list)
print len(background_estimator.obs_list())
for obs in  background_estimator.obs_list:
    print 1
    
    
background_estimator.run()
background_estimator1=ReflectedRegionsBackgroundEstimator(obs_list=fakerun1,on_region=on_region, exclusion_mask = exclusion_mask)
print background_estimator1.obs_list
for obs in  background_estimator1.obs_list:
    print 1
    
    
len(background_estimator1.obs_list)
len(background_estimator.obs_list)
a=[1]
for i in a:
    print i
    
a1=np.array(a)
a1
for i in a1:
    print i
    
inbuiltbkg
plt.plot(inbuiltbkg)
plt.show()
plt.plot(inbuiltbkg)
plt.show()
bkg_model = PowerLaw(index=bkg_index, amplitude=bkg_amplitude, reference=reference)
bkg_index = 2.5 * u.Unit('')
bkg_amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1')
reference = 1 * u.TeV

bkg_model = PowerLaw(index=bkg_index, amplitude=bkg_amplitude, reference=reference)
bkg_model
print bkg_model
sim = SpectrumSimulation(aeff=aeff,
                         edisp=edisp,
                         source_model=model,
                         livetime=livetime,
                         background_model=inbuiltbkg, alpha=alpha)
alpha=0.2
sim = SpectrumSimulation(aeff=aeff,
                         edisp=edisp,
                         source_model=model,
                         livetime=livetime,
                         background_model=inbuiltbkg, alpha=alpha)
index
sim.simulate_obs(seed=42, obs_id=0)
sim = SpectrumSimulation(aeff=aeff,
                         edisp=edisp,
                         source_model=model,
                         livetime=livetime,
                         background_model=none, alpha=alpha)
sim = SpectrumSimulation(aeff=aeff,
                         edisp=edisp,
                         source_model=model,
                         livetime=livetime,
                         background_model="none", alpha=alpha)
sim.simulate_obs(seed=42, obs_id=0)
sim = SpectrumSimulation(aeff=aeff, edisp=edisp, source_model=model, livetime=livetime)
sim.simulate_obs(seed=42, obs_id=0)
get_ipython().magic(u'pinfo SpectrumFit')
fit = SpectrumFit(obs_list=sim.obs, model=model, background_model=inbuiltbkg, stat='cash')
fit.run()
get_ipython().magic(u'pinfo SpectrumFit')
fit = SpectrumFit(obs_list=sim.obs, model=model, background_model=inbuiltbkg, stat='wstat')
bkg_model
fit = SpectrumFit(obs_list=sim.obs, model=model, background_model=bkg_model, stat='wstat')
fit = SpectrumFit(obs_list=sim.obs, model=model, background_model=background_estimator1, stat='wstat')
fit = SpectrumFit(obs_list=sim.obs, model=model, background_model=background_estimator, stat='wstat')
background_estimator1.run()
fit = SpectrumFit(obs_list=sim.obs, model=model, background_model=background_estimator1, stat='wstat')
background_estimator1.plot()
plt.show()
background_estimator.obs_list
print background_estimator.obs_list
print background_estimator1.obs_list
background_estimator.plot()
datastore.obs(obsid[0])
datastore.obs([obsid[0]])
mylist = datastore.obs(obsid[0])
mylist
mylist = datastore.obs((obsid[0]))
mylist
datastore.obs_list
get_ipython().magic(u'pinfo datastore.obs_list')
mylist = datastore.obs_list((obsid[0]))
mylist = datastore.obs_list(*(obsid[0]))
mylist = datastore.obs_list((obsid[0],))
mylist
mylist[0]
mylist[1]
sim
print sim.result
get_ipython().magic(u'pinfo np.random.choice')
nbkg=356
alpha=0.1
nOFF=356
nbkg=alpha*nOFF
nbkg = np.random.poisson(nbkg)
nbkg
np.choice(np.arange(nOFF),nbkg,replace=False)
np.rando;.choice(np.arange(nOFF),nbkg,replace=False)
np.random.choice(np.arange(nOFF),nbkg,replace=False)
idxON=np.random.choice(np.arange(nOFF),nbkg,replace=False)
get_ipython().magic(u'pinfo fakerun.bkg')
get_ipython().magic(u'pinfo fakerun')
sim.alpha
sim.off_vector
sim = SpectrumSimulation(aeff=aeff, edisp=edisp, source_model=model, livetime=livetime)
sim.simulate_obs(seed=42, obs_id=0)
sim.obs.peek()
sim.off_vector
sim.on_vector
sim.on_vector.energy
sim.on_vector.energy.bins
sim.on_vector.energy.hi
sim.on_vector.energy.lo
sim.on_vector
b1=fakerun.bkg
b1
print b1
print b1.offset=offset
inbuiltbkg
plt.plot(inbuiltbkg)
plt.show()
plt.plot(inbuiltbkg)
plt.show()
plt.loglog(inbuiltbkg['energy'],inbuiltbkg['value'])
plt.show()
get_ipython().magic(u'pinfo np.hypot')
background_estimator.finder
print background_estimator.finder
fakerun.check_observation()
fakerun
print fakerun
print fakerun.events
print fakerun.aeff
print fakerun.edisp
print fakerun.location
fakerun.data_store.info
print fakerun.data_store.info
print fakerun.muoneff
print fakerun.target_radec
print fakerun
crab
fakerun1.data
datastore
fakerun
fakerun1
fakerun.obs_id.view
fakerun.obs_id
fakerun.obs_id["Pointing pos"]
fakerun.obs_id.item
print fakerun.obs_id.item
print fakerun.pointing_radec
crab
SkyCoord.separation(fakerun.pointing_radec,crab)
SkyCoord.separation(fakerun1.pointing_radec,crab)
SkyCoord.separation(fakerun1[0].pointing_radec,crab)
SkyCoord.separation(fakerun1[1].pointing_radec,crab)
mylist
fakerun
fakerun1
mylist[0]
SkyCoord.separation(mylist[0].pointing_radec,crab)
sim.alpha
sim
print sim
len(fakerun1)
