# coding: utf-8
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
from gammapy.irf import EnergyDispersion, EnergyDispersion2D, EffectiveAreaTable, EffectiveAreaTable2D
from gammapy.irf import EnergyDependentMultiGaussPSF 
from gammapy.spectrum import SpectrumSimulation, SpectrumFit
from gammapy.spectrum.models import PowerLaw
import os
from gammapy.background.reflected import ReflectedRegionsBackgroundEstimator
from gammapy.spectrum import PHACountsSpectrum


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
livetime_src = 10.0 * u.hr
livetime = mylist[0].observation_live_time_duration
offset = SkyCoord.separation(mylist[0].pointing_radec,src)
lo_threshold = 0.1 * u.TeV
hi_threshold = 60 * u.TeV

# Define spectral model
index = 2.0 * u.Unit('')
amplitude = 2.5 * 1e-10 * u.Unit('cm-2 s-1 TeV-1')
reference = 1 * u.TeV
model = PowerLaw(index=index, amplitude=amplitude, reference=reference)

#do a user defined model too!





edisp=mylist[0].edisp.to_energy_dispersion(offset=offset)
aeff=mylist[0].aeff.to_effective_area_table(offset=offset)

aeff.lo_threshold = lo_threshold
aeff.hi_threshold = hi_threshold



#inbuiltbkg=fakerun.bkg.evaluate_at_offset(offset=offset)
#fig, axes = plt.subplots(1, 2, figsize=(12, 6))
#edisp.plot_matrix(ax=axes[0])
#aeff.plot(ax=axes[1])
#plt.show()

#aeff.lo_threshold = lo_threshold
#aeff.hi_threshold = hi_threshold

n_obs=1
seeds = np.arange(n_obs)
sim = SpectrumSimulation(aeff=aeff, edisp=edisp, source_model=model, livetime=livetime_src)
sim.run(seeds)

n_on=sim.result[0].total_stats.n_on

# now fitting different backgrounds.

EXCLUSION_FILE = '$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits'
allsky_mask = SkyImage.read(EXCLUSION_FILE)
exclusion_mask = allsky_mask.cutout(
    position=src,
    size=Angle('6 deg'),
)
on_region=CircleSkyRegion(src,0.11 * u.deg)
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=mylist, on_region=on_region, exclusion_mask = exclusion_mask)

background_estimator.run()
#background_estimator.plot()
#plt.show()

bkg_res=background_estimator.result[0]

a_off=bkg_res.a_off
a_on = bkg_res.a_on
nOFF= len(bkg_res.off_events.table)
alpha=float(a_on)/float(a_off)
nbkg=alpha*nOFF
nbkg = np.random.poisson(nbkg)
idxON=np.random.choice(np.arange(nOFF),nbkg,replace=False)
excess=n_on+nbkg - alpha*nOFF

print "excess=", excess

bkg_ev=bkg_res.off_events.table[idxON]
bkg_hist,edge=np.histogram(bkg_ev["ENERGY"],sim.on_vector.energy.bins)
sim.on_vector.data.data += bkg_hist * u.ct

off_events=bkg_res.off_events.table
off_events.remove_rows(idxON) # should these be removed ?
alpha_new=float(a_on)/float(a_off - 1)

off_counts,edge=np.histogram(off_events["ENERGY"],sim.on_vector.energy.bins)
off_vector = PHACountsSpectrum(energy_lo=sim.on_vector.energy.lo,
                            energy_hi=sim.on_vector.energy.hi,
                            data=off_counts,
                            backscal=1. / alpha_new,
                            is_bkg=True)


off_vector.livetime = livetime
sim.result[0].off_vector = off_vector
obs=sim.result[0]

#extraction = SpectrumExtraction(
#    obs_list=sim.result,
#    bkg_estimate=background_estimator.result,
#    containment_correction=False,
#)


fit = SpectrumFit(obs, model=model, stat='wstat')
fit.model.parameters['index'].value = 2
fit.run()

fit.est_errors()

print fit.result[0]
ax0, ax1 = fit.result[0].plot(figsize=(8,8))

plt.show() 

#for flux points--

# how to change livetime ?




