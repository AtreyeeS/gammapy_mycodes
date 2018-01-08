# coding: utf-8
#importing!

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

# Setup the logger
import logging
logging.basicConfig()
log = logging.getLogger('gammapy.spectrum')
log.setLevel(logging.WARNING)

datastore = DataStore.from_dir("$HESS_DATA") 
crab=SkyCoord.from_name("Crab")
sep=SkyCoord.separation(crab,datastore.obs_table.pointing_radec)
get_ipython().magic(u'pinfo datastore.obs_table.pointing_galactic')
sep1=SkyCoord.separation(crab,datastore.obs_table.pointing_galactic)
datastore.obs_table
datastore.obs_table["TSTART"]> 55562.0
datastore.obs_table["TSTART"]> 55562.0 && datastore.obs_table["TSTOP"]<2011-01-01
datastore.obs_table["TSTART"]> 55562.0 & datastore.obs_table["TSTOP"]<2011-01-01
(datastore.obs_table["TSTART"]> 55562.0) & (datastore.obs_table["TSTOP"] < 2011-01-01)
(datastore.obs_table["TSTART"]> 55562.0) & (datastore.obs_table["TSTOP"] < 55927)
time = (datastore.obs_table["TSTART"]> 55562.0) & (datastore.obs_table["TSTOP"] < 55927)
crablist=datastore.obs_list(obsid)
crabrun=(datastore.obs_table[sep<2.0*u.deg]) 
obsid=crabrun['OBS_ID'].data
crablist=datastore.obs_list(obsid)
crab2011=crablist[time]
crab2011=crablist(time)
crab2011=crablist[time]
time
len(time)
crablist=datastore.obs_list[time]
crablist=datastore.obs_list[time]
crablist=datastore[time]
crabrun=(datastore.obs_table[sep<2.0*u.deg]) 
obsid=crabrun['OBS_ID'].data
crablist=datastore.obs_list(obsid)
crablist
crablist[1]
print crablist[1]
(crablist["TSTART"]> 55562.0) 
(crablist[1]["TSTART"]> 55562.0) 
sel=(time) & (sep<2.0)
time
sep<2.0*u.deg
sel=(time) & (sep<2.0*u.deg)
sel
crabrun=(datastore.obs_table[sel])
obsid=crabrun['OBS_ID'].data
crablist=datastore.obs_list(obsid)
crablist
np.where(time)
np.where(sel)
get_ipython().magic(u'pinfo plt')
import numpy as np
import matplotlib.pyplot as plt
x = np.arange(0, 5, 0.1);
y = np.sin(x
)
plt.plot(x, y)
plt.show()
time = (datastore.obs_table["TSTART"]> 55562.0) & (datastore.obs_table["TSTOP"] < 56927)
sel=(time) & (sep<2.0*u.deg)
sel
np.where(sel)
crabrun=(datastore.obs_table[sel])
on_region=CircleSkyRegion(crab,0.15 * u.deg)
model = models.LogParabola(
    alpha = 2.3,
    beta = 0,
    amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference = 1 * u.TeV,
)

flux_point_binning = EnergyBounds.equal_log_spacing(0.7, 30, 5, u.TeV)

exclusion_mask = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
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
on_region=CircleSkyRegion(crab,0.15 * u.deg)
model = models.LogParabola(
    alpha = 2.3,
    beta = 0,
    amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference = 1 * u.TeV,
)

flux_point_binning = EnergyBounds.equal_log_spacing(0.7, 30, 5, u.TeV)

exclusion_mask = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
import logging
logging.basicConfig()
log = logging.getLogger('gammapy.spectrum')
log.setLevel(logging.WARNING)
from gammapy.spectrum import models
on_region=CircleSkyRegion(crab,0.15 * u.deg)
model = models.LogParabola(
    alpha = 2.3,
    beta = 0,
    amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference = 1 * u.TeV,
)

flux_point_binning = EnergyBounds.equal_log_spacing(0.7, 30, 5, u.TeV)

exclusion_mask = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
ana = SpectrumAnalysisIACT(
        observations=crabrun, config=config,)
from gammapy.data import DataStore
from gammapy.scripts import SpectrumAnalysisIACT
ana = SpectrumAnalysisIACT(
        observations=crabrun, config=config,)
config = dict(
    outdir = None,
    background = dict(
        on_region=on_region,
        exclusion_mask=exclusion_mask,
        min_distance = 0.1 * u.rad,
    ),
    extraction = dict(containment_correction=False),
    fit = dict(
        model=model,
        stat='wstat',
        forward_folded=True,
        fit_range = flux_point_binning[[0, -1]]
    ),
    fp_binning=flux_point_binning
)
ana = SpectrumAnalysisIACT(
        observations=crabrun, config=config,)
ana.run()
crabrun=(datastore.obs_table[sel])
obsid=crabrun['OBS_ID'].data
crablist=datastore.obs_list(obsid)
ana = SpectrumAnalysisIACT(
        observations=crablist,
        config=config,)

ana.run()
EXCLUSION_FILE = '$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits'

allsky_mask = SkyImage.read(EXCLUSION_FILE)
exclusion_mask = allsky_mask.cutout(
    position=on_region.center,
    size=Angle('6 deg'),
)
background_estimator = ReflectedRegionsBackgroundEstimator(
    obs_list=crablist,on_region=on_region,
    exclusion_mask = exclusion_mask)
background_estimator.run()
print(background_estimator.result[0])
plt.figure(figsize=(8,8))
background_estimator.plot()
plt.show()
plt.show()
plt.figure(figsize=(8,8))
background_estimator.plot()
plt.show()
stats = []
for obs, bkg in zip(crablist,background_estimator.result):
    stats.append(ObservationStats.from_obs(obs,bkg))
    
stats
stats[0]
stats[1]
print stats[1]
print stats[0]
len(stats)
print stats[44]
print stats[:]["Sigma"]
print stats[1]["Sigma"]
print stats[1].sigma
print stats[1]
print stats[0:].sigma
print stats[1:5].sigma
sig=[]
for i in range(len(stats)):
    sig.append(stats[i].sigma)
    
sig
obs_summary = ObservationSummary(stats)
obs_summary
obs_summary.sigma
obs_summary.bg_rate
fig = plt.figure(figsize=(10,6))
ax1=fig.add_subplot(121)
obs_summary.plot_excess_vs_livetime(ax=ax1)
ax2=fig.add_subplot(122)
obs_summary.plot_significance_vs_livetime(ax=ax2)
plt.plot()
plt.show()
e_reco = EnergyBounds.equal_log_spacing(0.1, 40, 40, unit='TeV')
e_true = EnergyBounds.equal_log_spacing(0.05, 100., 200, unit='TeV')
ANALYSIS_DIR = 'crab_analysis'
extraction = SpectrumExtraction(obs_list=crablist, bkg_estimate=background_estimator.result, containment_correction=False,)
extraction.run()
extraction.compute_energy_threshold(method_lo='area_max', area_percent_lo=10.0)

print(extraction.observations[0])
extraction.run(obs_list=obs_list, bkg_estimate=background_estimator.result, outdir=ANALYSIS_DIR)
extraction.run(obs_list=crablist, bkg_estimate=background_estimator.result, outdir=ANALYSIS_DIR)
extraction.write(obs_list=crablist, bkg_estimate=background_estimator.result, outdir=ANALYSIS_DIR)
extraction.write(crablist)
extraction.write(crablist, bkg_estimate=background_estimator.result, outdir=ANALYSIS_DIR)
extraction.write(crablist, background_estimator.result, outdir=ANALYSIS_DIR)
get_ipython().magic(u'pinfo extraction.write')
extraction.observations[0].peek()
model = PowerLaw(
    index=2 * u.Unit(''),
    amplitude=2e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference=1 * u.TeV,
)

joint_fit = SpectrumFit(obs_list=extraction.observations, model=model)
joint_fit.fit()
joint_fit.est_errors()
#fit.run(outdir = ANALYSIS_DIR)

joint_result = joint_fit.result
ax0, ax1 = joint_result[0].plot(figsize=(8,8))
ax0.set_ylim(0, 20)
print(joint_result[0])
plt.show()
ebounds = [0.3, 1.1, 3, 10.1, 30] * u.TeV

stacked_obs = extraction.observations.stack()

seg = SpectrumEnergyGroupMaker(obs=stacked_obs)
seg.compute_range_safe()
seg.compute_groups_fixed(ebounds=ebounds)

print(seg.groups)
fpe = FluxPointEstimator(
    obs=stacked_obs,
    groups=seg.groups,
    model=joint_result[0].model,
)
fpe.compute_points()
fpe.flux_points.plot()
fpe.flux_points.table
spectrum_result = SpectrumResult(
    points=fpe.flux_points,
    model=joint_result[0].model,
)
ax0, ax1 = spectrum_result.plot(
    energy_range=joint_fit.result[0].fit_range,
    energy_power=2, flux_unit='erg-1 cm-2 s-1',
    fig_kwargs=dict(figsize=(8,8)),
    point_kwargs=dict(color='navy')
)

ax0.set_xlim(0.4, 50)
plt.show()
stacked_obs = extraction.observations.stack()

stacked_fit = SpectrumFit(obs_list=stacked_obs, model=model)
stacked_fit.fit()
stacked_fit.est_errors()


stacked_result = stacked_fit.result
print(stacked_result[0])

stacked_table = stacked_result[0].to_table(format='.3g')
stacked_table['method'] = 'stacked'
joint_table = joint_result[0].to_table(format='.3g')
joint_table['method'] = 'joint'
total_table = vstack_table([stacked_table, joint_table])
print(total_table['method', 'index', 'index_err', 'amplitude', 'amplitude_err'])
