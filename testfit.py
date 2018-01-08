# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from gammapy.spectrum import models
from gammapy.utils.modeling import Parameter, ParameterList
pwl=models.PowerLaw()
pwl
print pwl
pwl.parameters.names
pwl.parameters.names[index]=2.3
pwl.parameters.names["index"]=2.3
pwl.parameters["index"]=2.3
pwl["index"]=2.3
pwl("index")=2.3
pwl("index")
pwl["index"]
pwl.parameters
print pwl.parameters
pwl.parameters["index"].value=2.3
pwl.parameters["index"].min=1.0
print pwl
pwl.parameters["index"].value=2.6
print pwl
pwl.parameters["index"]
print pwl.parameters["index"]
pwl.parameters["index"].min=1.0
pwl.parameters["index"].min
pwl.parameters["index"].max
pwl.parameters["index"].max=3.0
pwl.parameters["index"].max
print pwl.parameters["index"]
models.SpectralModel.__subclasses__()
energy_range = [0.1, 10] * u.TeV
pwl.plot(energy_range, energy_power=2, energy_unit='GeV')
plt.show()
errors = dict(
    index = 0.2 * u.Unit(''),
    amplitude = 0.1 * pwl.parameters['amplitude'].quantity
)
pwl.parameters.set_parameter_errors(errors)
print(pwl)
print(pwl.parameters.covariance)
print(pwl.parameters.error('index'))
ax = pwl.plot_error(energy_range, color='blue', alpha=0.2)
pwl.plot(energy_range, ax=ax, color='blue')
import uncertainties
print('gammapy:', gammapy.__version__)
import gammapy
print('gammapy:', gammapy.__version__)
import uncertainties
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
datastore = DataStore.from_dir($HESS_DATA)
datastore = DataStore.from_dir("$HESS_DATA")
datastore
datastore.obs
print datastore.obs
print datastore.obs_list
print datastore.obs_table
print datastore.DEFAULT_NAME
print datastore.DEFAULT_OBS_TABLE
from gammapy.data import DataStore, ObservationTable
ObservationTable.select_observations(datastore)
get_ipython().magic(u'pinfo ObservationTable.select_observations')
get_ipython().magic(u'pinfo datastore.obs')
datastore.obs(22222)
EXCLUSION_FILE = '$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits'
allsky_mask = SkyImage.read(EXCLUSION_FILE)
exclusion_mask = allsky_mask.cutout(
    position=on_region.center,
    size=Angle('6 deg'),
)
obs_id = data_store.obs_table['OBS_ID'].data
obs_id = datastore.obs_table['OBS_ID'].data
obs_id
print("Use observations {}".format(obs_id))
datastore.obs_table['OBS_ID'].data
datastore.obs_table['OBS_ID']
l1=[122091,122101]
datastore.obs_list(l1)
d1=datastore.obs_list(l1)
d1
d1.obs_table['OBS_ID']
d1
print d1
crab_pos = SkyCoord.from_name('crab')
crab_pos
datastore
datastore.table
datastore.Table
datastore.obs_table
datastore.obs_table["ZEN_PNT"]<20.0
from gammapy.data import ObservationGroups, ObservationGroupAxis
zenith = Angle([0, 30, 40, 50], 'deg')
ntels = [3, 4]
obs_groups = ObservationGroups([
    ObservationGroupAxis('ZENITH', zenith, fmt='edges'),
    ObservationGroupAxis('N_TELS', ntels, fmt='values'),
])
print(obs_groups.info)
obs1=obs_groups.apply(datastore.obs_table)
zenith = Angle([0, 30, 40, 50], 'deg')
ntels = [3, 4]
obs_groups = ObservationGroups([
    ObservationGroupAxis('ZEN_PNT', zenith, fmt='edges'),
    ObservationGroupAxis('N_TELS', ntels, fmt='values'),
])
obs_table

zenith = Angle([0, 30, 40, 50], 'deg')
ntels = [3, 4]
obs_groups = ObservationGroups([
    ObservationGroupAxis('ZEN_PNT', zenith, fmt='edges'),
    ObservationGroupAxis('N_TELS', ntels, fmt='values'),
])
obs1=obs_groups.apply(datastore.obs_table)
obs1
print(obs1)
print(datastore.obs_table)
obs_gr=obs_groups.get_group_of_observations(obs1,2)
obs_gr
datastore
datastore.obs_table
datastore.obs_table.Column
datastore.obs_table.field
crab=SkyCoord.from_name("Crab")
crab
crab_table=datastore.obs_table(SkyCoord.separation(crab))
crab_table=datastore.obs_table(SkyCoord.separation(crab,2.0))
get_ipython().magic(u'pinfo SkyCoord.separation')
SkyCoord.separation(datastore.obs_table,crab)
SkyCoord.separation(crab,datastore.obs_table)
geminga=SkyCoord.from_name("geminga")
SkyCoord.separation(geminga,crab)
get_ipython().magic(u'pinfo datastore.obs_table.group_by')
get_ipython().magic(u'pinfo datastore.obs_table')
datastore.obs_table.keys
get_ipython().magic(u'pinfo datastore.obs_table.keys')
get_ipython().magic(u'pinfo datastore.obs_table.pointing_radec')
datastore.obs_table.pointing_radec
SkyCoord.separation(crab,datastore.obs_table.pointing_radec))
SkyCoord.separation(crab,datastore.obs_table.pointing_radec)
sep=SkyCoord.separation(crab,datastore.obs_table.pointing_radec)
crab_table=datastore(datastore.obs_table[sep<2.0])
(datastore.obs_table[sep<2.0])
c=[sep<2.0]
c=sep<2.0
c=sep<2.0*u.deg
(datastore.obs_table[sep<2.0*u.deg])
crab_table=datastore(datastore.obs_table[sep<2.0*u.deg])
crabrun=(datastore.obs_table[sep<2.0*u.deg])
crabrun
crab
crab*u.gal
crab*u.alt
crab
crab
DATA_DIR = '$GAMMAPY_EXTRA/datasets/hess-crab4-hd-hap-prod2'
d1=DataStore.from_dir(DATA_DIR)
to=d1.obs_list([23523, 23526, 23559, 23592])
to
d1.obs_table
crab
crabrun
crabrun=(datastore.obs_table[sep<2.0*u.deg])
crabrun=(datastore[sep<2.0*u.deg])
on_region
crabrun=(datastore[sep<0.2*u.deg])
crabrun=(datastore.obs_table[sep<0.2*u.deg])
crabrun
crabrun=(datastore.obs_table[sep<2.0*u.deg])
crabrun
on_region=CircleSkyRegion(crab,0.15 * u.deg)
model = models.LogParabola(
    alpha = 2.3,
    beta = 0,
    amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference = 1 * u.TeV,
)

flux_point_binning = EnergyBounds.equal_log_spacing(0.7, 30, 5, u.TeV)

exclusion_mask = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
model = models.LogParabola(
    alpha = 2.3,
    beta = 0,
    amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference = 1 * u.TeV,
)

flux_point_binning = EnergyBounds.equal_log_spacing(0.7, 30, 5, u.TeV)

exclusion_mask = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits'


)
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
    observations=crabrun,
    config=config,
)
ana = SpectrumAnalysisIACT(
    observations=obs_list,
    config=config,
)

get_ipython().magic(u'matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt

from gammapy.data import DataStore
from gammapy.scripts import SpectrumAnalysisIACT

import numpy as np
import matplotlib.pyplot as plt

from gammapy.data import DataStore
from gammapy.scripts import SpectrumAnalysisIACT
ana = SpectrumAnalysisIACT(
    observations=obs_list,
    config=config,
)
ana = SpectrumAnalysisIACT(
    observations=crabrun,
    config=config,
)
ana.run
ana.run()
crabrun
crabrun.get_obs_idx
crabrun
get_ipython().magic(u'pinfo crabrun')
get_ipython().magic(u'pinfo datastore.obs_table')
get_ipython().magic(u'pinfo datastore.obs_list')
o1
d1
d1.obs_list([23523 23526 23559 23592]
)
d1.obs_list([23523 23526 23559 23592])
obsid=crabrun['OBS_ID'].data
obsid
crablist=datastore.obs_list(obsid)
ana = SpectrumAnalysisIACT(
    observations=crablist,
    config=config,
)
ana.run()
print(ana.fit.result[0])
ana.spectrum_result.plot(
    energy_range=ana.fit.fit_range,
    energy_power=2,
    flux_unit='erg-1 cm-2 s-1',
    fig_kwargs=dict(

))
ana.spectrum_result.plot(
    energy_range=ana.fit.fit_range,
    energy_power=2,
    flux_unit='erg-1 cm-2 s-1',
    fig_kwargs=dict(figsize = (8,8)),
)
plt.show()
plt.show()
plt.plot()
plt.show()
a=np.linspace(1,10)
b=a
plt.plot(a,b)
plt.show()
import matplotlib.pyplot as plt
ana.spectrum_result.plot(
    energy_range=ana.fit.fit_range,
    energy_power=2,
    flux_unit='erg-1 cm-2 s-1',
    fig_kwargs=dict(figsize = (8,8)),
)
plt.show()
plt.plot()
x = np.arange(0, 5, 0.1);
y = np.sin(x)
plt.plot(x, y)
plt.show()
allsky_mask = SkyImage.read(EXCLUSION_FILE)
exclusion_mask = allsky_mask.cutout(
    position=on_region.center,
    size=Angle('6 deg'),
)
background_estimator = ReflectedRegionsBackgroundEstimator(obs_list=crablist, 
on_region=on_region, exclusion_mask = exclusion_mask)
background_estimator.run()
print(background_estimator.result[0])
plt.figure(figsize=(8,8))
background_estimator.plot()
plt.plot()
plt.show()
get_ipython().magic(u'save testfit.py 1-194')
