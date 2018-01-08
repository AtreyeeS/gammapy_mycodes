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
from gammapy.spectrum import models
# Setup the logger
import logging
logging.basicConfig()
log = logging.getLogger('gammapy.spectrum')
log.setLevel(logging.WARNING)

datastore = DataStore.from_dir("$HESS_DATA") 
crab=SkyCoord.from_name("Crab")
sep=SkyCoord.separation(crab,datastore.obs_table.pointing_radec)
crabrun=(datastore.obs_table[sep<2.0*u.deg]) 
obsid=crabrun['OBS_ID'].data
crablist=datastore.obs_list(obsid)

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
        observations=crablist,
        config=config,)

ana.run()
