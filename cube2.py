import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from regions import CircleSkyRegion

from gammapy.extern.pathlib import Path
from gammapy.irf import EnergyDispersion
from gammapy.cube import SkyCube
from gammapy.cube.sherpa_ import (
    CombinedModel3DInt,
    CombinedModel3DIntConvolveEdisp,
    NormGauss2DInt,
)

from sherpa.models import PowLaw1D, TableModel
from sherpa.estmethods import Covariance
from sherpa.optmethods import NelderMead
from sherpa.stats import Cash
from sherpa.fit import Fit

cube=counts.to_sherpa_data3d(dstype='Data3DInt')

bkg = TableModel('bkg')
bkg.load(None, background.data.value.ravel())
bkg.ampl = 1
bkg.ampl.freeze()

exposure_3d=exposure
exposure_3d = SkyCube.read(cube_dir / 'exposure_cube.fits')
i_nan = np.where(np.isnan(exposure_3d.data))
exposure_3d.data[i_nan] = 0
# In order to have the exposure in cm2 s
exposure_3d.data = exposure_3d.data * u.Unit('m2 / cm2').to('')
psf_3d=mean_psf_cube

# Define a 2D gaussian for the spatial model
spatial_model = NormGauss2DInt('spatial-model')

# Define a power law for the spectral model
spectral_model = PowLaw1D('spectral-model')

# Combine spectral and spatial model
coord = counts.sky_image_ref.coordinates(mode="edges")
energies = counts.energies(mode='edges').to("TeV")
# Here the source model will be convolve by the psf:
# PSF * (source_model * exposure)
source_model = CombinedModel3DInt(
    coord=coord,
    energies=energies,
    use_psf=True,
    exposure=exposure_3d,
    psf=psf_3d,
    spatial_model=spatial_model,
    spectral_model=spectral_model,
)

center = SkyCoord(83.633083, 22.0145, unit="deg").galactic
source_model.gamma = 2.2
source_model.xpos = center.l.value
source_model.ypos = center.b.value
source_model.fwhm = 0.12
source_model.ampl = 1.0

# Define the model
flux_factor = 1e-11
model = bkg + flux_factor * source_model

# Fit to the counts Cube sherpa object
fit = Fit(
    data=cube,
    model=model,
    stat=Cash(),
    method=NelderMead(),
    estmethod=Covariance(),
)
fit_results = fit.fit()
print(fit_results.format())

err_est_results = fit.est_errors()
print(err_est_results.format())

#now for the mask....

exclusion_region = CircleSkyRegion(
    center=SkyCoord(83.60, 21.88, unit='deg'),
    radius=Angle(0.1, "deg"),
)

sky_mask_cube = counts.region_mask(exclusion_region)
sky_mask_cube.data = sky_mask_cube.data.astype(int)
index_region_selected_3d = np.where(sky_mask_cube.data.value == 1)

# Set the counts
cube = counts_3d.to_sherpa_data3d(dstype='Data3DInt')
cube.mask = sky_mask_cube.data.value.ravel()

# Set the bkg and select only the data points of the selected region
bkg = TableModel('bkg')
bkg.load(None, background.data.value[index_region_selected_3d].ravel())
bkg.ampl = 1
bkg.ampl.freeze()

# The model is evaluated on all the points then it is compared with the data only on the selected_region
source_model = CombinedModel3DInt(
    coord=coord,
    energies=energies,
    use_psf=True,
    exposure=exposure_3d,
    psf=psf_3d,
    spatial_model=spatial_model,
    spectral_model=spectral_model,
    select_region=True,
    index_selected_region=index_region_selected_3d,
)

# Set starting values
source_model.gamma = 2.2
source_model.xpos = center.l.value
source_model.ypos = center.b.value
source_model.fwhm = 0.12
source_model.ampl = 1.0

# Define the model
model = bkg + flux_factor * (source_model)

fit = Fit(
    data=cube,
    model=model,
    stat=Cash(),
    method=NelderMead(),
    estmethod=Covariance(),
)
fit_results = fit.fit()
print(fit_results.format())

err_est_results = fit.est_errors()

print(err_est_results.format())

#now with energy dispersion

# Set the counts
counts_3d = SkyCube.read(cube_dir / "counts_cube.fits")
cube = counts_3d.to_sherpa_data3d(dstype='Data3DInt')

# Set the bkg
bkg_3d = SkyCube.read(cube_dir / 'bkg_cube.fits')
bkg = TableModel('bkg')
bkg.load(None, bkg_3d.data.value.ravel())
bkg.ampl = 1
bkg.ampl.freeze()

# Set the exposure
exposure_3d = SkyCube.read(cube_dir / 'exposure_cube_etrue.fits')
i_nan = np.where(np.isnan(exposure_3d.data))
exposure_3d.data[i_nan] = 0
exposure_3d.data = exposure_3d.data * 1e4

# Set the mean psf model
psf_3d = SkyCube.read(cube_dir / 'psf_cube_etrue.fits')

# Load the mean rmf calculated for the 4 Crab runs
rmf = EnergyDispersion.read(cube_dir / 'rmf.fits')

# Setup combined spatial and spectral model
spatial_model = NormGauss2DInt('spatial-model')
spectral_model = PowLaw1D('spectral-model')

coord = counts_3d.sky_image_ref.coordinates(mode="edges")
energies = counts_3d.energies(mode='edges').to("TeV")
source_model = CombinedModel3DIntConvolveEdisp(
    coord=coord,
    energies=energies,
    use_psf=True,
    exposure=exposure_3d,
    psf=psf_3d,
    spatial_model=spatial_model,
    spectral_model=spectral_model,
    edisp=rmf.data.data,
)

# Set starting values
center = SkyCoord(83.633083, 22.0145, unit="deg").galactic
source_model.gamma = 2.2
source_model.xpos = center.l.value
source_model.ypos = center.b.value
source_model.fwhm = 0.12
source_model.ampl = 1.0

# Define the model
model = bkg + flux_factor * source_model

fit = Fit(
    data=cube,
    model=model,
    stat=Cash(),
    method=NelderMead(),
    estmethod=Covariance(),
)
fit_results = fit.fit()
print(fit_results.format())

err_est_results = fit.est_errors()
print(err_est_results.format())


