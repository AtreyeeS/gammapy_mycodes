import matplotlib.pyplot as plt
import logging
logging.getLogger().setLevel(logging.ERROR)
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.visualization import simple_norm

from gammapy.data import DataStore
from gammapy.image import SkyImage
from gammapy.cube import SkyCube, StackedObsCubeMaker

# define WCS specification
WCS_SPEC = {'nxpix': 50,
            'nypix': 50,
            'binsz': 0.05,
            'xref': 83.63,
            'yref': 22.01,
            'proj': 'TAN',
            'coordsys': 'CEL'}

# define reconstructed energy binning
ENERGY_SPEC = {'mode': 'edges',
               'enumbins': 5,
               'emin': 0.5,
               'emax': 40,
               'eunit': 'TeV'}

REF_CUBE = SkyCube.empty(emin=0.5, emax=40.0, enumbins=5, **WCS_SPEC)
# setting up the data store
data_store = DataStore.from_dir("$GAMMAPY_EXTRA/test_datasets/cube/data")

# temporary fix for load psftable for one of the run that is not implemented yet...
data_store.hdu_table.remove_row(14)
# read in TeVCat exclusion mask
exclusion_mask = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')

# reproject exclusion mask to reference cube
exclusion_mask = exclusion_mask.reproject(reference=REF_CUBE.sky_image_ref, order='nearest-neighbor')
exclusion_mask.show()

#Select the offset band on which you want to select the events in the FOV of each observation
offset_band = Angle([0, 2.49], 'deg')

# instanciate StackedObsCubeMaker
cube_maker = StackedObsCubeMaker(
    empty_cube_images=REF_CUBE,
    empty_exposure_cube=REF_CUBE,
    offset_band=offset_band,
    data_store=data_store,
    obs_table=data_store.obs_table,
    exclusion_mask=exclusion_mask,
    save_bkg_scale=True,
)

# run cube maker
cube_maker.make_cubes(make_background_image=True, radius=4.)

# get counts cube
counts = cube_maker.counts_cube

energies = counts.energies(mode="edges")
print("Bin edges of the count cubes: \n{0}".format(energies))

# show counts image with idx = 0
idx = 0
counts.sky_image_idx(idx=idx).show()
print('Energy range: {emin:.2f} to {emax:.2f}'.format(emin=energies[idx], emax=energies[idx + 1]))

#background
background = cube_maker.bkg_cube

idx = 0
background.sky_image_idx(idx=idx).show()
print('Energy range: {emin:.2f} to {emax:.2f}'.format(emin=energies[idx], emax=energies[idx + 1]))

#excess
excess = cube_maker.excess_cube

idx = 0
excess.sky_image_idx(idx=idx).show()
print('Energy range: {emin:.2f} to {emax:.2f}'.format(emin=energies[idx], emax=energies[idx + 1]))

significance = cube_maker.significance_cube

idx = 0
significance.sky_image_idx(idx=idx).show()
print('Energy range: {emin:.2f} to {emax:.2f}'.format(emin=energies[idx], emax=energies[idx + 1]))

exposure = cube_maker.exposure_cube

idx = 0
exposure.sky_image_idx(idx=idx).show()
print('Energy range: {emin:.2f} to {emax:.2f}'.format(emin=energies[idx], emax=energies[idx + 1]))

#to make a mean psf cube
with np.errstate(divide='ignore', invalid='ignore'):
    mean_psf_cube = cube_maker.make_mean_psf_cube(REF_CUBE)

idx = 3
psf_image = mean_psf_cube.sky_image_idx(idx=idx)
norm = simple_norm(psf_image.data, stretch='log')
psf_image.show(norm=norm) #NOT WORKING
print('Energy range: {emin:.2f} to {emax:.2f}'.format(emin=energies[idx], emax=energies[idx + 1]))

# repeat for true...
REF_CUBE_TRUE = SkyCube.empty(emin=0.1, emax=100, enumbins=20, **WCS_SPEC)

with np.errstate(divide='ignore', invalid='ignore'):
    mean_psf_cube = cube_maker.make_mean_psf_cube(REF_CUBE_TRUE)

e_reco = cube_maker.counts_cube.energies(mode="edges")
print("Energy bin edges of the count cube: \n{}\n".format(e_reco))

e_true = cube_maker.exposure_cube.energies(mode="edges")
print("Energy bin edges of the exposure cube: \n{}".format(e_true))

# get reference energies
obs_list = data_store.obs_list(data_store.obs_table['OBS_ID'])
position = REF_CUBE.sky_image_ref.center
e_true = REF_CUBE_TRUE.energies('edges')
e_reco = REF_CUBE.energies('edges')

with np.errstate(divide='ignore', invalid='ignore'):
    mean_rmf = obs_list.make_mean_edisp(position=position, e_true=e_true, e_reco=e_reco)

mean_rmf.plot_matrix()

