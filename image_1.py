import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
from astropy.visualization import simple_norm

from gammapy.data import DataStore
from gammapy.image import SkyImage
from gammapy.scripts import StackedObsImageMaker
from gammapy.utils.energy import Energy
import copy
#choose the obsevartion
name="Crab"
datastore = DataStore.from_dir("$HESS_DATA") 
src=SkyCoord.from_name(name)
sep=SkyCoord.separation(src,datastore.obs_table.pointing_radec)
srcruns=(datastore.obs_table[sep<2.5*u.deg]) 
obsid=srcruns['OBS_ID'].data
#mylist=datastore.obs_list((obsid[0],))
mylist=datastore.obs_list(obsid[0:10])
obs_table=Table()
obs_table['OBS_ID']=obsid
# Define sky image
ref_image = SkyImage.empty(name="pks1",
    nxpix=400, nypix=400, binsz=0.02,
    xref=src.ra.deg, yref=src.dec.deg,
    proj='TAN', coordsys='CEL',
)

# Define energy band
energy_band = Energy([0.1, 10], 'TeV')

# Define maximum field of view offset cut
offset_band = Angle([0, 2.49], 'deg')

# Define exclusion mask (known gamma-ray sources)
# This is used in the background model image estimation
exclusion_mask = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
exclusion_mask = exclusion_mask.reproject(reference=ref_image)
mask1=copy.copy(exclusion_mask)
mask1.data=np.invert(exclusion_mask.data.astype(bool))

image_maker = StackedObsImageMaker(
    empty_image=ref_image,
    energy_band=energy_band,
    offset_band=offset_band,
    data_store=datastore,
    obs_table=obs_table,
    exclusion_mask=mask1
)

image_maker.make_images(
    make_background_image=True,
    for_integral_flux=True,
    radius=4,
)

counts_image0 = image_maker0.images['counts']
norm = simple_norm(counts_image0.data, stretch='sqrt', min_cut=0, max_cut=0.9)
counts_image0.smooth(radius=0.08 * u.deg).plot(norm=norm, add_cbar=True)
plt.show()

background_image = image_maker.images['bkg']
norm = simple_norm(background_image.data, stretch='sqrt', min_cut=0, max_cut=0.2)
background_image.plot(norm=norm, add_cbar=True)
plt.show()

excess_image = image_maker.images['excess']
norm = simple_norm(excess_image.data, stretch='sqrt', min_cut=0, max_cut=0.9)
excess_image.smooth(radius=0.08 * u.deg).plot(norm=norm,add_cbar=True)
plt.show()

image_maker.images["significance"].plot(add_cbar=True)
plt.show()

image_maker.images["exposure"].plot(add_cbar=True)

plt.show()
