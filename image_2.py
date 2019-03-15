
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.convolution import Ring2DKernel, Tophat2DKernel
from astropy.visualization import simple_norm

from gammapy.data import DataStore
from gammapy.image import SkyImage, SkyImageList
from gammapy.detect import KernelBackgroundEstimator as KBE

#choose the obsevartion
name="PKS 2155-304"
datastore = DataStore.from_dir("$HESS_DATA") 
src=SkyCoord.from_name(name)
sep=SkyCoord.separation(src,datastore.obs_table.pointing_radec)
srcruns=(datastore.obs_table[sep<2.0*u.deg]) 
obsid=srcruns['OBS_ID'].data
#mylist=datastore.obs_list((obsid[0],))
mylist=datastore.obs_list(obsid[0:5])
obs_table=Table()
obs_table['OBS_ID']=obsid[0:5]


ref_image = SkyImage.empty(
    nxpix=400, nypix=400, binsz=0.02,
    xref=src.ra.deg, yref=src.dec.deg,
    coordsys='CEL', proj='TAN',
)

# Make a counts image for a single observation
events = datastore.obs(obs_id=obsid[0]).events
counts_image = SkyImage.empty_like(ref_image)
counts_image.fill_events(events)


norm = simple_norm(counts_image.data, stretch='sqrt', min_cut=0, max_cut=0.3)
counts_image.smooth(radius=0.1 * u.deg).plot(norm=norm, add_cbar=True)

#kernel background???

source_kernel = Tophat2DKernel(radius=5)
source_kernel.normalize(mode='peak')
source_kernel = source_kernel.array

background_kernel = Ring2DKernel(radius_in=20, width=10)
background_kernel.normalize(mode='peak')
background_kernel = background_kernel.array
plt.imshow(source_kernel, interpolation='nearest', cmap='gray')
plt.colorbar()
plt.grid('off')
plt.imshow(background_kernel, interpolation='nearest', cmap='gray')
plt.colorbar()
plt.grid('off')
obs_ids=obsid[0:3]
counts_image2 = SkyImage.empty_like(ref_image)
for obs_id in obs_ids:
    events = datastore.obs(obs_id=obs_id).events
    counts_image2.fill_events(events)

# To use the `KernelBackgroundEstimator` you first have to set
# up a source and background kernel and put the counts image input
# into a container `SkyImageList` class.
images = SkyImageList()
images['counts'] = counts_image2

kbe = KBE(
    kernel_src=source_kernel,
    kernel_bkg=background_kernel,
    significance_threshold=5,
    mask_dilation_radius=0.06 * u.deg,
)
# This takes about 10 seconds on my machine
result = kbe.run(images)

background_image = result['background']
norm = simple_norm(background_image.data, stretch='sqrt', min_cut=0, max_cut=0.5)
background_image.plot(norm=norm, add_cbar=True)

significance_image = result['significance']
significance_image.plot(add_cbar=True, vmin=-3, vmax=20)

import regions 
exclusion_region = regions.CircleSkyRegion(src,0.3*u.deg)
mask = counts_image.region_mask(exclusion_region)

#background simulations...


