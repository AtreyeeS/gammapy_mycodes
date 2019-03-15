import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
from astropy.visualization import simple_norm

from gammapy.data import DataStore
from gammapy.image import SkyImage
from gammapy.utils.energy import Energy,  EnergyBounds
import regions

import scipy.stats as ss
import copy
#get the observations
name="PKS 2155-304"
#name="Mkn 421"
#name="1ES 0347-121"
datastore = DataStore.from_dir("$HESS_DATA") 
src=SkyCoord.from_name(name)
sep=SkyCoord.separation(src,datastore.obs_table.pointing_radec)
srcruns=(datastore.obs_table[sep<1.5*u.deg]) 
obsid=srcruns['OBS_ID'].data
good_obs=np.loadtxt("PKS2155_PA_201801.list.txt",usecols=(0))
obsid1=np.intersect1d(good_obs,obsid)
mylist=datastore.obs_list(obsid1)

ref_image = SkyImage.empty(
    nxpix=400, nypix=400, binsz=0.02,
    xref=src.ra.deg, yref=src.dec.deg,
    coordsys='CEL', proj='TAN',
)

#exclusion_mask_tevcat = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
#exclusion_mask_tevcat = exclusion_mask_tevcat.reproject(reference=ref_image)

energy_band = Energy([0.5, 1.0], 'TeV')
#energy_band = Energy([0.5, 50.0], 'TeV')
offset_band = Angle([0.0, 1.5], 'deg')
backnorm=[]
Ncounts=[]
list_zen=filter(lambda o1: o1.pointing_zen.value<34.3, mylist)
N=len(list_zen)
print N

for i in range N:
    obs=list_zen[i]

    o1=obs
    o1.events = events.select_offset(self.offset_band)


