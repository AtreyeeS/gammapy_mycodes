import matplotlib.pyplot as plt
import shutil
import os
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from gammapy.extern.pathlib import Path
from gammapy.utils.energy import EnergyBounds
from gammapy.utils.nddata import sqrt_space
from gammapy.data import DataStore, ObservationGroupAxis, ObservationGroups
from gammapy.background import EnergyOffsetBackgroundModel
from gammapy.background import OffDataBackgroundMaker

#create a directory
os.mkdir("background")

#observation list
name="PKS 2155-304"
name="Crab"
datastore = DataStore.from_dir("$HESS_DATA") 
src=SkyCoord.from_name(name)
sep=SkyCoord.separation(src,datastore.obs_table.pointing_radec)
srcruns=(datastore.obs_table[sep<2.0*u.deg]) 
obsid=srcruns['OBS_ID'].data
mylist=datastore.obs_list(obsid[:30])

# Define the grouping
zenith_bins=np.linspace(0,90,6)
axes = [ObservationGroupAxis('ZEN_PNT', zenith_bins, fmt='edges')]

# Create the ObservationGroups object
obs_groups = ObservationGroups(axes)

# write it to file
filename = str(scratch_dir / 'group-def.fits')
obs_groups.obs_groups_table.write(filename, overwrite=True)

obs_table_with_group_id = obs_groups.apply(srcruns[0:30])
