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
from gammapy.catalog import SourceCatalogGammaCat
#create a directory

outdir="background_files1"
os.mkdir(outdir)

#observation list
datastore = DataStore.from_dir("$HESS_DATA") 

#runs above 5 deg lat

radec=datastore.obs_table.pointing_radec
lat=radec.galactic.b.value
agnrun=datastore.obs_table[np.abs(lat)>5.0]
obsid=agnrun['OBS_ID'].data

#remove LMC and SN1006

SN1006=SkyCoord.from_name("SN 1006")
sep=SkyCoord.separation(SN1006,radec)
SNrun=datastore.obs_table[sep<2.0*u.deg]
SNid=SNrun['OBS_ID'].data

LMC=SkyCoord.from_name("LMC")
sep=SkyCoord.separation(LMC,radec)
LMCrun=datastore.obs_table[sep<2.0*u.deg]
LMCid=LMCrun['OBS_ID'].data

#get required runs

agnid=list(set(obsid)-set(LMCid)-set(SNid))
mylist=datastore.obs_list(agnid)
zen_ang=[o1.pointing_zen.value for o1 in mylist] 


# Define the grouping
nbins=10
zenith_bins=[0,10,20,30,40,50,90]
#zenith_bins=[min(zen_ang), 10.0, 20.0, 30.0, 45.0, 60.0, max(zen_ang)]
zenith_bins=zenith_bins*u.deg
axes = [ObservationGroupAxis('ZEN_PNT', zenith_bins, fmt='edges')]

# Create the ObservationGroups object
obs_groups = ObservationGroups(axes)

# write it to file
filename = str(outdir + "/group-def.fits")
obs_groups.obs_groups_table.write(filename, overwrite=True)

obs_table_with_group_id = obs_groups.apply(datastore.obs_table.select_obs_id(agnid))
#gammacat exclusion mask

fil_gammacat="/Users/asinha/Gammapy-dev/gammapy-extra/datasets/catalogs/gammacat/gammacat.fits.gz"
cat = SourceCatalogGammaCat(filename=fil_gammacat)
exclusion_table = cat.table.copy()
exclusion_table.rename_column('ra', 'RA')
exclusion_table.rename_column('dec', 'DEC')
radius = exclusion_table['morph_sigma'].data
radius[np.isnan(radius)] = 0.3
exclusion_table['Radius'] = radius * u.deg
exclusion_table = Table(exclusion_table)

#now run the bgmaker
bgmaker = OffDataBackgroundMaker(
    data_store=datastore,
    outdir=outdir,
    run_list=None,
    obs_table=obs_table_with_group_id,
    ntot_group=obs_groups.n_groups,
    excluded_sources=exclusion_table,
)

# Define the energy and offset binning to use
ebounds = EnergyBounds.equal_log_spacing(0.1, 100, 15, 'TeV')
#offset = sqrt_space(start=0, stop=2.5, num=20) * u.deg
offset=np.linspace(0,2.5,20) * u.deg

# Make the model (i.e. stack counts and livetime)
bgmaker.make_model("2D", ebounds=ebounds, offset=offset)

# Smooth the model
bgmaker.smooth_models("2D")
# Write the model to disk
bgmaker.save_models("2D")
bgmaker.save_models(modeltype="2D", smooth=True)

#now copy the background files as bkg into the source runs
data_dir="data_new"
shutil.move(outdir, data_dir) 
datastore= DataStore.from_dir("$HESS_DATA") 
datastore.copy_obs(datastore.obs_table,data_dir) 

group_filename = data_dir + '/background/group-def.fits'  
data_store = DataStore.from_dir(data_dir)                                    
hdu_index_table = bgmaker.make_total_index_table(                            
    data_store=data_store,
    modeltype='2D',
    out_dir_background_model=outdir,
    filename_obs_group_table=str(group_filename),
    smooth=False,
)
filename = data_dir + '/hdu-index.fits.gz'
hdu_index_table.write(str(filename), overwrite=True) 

