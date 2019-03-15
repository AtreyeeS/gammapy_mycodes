import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
from astropy.visualization import simple_norm

from gammapy.data import DataStore
from gammapy.image import SkyImage
from gammapy.utils.energy import Energy
import regions

import scipy.stats as ss
import copy

from scaling_compare import *

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

exclusion_mask_tevcat = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
exclusion_mask_tevcat = exclusion_mask_tevcat.reproject(reference=ref_image)

energy_band = Energy([0.5, 1.0], 'TeV')
#energy_band = Energy([0.5, 50.0], 'TeV')
offset_band = Angle([0.0, 1.5], 'deg')
backnorm=[]
Ncounts=[]
list_zen=filter(lambda o1: o1.pointing_zen.value<34.3, mylist)
N=len(list_zen)
#print N

N=1

for i in range(N):
#    if i%10 ==0:
#        print "----running observation -- ",i
    
    obs=list_zen[i]

    events = datastore.obs(obs.obs_id).events
    counts_image = SkyImage.empty_like(ref_image)
    counts_image.fill_events(events)

    exclusion_region = regions.CircleSkyRegion(src,0.3*u.deg)
    mask = counts_image.region_mask(exclusion_region)

    mask1=copy.copy(mask)
    mask1.data=np.invert(mask.data)

    image_maker = SingleObsImageMaker(
        obs=obs,
        empty_image=ref_image,
        energy_band=energy_band,
        offset_band=offset_band,
        exclusion_mask=mask1,
    )

    image_maker.counts_image()
    counts_image = image_maker.images['counts']
    #norm = simple_norm(counts_image.data, stretch='sqrt', min_cut=0, max_cut=0.9)
    #counts_image.smooth(radius=0.08 * u.deg).plot(norm=norm, add_cbar=True)
    #plt.show()
    

    image_maker.bkg_image()
    background_image = image_maker.images['bkg']
    #norm = simple_norm(background_image.data, stretch='sqrt', min_cut=0, max_cut=0.2)
    #background_image.plot(norm=norm, add_cbar=True)
    #plt.show()

    backnorm.append(image_maker.table_bkg_scale['bkg_scale'].data[0])
    Ncounts.append(image_maker.table_bkg_scale['N_counts'].data[0])
   
 
print "exited from loop...now plot..." 

print Ncounts, backnorm

"""
plt.hist(backnorm)
plt.title("Energy band: "+str(energy_band))
plt.ylabel("frequency of occurance")
plt.xlabel("Backscale factor")
plt.show()

#plt.plot(Ncounts,backnorm,"b.")
#plt.xlabel("Background counts")
#plt.show()

time=[o1.observation_live_time_duration.value for o1 in list_zen]    
#plt.plot(time[:-1],backnorm,"r+")
#plt.xlabel("Observation Duration")
#plt.show()

zen_ang=[o1.pointing_zen.value for o1 in list_zen]    
#plt.plot(zen_ang[:-1],backnorm,"r+")
#plt.xlabel("Zenith angle")
#plt.show()

muonef=[o1.muoneff for o1 in list_zen]    
#plt.plot(muonef[:-1],backnorm,"r+")
#plt.xlabel("Muon effectivity")
#plt.show()


bkg_counts = np.divide(Ncounts,backnorm)

plt.plot(zen_ang[:-1],np.divide(Ncounts,time[:-1]),"r.",label="Observed")
plt.plot(zen_ang[:-1],np.divide(bkg_counts,time[:-1]),"b.",label="Model predicted")
plt.xlabel("Zenith angle")
plt.ylabel("Count rate [/s]")
plt.title("Energy band: "+str(energy_band))
plt.legend()
plt.show()

plt.plot(muonef[:-1],np.divide(Ncounts,time[:-1]),"r.",label="Observed")
plt.plot(muonef[:-1],np.divide(bkg_counts,time[:-1]),"b.",label="Model predicted")
plt.xlabel("Muon effectivity")
plt.ylabel("Count rate [/s]")
plt.title("Energy band: "+str(energy_band))
plt.legend()
plt.show()

dev = np.absolute(np.subtract(1.0,backnorm))
"""
