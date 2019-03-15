import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
from astropy.visualization import simple_norm

from gammapy.data import DataStore
from gammapy.image import SkyImage
from gammapy.scripts import SingleObsImageMaker
from gammapy.utils.energy import Energy
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
#N=len(mylist)

ref_image = SkyImage.empty(
    nxpix=400, nypix=400, binsz=0.02,
    xref=src.ra.deg, yref=src.dec.deg,
    coordsys='CEL', proj='TAN',
)

exclusion_mask_tevcat = SkyImage.read('$GAMMAPY_EXTRA/datasets/exclusion_masks/tevcat_exclusion.fits')
exclusion_mask_tevcat = exclusion_mask_tevcat.reproject(reference=ref_image)

energy_band = Energy([0.52480746, 1.58489319], 'TeV')
energy_band = Energy([0.5, 1.5], 'TeV')
of1=Angle([0.0, 0.50251891], 'deg')
of2=Angle([0.50251891, 1.20499801], 'deg')
of3=Angle([1.20499801, 2.01007563], 'deg')
of4=Angle([2.01007563, 2.48734169], 'deg')
#of1=Angle([0.0, 1.2], 'deg')
#of2=Angle([0.50251891, 1.2], 'deg')
offset_band = [of1,of2,of3,of4]
backnorm=[]
Ncounts=[]



list_zen=filter(lambda o1: o1.pointing_zen.value<34.3, mylist)

time=[o1.observation_live_time_duration.value for o1 in list_zen]    
zen_ang=[o1.pointing_zen.value for o1 in list_zen]    
muonef=[o1.muoneff for o1 in list_zen]    


N=len(list_zen)
print N

for i in range(N):
    if i%10 ==0:
        print "----running observation -- ",i
    
    obs=list_zen[i]

    events = datastore.obs(obs.obs_id).events
    counts_image = SkyImage.empty_like(ref_image)
    counts_image.fill_events(events)

    exclusion_region = regions.CircleSkyRegion(src,0.3*u.deg)
    mask = counts_image.region_mask(exclusion_region)

    mask1=copy.copy(mask)
    mask1.data=np.invert(mask.data)
    norm_off=[]
    Nc=[]
    for j in range(len(offset_band)):

        image_maker = SingleObsImageMaker(
            obs=obs,
            empty_image=ref_image,
            energy_band=energy_band,
            offset_band=offset_band[j],
           exclusion_mask=mask1,
        )

        image_maker.counts_image()
        counts_image = image_maker.images['counts']
    

        image_maker.bkg_image()
        background_image = image_maker.images['bkg']
        #norm = simple_norm(background_image.data, stretch='sqrt', min_cut=0, max_cut=0.2)
        #background_image.plot(norm=norm, add_cbar=True)
        #plt.show()

        norm_off.append(image_maker.table_bkg_scale['bkg_scale'].data[0])
        Nc.append(image_maker.table_bkg_scale['N_counts'].data[0])

    backnorm.append(norm_off)
    Ncounts.append(Nc)
   
   
print "exited from loop...now plot..." 

#plt.hist(backnorm)
#plt.title("Energy band: "+str(energy_band))
#plt.ylabel("frequency of occurance")
#plt.xlabel("Backscale factor")
#plt.show()

#plt.plot(Ncounts,backnorm,"b.")
#plt.xlabel("Background counts")
#plt.show()

#time=[o1.observation_live_time_duration.value for o1 in mylist]    
#plt.plot(time[:-1],backnorm,"r+")
#plt.xlabel("Observation Duration")
#plt.show()

#zen_ang=[o1.pointing_zen.value for o1 in mylist]    
#plt.plot(zen_ang[:-1],backnorm,"r+")
#plt.xlabel("Zenith angle")
#plt.show()

#muonef=[o1.muoneff for o1 in mylist]    
#plt.plot(muonef[:-1],backnorm,"r+")
#plt.xlabel("Muon effectivity")
#plt.show()



bkg_counts = np.divide(Ncounts,backnorm)

#plt.plot(zen_ang,np.divide(Ncounts,time),"r.",label="Observed")
#plt.plot(zen_ang,np.divide(bkg_counts,time),"b.",label="Model predicted")
#plt.xlabel("Zenith angle")
#plt.ylabel("Count rate [/s]")
#plt.title("Energy band: "+str(energy_band))
#plt.legend()
#plt.show()

#plt.plot(muonef,np.divide(Ncounts,time),"r.",label="Observed")
#plt.plot(muonef,np.divide(bkg_counts,time),"b.",label="Model predicted")
#plt.xlabel("Muon effectivity")
#plt.ylabel("Count rate [/s]")
#plt.title("Energy band: "+str(energy_band))
#plt.legend()
#plt.show()

#dev = np.absolute(np.subtract(1.0,backnorm))


Nc0=[row[0] for row in Ncounts]
Nc1=[row[1] for row in Ncounts]
Nc2=[row[2] for row in Ncounts]
Nc3=[row[3] for row in Ncounts]

b0=[row[0] for row in bkg_counts]
b1=[row[1] for row in bkg_counts]
b2=[row[2] for row in bkg_counts]
b3=[row[3] for row in bkg_counts]

bn0=[row[0] for row in backnorm]
bn1=[row[1] for row in backnorm]
bn2=[row[2] for row in backnorm]
bn3=[row[3] for row in backnorm]


plt.hist(bn0,alpha=0.5,label=str(of1))
plt.hist(bn1,alpha=0.5,label=str(of2))
plt.hist(bn2,alpha=0.5,label=str(of3))
plt.hist(bn3,alpha=0.5,label=str(of4))
plt.legend()
plt.show()


zenc=10

bn0_zc=[]
bn1_zc=[]
bn2_zc=[]

for b,o1 in zip(bn0,list_zen):
    if o1.pointing_zen.value<zenc:
        bn0_zc.append(b)

for b,o1 in zip(bn1,list_zen):
    if o1.pointing_zen.value<zenc:
        bn1_zc.append(b)


for b,o1 in zip(bn2,list_zen):
    if o1.pointing_zen.value<zenc:
        bn2_zc.append(b)



plt.hist(bn0_zc,alpha=0.5,label=str(of1))
plt.hist(bn1_zc,alpha=0.5,label=str(of2))
plt.hist(bn2_zc,alpha=0.5,label=str(of3))
plt.legend()
plt.show()


N0_zc=[]
N1_zc=[]
N2_zc=[]


for b,o1 in zip(Nc0,list_zen):
    if o1.pointing_zen.value<zenc:
        N0_zc.append(b)

for b,o1 in zip(Nc1,list_zen):
    if o1.pointing_zen.value<zenc:
        N1_zc.append(b)


for b,o1 in zip(Nc2,list_zen):
    if o1.pointing_zen.value<zenc:
        N2_zc.append(b)


b0_zc=[]
b1_zc=[]
b2_zc=[]
b3_zc=[]


for b,o1 in zip(b0,list_zen):
    if o1.pointing_zen.value<zenc:
        b0_zc.append(b)

for b,o1 in zip(b1,list_zen):
    if o1.pointing_zen.value<zenc:
        b1_zc.append(b)


for b,o1 in zip(b2,list_zen):
    if o1.pointing_zen.value<zenc:
        b2_zc.append(b)


list_zen1=filter(lambda o1: o1.pointing_zen.value<zenc, mylist)
mul=[o1.muoneff for o1 in list_zen1]

list_zen1=filter(lambda o1: o1.pointing_zen.value<zenc, mylist)
time1=[o1.observation_live_time_duration.value for o1 in list_zen1]

zen=[o1.pointing_zen.value for o1 in list_zen1]

plt.plot(zen_ang,np.divide(Nc0,time),"r.",label=of1*u.deg)
#plt.plot(zen_ang,np.divide(b0,time),"r+")
plt.plot(zen_ang,np.divide(Nc1,time),"g.",label=of2*u.deg)
#plt.plot(zen_ang,np.divide(b1,time),"g+")
plt.plot(zen_ang,np.divide(Nc2,time),"b.",label=(of3*u.deg))
#plt.plot(zen_ang,np.divide(b2,time),"b+")
plt.plot(zen_ang,np.divide(Nc3,time),"c.",label=(of4*u.deg))
#plt.plot(zen_ang,np.divide(b3,time),"c+")
plt.legend()
plt.show()


plt.plot(zen,np.divide(N0_zc,time1),"r.")

mu1=0.773
mu2=0.774


list_mu1=filter(lambda o1: o1.muoneff<mu2, list_zen1)
list_mu1=filter(lambda o1: o1.muoneff>mu1, list_mu1)

Nm0=[]

for b,o1 in zip(N0_zc,list_zen1):
    if (o1.muoneff>mu1 and o1.muoneff<mu2):
        Nm0.append(b)

tm=[o1.observation_live_time_duration.value for o1 in list_mu1]
zm=[o1.pointing_zen.value for o1 in list_mu1]

plt.plot(zen,np.divide(N0_zc,time1),"r.",label=of1*u.deg)
plt.plot(zm,np.divide(Nm0,tm),"b+")
plt.show()

