#from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from numpy.testing import assert_allclose
import gammapy
import numpy as np
import astropy
import regions
import sherpa

from gammapy.utils.testing import requires_dependency
from gammapy.irf import EnergyDispersion, EffectiveAreaTable, EffectiveAreaTable2D
from gammapy.spectrum import SpectrumExtraction, SpectrumSimulation, SpectrumFit, SpectrumObservationList, SpectrumObservationStacker
from gammapy.spectrum.models import PowerLaw, SpectralModel
from gammapy.utils.modeling import Parameter, ParameterList
from gammapy.utils.energy import EnergyBounds

from astropy.coordinates import SkyCoord, Angle
from astropy.table import vstack as vstack_table
from regions import CircleSkyRegion
from gammapy.data import DataStore, ObservationList


#factor: between 0.5 and 1 -> shift the acceptance area (an example)

plt.ion()

 #======================================
 
#def Test_stacker(alpha1,alpha2,time1,time2,factor):

alpha1=0.2
alpha2=0.2
time1=1
time2=1
factor=1
    
e_true = np.logspace(-2, 2.5, 109) * u.TeV
e_reco = np.logspace(-2,2, 79) * u.TeV

edisp1 = EnergyDispersion.from_gauss(e_true=e_true, e_reco=e_reco, sigma=0.2, bias=0,)
edisp2 = EnergyDispersion.from_gauss(e_true=e_true, e_reco=e_reco, sigma=0.2, bias=0,)

#Aeff
ee = EnergyBounds(np.logspace(-2,2.5,109)*u.TeV)
p1=6.85e9
p2=0.0891
p3=5.0e5

f = lambda x : p1 * (x/u.MeV)**(-p2) * np.exp(((-p3)*u.MeV)/x)
value=f(ee.log_centers.to('MeV'))
data = value * u.cm ** 2
aeff1=EffectiveAreaTable(ee.lower_bounds,ee.upper_bounds,data)

f2 = lambda x : p1 * (x/u.MeV)**(-p2) * np.exp(((-p3*factor)*u.MeV)/x)
value2=f2(ee.log_centers.to('MeV'))
data2 = value2 * u.cm ** 2
aeff2=EffectiveAreaTable(ee.lower_bounds,ee.upper_bounds,data2)

#======================================

#DEFINE MODEL
   
#Model for the source 
index = 2.1 * u.Unit('')
amplitude = 3.5e-12 * u.Unit('cm-2 s-1 TeV-1')
reference = 1 * u.TeV
pwl = PowerLaw(index=index, amplitude=amplitude, reference=reference)
print(pwl)

# Bkg :
bkg_index = 3 * u.Unit('')
bkg_amplitude = 4.31e-13 * u.Unit('cm-2 s-1 TeV-1')
reference = 1 * u.TeV
bkg_model = PowerLaw(index=bkg_index, amplitude=bkg_amplitude, reference=reference)

 
   
#======================================

#SIMULATE SPECTRA

livetime1 = time1 * u.h
livetime2= time2 * u.h

n_obs = 100
seeds = np.arange(n_obs)
sim1 = SpectrumSimulation(aeff=aeff1,
                          edisp=edisp1,
                          source_model=pwl,
                          livetime=livetime1,
                          background_model=bkg_model,
                          alpha=alpha1)
    
sim1.run(seeds)
#print(sim1.result)
    
#sim2 = SpectrumSimulation(aeff=aeff2,
#                          edisp=edisp2,
#                          source_model=pwl,
#                          livetime=livetime2,
#                          background_model=bkg_model,
#                          alpha=alpha2)
    
#sim2.run(seeds)
#print(sim2.result)

sim2=sim1

Indiv_best_fit_index = []
best_fit_index = []

i=0
 
while i<n_obs: 
    sim_result=(sim1.result[i],sim2.result[i])
    #FIT SPECTRA:
    pwl.parameters['index'].parmax = 10
        
    for obs in sim_result:
        fit = SpectrumFit(obs, pwl.copy(), stat='wstat')
        fit.model.parameters['index'].value = 2
        fit.fit()
        fit.est_errors()
        Indiv_best_fit_index.append(fit.result[0].model.parameters['index'].value)      
        #print(fit.result[0])
        #print(' Indiv_best_fit_index = ', Indiv_best_fit_index)
        
    #i+=1    
    #STACK SPECTRA 2 by 2
    # Add the two spectra
    obs_stacker = SpectrumObservationStacker(sim_result)
    #print('sim_result is the list = ',sim_result)   
    obs_stacker.run()
    
    fit_stacked = SpectrumFit(obs_stacker.stacked_obs, pwl.copy(), stat='wstat')
    fit_stacked.model.parameters['index'].value = 2
    fit_stacked.fit()
    fit_stacked.est_errors()
    best_fit_index.append(fit_stacked.result[0].model.parameters['index'].value)
    #print('  best_fit_index =',best_fit_index)
    print(i)
    i+=1


#ax0, ax1  =  Indiv_best_fit_index.plot(figsize=(8,8))
#ax0.set_ylim(0, 20)

#ax0, ax1  =  best_fit_index.plot(figsize=(8,8))
#ax0.set_ylim(0, 20)

fix, axes = plt.subplots(1,2, figsize=(12, 4))
axes[0].hist(Indiv_best_fit_index)
axes[0].set_xlabel('Index for Indiv spectra')
axes[1].hist(best_fit_index)
axes[1].set_xlabel('Index for stacked spectra')
print("Indiv = ",np.mean(Indiv_best_fit_index),np.std(Indiv_best_fit_index))
print("Stacked = ",np.mean(best_fit_index),np.std(best_fit_index))

