# start by simulating crab spectrum!

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from gammapy.irf import EnergyDispersion, EnergyDispersion2D, EffectiveAreaTable, EffectiveAreaTable2D
from gammapy.irf import EnergyDependentMultiGaussPSF 
from gammapy.spectrum import SpectrumSimulation, SpectrumFit
from gammapy.spectrum.models import PowerLaw

import os


#now reading IRFs from one particular observation
DIR = "/Users/asinha/HESS_data/ash_stereo/run054200-054399/run054213/"
files=os.listdir(DIR)
for afile in files:
    if afile[0:4]=="aeff":
        f_ar=str(DIR+afile)
    if afile[0:5]=="edisp":
        f_edp=str(DIR+afile)
    if afile[0:3]=="psf":
        f_psf=str(DIR+afile)

aeff = EffectiveAreaTable2D.read(f_ar, hdu='EFFECTIVE AREA')
edisp = EnergyDispersion2D.read(f_edp)
#psf = EnergyDependentMultiGaussPSF.read(f_psf, hdu='POINT SPREAD FUNCTION')

# Define obs parameters
livetime = 1.0 * u.hr
offset = 0.01 * u.deg
lo_threshold = 0.1 * u.TeV
hi_threshold = 60 * u.TeV

# Define spectral model
index = 2.0 * u.Unit('')
amplitude = 2.5 * 1e-10 * u.Unit('cm-2 s-1 TeV-1')
reference = 1 * u.TeV
model = PowerLaw(index=index, amplitude=amplitude, reference=reference)



edisp=edisp.to_energy_dispersion(offset=offset)
aeff = aeff.to_effective_area_table(offset=offset)

#fig, axes = plt.subplots(1, 2, figsize=(12, 6))
#edisp.plot_matrix(ax=axes[0])
#aeff.plot(ax=axes[1])
#plt.show()

aeff.lo_threshold = lo_threshold
aeff.hi_threshold = hi_threshold
sim = SpectrumSimulation(aeff=aeff, edisp=edisp, source_model=model, livetime=livetime)
sim.simulate_obs(seed=42, obs_id=0)

sim.obs.peek()
print sim.obs
plt.show()

fit = SpectrumFit(obs_list=sim.obs, model=model, stat='cash')
fit.run()
result = fit.result[0]
print result


energy_range = [0.1, 100] * u.TeV
model.plot(energy_range=energy_range, energy_power=2)
result.model.plot(energy_range=energy_range, energy_power=2)
result.model.plot_error(energy_range=energy_range, energy_power=2)

plt.show()

#now importing a background

#from a background model
bkg_index = 2.5 * u.Unit('')
bkg_amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1')
reference = 1 * u.TeV

bkg_model = PowerLaw(index=bkg_index, amplitude=bkg_amplitude, reference=reference)
alpha = 0.2

n_obs = 1
seeds = np.arange(n_obs)

sim = SpectrumSimulation(aeff=aeff,
                         edisp=edisp,
                         source_model=model,
                         livetime=livetime,
                         background_model=bkg_model,
                         alpha=alpha)

sim.run(seeds)
print(sim.result)
print(sim.result[0])

n_on = [obs.total_stats.n_on for obs in sim.result]
n_off = [obs.total_stats.n_off for obs in sim.result]
excess = [obs.total_stats.excess for obs in sim.result]
fix, axes = plt.subplots(1,3, figsize=(12, 4))
axes[0].hist(n_on)
axes[0].set_xlabel('n_on')
axes[1].hist(n_off)
axes[1].set_xlabel('n_off')
axes[2].hist(excess)
axes[2].set_xlabel('excess')


best_fit_index = []

for obs in sim.result:
    fit = SpectrumFit(obs, model=model, stat='wstat')
    fit.model.parameters['index'].value = 2
    fit.run()
    best_fit_index.append(fit.result[0].model.parameters['index'].value)

plt.hist(best_fit_index)
print('best_fit_index:', best_fit_index)
plt.show()
