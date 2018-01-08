#simulating a simple power law from the irfs... following the notebook cta_simulation.ipynb

import matplotlib.pyplot as plt
from gammapy.catalog import SourceCatalog3FHL
redshift = source.data['Redshift']
src_spectral_model = source.spectral_model

import astropy.units as u
src_spectral_model.plot(energy_range=[10 * u.GeV, 2 *u.TeV])


from gammapy.scripts import CTAPerf

from gammapy.scripts import CTAPerf
# South array optimisation for faint source 
filename = '$GAMMAPY_EXTRA/datasets/cta/perf_prod2/point_like_non_smoothed/South_50h.fits.gz'
cta_perf = CTAPerf.read(filename)
cta_perf.peek()

prod_dir = '$GAMMAPY_EXTRA/datasets/cta/perf_prod2/point_like_non_smoothed/'
opti = ['0.5h', '5h', '50h']
site = ['South', 'North']
cta_perf_list = []  # will be filled with different performance
labels = []  # will be filled with different performance labels for the legend
for isite in site: 
    for iopti in opti:
        filename = prod_dir + '/' + isite + '_' + iopti + '.fits.gz'
        cta_perf = CTAPerf.read(filename)
        cta_perf_list.append(cta_perf)
        labels.append(isite + ' (' + iopti + ')')

CTAPerf.superpose_perf(cta_perf_list, labels)

spectral_model=src_spectral_model


from gammapy.scripts.cta_utils import Target
target = Target(
    name=source.data['Source_Name'],  # from the 3FGL catalogue source class
    model=source.spectral_model,  # defined above
)

from gammapy.scripts.cta_utils import ObservationParameters
alpha = 0.2 * u.Unit('')  # normalisation between ON and OFF regions
livetime = 5. * u.h
# energy range used for statistics (excess, significance, etc.)
emin, emax = 0.05 * u.TeV, 5 * u.TeV
params = ObservationParameters(
    alpha=alpha, livetime=livetime,
    emin=emin, emax=emax,
)

from gammapy.scripts import CTAPerf
# PKS 2155-304 is 10 % of Crab at 1 TeV ==> intermediate source
filename = '$GAMMAPY_EXTRA/datasets/cta/perf_prod2/point_like_non_smoothed/South_5h.fits.gz'
perf = CTAPerf.read(filename)

from gammapy.scripts.cta_utils import CTAObservationSimulation

simu = CTAObservationSimulation.simulate_obs(
    perf=perf, target=target, obs_param=params)

# print simulation results
print(simu)

CTAObservationSimulation.plot_simu(simu, target)

stats = simu.total_stats_safe_range
stats_dict = stats.to_dict()
print('excess: {}'.format(stats_dict['excess']))
print('sigma: {:.1f}'.format(stats_dict['sigma']))


