# coding: utf-8
from gammapy.data import EventList
filename = '$GAMMAPY_EXTRA/test_datasets/unbundled/hess/run_0023037_hard_eventlist.fits.gz'
events=EventList.read(filename)
events
filename = '$GAMMAPY_EXTRA/datasets/fermi_2fhl/2fhl_events.fits.gz'
events.observation_live_time_duration
print(events.observation_live_time_duration)
print(events.observation_dead_time_fraction)
print(events.observation_time_duration)
events.observation_time_duration * events.observation_dead_time_fraction
1577.0 s - 56.39856698736619 s
1577.0  - 56.39856698736619 
1.0-0.0357632003725
0.9642367996275 * events.observation_time_duration
