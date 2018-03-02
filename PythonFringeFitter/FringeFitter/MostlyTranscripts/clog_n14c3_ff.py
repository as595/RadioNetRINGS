import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils
import sys

# msname='N14C3.MS'

# date = '2014-10-22'
# times = [('13:18:00','13:20:00')]
# timeqs = utils.mktimeqs(date, times)
msname = '/scratch/small/NME_challenge/N14C3_scan2_2.MS'
timeqs = ['True']



ref_station_name = 'Jb' 
polind = 0 # FIXME

solint = 120
threshold = 0.0

swids = utils.get_spectral_windows(msname)
all_delays = []
all_phases = []
all_rates = []

for timeq, swid in itertools.product(timeqs, swids):
    # ref_station = utils.station_by_name(msname, ref_station_name)
    # stations2 = utils.get_stations(msname, timeq)
    ref_station = 1
    # stations2 = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11]
    stations2 = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11]
    delays, phases, rates, sigs = [], [], [], []
    for pol_id in [0,3]:
        dels, phs, rs, sig = fringer.fit_fringe_lm(msname, stations2, ref_station,
                                                   swid, polind, pol_id, timeq,
                                                   solint=solint, datacol="CORRECTED_DATA",
                                                   threshold=threshold, ref_weights=True)
        delays.append(dels)
        phases.append(phs)
        rates.append(rs)
        sigs.append(sig)
    all_delays.append(delays)
    all_phases.append(phases)
    all_rates.append(rates)
    print "#\n#{}, IF={}\n#".format(timeq, swid)
    station_map = utils.get_station_map(msname)
    utils.print_res(station_map, stations2, phases, delays, rates)

all_delays = np.array(all_delays)
all_rates = np.array(all_rates)
all_phases = np.array(all_phases)


sm = utils.get_station_map(msname, 'NAME')
ism = utils.invert_map(sm)
stationt = ",".join([ism[s] for s in stations2])
spwt = ",".join(str(i) for i in swids)
            
l_delay = utils.flattenificate(-1e9*all_delays)
l_phase = utils.flattenificate(-1*180/np.pi*all_phases)

gencal(vis=msname, caltable='huh.G', caltype='sbd', spw=spwt,
       antenna=stationt, pol='R,L', parameter = l_delay)

gencal(vis=msname, caltable='huh_ph.G', caltype='ph', spw=spwt,
       antenna=stationt, pol='R,L', parameter = l_phase)


sys.exit(0)

# ------ Don't do this:
plotms(vis=msname, xaxis='frequency', yaxis='phase',
       scan='2', averagedata=True,
       avgtime='600s', # don't know how to say "long enough"
       iteraxis='baseline', correlation="ll", antenna="EF&NT", ydatacolumn='corrected')

plotms(vis=msname, xaxis='frequency', yaxis='phase',
       scan='2', averagedata=True,
       avgtime='600s', # don't know how to say "long enough"
       iteraxis='baseline', correlation="ll", antenna="EF&NT", ydatacolumn='corrected')

import refonly
refonly.all_cals_to_ref(msname, stations2, swids, polind, pol_id, timeq, solint=300)


plotms(vis=msname, xaxis='frequency', yaxis='phase',
       scan='2', iteraxis='time', correlation="ll", antenna="EF&SV", ydatacolumn='corrected')

plotms(vis=msname, xaxis='time', yaxis='phase',
       averagedata=True, avgchannel='32', 
       correlation="ll", antenna="EF&SV", ydatacolumn='corrected')
