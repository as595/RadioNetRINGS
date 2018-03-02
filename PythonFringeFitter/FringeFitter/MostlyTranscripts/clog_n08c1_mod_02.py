import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils


msname='N08C1_mod_fringready.MS'
ms.open(msname)
# ms.summary()
# Better:
# listobs(msname)

date = '2008-03-10'
times = [('17:06:01.0', '17:09:00.0'),
         ('17:09:00.0', '17:11:00.0')]

timeqs = utils.mktimeqs(date, times)


station_map = utils.get_station_map(msname)
station_names = sorted(station_map.keys())
stations2 = sorted(map(station_map.get, station_names))
stations2 = [0, 1, 2, 3, 4, 5, 6, 7, 9]




#ref_station_name = '1' # Effelsberg
ref_station_name = '3' 
ref_station2 = station_map[ref_station_name]
ref_s_ind2 = stations2.index(ref_station2)

polind = 0

for timeq, swid in itertools.product(timeqs[:1], range(4)[:1]):
    timeq2 = ffd.actual_timerangeq(msname, timeq)
    delays, phases, rates, sigs = [], [], [], []
    for pol_id in [0,3]:
        dels, phs, rs, sig = fringer.fit_fringe_lm(msname, stations2, ref_station2,
                                                   swid, polind, pol_id, timeq2,
                                                   solint=300)
        delays.append(dels)
        phases.append(phs)
        rates.append(rs)
        sigs.append(sig)
    print "#\n#{}, IF={}\n#".format(timeq, swid)
    utils.print_res(station_map, stations2, delays, rates)


# for pol r, presumably.
aips_phases = [-1.52, -1.65, 0.0, -0.31, -2.78, 1.71, -2.05, 0.85, 0.24, 0.99]
