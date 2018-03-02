import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils


msname='EVN_ST2_NT600_5AM.MS'
ms.open(msname)
# ms.summary()
# Better:
# listobs(msname)

date = '2007-06-11'
times = [('05:05:00', '05:10:00')]

timeqs = utils.mktimeqs(date, times)


station_map = utils.get_station_map(msname)
station_names = sorted(station_map.keys())
stations2 = sorted(map(station_map.get, station_names))


ref_station_name = '1' # Effelsberg
ref_station2 = station_map[ref_station_name]
ref_s_ind2 = stations2.index(ref_station2)
polind = 0

for timeq, swid in itertools.product(timeqs, range(1)):
    timeq2 = ffd.actual_timerangeq(msname, timeq)
    delays, phases, rates, sigs = [], [], [], []
    for pol_id in [0,3]:
        dels, phs, rs, sig = fringer.fit_fringe_lm(msname, stations2, ref_s_ind2,
                                                   swid, polind, pol_id, timeq2,
                                                   solint=300)
        delays.append(dels)
        phases.append(phs)
        rates.append(rs)
        sigs.append(sig)
    print "#\n#{}, IF={}\n#".format(timeq, swid)
    # utils.print_res(station_map, stations2, delays, rates, sigs=sigs)


