import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils


msname='N08C1_4C39.25.MS'
ms.open(msname)
ms.summary()
# Better:
listobs(msname)

date = '2008-03-10'
times = [('17:06:01.0', '17:09:00.0')]
mysolint = 300 # seconds
myrefant = '1' # string for 1-based refant number, 1 = EF


timeqs = ffd.mktimeqs(date, times)


station_map = utils.get_station_map(msname)
station_names = sorted(station_map.keys())
# stations2 = sorted(map(station_map.get, station_names))
# station 8 missing; 3 misses a bunch of stuff.
stations2 = [0, 1, 2, 3, 4, 5, 6, 7, 9, 10]


ref_station_name = myrefant # '1' is Effelsberg
ref_station2 = station_map[ref_station_name]
ref_s_ind2 = stations2.index(ref_station2)
swid, polind, = (0, 0)

#stations2 = [0,1,2,3,4,5,6,7,9,10]
stations2 = [0,1,2,3,4,5,6,7,9,10]


for timeq, swid in itertools.product(timeqs, range(4)):
    timeq2 = ffd.actual_timerangeq(msname, timeq)
    delays, phases, rates, sigs = [], [], [], []
    for pol_id in [0,3]:
        print "# Doing a thing", swid, pol_id
        dels, phs, rs, sig = fringer.fit_fringe_lm(msname, stations2, ref_s_ind2,
                                                   swid, polind, pol_id, timeq2, solint=mysolint, 
                                                   threshold=1000)
        delays.append(dels)
        phases.append(phs)
        rates.append(rs)
        sigs.append(sig)
    print "#\n#{}\n# IF {}\n#".format(timeq, swid+1)
    utils.print_res(station_map, stations2,phases, delays, rates)


