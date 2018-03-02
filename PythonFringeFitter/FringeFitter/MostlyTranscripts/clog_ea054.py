import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils

msname = "ea054_continuum.MS"

if False:
    flagdata(msname, spw='*:0~1;30~31', mode='manual', flagbackup=True) 

ms.open(msname)
# Better:
#listobs(msname)

date = '2013-10-30'
times = [('23:11:10', '23:21:10')]

timeqs = utils.mktimeqs(date, times)

timeqs = ["SCAN_NUMBER=43"]
station_map = utils.get_station_map(msname)
# station_names = sorted(station_map.keys())
#['BD', 'EF', 'JB', 'MC', 'NT', 'ON', 'SV', 'TR', 'UR', 'WB', 'ZC']

# station_names = ['BD', 'EF', 'JB', 'MC', 'ON', 'SV', 'TR', 'UR', 'WB', 'ZC']
stations2 = sorted(utils.actual_stations(msname, timeqs[0]).keys())
ism = utils.invert_map(station_map)
station_names = [ism[s] for s in stations2]


ref_station_name = 'EF' 
ref_station2 = station_map[ref_station_name]
ref_s_ind2 = stations2.index(ref_station2)
polind = 0
swids = range(8)

delays, phases, rates, sigs = [], [], [], []
for timeq, swid in itertools.product(timeqs, swids):
    print "doing a thing", swid
    timeq2 = utils.actual_timerangeq(msname, timeq)
    for pol_id in [0,1]:
        dels, phs, rs, sig = fringer.fit_fringe_lm(msname, stations2, ref_s_ind2,
                                                   swid, polind, pol_id, timeq2,
                                                   solint=300, threshold=1000)
        delays.append(dels)
        phases.append(phs)
        rates.append(rs)
        sigs.append(sig)
    print "#\n#{}, IF={}\n#".format(timeq, swid)
    # utils.print_res(station_map, stations2, delays, rates)


stationt = ",".join(map(str, stations2))
spwt = ",".join(map(str, swids))

darr = np.array(delays)
darr_r = -1e9*darr[0::2, :]
darr_l = -1e9*darr[1::2, :]

plist_r = np.ndarray.flatten(darr_r)
plist_l = np.ndarray.flatten(darr_l)

pharr = np.array(phases)
phlist_r = np.ndarray.flatten(-180*pharr[0::2, :]/np.pi)
phlist_l = np.ndarray.flatten(-180*pharr[1::2, :]/np.pi)


gencal(vis=msname, caltable='man_pcal_r.G',
       caltype='sbd',
       spw=spwt,
       antenna=stationt,
       pol='R',
       parameter = plist_r)

gencal(vis=msname, caltable='man_pcal_l.G',
       caltype='sbd',
       spw=spwt,
       antenna=stationt,
       pol='L',
       parameter = plist_l)

gencal(vis=msname,
       caltable='man_ph_r.G',
       caltype='ph',
       spw=spwt,
       antenna=stationt,
       pol='R',
       parameter = phlist_r)

gencal(vis=msname,
       caltable='man_ph_l.G',
       caltype='ph',
       spw=spwt,
       antenna=stationt,
       pol='L',
       parameter = phlist_l)





