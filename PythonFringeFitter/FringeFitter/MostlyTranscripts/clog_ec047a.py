import casa
import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils

msname = globals()['msname']

if False:
    flagdata(msname, spw='*:0~1;30~31', mode='manual', flagbackup=True) 

ms.open(msname)

timeqs = ["SCAN_NUMBER=2"]
station_map = utils.get_station_map(msname)
stations2 = sorted(ffd.actual_stations(msname, timeqs[0]).keys())
ism = utils.invert_map(station_map)
station_names = [ism[s] for s in stations2]


ref_station_name = 'EF' 
ref_station2 = station_map[ref_station_name]
ref_s_ind2 = stations2.index(ref_station2)
polind = 0
swids = range(8)

shape = (2, len(stations2), len(swids))
delays = np.zeros(shape, np.float, order='F') # Looks like we don't need 'F'.
phases = np.zeros(shape, np.float, order='F')
rates = np.zeros(shape, np.float, order='F')
sigs = []
for timeq, swid in itertools.product(timeqs, swids):
    casalog.post("doing a thing {}".format(swid))
    timeq2 = ffd.actual_timerangeq(msname, timeq)
    for pol_id in [0,1]:
        dels, phs, rs, sig = fringer.fit_fringe_lm(msname, stations2, ref_s_ind2,
                                                   swid, polind, pol_id, timeq2,
                                                   datacol="CORRECTED_DATA",
                                                   solint=300, threshold=1000)
        delays[pol_id, :, swid] = dels
        phases[pol_id, :, swid] = phs
        rates[pol_id, :, swid] = rs
        sigs.append(sig)
    print "#\n#{}, IF={}\n#".format(timeq, swid)
    # utils.print_res(station_map, stations2, delays, rates)


stationt = ",".join(map(str, stations2))
spwt = ",".join(map(str, swids))


darr = -delays*1e9
pharr = -180*phases/np.pi

# casa.gencal(vis=msname, caltable='man_pcal_single.G',
#             caltype='sbd',
#             spw=spwt,
#             antenna=stationt,
#             pol='R,L',
#             parameter = np.ndarray.flatten(darr, order='F'))


# casa.gencal(vis=msname,
#             caltable='man_ph.G',
#             caltype='ph',
#             spw=spwt,
#             antenna=stationt,
#             pol='R,L',
#             parameter = np.ndarray.flatten(pharr, order='F'))



