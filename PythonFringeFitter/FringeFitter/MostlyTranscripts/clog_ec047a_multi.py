import casa
from taskinit import casalog
import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils

# NB: this assumes manual phase cal has been done first.

msname = globals()['msname']

if False:
    flagdata(msname, spw='*:0~1;30~31', mode='manual', flagbackup=True) 

ms.open(msname)

timeqs = ["SCAN_NUMBER=44"]
station_map = utils.get_station_map(msname)
stations2 = sorted(ffd.actual_stations(msname, timeqs[0]).keys())
ism = utils.invert_map(station_map)
station_names = [ism[s] for s in stations2]


ref_station_name = 'EF' 
ref_station2 = station_map[ref_station_name]
ref_s_ind2 = stations2.index(ref_station2)
polind = 0
swids = range(8)


shape = (2, len(stations2))
delays = np.zeros(shape, np.float, order='F') # Looks like we don't need 'F'.
phases = np.zeros(shape, np.float, order='F')
rates = np.zeros(shape, np.float, order='F')
sigs = []
for timeq in timeqs:
    timeq2 = ffd.actual_timerangeq(msname, timeq)
    for pol_id in [0,1]:
        casalog.post("Getting data")
        anffd = ffd.FFData.make_FFD_multiband(msname, stations2, polind, pol_id, timeq2,
                                            datacol="CORRECTED_DATA", solint=500)
        casalog.post("Fitting fringes")
        dels, phs, rs, sig = fringer.fit_fringe_ffd(anffd, ref_station2, stations2)
        delays[pol_id, :] = dels
        phases[pol_id, :] = phs
        rates[pol_id, :] = rs
        sigs.append(sig)
        # utils.print_res(station_map, stations2, delays, rates)


stationt = ",".join(map(str, stations2))
spwt = ",".join(map(str, swids))

darr = -delays*1e9
pharr = -180*phases/np.pi

phoff = 360* ((anffd.freqs[0] * delays) % 1.0)

# We need the 'F' in the flattening.
# casa.gencal(vis=msname, caltable='multi_delay.G',
#             caltype='mbd', pol='R,L',
#             antenna=stationt,
#             parameter = np.ndarray.flatten(darr, order='F'))

# casa.gencal(vis=msname,
#             caltable='multi_phase.G',
#             caltype='ph', pol='R,L',
#             antenna=stationt,
#             parameter = np.ndarray.flatten(pharr, order='F'))


# casa.gencal(vis=msname,
#             caltable='multi_offset.G',
#             caltype='ph', pol='R,L',
#             antenna=stationt,
#             parameter = np.ndarray.flatten(phoff, order='F'))


# apply ALL of the calibrations:
# casa.applycal(vis=msname, gaintable=["ec047a.gc", "ec047a.tsys"] + glob.glob("man*.G") + glob.glob("multi_[lr].G") + glob.glob("multi_ph*G") + glob.glob("offset*.G"), scan='44', parang=True)

# To see results:
# casa.plotms(vis=msname, xaxis="frequency", yaxis="phase", spw="0~7", antenna="EF", averagedata=T, avgtime="5000", correlation="RR", ydatacolumn="corrected", scan="44")
