import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils


msname='N14C3.MS'
msname = '/scratch/small/N14C3/NME_challenge/N14C3_scan2_2.MS'

timeqs = ['SCAN_NUMBER=2']

station_map = utils.get_station_map(msname)
station_names =  [t[1] for t in sorted([(i,n) for (n, i) in station_map.items()])]
# stations2 = sorted(map(station_map.get, station_names))
stations2 = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11]
n_stations = len(stations2)

# ref_station_name = '1' # Effelsberg
ref_station_name = 'EF' 
ref_station = station_map[ref_station_name]

polind = 0
solint = 600
threshold = 0.0

all_delays = []
all_phases = []
all_rates = []
for timeq, swid in itertools.product(timeqs, range(8)):
    timeq2 = ffd.actual_timerangeq(msname, timeq)
    delays, phases, rates, sigs = [], [], [], []
    for pol_id in [0,3]:
        anffd = ffd.FFData.make_FFD(msname, stations2, swid, polind, pol_id, timeq2,
                                   solint=solint)
        anffd.fft(pad=32)
        bl_ref_station = 0 # Ef
        anffd.weights = utils.make_ref_weights(anffd, stations2, bl_ref_station)
        params2 = anffd.get_params(ref_station=ref_station)
        if threshold > 0.0:
            anffd.mask_fringes_below_threshold(threshold)
        e_ref_station = anffd.get_station_index(ref_station)
        v_dash, refs = lsqrs.remove_ref(params2, e_ref_station)
        args=(anffd.tgrid0, anffd.fgrid0,
              len(stations2), anffd.data, anffd.weights, e_ref_station)
        res_t = scipy.optimize.leastsq(lsqrs.vector_s3_test, v_dash,
                                       full_output=1, args=args,
                                       col_deriv=True, Dfun=lsqrs.matrix_j_s3)
        sol_out = list(res_t[0])
        pout2d = lsqrs.restore_ref(sol_out, refs, e_ref_station)
        dels = lsqrs.get_delays(pout2d)
        phs = lsqrs.get_phases(pout2d)
        rs = lsqrs.get_rates(pout2d)/anffd.ref_freq
        delays.append(dels)
        phases.append(phs)
        rates.append(rs)
    all_delays.append(delays)
    all_phases.append(phases)
    all_rates.append(rates)
    print "#\n#{}, IF={}\n#".format(timeq, swid)
    utils.print_res(station_map, stations2, phases, delays, rates)

all_delays = np.array(all_delays)
all_rates = np.array(all_rates)
all_phases = np.array(all_phases)



sm = utils.get_station_map(msname, 'NAME')
ism = utils.invert_map(sm)
stationt = ",".join([ism[s] for s in stations2])
spwt = ",".join(str(i) for i in swids)
freqs=utils.get_reference_freqs(msname)
l_delay = utils.flattenificate(-1e9*all_delays)
l_phase = utils.flattenificate(-1*180/np.pi*all_phases)

if False:
    gencal(vis=msname, caltable='huh.G', caltype='sbd',
           spw=spwt, antenna=stationt,   pol='R,L',
           parameter = l_delay)
    gencal(vis=msname, caltable='huh_ph.G', caltype='ph',
           spw=spwt, antenna=stationt,   pol='R,L',
           parameter = l_phase)



