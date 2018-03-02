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
n_stations = len(stations2)

# ref_station_name = '1' # Effelsberg
ref_station_name = '3' 
ref_station2 = station_map[ref_station_name]
ref_s_ind2 = stations2.index(ref_station2)


polind = 0
ref_station = ref_station2
solint = 300
threshold = 0.0

def make_ref_weights(stations2, bl_ref_station):
    weights = np.zeros((n_stations, n_stations, 1, 1))
    for i0 in stations2:
        i = anffd.get_station_index(i0)
        if i==bl_ref_station: continue
        weights[i, bl_ref_station] = 1.0
        weights[bl_ref_station, i] = 1.0
    return weights
        
for timeq, swid in itertools.product(timeqs[:1], range(4)):
    timeq2 = ffd.actual_timerangeq(msname, timeq)
    delays, phases, rates, sigs = [], [], [], []
    for pol_id in [0,3]:
        anffd = ffd.FFData.make_FFD(msname, stations2, swid, polind, pol_id, timeq2,
                                   solint=solint)
        anffd.fft(pad=32)
        bl_ref_station = 0 # Ef
        weights = make_ref_weights(stations2, bl_ref_station)
        # weights = np.ones((n_stations, n_stations, 1, 1))
        # weights = utils.total_perspective_vortex(msname, stations2, timeq2)
        anffd.weights = weights
        params2 = anffd.get_params(ref_station=ref_station)
        if threshold > 0.0:
            anffd.mask_fringes_below_threshold(threshold)
        e_ref_station = anffd.get_station_index(ref_station)
        r_scale, tau_scale = 1e9, 1e5
        params2d = lsqrs.scale_params(params2, [1, r_scale, tau_scale])
        v_dash, refs = lsqrs.remove_ref(params2d, e_ref_station)
        args=(anffd.tgrid0/r_scale, anffd.fgrid0/tau_scale,
              len(stations2), anffd.data, anffd.weights, e_ref_station)
        res_t = scipy.optimize.leastsq(lsqrs.vector_s3_test, v_dash,
                                       full_output=1, args=args,
                                       col_deriv=True, Dfun=lsqrs.matrix_j_s3)
        sol_out = list(res_t[0])
        # http://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-\
        # using-the-optimize-leastsq-method-i
        # all this to get sigma_p
        sigma_p = lsqrs.get_phases(sol_out)
        sol = lsqrs.restore_ref(sol_out, refs, e_ref_station)
        pout2d = lsqrs.scale_params(sol, [1, 1/r_scale, 1/tau_scale])

        dels = lsqrs.get_delays(pout2d)
        phs = lsqrs.get_phases(pout2d)
        rs = lsqrs.get_rates(pout2d)/anffd.ref_freq
        delays.append(dels)
        phases.append(phs)
        rates.append(rs)
        sigs.append(sig)
    print "#\n#{}, IF={}\n#".format(timeq, swid)
    utils.print_res(station_map, stations2, delays, rates)


# inlining
def fit_fringe_lm(msname, stations2, ref_station, swid, polind, pol_id, timeq,
                  solint=None, threshold=0.0):
    anffd = ffd.FFData.make_FFD(msname, stations2, swid, polind, pol_id, timeq, solint=solint)
    anffd.fft(pad=32)
    params2 = anffd.get_params(ref_station=ref_station)
    if threshold > 0.0:
        anffd.mask_fringes_below_threshold(threshold)
    e_ref_station = anffd.get_station_index(ref_station)
    r_scale, tau_scale = 1e9, 1e5
    params2d = lsqrs.scale_params(params2, [1, r_scale, tau_scale])
    v_dash, refs = lsqrs.remove_ref(params2d, e_ref_station)
    args=(anffd.tgrid0/r_scale, anffd.fgrid0/tau_scale,
          len(stations2), anffd.data, anffd.weights, e_ref_station)
    res_t = scipy.optimize.leastsq(lsqrs.vector_s3_test, v_dash,
                                   full_output=1, args=args,
                                   col_deriv=True, Dfun=lsqrs.matrix_j_s3)
    sol_out = list(res_t[0])
    # http://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-\
    # using-the-optimize-leastsq-method-i
    # all this to get sigma_p
    if False:
        cov_x = res_t[1]
        v = (sol_out,) + args
        dv = lsqrs.vector_s3_test(*v)
        print dv
        s_sq = ((dv)**2).sum()/(len(dv)-len(sol_out))
        pcov = cov_x * s_sq
        sigma_p0 = np.sqrt(np.diagonal(pcov))
        sigma_p1 = lsqrs.scale_params(sigma_p0,  [1, 1/r_scale/anffd.ref_freq, 1/tau_scale])
        sigma_p = lsqrs.restore_ref(list(sigma_p1), [0,0,0], e_ref_station)
        # end sigma_p block
    else:
        sigma_p = lsqrs.get_phases(sol_out)
    sol = lsqrs.restore_ref(sol_out, refs, e_ref_station)
    pout2d = lsqrs.scale_params(sol, [1, 1/r_scale, 1/tau_scale])

    dels = lsqrs.get_delays(pout2d)
    phs = lsqrs.get_phases(pout2d)
    rs = lsqrs.get_rates(pout2d)/anffd.ref_freq
    return dels, phs, rs, sigma_p

import sys
sys.exit(0)

gencal(vis=msname,
       caltable='man_pcal_r_3.G',
       caltype='sbd',
       spw=spwt,
       antenna=stationt,
       pol='R',
       parameter = -plist)

applycal(vis=msname, gaintable="man_pcal_r_3.G", field="*")
stationt = ",".join(map(str, stations2))
spwt = ",".join(map(str, range(8)))
darr = np.array(delays)
plist = list(1e9*darr.transpose())

gencal(vis=msname,
       caltable='man_ph_r_min.G',
       caltype='ph', # or ph or mbd
       spw=spwt,
       antenna=stationt,
       pol='R',
       parameter =180/np.pi*phlist)
