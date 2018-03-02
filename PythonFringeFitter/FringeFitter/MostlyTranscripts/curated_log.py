import numpy as np, operator
import lsqrs, full, ffd

msname = '/jop88_1/small/Datas/EY015D/EY015D.ms'
# msname = '/jop88_1/small/Datas/N14C1/N14C1.ms'
swid, pol_id, polind, = (1, 0, 0)
timeq = 'TIME > MJD("2012-10-24-21:29:10") and TIME < MJD("2012-10-24-21:30:00")'
stations2 = [0, 1, 4]
ref_station2 = 0
ref_s_ind2 = stations2.index(ref_station2)
#
anffd = ffd.FFData.make_FFD(msname, [0, 1, 4], swid, polind, pol_id, timeq, startrow=0)
anffd.fft()
params2 = anffd.get_params()
weights2 = np.ones(anffd.data.shape)
r_scale, tau_scale = 1e4, 1e7
params2d = lsqrs.scale_params(params2, [1, r_scale, tau_scale])
pout2 = lsqrs.minimise_wrapper(lsqrs.fun_s3, params2d,
                               args=(anffd.tgrid0/r_scale, anffd.fgrid0/tau_scale,
                                     len(stations2), anffd.data, anffd.weights, ref_s_ind2),
                               method='TNC', jac=lsqrs.jac_s3, options={'maxiter':200})[0]
pout2d = lsqrs.scale_params(pout2, [1, 1/r_scale, 1/tau_scale])
anffd.plot_fringes(pout2d)

