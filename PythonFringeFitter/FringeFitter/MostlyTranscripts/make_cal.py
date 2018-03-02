# flagmanager(vis=msname, mode='list')
def get_nchans(vis):
    tb.open(vis)
    q = 'select distinct NUM_CHAN from {}::SPECTRAL_WINDOW'.format(vis)
    (nchans,) = tuple(tb.taql(q).getcol('NUM_CHAN'))
    return nchans

# Flagging command:
def flag_edges(vis, width=2):
    nchans = get_nchans(vis)
    s = '*:{}~{};{}~{}'.format(0, width, 
                               nchans-width-1, nchans-1)
    flagdata(msname, spw=s, mode='manual', flagbackup=True))
    return None

## FINALs

delays = []
for swid in range(8):
    anffd = ffd.FFData.make_FFD(msname, stations2, swid, polind, 
                               pol_id, timeq, startrow=0)
    anffd.fft()
    params2 = anffd.get_params()
    r_scale, tau_scale = 1e4, 1e7
    params2d = lsqrs.scale_params(params2, [1, r_scale, tau_scale])
    pout2 = lsqrs.minimise_wrapper(
                  lsqrs.fun_s3, params2d,
                  args=(anffd.tgrid0/r_scale, anffd.fgrid0/tau_scale,
                        len(stations2), anffd.data, anffd.weights, ref_s_ind2),
                   method='TNC', jac=lsqrs.jac_s3, options={'maxiter':200})[0]
    pout2d = lsqrs.scale_params(pout2, [1, 1/r_scale, 1/tau_scale])
    dels = lsqrs.get_delays(pout2d)
    delays.append(dels)
    
plist = np.ndarray.flatten(1e9*darr)
spwt = ",".join(map(str, range(8)))

gencal(vis=msname, caltable='man_pcal_r_3.G',
       caltype='sbd',
       spw=spwt,
       antenna=stationt,
       pol='R',
       parameter = -plist)


gencal(vis=msname,
       caltable='man_ph_r_min.G',
       caltype='ph',
       spw=spwt,
       antenna=stationt,
       pol='R',
       parameter =180/np.pi*phlist)

# applycal(vis=msname, gaintable="man_pcal_r_3.G", field="*")
# plotms(vis=msname, xaxis='frequency', yaxis='phase',
#        scan='85', spw='*', averagedata=True,
#        avgtime='600s', # don't know how to say "long enough"
#        iteraxis='baseline', correlation="rr", antenna="EF&ON")
