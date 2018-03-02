import clog_classified
import tasks


casalog.filter(level="DEBUG")
anmsname = 'n14c2.ms'
fj_table = "sbd_gah.fj"

# Set corrected data to original
## Single bands or "manual phase cals
tasks.flagmanager(anmsname, versionname='applycal_10', mode='restore')
casalog.post("Reset CORRECTED_DATA to original values.")

# Set the boolean flags to do single-band delays, multi-band delays or both.
if False:
    tasks.applycal(anmsname, gaintable=['n14c2.gcal', 'n14c2.tsys'], parang=True)
    ff = clog_classified.FringeFitter(anmsname, fj_table, ref_antenna_name='EF', scans=[5],
                                      # antennas=None,
                                      # antennas = [0,4,6,7,8,10],
                                      snr_threshold=30.0,
                                      solint=60.0,
                                      solsub=2)
    ff.run()
    new_antenna_list = [s for s in ff.antennas2 if s not in ff.bad_antennas]
tasks.applycal(anmsname, gaintable=['n14c2.gcal', 'n14c2.tsys', fj_table], parang=True)


if True:
    new_antenna_list = [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15]
    ## Compute multiband calibration
    mb_table = 'mbd_meh.fj'
    ff2 = clog_classified.MultiBandFringeFitter(anmsname, mb_table, ref_antenna_name='EF',
                                                scans=[29], antennas=new_antenna_list,
                                                snr_threshold=30.0, solint=60.0, solsub=2)
    ff2.run()
    ## Then apply to scan to check flatness
    tasks.applycal(anmsname, gaintable=['n14c2.gcal', 'n14c2.tsys', fj_table, mb_table], scan='29', parang=True)
