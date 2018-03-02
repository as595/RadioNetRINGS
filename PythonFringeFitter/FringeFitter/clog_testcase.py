# Casalog filter levels include:
# DEBUG1 DEBUG2 ERROR WARN INFO INFO1 INFO2 INFO3 INFO4 INFO5 DEBUG DEBUG1 DEBUG2 INFO
import clog_classified

casalog.filter(level="DEBUG")
# anmsname = globals()['msname']
anmsname = '/scratch/small/N14C1/n14c1.ms.calibrated'
ff = clog_classified.FringeFitter(anmsname, ref_antenna_name='EF', scans=[3])
ff.run()
ff.write_tables('man_pcal_single.G', 'man_ph.G', ff.phases, ff.delays, ff.rates)

