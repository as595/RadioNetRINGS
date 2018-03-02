# Casalog filter levels include:
# DEBUG1 DEBUG2 ERROR WARN INFO INFO1 INFO2 INFO3 INFO4 INFO5 DEBUG DEBUG1 DEBUG2 INFO
import clog_classified

casalog.filter(level="DEBUG")
# anmsname = globals()['msname']
anmsname = '/scratch/small/N14C1/n14c1.ms.calibrated'
ff2 = clog_classified.MultiBandFringeFitter(anmsname, 'mbd.fj', ref_antenna_name='EF', scans=[3])
ff2.run()

