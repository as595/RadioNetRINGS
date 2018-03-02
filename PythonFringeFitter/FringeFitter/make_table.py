# to use this code on casadev, do the following:
# 
# $ . /home/casa/src/casa/casainit.sh
# $ export PYTHONPATH=/home/casa/src/casa/linux64/python/2.7
#
# After doing that, you can just run the script using the normal python
# interpreter.
# 
# Cheers,
# 
# Mark

import sys
sys.path.insert(1, '/home/casa/src/casa/linux64/python/2.7')

from __casac__.calibrater import *
from __casac__.table import *
import numpy as np

def make_table(msname, caltable): # both strings
    cb = calibrater()
    cb.open(msname, addcorr=False, addmodel=False)
    cb.createcaltable(caltable, 'Real', 'Fringe Jones', True)
    cb.close()
    
def add_row(ctname, i, time, interval, antenna, field, scan, obsid, swid, param):
    # (6,1) array of Phase (radians), delay (currently nanoseconds; will
    # likely change to seconds), rates (s/s; i.e., *sensibly*
    # dimensionless, as opposed to picoseconds per second which occurs
    # elsewhere).
    tb = table()
    tb.open(ctname, nomodify=False)
    tb.addrows(1) # nrows to add
    # The following parameters are stubbed for now,
    paramerr = -np.ones(shape=(6,1), dtype='float32')
    flag = np.zeros(shape=(6,1), dtype='bool')
    snr = np.ones(shape=(6,1), dtype='float32')
    weight = np.ones(shape=(6,1), dtype='float32')
    [tb.putcell(k, i, v) for (k, v) in [('TIME', time),
                                        ('INTERVAL', interval),
                                        ('ANTENNA1', antenna),
                                        ('ANTENNA2', -1),
                                        ('FIELD_ID', field),
                                        ('SCAN_NUMBER', scan),
                                        ('OBSERVATION_ID', obsid),
                                        ('SPECTRAL_WINDOW_ID', swid),
                                        ('FPARAM', param),
                                        ('PARAMERR', paramerr),
                                        ('FLAG', flag),
                                        ('SNR', snr),
                                        ('WEIGHT', weight)]]
    tb.close()

# ms = "test.ms"
# ct = "test.fringecal"
# make_table(ms, ct)


# # s = '2016-07-06T13:55:10'
# # tu = qa.quantity(s) # => {'unit': 'd', 'value': 57575.579976851855}
# # qa.time(tu,form="fits") # => ['2016-07-06T13:55:10']


# time = 4920700000.0 # Casa uses seconds after something.
# interval = 600
# antenna = 1
# spwid = 3
# field = 1
# scan = 0
# obsid = 0

# param = np.zeros(shape=(6,1), dtype='float32')

