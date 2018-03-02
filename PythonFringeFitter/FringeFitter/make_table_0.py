from __casac__.calibrater import *
from __casac__.table import *
import numpy as np

ms = "test.ms"
ct = "test.fringecal"

cb = calibrater()

cb.open(ms, addcorr=False, addmodel=False)
cb.createcaltable(ct, 'Real', 'Fringe Jones', True)
cb.close()

tb = table()

time = 4920700000.0
interval = 600
antenna = 1
field = 1
scan = 0
obsid = 0
spwid = 3

param = np.zeros(shape=(6,1), dtype='float32')
paramerr = -np.ones(shape=(6,1), dtype='float32')
flag = np.zeros(shape=(6,1), dtype='bool')
snr = np.ones(shape=(6,1), dtype='float32')
weight = np.ones(shape=(6,1), dtype='float32')

tb.open(ct, nomodify=False)
tb.addrows(1)
tb.putcell('TIME', 0, time)
tb.putcell('INTERVAL', 0, interval)
tb.putcell('ANTENNA1', 0, antenna)
tb.putcell('ANTENNA2', 0, -1)
tb.putcell('FIELD_ID', 0, field)
tb.putcell('SCAN_NUMBER', 0, scan)
tb.putcell('OBSERVATION_ID', 0, obsid)
tb.putcell('SPECTRAL_WINDOW_ID', 0, spwid)
tb.putcell('FPARAM', 0, param)
tb.putcell('PARAMERR', 0, paramerr)
tb.putcell('FLAG', 0, flag)
tb.putcell('SNR', 0, snr)
tb.putcell('WEIGHT', 0, weight)
tb.close()
