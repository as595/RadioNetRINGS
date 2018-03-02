import casa
from taskinit import casalog
import numpy as np, scipy, itertools
import glob, lsqrs, ffd, fringer, utils
import make_table

# NB: this assumes manual phase cal has been done first.

# Don't forget to Source The Magic:
#. /home/casa/src/casa/casainit.sh
# and use the blessed version of Casa
#/home/casa/src/casa/linux64/bin/casa

# You will need to set the path to Kettenis's magic version:
# import sys
# sys.path.insert(1, '/home/casa/src/casa/linux64/python/2.7')
# from __casac__.calibrater import *

def fit_multiband_fringe(msname, scan_number, ctname):
    ms.open(msname)
    timeqs = ["SCAN_NUMBER={}".format(scan_number)]
    station_map = utils.get_station_map(msname)
    stations2 = sorted(ffd.actual_stations(msname, timeqs[0]).keys())
    ism = utils.invert_map(station_map)
    station_names = [ism[s] for s in stations2]

    ref_station_name = 'EF' 
    ref_station2 = station_map[ref_station_name]
    ref_s_ind2 = stations2.index(ref_station2)
    polind = 0
    swids = range(8)
    reffreqs = utils.get_min_freqs(msname)
    minfreq = utils.get_min_freqs(msname)[0]

    make_table.make_table(msname, ctname)

    shape = (2, len(stations2))
    delays = np.zeros(shape, np.float)
    phases = np.zeros(shape, np.float)
    rates = np.zeros(shape, np.float)
    sigs = []

    rowcount = 0
    for timeq in timeqs[:1]:
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

        obsid, field, scan = [ffd.distinct_thing(msname, timeq, col)
                              for col in ['OBSERVATION_ID', 'FIELD_ID', 'SCAN_NUMBER']]
        darr = -delays*1e9
        pharr = -phases # radians!
        rarr = -rates
        for i,s in enumerate(stations2):
            antenna = s
            assert (anffd.get_station_index(s) == i)
            time = anffd.get_ref_time()
            # time = anffd.times[0]
            interval = anffd.get_interval()
            midfreqs = utils.get_mid_freqs(msname)
            for swid in swids:
                # Df = (reffreqs[swid]-minfreq)
                Df = (midfreqs[swid]-minfreq)
                phase_offsets =  utils.turns_to_radians(Df * darr/1e9 +
                                                      interval/2 * rarr)

                ph = pharr + phase_offsets

                param = np.zeros(shape=(6,1), dtype='float32')
                param[:, 0] = [ph[0, i], darr[0, i], rates[0,i], 
                               ph[1, i], darr[1, i], rates[1,i] ]
                make_table.add_row(ctname, rowcount, time, interval, antenna, 
                                  field, scan, obsid, swid, param)
                rowcount += 1
            



# apply preliminary calibrations
#
# casa.applycal(vis=msname, gaintable=["ec047a.gc", "ec047a.tsys", "man_pcal_single.G", "man_ph.G"], scan='44', parang=True)

#
# To see results, first apply with new calibration table:
# sscan_number = str(scan_number)
# casa.applycal(vis=msname, gaintable=['ec047a.gc', 'ec047a.tsys', 'man_pcal_single.G',  'man_ph.G', ctname],
#               interp=['', '', '', '', 'linearPR'], scan=sscan_number, parang=True)
#
# And then plot like this: 
#
# casa.plotms(vis=msname, xaxis="frequency", yaxis="phase", spw="0~7", antenna="EF", averagedata=T, avgtime="5000", correlation="RR", ydatacolumn="corrected", scan="44", iteraxis="baseline")

scan_number = 36
ctname = "multi-{}.fj".format(scan_number)
msname = globals()['msname']


