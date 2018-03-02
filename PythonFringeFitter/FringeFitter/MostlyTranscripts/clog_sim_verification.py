import os, glob
import lsqrs, ffd, fringer, utils

timeqs = ['TIME > MJD("2013-06-03-20:10:00") and TIME < MJD("2013-06-03-20:15:00")']


def process_ms(msname, outfn):
    station_map = utils.get_station_map(msname)
    station_names = sorted(station_map.keys())
    stations2 = sorted(map(station_map.get, station_names))
    ref_station_name = '1' # ANTENNA table has ID, Name and 'Station' fields'; here we want the last of these
    ref_station2 = station_map[ref_station_name]
    ref_s_ind2 = stations2.index(ref_station2)
    swid, pol_id, polind, = (0, 0, 0)

    for timeq in timeqs:
        delays, phases, rates, sigs = [], [], [], []
        for pol_id in [0,3]:
            dels, phs, rs, sig = fringer.fit_fringe_lm(msname, stations2, ref_s_ind2,
                                                       swid, polind, pol_id, timeq, solint=300)
            delays.append(dels)
            phases.append(phs)
            rates.append(rs)
            sigs.append(sig)
        outf = file(outfn, 'w')
        print >>outf, "#\n#{}\n#".format(timeq)
        utils.print_res(station_map, stations2, phases, delays, rates, sigs=sigs, outf=outf)


msname = globals()['msname']
outfn = os.path.splitext(msname)[0] + '-1-old.out' 
process_ms(msname, outfn)

    
