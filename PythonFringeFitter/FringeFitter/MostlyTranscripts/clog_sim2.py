import glob, lsqrs, ffd, fringer, utils

timeq = 'TIME > MJD("2013-06-03-20:30:00") and TIME < MJD("2013-06-03-20:35:00")'
timeq = 'TIME > MJD("2013-06-03-20:35:00") and TIME < MJD("2013-06-03-20:40:00")'

msname = "sim2_data.ms"

station_map = utils.get_station_map(msname)
station_names = sorted(station_map.keys())
stations2 = sorted(map(station_map.get, station_names))

ref_station_name = '7'
ref_station2 = station_map[ref_station_name]
ref_s_ind2 = stations2.index(ref_station2)
swid, pol_id, polind, = (0, 0, 0)


delays, phases, rates = [], [], []
for pol_id in [0,3]:
    dels, phs, rs = fringer.fit_fringe_lm(msname, stations2, ref_s_ind2, swid, polind, pol_id, timeq, solint=300)
    delays.append(dels)
    phases.append(phs)
    rates.append(rs)

def print_res(stations2, delays, rates):
    fmts = " ".join(["{:3>}"] + 4*["{:>20}"])
    print fmts.format('', "Delay1", "Delay2", "Rate1", "Rate2")
    for i, s in enumerate(stations2):
        d1, r1 = delays[0][i], rates[0][i]
        d2, r2 = delays[1][i], rates[1][i]
        print fmts.format(s, d1, d2, r1, r2)

    
