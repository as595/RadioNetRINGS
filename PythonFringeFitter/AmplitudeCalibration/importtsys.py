import os
import sys
import key
import math
import tempfile
import time
import StringIO
import scipy
import numpy as np

antab = 'n14c1.antab'
vis = 'n14c1.ms'

os.environ['TZ'] = 'UTC'
time.tzset()

columnnames = [
    "ANTENNA_ID",
    "ARRAY_ID",
    "FEED_ID",
    "INTERVAL",
    "SPECTRAL_WINDOW_ID",
    "TIME",
    "NUM_RECEPTORS",
    "TSYS"
]

datatypes = [
    "I",
    "I",
    "I",
    "D",
    "I",
    "D",
    "I",
    "R2"
]

tsys_values = {}
tsys_times = {}

def update_map(spws, spwmap, index):
    idx = 0
    for labels in index:
        for label in labels.split('|'):
            pol = label[0]
            rng = label[1:].split(':')
            if pol != 'X':
                if len(rng) == 1:
                    rng.append(rng[0])
                    pass
                rng = [int(x) - 1 for x in rng]
                for spw in range(rng[0], rng[1] + 1):
                    if not spw in spws:
                        spws.append(spw)
                        pass
                    spwmap[(pol, spw)] = idx
                    continue
                pass
            continue
        idx += 1
        continue
    spws = sorted(spws)
    return

def process_values(infp, outfp, keys, vis):
    tb.open(vis)
    secs = tb.getcell('TIME') - (40587 * 86400)
    year = time.gmtime(secs).tm_year
    tb.close()
    tb.open(vis + '/ANTENNA')
    namelist = tb.getcol("NAME").tolist()
    antenna = namelist.index(keys[0][1][0])
    tb.close()
    keys = dict(keys[0])
    scan = 0
    pols = ['R', 'L']
    spws = []
    spwmap = {}
    update_map(spws, spwmap, keys['INDEX'])
    if 'INDEX2' in keys:
        update_map(spws, spwmap, keys['INDEX2'])
        pass
    timeoff = 0
    if 'TIMEOFF' in keys:
        timeoff = float(keys['TIMEOFF'])
    for line in infp:
        if line.startswith('!'):
            continue
        fields = line.split()
        if len(fields) > 1:
            tm_year = year
            tm_yday = int(fields[0])
            tm_hour = int(fields[1].split(':')[0])
            tm_min = math.modf(float(fields[1].split(':')[1]))
            tm_sec = int(60 * tm_min[0])
            tm_min = int(tm_min[1])
            t = "%dy%03dd%02dh%02dm%02ds" % \
                (tm_year, tm_yday, tm_hour, tm_min, tm_sec)
            t = time.strptime(t, "%Yy%jd%Hh%Mm%Ss")
            secs = time.mktime(t) + timeoff
            values = fields[2:]
            secs = secs + (40587.0 * 86400)
            if secs <= scan_times[-1][1]:
                while secs > scan_times[scan][1]:
                    scan += 1
                    continue
                for pol in pols:
                    for spw in spws:
                        idx = (antenna, scan, spw)
                        if not idx in tsys_values:
                            tsys_values[idx] = {}
                            tsys_times[idx] = {}
                            pass
                        if not pol in tsys_values[idx]:
                            tsys_values[idx][pol] = []
                            tsys_times[idx][pol] = []
                            pass
                        value = float(values[spwmap[(pol, spw)]])
                        if value > 0:
                            tsys_values[idx][pol].append(value)
                            tsys_times[idx][pol].append(secs)
                            pass
                        continue
                    continue
                pass
            pass
        if line.strip().endswith('/'):
            break
        continue
    return

def write_values(outfp):
    pols = ['R', 'L']

    keys = tsys_values.keys()
    for idx in sorted(keys):
        antenna = idx[0]
        scan = idx[1]
        spw = idx[2]
        x = tsys_times[idx]
        y = tsys_values[idx]
        for pol in pols:
            if len(y[pol]) == 0:
                x[pol] = [(scan_times[scan][0] + scan_times[scan][1]) / 2]
                y[pol] = [-1.0]
                pass
            continue

        secs = scan_times[scan][0]
        while secs <= scan_times[scan][1]:
            # ANTENNA_ID
            print >> outfp, antenna,
            # ARRAY_ID
            print >> outfp, 0,
            # FEED_ID
            print >> outfp, 0,
            # INTERVAL
            print >> outfp, 0,
            # SPECTRAL_WINDOW_ID
            print >> outfp, spw,
            # TIME
            print >> outfp, secs,
            # NUM_RECEPTORS
            print >> outfp, 2,
            # TSYS
            for pol in pols:
                print >> outfp, scipy.interp(secs, x[pol], y[pol]),
                continue
            print >> outfp
            secs += 30
            continue
        continue
    return

keys = StringIO.StringIO()
section = 0

ms.open(vis)
scans = ms.getscansummary()
ms.close()

scan_times = []
for scan in scans:
    integration_time = scans[scan]['0']['IntegrationTime']
    start = scans[scan]['0']['BeginTime'] * 86400 - integration_time
    end = scans[scan]['0']['EndTime'] * 86400 + integration_time
    scan_times.append([start, end])
scan_times = sorted(scan_times)

outfp = tempfile.NamedTemporaryFile('w')

fp = open(antab, 'r')
for line in fp:
    if line.startswith('!'):
        continue
    keys.write(line)
    if line.strip().endswith('/'):
        keys.seek(0)
        tsys = key.read_keyfile(keys)
        if tsys[0][0][0] == 'TSYS':
            process_values(fp, outfp, tsys, vis)
            pass
        keys = StringIO.StringIO()
        continue
    continue

write_values(outfp)

outfp.flush()

tb.open(vis)
unit = tb.getcolkeyword("TIME", "QuantumUnits")
meas = tb.getcolkeyword("TIME", "MEASINFO")
tb.close()

tb.fromascii(tablename=vis + '/SYSCAL', asciifile=outfp.name, sep=' ',
             columnnames=columnnames, datatypes=datatypes)
tb.open(vis + '/SYSCAL', nomodify=False)
tb.putcolkeyword("INTERVAL", "QuantumUnits", unit)
tb.putcolkeyword("TIME", "QuantumUnits", unit)
tb.putcolkeyword("TIME", "MEASINFO", meas)
tb.putcolkeyword("TSYS", "QuantumUnits", "K")
tb.close()

tb.open(vis, nomodify=False)
tb.putkeyword('SYSCAL', 'Table: ' + vis + '/SYSCAL')
tb.close()

outfp.close()
fp.close()
