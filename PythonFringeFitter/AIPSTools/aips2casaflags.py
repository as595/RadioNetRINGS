#! /usr/bin/env ParselTongue

import datetime
import sys

import AIPS
from AIPSData import AIPSUVData

AIPS.userno = int(sys.argv[1])

name = sys.argv[2]
klass = sys.argv[3]
disk = int(sys.argv[4])
seqno = int(sys.argv[5])
uvdata = AIPSUVData(name, klass, disk, seqno)
date_obs = uvdata.header['date_obs']
date_obs = datetime.datetime.strptime(date_obs, '%Y-%m-%d')
nchan = uvdata.header.naxis[uvdata.header.ctype.index('FREQ')]
f = open(sys.argv[6], "w")

fgtable = uvdata.table('FG', 0)
for row in fgtable:
    cmd = ''
    if row.time_range[0] != 0 or row.time_range[0] != 0:
        date0 = date_obs + datetime.timedelta(row.time_range[0])
        timerang0 = date0.strftime("%Y/%m/%d/%H:%M:%S")
        date1 = date_obs + datetime.timedelta(row.time_range[1])
        timerang1 = date1.strftime("%Y/%m/%d/%H:%M:%S")
        cmd = cmd + "timerange='%s~%s'" % (timerang0, timerang1)
        pass
    if row.ants[0] > 0:
        if cmd:
            cmd = cmd + ' '
            pass
        cmd = cmd + "antenna='%d" % (row.ants[0] - 1)
        if row.ants[1] > 0:
            cmd = cmd + "&%d" % (row.ants[1] - 1)
            pass
        cmd = cmd + "'"
        pass
    if row.source > 0:
        if cmd:
            cmd = cmd + ' '
            pass
        cmd = cmd + "field='%d'" % (row.source - 1)
        pass
    spw = ''
    if row.ifs[0] > 0:
        spw = spw + "%d" % (row.ifs[0] - 1)
        if row.ifs[1] > 0:
            spw = spw + "~%d" % (row.ifs[1] - 1)
            pass
        pass
    if row.chans[0] > 0:
        if spw:
            spw = spw + ":"
            pass
        spw = spw + "%d" % (row.chans[0] - 1)
        if row.chans[1] == 0:
            spw = spw + "~%d" % (nchan - 1)
        elif row.chans[1] > row.chans[0]:
            spw = spw + "~%d" % (row.chans[1] - 1)
            pass
        pass
    if spw:
        cmd = cmd + " spw='%s'" % spw
        pass
    f.write(cmd + "\n")
    continue

f.close()
