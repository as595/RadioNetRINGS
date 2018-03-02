#!/usr/bin/env ParselTongue
#
# Tool to convert AIPS Calibration tables to Measurement Set Calibration tables
# Stephen Bourke, JIVE
# Version 2010.11.29

import AIPS
from AIPSTask import AIPSTask
from AIPSData import AIPSUVData
import Wizardry.AIPSData as wiz
from pyrap.tables import table
import sys
import shutil
import datetime
import calendar
import numpy
import os.path
from __casac__.calibrater import *

AIPS.userno = 666

class progress:
	__doc__ = 'Print out progress as dots and percentages'

	def __init__(self, done_value=100.0, dot_value=2.51):
		self.done_value = done_value
		self.progress_next = 0
		self.dot_value = dot_value
		self.progress_dot = dot_value

	def update(self, value):
		progress_curr = float(value) / self.done_value * 100
		if progress_curr > self.progress_next:
			sys.stdout.write('%d%%' % self.progress_next)
			sys.stdout.flush()
			self.progress_dot = self.progress_next + self.dot_value
			self.progress_next += 10
		if progress_curr > self.progress_dot:
			sys.stdout.write('.')
			sys.stdout.flush()
			self.progress_dot += self.dot_value

	def done(self):
		sys.stdout.write('100%\n')
		sys.stdout.flush()

cb = calibrater()

mod_jul_day = 3506716800.0

# Mode can be BP (bandpass), SN (solution), or CL (calibration)
if '-b' in sys.argv:
	CAL_MODE = 'BP'
	del sys.argv[sys.argv.index('-b')]
elif '-s' in sys.argv:
	CAL_MODE = 'SN'
	del sys.argv[sys.argv.index('-s')]
else:
	CAL_MODE = 'CL'
DOBAND = CAL_MODE == 'BP'
DOCALIB = not DOBAND
	
if len(sys.argv) != 5 and len(sys.argv) != 6:
	print >> sys.stderr, 'Usage: %s <file> <inver> <ms> <outcal> [-b|-s] [antfile]' % sys.argv[0]
	sys.exit()
inver = int(sys.argv[2])
ms_name, cal_base = sys.argv[3:5]
try:
	antfile = sys.argv[6]
except IndexError:
	antfile = None

path = sys.argv[1]
if not path.startswith("/") and not path.startswith("."):
	path = "./" + path

uv = AIPSUVData('CASA', 'TASAV', 1, 1)

fitld = AIPSTask('fitld')
fitld.datain = path
fitld.outdata = uv
fitld.go()

# Select dataset and table
uv = wiz.AIPSUVData(uv)
obs_date = datetime.datetime(*[int(i) for i in uv.header.date_obs.split('-')])
print AIPSUVData(uv), '%s%d' % (CAL_MODE, inver), '->', ms_name, cal_base

cal_aips = uv.table(CAL_MODE, inver)
pols = cal_aips.keywords['NO_POL']
ifs = cal_aips.keywords['NO_IF']	# AIPS IF / Casa SPW
if DOBAND:
	chans = cal_aips.keywords['NO_CHAN']
else:
	chans = 1

# User specified mapping from AIPS ants to MS ant ids
ant_mapping = {}
if antfile:
	inf = open(antfile)
	for line in inf:
		aips_ant_no, ms_ant_id = [int(x) for x in line.split()]
		ant_mapping[aips_ant_no] = ms_ant_id
	inf.close()

def ant_aips2ms(ant_num):
	if ant_mapping:
		return ant_mapping[ant_num]
	else:
		return ant_num - 1

prgs = progress(len(cal_aips))
sys.stdout.write('Reading %s table\n' % CAL_MODE)
sys.stdout.flush()

# Holds all table data as a multi-level dict indexed by time, antenna,
# IF, mode. Where mode is one of 'common', 'sbdelay', 'mbdelay'. Format
# is the same as in the final Table. Casa requires a fairly regular
# table with rows present for all antennas, and IFs. 
cals = {}
	  
for i in range(len(cal_aips)):
	prgs.update(i) # Progress bar
	sol_time = obs_date + datetime.timedelta(float(cal_aips[i]['time']))
	time_unix = calendar.timegm(sol_time.timetuple())
	time_casa = time_unix + mod_jul_day
	if time_casa not in cals:
		cals[time_casa] = {}
	if DOBAND:
		time_interval = cal_aips[i]['interval']
	else:
		time_interval = cal_aips[i]['time_interval']
	time_interval *= 24 * 60 * 60 # Convert days to seconds
	field_id = cal_aips[i]['source_id'] - 1
	try:
		if DOBAND:
			antenna1 = ant_aips2ms(cal_aips[i]['antenna'])
		else:
			antenna1 = ant_aips2ms(cal_aips[i]['antenna_no'])
	except KeyError:
		# Ignore antennas not present in MS
		continue
	if antenna1 not in cals[time_casa]:
		cals[time_casa][antenna1] = {}
	for iif in range(ifs):
		if not iif in cals[time_casa][antenna1]:
			cals[time_casa][antenna1][iif] = {}
		row = iif * len(cal_aips) + i
		gain = numpy.ones(shape=(chans,pols), dtype='complex64')
		mbdelay = numpy.ones(shape=(chans,pols), dtype='float32')
		sbdelay = numpy.zeros(shape=(chans,3*pols), dtype='float32')
		snr = numpy.ones(shape=(chans,pols), dtype='float32')
		sbdsnr = numpy.ones(shape=(chans,3*pols), dtype='float32')
		for p in range(pols):
			if DOBAND:
				gain[:,p].real = cal_aips[i]['real_%d' % (p+1)][iif*chans:(iif+1)*chans]
				gain[:,p].imag = cal_aips[i]['imag_%d' % (p+1)][iif*chans:(iif+1)*chans]
			if ifs == 1:
				if DOCALIB:
					gain[:,p] = cal_aips[i]['real%d' % (p+1)] + 1j * cal_aips[i]['imag%d' % (p+1)]
					sbdelay[:,3*p+0] = numpy.angle(gain[:,p])
					sbdelay[:,3*p+1] = cal_aips[i]['delay_%d' % (p+1)] * 1e9 # convert to nano-seconds
					sbdelay[:,3*p+2] = cal_aips[i]['rate_%d' % (p+1)][iif]
				snr[:,p] = cal_aips[i]['weight_%d' % (p+1)]
				sbdsnr[:,3*p+1] = cal_aips[i]['weight_%d' % (p+1)]
			else:
				if DOCALIB:
					gain[:,p] = cal_aips[i]['real%d' % (p+1)][iif] + 1j * cal_aips[i]['imag%d' % (p+1)][iif]
					sbdelay[:,3*p+0] = numpy.angle(gain[:,p])
					sbdelay[:,3*p+1] = cal_aips[i]['delay_%d' % (p+1)][iif] * 1e9 # convert to nano-seconds
					sbdelay[:,3*p+2] = cal_aips[i]['rate_%d' % (p+1)][iif]
				snr[:,p] = cal_aips[i]['weight_%d' % (p+1)][iif]
				sbdsnr[:,3*p+1] = cal_aips[i]['weight_%d' % (p+1)][iif]
			if DOCALIB:
				mbdelay[:,p] = cal_aips[i]['mbdelay%d' % (p+1)] * 1e9 # convert to nano-seconds
		zero_wt = numpy.logical_not(numpy.array(snr, dtype='bool'))
		# AIPS seems to use 3140 to indicate a bad solution
		flag_val = numpy.logical_and(numpy.array(gain.real, dtype='int') == 3140,
			   numpy.array(gain.imag, dtype='int') == 3140)
		flag = numpy.logical_or(zero_wt, flag_val)
		sbdflag = numpy.zeros(shape=(chans,3*pols), dtype='bool')
		cals[time_casa][antenna1][iif]['common'] = {
			   'TIME': time_casa, 'INTERVAL': time_interval,
			   'SNR': snr, 'ANTENNA1': antenna1, 'ANTENNA2': -1,
			   'SPECTRAL_WINDOW_ID': iif, 'FIELD_ID': field_id}
		cals[time_casa][antenna1][iif]['gain'] = {
			   'CPARAM': gain, 'FLAG': flag, 'SNR': snr}
		cals[time_casa][antenna1][iif]['mbdelay'] = {
			   'FPARAM': mbdelay, 'FLAG': flag, 'SNR': snr}
		cals[time_casa][antenna1][iif]['sbdelay'] = {
			   'FPARAM': sbdelay, 'FLAG': sbdflag, 'SNR': sbdsnr}
prgs.done()
uv.zap()

if DOBAND:
	cal_type_list = ['bcal']
else:
	cal_type_list = ['gcal', 'fringecal']

for cal_type in cal_type_list:
	cal_table = cal_base + '.' + cal_type
	if cal_type == 'bcal':
		cal_type = 'B Jones'
		par_type = 'Complex'
		singlechan = False
	elif cal_type == 'gcal':
		cal_type = 'G Jones'
		par_type = 'Complex'
		singlechan = True
	elif cal_type == 'fringecal':
		cal_type = 'Fringe Jones'
		par_type = 'Real'
		singlechan = True
	cb.open(ms_name, addcorr=False, addmodel=False)
	cb.createcaltable(cal_table, par_type, cal_type, singlechan)
	cb.close()

# Make the CAL_MAIN table
defs = {'PARAMERR': -numpy.ones(shape=(1,pols), dtype='float32'),
	'WEIGHT': numpy.ones(shape=(1,pols), dtype='float32')}
sbddefs = {'PARAMERR': -numpy.ones(shape=(1,3*pols), dtype='float32'),
	   'WEIGHT': numpy.ones(shape=(1,3*pols), dtype='float32')}

defs_bad = {'FLAG': numpy.ones(shape=(1,pols), dtype='bool'),
	    'CPARAM': numpy.ones(shape=(1,pols), dtype='complex64'),
	    'FPARAM': numpy.ones(shape=(1,pols), dtype='float32'),
	    'SNR': numpy.zeros(shape=(1,pols), dtype='float32')}
sbddefs_bad = {'FLAG': numpy.ones(shape=(1,3*pols), dtype='bool'),
	       'CPARAM': numpy.ones(shape=(1,3*pols), dtype='complex64'),
	       'FPARAM': numpy.ones(shape=(1,3*pols), dtype='float32'),
	       'SNR': numpy.zeros(shape=(1,3*pols), dtype='float32')}

ms_ant = table(ms_name + '/ANTENNA', ack=False)
ms_num_ants = len(ms_ant)
ms_ant.close()

# Open table and create required number of rows
if DOBAND:
	bcal = table(cal_base + '.bcal', ack=False, readonly=False)
	bcal.addrows(ms_num_ants * len(cals) * ifs)
	tb_cals = [bcal]
else:
	gcal = table(cal_base + '.gcal', ack=False, readonly=False)
	gcal.addrows(ms_num_ants * len(cals) * ifs)
	sbdcal = table(cal_base + '.fringecal', ack=False, readonly=False)
	sbdcal.addrows(ms_num_ants * len(cals) * ifs)
	tb_cals = [gcal, sbdcal]

def update_row(tb, row, *vals):
	for val in vals:
		tb[row] = val
	tb[row].update()

timestamps = cals.keys()
timestamps.sort()
i = 0
prgs = progress(len(tb_cals[0]))
print 'Writing cal tables'
for iif in range(ifs):
	for ts in timestamps:
		for ant in range(ms_num_ants):
			try:
				c_common = cals[ts][ant][iif]['common']
				c_gain = cals[ts][ant][iif]['gain']
				if DOCALIB:
					c_mbdelay = cals[ts][ant][iif]['mbdelay']
					c_sbdelay = cals[ts][ant][iif]['sbdelay']
				defs2 = {}
				sbddefs2 = {}
			except KeyError:
				# If a solution does not exist for this ant/IF
				# create an standard one and mark it as bad.
				# Use another solution as a base.
				dummy_ant = cals[ts].keys()[0]
				dummy_if = cals[ts][dummy_ant].keys()[0]
				c_common = cals[ts][dummy_ant][dummy_if]['common']
				c_gain = cals[ts][dummy_ant][dummy_if]['gain']
				if DOCALIB:
					c_mbdelay = cals[ts][dummy_ant][dummy_if]['mbdelay']
					c_sbdelay = cals[ts][dummy_ant][dummy_if]['sbdelay']
				c_common['ANTENNA1'] = ant
				c_common['SPECTRAL_WINDOW_ID'] = iif
				defs2 = defs_bad
				sbddefs2 = sbddefs_bad
			if DOBAND:
				update_row(tb_cals[0], i, c_common, c_gain, defs, defs2)
			else:
				update_row(tb_cals[0], i, c_common, c_gain, defs, defs2)
				update_row(tb_cals[1], i, c_common, c_sbdelay, sbddefs, sbddefs2)
			i += 1
			prgs.update(i)
prgs.done()
