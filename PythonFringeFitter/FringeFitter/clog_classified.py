import sys
import numpy as np, itertools
import ffd, fringer, utils
from taskinit import casalog
import make_table


class UnhandledCase(Exception):
    def __init__(self, s):
        Exception.__init__(s)

class FringeFitter(object):
    def __init__(self, msname, fj_name, ref_antenna_name=None, scans=None,
                 threshold=None, snr_threshold=None, antennas=None,
                 spectral_windows=None,
                 solint=None, solmin=0, solsub=1, dofloat=True):
        self.msname = msname
        if scans is None:
            raise UnhandledCase("No scans selected!")
        else:
            self.scans = scans
        self.threshold_method = None
        self.solint = solint
        self.solmin = solmin
        self.solsub = solsub
        self.dofloat = dofloat
        if snr_threshold is None:
            self.snr_threshold = 5.0
        else:
            self.snr_threshold = snr_threshold
            self.threshhold_method = 'snr'
        if threshold is None:
            self.threshold = 2000
        else:
            self.threshold = threshold
            # SNR takes precedence if both thresholds are set.
            if self.threshold_method is None:
                self.threshhold_method = 'raw'
        if self.threshold_method is None:
            self.threshold_method = 'snr'
        #
        if spectral_windows is None:
            self.spectral_windows = utils.get_spectral_windows(self.msname)
        else:
            self.spectral_windows = spectral_windows
        self.timeqs = self.make_time_qs_from_scans(scans)
        self.antenna_map = utils.get_antenna_map(self.msname)
        self.ism = utils.invert_map(self.antenna_map)
        if antennas is None:
            # FIXME! Should do this per scan!
            self.antennas2 = sorted(
                ffd.actual_antennas(self.msname, self.timeqs[0]).keys())
        else:
            self.antennas2 = antennas
        self.antenna_names = [self.ism[s] for s in self.antennas2]
        if ref_antenna_name is None:
            self.ref_antenna = self.antennas2[0]
            self.ref_antenna_name = self.ism[self.ref_antenna]
            casalog.post("No reference antenna selected; using {} ({})"
                        "".format(self.ref_antenna_name, self.ref_antenna))
        else:
            self.ref_antenna_name = ref_antenna_name
            self.ref_antenna = self.antenna_map[self.ref_antenna_name]
        pol_ids = utils.get_polarization_ids(msname)
        if len(pol_ids) > 1:
            raise UnhandledCase("Non-unique polarisation id")
        else:
            self.pol_id = pol_ids[0]
        n_pols = utils.get_n_polarizations(msname, self.pol_id)
        if n_pols == 4:
            self.polinds = [0, 3]
        elif n_pols == 2:
            self.polinds = [0]
        else:
            raise UnhandledCase("Can't do {} polarizations".format(n_pols))
        self.bad_antennas = set()
        self.fj_name = fj_name
        make_table.make_table(self.msname, self.fj_name)
        self.rowcount = 0
        #
        shape = (len(self.polinds), len(self.antennas2))
        # Looks like we don't need 'F'? Leave it in for now.
        # FIXME: find out one way or another.
        self.delays = np.zeros(shape, np.float, order='F')
        self.phases = np.zeros(shape, np.float, order='F')
        self.rates = np.zeros(shape, np.float, order='F')
        self.flags = np.zeros(shape, np.bool, order='F')
        self.sigs = []
    def make_FFD(self, timeq, swid, polind, solint=300):
        timeq2 = ffd.actual_timerangeq(self.msname, timeq)
        casalog.post("Processing spectral window {} for data selection {}"
                     "".format(swid, timeq), "DEBUG")
        anffd = ffd.FFData.make_FFD(self.msname, self.antennas2, swid,
                                    self.pol_id, polind, timeq2,
                                    datacol="CORRECTED_DATA", solint=solint)
        return anffd
    def run(self):
        for scan in self.scans:
            if self.solint is None:
                solint = utils.get_scan_length(self.msname, scan)
            else:
                solint = self.solint
            scanq = self.make_time_q_from_scan(scan)
            solintqs = ffd.divide_up_timerange(self.msname, scanq, solint,
                                               self.solmin, self.solsub, self.dofloat)
            for timeq, swid in itertools.product(solintqs, self.spectral_windows):
                for pi, pol_ind in enumerate(self.polinds):
                    self.anffd = self.make_FFD(timeq, swid, pol_ind, solint)
                    t = fringer.fit_fringe_ffd(
                        self.anffd, self.ref_antenna, self.antennas2,
                        threshold=self.threshold,
                        threshold_method=self.threshold_method,
                        snr_threshold=self.snr_threshold)
                    flags, dels, phs, rs, sig = t
                    self.delays[pi, :] = dels
                    self.phases[pi, :] = phs
                    self.rates[pi, :] = rs
                    self.flags[pi, :] = flags
                    self.sigs.append(sig)
                    # FIXME: Both here and in fringer.fit_fringe_ffd we
                    # have explicit code to handle the choice of methods.
                    # This smells bad.
                    if self.threshold_method == 'snr':
                        bad_antennas = self.anffd.get_p_antennas_below_snr_threshold(self.ref_antenna, self.snr_threshold)
                    elif self.threshold_method == 'raw':
                        bad_antennas = self.anffd.get_p_antennas_below_raw_threshold(self.ref_antenna, self.snr_threshold)
                    else:
                        bad_antennas = []
                    self.bad_antennas |= set(bad_antennas)
                self.write_table(timeq, self.flags, swid, self.phases, self.delays, self.rates)
    def getBadAntennas(self):
        return self.bad_antennas()
    def make_time_q_from_scan(self, scan):
        return "SCAN_NUMBER={}".format(scan)
    def make_time_qs_from_scans(self, scans):
        return [self.make_time_q_from_scan(s) for s in scans]
    def write_table(self, timeq, flags, swid, phases, delays, rates):
        """Single-band version."""
        timeq2 = timeq + " AND (ANTENNA1 = {} OR ANTENNA2 = {})".format(self.ref_antenna, self.ref_antenna)
        obsid, field, scan = [ffd.distinct_thing(self.msname, timeq2, col)
                              for col in ['OBSERVATION_ID', 'FIELD_ID', 'SCAN_NUMBER']]
        anffd = self.anffd
        darr = -delays*1e9
        pharr = -phases # Now radians!
        rarr = -rates
        # We don't write the rates!
        for i, s in enumerate(self.antennas2):
            antenna = s
            casalog.post("Writing row {} for antenna {}".format(i, antenna))
            assert anffd.get_antenna_index(s) == i
            # time = anffd.get_ref_time()
            time = anffd.get_min_time()
            interval = anffd.get_interval()
            ref_freq = anffd.ref_freq
            ph = pharr + 0.5*interval*2*np.pi*rarr*ref_freq
            if (flags[0, i] or flags[1, i]):
                casalog.post("Skipping station {} swid {}"
                             "".format(s, swid))
                continue
            casalog.post("Writing station {} swid {}"
                         "".format(s, swid))
            # tables require (6,1) arrays.
            v = [ph[0, i], darr[0, i], rarr[0, i],
                 ph[1, i], darr[1, i], rarr[1, i]],
            param = np.array(v, dtype='float32', ndmin=2).transpose()
            make_table.add_row(self.fj_name, self.rowcount, time + interval/2, interval, antenna,
                               field, scan, obsid, swid, param)
            self.rowcount += 1
        return None

class MultiBandFringeFitter(FringeFitter):
    def __init__(self, msname, ctname, ref_antenna_name=None, scans=None,
                 snr_threshold=None, threshold=None, antennas=None, pad=2,
                 solint=None, solmin=0, solsub=1, dofloat=True):
        FringeFitter.__init__(self, msname, ctname, ref_antenna_name=ref_antenna_name, scans=scans,
                              snr_threshold=snr_threshold, threshold=threshold, antennas=antennas,
                              solint=solint, solmin=solmin, solsub=solsub, dofloat=dofloat)
        # Just for multiband case:
        self.ctname = ctname
        self.reffreqs = utils.get_min_freqs(msname)
        self.minfreq = utils.get_min_freqs(msname)[0]
        self.pad = pad
        # ctname
        make_table.make_table(self.msname, self.ctname)
        self.rowcount = 0 # In the table.
    def run(self):
        shape = (2, len(self.antennas2))
        flags = np.zeros(shape, np.bool)
        delays = np.zeros(shape, np.float)
        phases = np.zeros(shape, np.float)
        rates = np.zeros(shape, np.float)
        sigs = []
        ref_freqs = utils.get_reference_freqs(self.msname)
        ref_freq_diffs = ref_freqs - ref_freqs[0]
        for scan in self.scans:
            if self.solint is None:
                solint = utils.get_scan_length(self.msname, scan)
            else:
                solint = self.solint
            scanq = self.make_time_q_from_scan(scan)
            solintqs = ffd.divide_up_timerange(self.msname, scanq, solint,
                                               self.solmin, self.solsub, self.dofloat)
            for timeq in solintqs:
                timeq2 = ffd.actual_timerangeq(self.msname, timeq)
                for pi, polind in enumerate(self.polinds):
                    casalog.post("Getting data")
                    anffd = ffd.FFData.make_FFD_multiband(self.msname, self.antennas2, self.pol_id,
                                                          polind, timeq2, datacol="CORRECTED_DATA",
                                                          solint=self.solint)
                    self.anffd = anffd # so developers can poke it offline.
                    casalog.post("Fitting fringes")
                    self.t = fringer.fit_fringe_ffd(anffd, self.ref_antenna, self.antennas2, pad=self.pad)
                    fs, dels, phs, rs, sig = self.t
                    flags[pi, :] = fs
                    delays[pi, :] = dels
                    phases[pi, :] = phs
                    rates[pi, :] = rs
                    sigs.append(sig)
                for swid in self.spectral_windows:
                    diffs = 2*np.pi*((ref_freq_diffs[swid]*delays) % 1.0)
                    ph = phases + diffs
                    # casalog.post("phases {} diffs {}".format(phases, diffs))
                    self.write_table(anffd, timeq2, flags, swid, delays, ph, rates, sigs)
    def write_table(self, anffd, timeq, flags, swid, delays, phases, rates, sigs):
        """Multi-band version."""
        timeq2 = timeq + " AND (ANTENNA1 = {} OR ANTENNA2 = {})".format(self.ref_antenna, self.ref_antenna)
        print >>sys.stderr, timeq2
        obsid, field, scan = [ffd.distinct_thing(self.msname, timeq2, col)
                              for col in ['OBSERVATION_ID', 'FIELD_ID', 'SCAN_NUMBER']]
        darr = -delays*1e9
        pharr = -phases # radians!
        rarr = -rates
        for i, s in enumerate(self.antennas2):
            antenna = s
            assert anffd.get_antenna_index(s) == i
            if (flags[0, i] or flags[1, i]):
                continue
            # time = anffd.get_ref_time()
            time = anffd.get_min_time()
            interval = anffd.get_interval()
            # tables require (6,1) arrays.
            param = np.array([pharr[0, i], darr[0, i], rarr[0, i],
                              pharr[1, i], darr[1, i], rarr[1, i]],
                             dtype='float32', ndmin=2).transpose()
            make_table.add_row(self.ctname, self.rowcount, time, interval, antenna,
                               field, scan, obsid, swid, param)
            self.rowcount += 1
        return None

