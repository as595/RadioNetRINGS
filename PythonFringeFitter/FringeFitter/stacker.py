tb.open("EY015D.ms")
# Scans are geb0rked by round-tripping to IDI, so we use time ranges
t = tb.query('TIME > MJD("2012-10-24-21:26:20") and TIME < MJD("2012-10-24-21:27:40")')


ants = set(t.getcol('ANTENNA1'))

ref = 0
k = 2

q = 'TIME > MJD("2012-10-24-21:26:20") and TIME < MJD("2012-10-24-21:27:40") AND ANTENNA1={} AND ANTENNA2={} AND DATA_DESC_ID=0'

data = sum_stack2(tb, q, ref, k)
f = fft.fftshift(fft.fft2(data))

imax, jmax = fringer.index_of_max(abs(f))

rms = fringer.rms_window(m, (imax, jmax), (3,3))

data2 = fringer.get_stackable_baseline(tb, q, ref, k)

phi_0_4 = fringer.get_stackable_baseline(tb, q, 0, 4, startrow=8) # now unitized for your pleasure
