import math, numpy as np, numpy.ma as ma
import operator
import __casac__

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

tb = __casac__.table.table()


def extract_masked_phases(anmsname, ant2, scan_number):
    correlation = 0
    tb.open(anmsname)
    t = tb.query("DATA_DESC_ID=0 AND ANTENNA1=0 AND ANTENNA2={} AND SCAN_NUMBER={}"
                 "".format(ant2, scan_number))
    data = t.getcol('CORRECTED_DATA')[correlation]
    f = t.getcol("FLAG")[correlation]
    fr = t.getcol("FLAG_ROW")
    times = t.getcol('TIME')

    flags = operator.or_(f, fr)
    # ? ma.masked_invalid()
    phi = ma.masked_array(180/np.pi*np.angle(data), mask=flags)
    return phi, times

def plot3d(anmsname, ant2, scan_number):
    ref_ant = 0
    phi, times = extract_masked_phases(anmsname, ant2, scan_number)
    fig = plt.figure()
    plt.title('Baseline {}-{}'.format(ref_ant, ant2))
    ax = fig.add_subplot(211, projection='3d')
    correlation = 0
    nf, nt = phi.shape
    times0 = times - times[0]
    T0s, Freqs = np.meshgrid(times0, range(nf))
    # p = phi.filled(0)
    p = phi
    # ax.scatter
    # ax.plot_surface(Freqs[16:48:2, :], T0s[16:48:2, :], p[16:48:2, :])
    skip = 2
    ax.plot_wireframe(Freqs[::skip, ::skip], T0s[::skip, ::skip], p[::skip, ::skip])
   # 
    ax2 = fig.add_subplot(212)
    phi_ave = np.mean(phi, axis=1)
    ax2.plot(phi_ave.filled(0))
    

def plot_mean(anmsname, ant2, scan_number):
    phi = extract_masked_phases(anmsname, ant2, scan_number)
    phi_ave = np.mean(phi+180.0, axis=1)-180.0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(phi_ave.filled(0))
