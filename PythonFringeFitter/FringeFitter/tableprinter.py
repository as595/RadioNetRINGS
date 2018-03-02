from casac import casac
import utils

tb = casac.table()

def formatTable(tname):
    tb.open(tname)
    times = [utils.render_time(t) for t in tb.getcol('TIME')]
    ants = tb.getcol('ANTENNA1')
    swids = tb.getcol('SPECTRAL_WINDOW_ID')
    flags = list(tb.getcol('FLAG').squeeze().transpose().all(axis=1))
    params = list(tb.getcol("FPARAM").squeeze())
    fs = "{:20} {:3} {:2} " + 6*"{:>15.8g} "
    s = "\n".join(fs.format(*r[1:]) for r in zip(flags, times, ants, swids, *params) if not r[0])
    tb.close()
    return s
