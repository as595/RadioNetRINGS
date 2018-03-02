import numpy as np

def gen_model(v, Fs):
    [k, tau] = v
    res = np.exp(1j*(k + 2*np.pi*tau*Fs))
    return res

def fun(v, Fs, data):
    delta = gen_model(v, Fs) - data
    return np.sum(np.ndarray.flatten(np.abs(delta)**2))
                
def jac(v, Fs, data):
    [k, tau] = v
    D0 = data * np.exp(-1j*(k + 2*np.pi*tau*Fs))
    return np.array([-2*np.sum(np.ndarray.flatten(np.imag(D0))),
                     -2*np.sum(np.ndarray.flatten(2*np.pi*Fs*np.imag(D0)))
                 ])


#scipy.optimize.minimize(ls1d.fun, [0,0], args=(ffd.fgrid0[:,8], fdata[:, 8]),
#                        method='Newton-CG', jac = ls1d.jac) 
#pout3 = lsqrs.minimise_wrapper(lsqrs.fun_s3, list(pout2),
#                               args=(ffd.tgrid0, ffd.fgrid0, len(antennas2), data2, weights2, ref_s_ind2),
#                               method='Newton-CG', jac=lsqrs.jac_s3,
#                               options={'scale':my_scales2, 'xtol':1e-12}
#                                    )[0]
