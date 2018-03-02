def scale_params(params, scales):
    # if we have an array we can't use its type as a constructor;
    # for lists, tuple and other iterables we can.
    t0 = type(params)
    t = {np.ndarray: np.array}.get(t0, t0)
    params2 = list(params)[:]
    for i in range(3):
        params2[i::3] = [x*scales[i] for x in params2[i::3]]
    return t(params2)

def get_phases(params): return params[0::3]
def get_rates(params):  return params[1::3]
def get_delays(params): return params[2::3]

def remove_ref(v, ref_ant):
    n = len(v)/3
    v_dash = v[:] # make a copy to mutate
    tau_ref = v_dash.pop(3*ref_ant+2)
    r_ref = v_dash.pop(3*ref_ant+1)
    psi_ref = v_dash.pop(3*ref_ant)
    for i in range(n-1):
        v_dash[3*i] -= psi_ref
        v_dash[3*i+1] -= r_ref
        v_dash[3*i+2] -= tau_ref
    return v_dash, (psi_ref, r_ref, tau_ref)

def restore_ref(params0, refs, ref_ant):
    n = len(params0)/3+1
    # Make a copy since it might be array rather than list
    sol = params0[:] 
    (psi_ref, r_ref, tau_ref) = refs
    sol.insert(3*ref_ant, 0.0)
    sol.insert(3*ref_ant+1, 0.0)
    sol.insert(3*ref_ant+2, 0.0)
    for i in range(n):
        sol[3*i+0] += psi_ref
        sol[3*i+1] += r_ref
        sol[3*i+2] += tau_ref
    return sol

def remove_non_ref(v, ant):
    v_dash = v[:] # make a copy to mutate
    tau = v_dash.pop(3*ant+2)
    r = v_dash.pop(3*ant+1)
    psi = v_dash.pop(3*ant)
    return v_dash, (psi, r, tau)

def restore_non_ref(params0, stored_params, ant):
    # We put values for phase, delay and rate for a given antenna
    # contiguously in an array.
    (psi, r, tau) = stored_params
    # Make a copy since it might be array rather than list
    # also out of politeness
    params = params0[:] 
    params.insert(3*ant, psi)
    params.insert(3*ant+1, r)
    params.insert(3*ant+2, tau)
    return params

def get_antenna_parameters(params, i):
    return tuple(params[3*i:3*(i+1)])

def add_ref_parameters(params0, rp, ref_ant):
    params = params0[:]
    n = len(params)/3
    (psi_ref, r_ref, tau_ref) = rp
    for i in range(n):
        if i == ref_ant:
            continue
        else:
            params[3*i+0] += psi_ref
            params[3*i+1] += r_ref
            params[3*i+2] += tau_ref
    return params

def remove_antennas(v, ants0):
    ants = reversed(sorted(ants0))
    vals = []
    for a in ants:
        v, p = remove_non_ref(v, a)
        # we list the results backwards because we reversed the list.
        vals.insert(0, p)
    return v, vals

def restore_antennas(params, frozen_params, ants0):
    # params is 3*number of initial array
    n_ants = len(ants0) + len(params)/3
    flags = [False for i in range(n_ants)]
    ants = sorted(ants0)
    for a, p in zip(ants, frozen_params):
        flags[a] = True
        params = restore_non_ref(params, p, a)
    return flags, params

