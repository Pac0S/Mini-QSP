from scipy.integrate import solve_ivp


def rhs(t, y, params):
    Ag, Ac, B = y
    
    ka = params["ka"]
    CL = params["CL"]
    V = params["V"]
    Imax = params["Imax"]
    IC50 = params["IC50"]
    kprod = params["kprod"]
    kdeg = params["kdeg"]
    
    C = Ac/V
    inhibition = Imax * C / (IC50 + C)
    
    dAg_dt = -ka * Ag
    dAc_dt = ka * Ag - CL * C
    
    dB_dt = kprod * (1-inhibition) - kdeg * B
    
    return [dAg_dt, dAc_dt, dB_dt]


def simulate(params, times):
    B0 = params["kprod"] / params["kdeg"]
    dose = params["dose"]
    
    #Initial state is Ag, Ac and B at t = 0
    initial_state = [dose, 0.0, B0]
    
    result = solve_ivp(
        fun=rhs,
        t_span=(times[0], times[-1]),
        y0=initial_state,
        t_eval=times,
        args=(params,),)
    
    if not result.success:
        raise RuntimeError(result.message)
        
    return result
    
def set_parameter(params, param, value):
    if param not in params:
        raise KeyError(f"Paramètre inconnu : {param}")

    modified_params = params.copy()
    modified_params[param] = value
    return modified_params

def set_parameters(in_params, var_params, var_values):
    modified_params = in_params.copy()
    for param, value in zip(var_params, var_values, strict = True):
        if param not in in_params:
            raise KeyError(f"Paramètre inconnu : {param}")
        modified_params[param] = value
    return modified_params

def residuals_ic50(x, reference_params, t_obs, B_obs):
    ic50_candidate = float(x[0])

    candidate_params = reference_params.copy()
    candidate_params["IC50"] = ic50_candidate

    result = simulate(params=candidate_params, times=t_obs)
    if not result.success:
        raise RuntimeError(result.message)

    B_pred = result.y[2]
    return B_pred - B_obs