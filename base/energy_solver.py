"""Bounded coupled temperature solve used when scalar updates stall.

The temperature-only fallback fixes L; the joint rescue also solves inverse L.
Sensible heat is an additional unknown, so convective gustiness and resistances
remain deterministic during finite-difference Jacobian evaluation.
"""
import numpy as np
from scipy.optimize import least_squares

from .eb import CP, air_density


TEMPERATURES = {'soilsunlit':'Tss','soilshaded':'Tsh',
                'streetsunlit':'Tss','streetshaded':'Tsh',
                'leafsunlit':'Tls','leafshaded':'Tlh',
                'roofsunlit':'Troof_sunlit','roofshaded':'Troof_shaded',
                'wallsunlit':'Twall_sunlit','wallshaded':'Twall_shaded'}


def solve_temperatures(state, evaluate, residuals, max_iter):
    """Mutate state to a bounded trust-region solution; return evaluations."""
    names = tuple(residuals())
    keys = tuple(TEMPERATURES[name] for name in names)
    start = np.array([min(max(float(state[key]),240.000001),379.999999) for key in keys]
                     + [float(state.get('scene_sensible_heat_W_m2',0.0))])
    lower = np.r_[np.full(len(keys),240.),-1e6]
    upper = np.r_[np.full(len(keys),380.),1e6]

    def objective(values):
        state.update(zip(keys,values[:-1]))
        state['scene_sensible_heat_W_m2'] = float(values[-1])
        evaluate(rebalance_state=False)
        imbalance = residuals()
        return np.array([imbalance[name] for name in names] +
                        [state['scene_sensible_heat_W_m2']-values[-1]])

    fit = least_squares(objective,start,bounds=(lower,upper),
                        x_scale=np.r_[np.full(len(keys),10.),100.],
                        max_nfev=max(1,int(max_iter)),ftol=1e-11,xtol=1e-11,gtol=1e-9)
    objective(fit.x)
    state['temperature_solver'] = 'bounded_trust_region_fallback'
    state['temperature_solver_success'] = bool(fit.success)
    state['temperature_solver_max_residual_W_m2'] = float(np.max(np.abs(fit.fun)))
    return int(fit.nfev)


def solve_coupled(state, meteo, evaluate, residuals, max_iter):
    """Joint T, gustiness H and inverse-L rescue for branch-changing periods."""
    names = tuple(residuals())
    keys = tuple(TEMPERATURES[name] for name in names)
    inverse = 1.0/state['L'] if np.isfinite(state['L']) and state['L'] else 0.0
    start = np.r_[[min(max(float(state[key]),240.000001),379.999999) for key in keys],
                  state.get('scene_sensible_heat_W_m2',0.0),inverse]
    lower = np.r_[np.full(len(keys),240.),-1e6,-np.inf]
    upper = np.r_[np.full(len(keys),380.),1e6,np.inf]
    coefficient = air_density(meteo)*CP*(float(meteo['air_temperature_c'])+273.15)/(0.4*9.81)

    def objective(values):
        state.update(zip(keys,values[:-2]))
        state['scene_sensible_heat_W_m2'] = float(values[-2])
        state['L'] = 1.0/values[-1] if abs(values[-1])>1e-14 else float('inf')
        evaluate(rebalance_state=False)
        h = state['scene_sensible_heat_W_m2']
        imbalance = residuals()
        return np.array([imbalance[name] for name in names] +
                        [h-values[-2],h+coefficient*state['ustar']**3*values[-1]])

    fit = least_squares(objective,start,bounds=(lower,upper),
                        x_scale=np.r_[np.full(len(keys),10.),100.,max(abs(inverse),0.1)],
                        max_nfev=max(1,int(max_iter)),ftol=1e-12,xtol=1e-12,gtol=1e-10)
    objective(fit.x)
    state['temperature_solver']='joint_temperature_gustiness_obukhov_trust_region'
    state['temperature_solver_success']=bool(fit.success)
    state['temperature_solver_max_residual_W_m2']=float(np.max(np.abs(fit.fun)))
    return int(fit.nfev)
