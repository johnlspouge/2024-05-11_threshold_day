#!/usr/bin/env python
"""
Epidemic Network
"""
from scipy.integrate import odeint
from jls_epidemic_compartment import Compartment 

def check(compartment2total):
    abbrs = ('I', 'R', 'S')
    assert sorted((c.name() for c,t in compartment2total.items())) == abbrs
    for c in compartment2total.values():
        assert c >= 0.0
    
def sir_odes(compartment2total, t_s, R0, mean_recovery_time):
    # The grid for ODE solution is, e.g., t_s = np.linspace(start, stop, num).
    check(compartment2total)
    # Normalize to probability.
    c2t = compartment2total
    total = sum(c2t.values())
    y0 = [x/total for x in [c2t[Compartment.Susceptible], c2t[Compartment.Infected], c2t[Compartment.Recovered]]] 
    # Contact rate, beta, and mean recovery rate, gamma.
    beta = R0/mean_recovery_time
    gamma = 1.0/mean_recovery_time
    # differential equations for the SIR model
    def deriv(y, t, beta, gamma):
        s, i, r = y
        ds_dt = -beta * s * i
        di_dt = beta * s * i - gamma * i
        dr_dt = gamma * i
        return ds_dt, di_dt, dr_dt
    # Integrates the SIR equations over the time grid, t.
    ys = odeint(deriv, y0, t_s, args=(beta, gamma))
    s_s, i_s, r_s = ys.T
    assert len(s_s) == len(i_s) == len(r_s) == len(t_s)
    return  s_s, i_s, r_s, t_s

def main(): 
    pass

if __name__ == "__main__":
    main()
