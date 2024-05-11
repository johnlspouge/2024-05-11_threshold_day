#!/usr/bin/env python
"""
Provides epidemic models with ODEs.

https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/
https://numpy.org/doc/stable/reference/generated/numpy.linspace.html # See documentation for parametrization of epidemic time.
"""
from math import isclose
import numpy as np
from scipy.integrate import odeint

from abc import ABC

class EpidemicODEModel(ABC):
    def __init__(self,name2arg):
        self.name2arg = name2arg
        EpidemicODEModel.check(name2arg)
        self.ts = name2arg['ts']
        super().__init__()
        # Derived classes must define the following attributes.
        #   self.variables = None # abbreviations for variables
        #   self.y0 = y0
    # Returns model name.
    def name(self):
        return type(self).__name__
    # Returns total of initial values.
    def n0(self):
        return sum(self.y0)
    # Returns derivatives of variables.
    @staticmethod
    def derivative(y, t, args):
        pass
    # Returns solutions.
    def solve(self):
        # Normalizes the initial values.
        self.y0 = self.name2arg['y0']
        y0 = [y/self.n0() for y in self.y0]
        # Integrates the SIR equations over the time grid, t.
        self.ys = odeint(self.derivative, y0, self.ts, args=tuple(self.name2arg['params']))
        solutions = {'t':list(self.ts)}
        for i,variable in enumerate(self.variables):
            solutions[variable] = list(self.ys[:,i])
        return solutions
    def assert_solve(self, rel_tol=1.0e-06, abs_tol=1.0e-06):
        for pts in self.ys:
            assert len(self.ts) == len(pts)
        assert isclose(self.n0,self.n,rel_tol,abs_tol)
    @staticmethod
    # Checks shared keys.    
    def check(name2arg):
        # times
        assert isinstance(name2arg['ts'], list) 
        assert len(name2arg['ts']) == 3
        for t in name2arg['ts']:
            assert isinstance(t, (float,int))
        start,stop,num = name2arg['ts']
        num = int(num)
        assert start < stop and 1 < num
        name2arg['ts'] = np.linspace(start, stop, num)
        # non-pharmaceutical intervention
        if name2arg['nu'] is None:
            name2arg['nu'] = 1.0
        elif isinstance(name2arg['nu'], (int, float)):
            assert 0.0 < name2arg['nu'] <= 1.0
            name2arg['nu'] = float(name2arg['nu'])
        else:
            raise ValueError('nu is None or isinstance(nu, (int, float)) and 0.0 < nu <= 1.0')
        name2arg['params'] = [name2arg['nu']]
class SIR(EpidemicODEModel):
    def __init__(self,name2arg):
        super().__init__(name2arg)
        SIR.check(name2arg)
        self.variables = ['s','i','r']
    @staticmethod
    def derivative(y, t, nu, beta, gamma, sigma):
        # Gets variables.
        s, i, r = y
        # ODEs
        ds_dt = -s * nu * beta * i             
        di_dt =  s * nu * beta * i - gamma * i + sigma * r * nu * beta * i 
        dr_dt =                      gamma * i - sigma * r * nu * beta * i
        return ds_dt, di_dt, dr_dt
    @staticmethod
    # Checks keys.    
    def check(name2arg, permit_0_infection=False):
        # initial variables
        name2arg['y0'] = []
        # susceptible
        assert isinstance(name2arg['susceptible'], list) and len(name2arg['susceptible']) == 2
        
        total, beta = name2arg['susceptible']
        assert isinstance(total, (float, int)) and 0 < total
        assert isinstance(beta, (float, int)) and 0.0 < beta

        name2arg['y0'].append(total) # 's'
        name2arg['beta'] = beta
        # infectious
        assert isinstance(name2arg['infectious'], list) and len(name2arg['infectious']) == 2
        
        total, gamma = name2arg['infectious']
        assert isinstance(total, (float, int)) and 0 <= total
        assert permit_0_infection or 0 < total
        assert isinstance(gamma, (float, int)) and 0.0 < gamma

        name2arg['y0'].append(total) # 'i'
        name2arg['gamma'] = gamma
        # recovered
        assert isinstance(name2arg['recovered'], list)
        if len(name2arg['recovered']) == 1:
            total = name2arg['recovered'][0]
            sigma = 0.0
        elif len(name2arg['recovered']) == 2:
            total, sigma = name2arg['recovered']
        assert isinstance(total, (float, int)) and 0.0 <= total
        assert isinstance(sigma, (float, int)) and 0.0 <= sigma

        name2arg['sigma'] = sigma
        name2arg['y0'].append(total) # 'r'
        name2arg['params'].extend([beta,gamma,sigma])
        assert len(name2arg['y0']) == 3

class SEIR(EpidemicODEModel):
    def __init__(self,name2arg):
        super().__init__(name2arg)
        SEIR.check(name2arg)
        self.variables = ['s','e','i','r']
    @staticmethod
    def derivative(y, t, nu, beta, eps, gamma, sigma):
        # Gets variables.
        s, e, i, r = y
        # ODEs
        ds_dt = -s * nu * beta * i
        de_dt =  s * nu * beta * i - eps * e             + sigma * r * nu * beta * i
        di_dt =                      eps * e - gamma * i
        dr_dt =                                gamma * i - sigma * r * nu * beta * i
        return ds_dt, de_dt, di_dt, dr_dt
    @staticmethod
    # Checks keys.    
    def check(name2arg):
        SIR.check(name2arg,permit_0_infection=True)
        assert isinstance(name2arg['exposed'], list) and len(name2arg['exposed']) == 2
        total_e, eps = name2arg['exposed']
            
        assert isinstance(total_e, (float, int)) and 0 <= total_e
        assert isinstance(eps, (float, int)) and 0.0 < eps
        
        total_i = name2arg['infectious'][0]
        assert 0 < total_e+total_i

        name2arg['eps'] = eps
        name2arg['y0'].insert(1,total_e)
        name2arg['params'].insert(2,eps)
        assert len(name2arg['y0']) == 4

class SIR_TwoType(EpidemicODEModel):
    def __init__(self,name2arg):
        super().__init__(name2arg)
        SIR_TwoType.check(name2arg)
        self.variables = ['s','i','r','i1','r1']
    @staticmethod
    def derivative(y, t, nu, beta, gamma, sigma):
        # Gets variables.
        i = [0]*2
        r = [0]*2
        di_dt = [0]*2
        dr_dt = [0]*2
        s, i[0], r[0], i[1], r[1] = y
        # ODEs
        ds_dt = 0.0
        for j in range(2):
            ds_dt +=  -s * nu * beta[j] * i[j]
            di_dt[j] = s * nu * beta[j] * i[j]
            for j0 in range(2):
                  di_dt[j] +=  r[j0] * nu * beta[j] * sigma[j][j0]
            di_dt[j] += -gamma[j] * i[j]
            dr_dt[j] =   gamma[j] * i[j]
            for j0 in range(2):
                  dr_dt[j] += -r[j0] * nu * beta[j] * sigma[j][j0]
        return ds_dt, di_dt[0], dr_dt[0], di_dt[1], dr_dt[1]
    @staticmethod
    # Checks keys.    
    def check(name2arg, permit_0_infection=False):
        # ancillary variables
        total = [0]*2
        beta = [0.0]*2
        gamma = [0.0]*2
        sigma = [[0.0]*2,[0.0]*2]
        # initial variables
        name2arg['y0'] = []
        y0 = []
        # susceptible
        assert isinstance(name2arg['susceptible'], list) and len(name2arg['susceptible']) == 3
        
        total[0], beta[0], beta[1] = name2arg['susceptible']
        assert isinstance(total[0], (float, int)) and 0 < total[0]
        assert isinstance(beta[0], float) and 0.0 < beta[0]
        assert isinstance(beta[1], float) and 0.0 < beta[1]

        name2arg['y0'].append(total[0]) # 's'
        name2arg['beta'] = beta
        # infectious
        assert isinstance(name2arg['infectious'], list) and len(name2arg['infectious']) == 4
        
        total[0], gamma[0], total[1], gamma[1] = name2arg['infectious']

        assert isinstance(total[0], (float, int)) and 0 <= total[0]
        assert permit_0_infection or 0 < total[0]
        assert isinstance(gamma[0], float) and 0.0 < gamma[0]
        assert isinstance(total[1], (float, int)) and 0 <= total[1]
        assert permit_0_infection or 0 < total[1]
        assert isinstance(gamma[1], float) and 0.0 < gamma[1]

        name2arg['y0'].append(total[0]) # 'i'
        y0.append(total[1]) # 'i1'
        # recovered
        assert isinstance(name2arg['recovered'], list)
        if len(name2arg['recovered']) == 6:
            total[0], sigma[0][0], sigma[0][1], total[1], sigma[1][0], sigma[1][1] = name2arg['recovered']
        elif len(name2arg['recovered']) == 2:
            total[0], total[1] = name2arg['recovered']
            sigma[0][0] = sigma[0][1] = sigma[1][0] = sigma[1][1] = 0.0
        else:
            raise ValueError('With a variant, len(argument.recovered) == 6 or len(argument.recovered) == 2 (no cross-infection)')
    
        assert isinstance(total[0], (float, int)) and 0 <= total[0]
        assert isinstance(total[1], (float, int)) and 0 <= total[1]
        assert isinstance(sigma[0][0], (float, int)) and 0.0 <= sigma[0][0]
        assert isinstance(sigma[0][1], (float, int)) and 0.0 <= sigma[0][1]
        assert isinstance(sigma[1][0], (float, int)) and 0.0 <= sigma[1][0]
        assert isinstance(sigma[1][1], (float, int)) and 0.0 <= sigma[1][1]

        name2arg['sigma'] = sigma
        name2arg['y0'].append(total[0]) # 'i'
        y0.append(total[1]) # 'r1'
        name2arg['y0'].extend(y0) # 'i1,r1'
        assert len(name2arg['y0']) == 5
        name2arg['params'].extend([beta,gamma,sigma])

class SEIR_TwoType(EpidemicODEModel):
    def __init__(self,name2arg):
        super().__init__(name2arg)
        SEIR_TwoType.check(name2arg)
        self.variables = ['s','e','i','r','e1','i1','r1']
    @staticmethod
    def derivative(y, t, nu, beta, eps, gamma, sigma):
        # Gets variables.
        e = [0]*2
        i = [0]*2
        r = [0]*2
        de_dt = [0]*2
        di_dt = [0]*2
        dr_dt = [0]*2
        s, e[0], i[0], r[0], e[1], i[1], r[1] = y
        # ODEs
        ds_dt = 0.0
        for j in range(2):
            ds_dt +=  -s * nu * beta[j] * i[j]
            de_dt[j] = s * nu * beta[j] * i[j]
            for j0 in range(2):
                  de_dt[j] +=  r[j0] * nu * beta[j] * sigma[j][j0]
            de_dt[j] += -eps[j] * e[j]
            di_dt[j] =   eps[j] * e[j] 
            di_dt[j] += -gamma[j] * i[j]
            dr_dt[j] =   gamma[j] * i[j]
            for j0 in range(2):
                  dr_dt[j] += -r[j0] * nu * beta[j] * sigma[j][j0]
        return ds_dt, de_dt[0], di_dt[0], dr_dt[0], de_dt[1], di_dt[1], dr_dt[1]
    @staticmethod
    # Checks keys.    
    def check(name2arg):
        SIR_TwoType.check(name2arg, permit_0_infection=True)

        total = [0]*2
        eps = [0.0]*2
        assert isinstance(name2arg['exposed'], list) and len(name2arg['exposed']) == 4
        total[0], eps[0], total[1], eps[1] = name2arg['exposed']
            
        assert isinstance(total[0], (float, int)) and 0 <= total[0]
        assert isinstance(eps[0], (float, int)) and 0.0 < eps[0]
        assert isinstance(total[1], (float, int)) and 0 <= total[1]
        assert isinstance(eps[1], (float, int)) and 0.0 < eps[1]

        name2arg['eps'] = eps
        name2arg['y0'].insert(1, total[0])
        name2arg['y0'].insert(4, total[1])
        assert len(name2arg['y0']) == 7
        name2arg['params'].insert(2, eps)

def epidemic_model_factory(name2arg):
    if not name2arg['is_variant']:
        if name2arg['exposed'] is None:
            return SIR(name2arg)
        else:
            return SEIR(name2arg)
    else:
        if name2arg['exposed'] is None:
            return SIR_TwoType(name2arg)
        else:
            return SEIR_TwoType(name2arg)
