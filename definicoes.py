import numpy as np
import scipy as sp
from constants import *


def Ep(p,M):
    return np.sqrt(p**2 + M**2)

def U(phi, phib, T):
    return (T**4)*(-0.5*(a0 + a1*(To/T) + a2*(To/T)**2 + a3*(To/T)**3)*phib*phi - (0.16667*b3)*(phi**3 + phib**3) + (0.25*b4)*(phib*phi)**2)

def zminus(phi, phib, u, T, p, M):
    return np.log(1 + 3*(phi + phib*np.exp((-1/T)*(Ep(p,M) - u)))*(np.exp((-1/T)*(Ep(p,M) - u))) + np.exp((-3/T)*(Ep(p,M) - u)))

def zplus(phi, phib, u, T, p, M):
    return np.log(1 + 3*(phib + phi*np.exp((-1/T)*(Ep(p,M) + u)))*(np.exp((-1/T)*(Ep(p,M) + u))) + np.exp((-3/T)*(Ep(p,M) + u)))

def potencial(phi, phib, u, T, M): ###Não vou precisar mais do potencial, eu acho. Como estou definindo todas as derivadas, este fica implícito nelas.
    a = sp.integrate.quad(lambda k: (k**2)*(zminus(phi, phib, u, T, k, M) + zplus(phi, phib, u, T, k, M)), 0, np.inf)[0]
    b = sp.integrate.quad(lambda k: (k**2)*Ep(k,M), 0, L)[0]
    return ((M-m)**2)/(2*G) - (1/(pi2))*(Nf*T*a + 3*Nf*b) + U(phi, phib, T)

##DERIVADAS

def dEp(p,M):
    return M/Ep(p,M)

def dUphi(phi,phib,T):
    return (T**4)*(-0.5*(a0 + a1*(To/T) + a2*(To/T)**2 + a3*(To/T)**3)*phib -0.5*b3*phi**2 + 0.5*b4*(phib**2)*phi)

def dUphib(phi,phib,T):
    return (T**4)*(-0.5*(a0 + a1*(To/T) + a2*(To/T)**2 + a3*(To/T)**3)*phi -0.5*b3*phib**2 + 0.5*b4*phib*(phi**2))

def dZplusphi(phi, phib, u, T, p, M):
    a = 3*np.exp(-2*(Ep(p,M) + u)/T)
    b = 1 + 3*(phib + phi*np.exp(-(Ep(p,M) + u)/T))*(np.exp(-(Ep(p,M) + u)/T)) + np.exp(-3*(Ep(p,M) + u)/T)
    return a/b

def dZplusphib(phi, phib, u, T, p, M):
    a = 3*np.exp(-(Ep(p,M) + u)/T)
    b = 1 + 3*(phib + phi*np.exp(-(Ep(p,M) + u)/T))*(np.exp(-(Ep(p,M) + u)/T)) + np.exp(-3*(Ep(p,M) + u)/T)
    return a/b

def dZminusphi(phi, phib, u, T, p, M):
    a = 3*np.exp(-(Ep(p,M) - u)/T)
    b = 1 + 3*(phi + phib*np.exp(-(Ep(p,M) - u)/T))*(np.exp(-(Ep(p,M) - u)/T)) + np.exp(-3*(Ep(p,M) - u)/T)
    return a/b

def dZminusphib(phi, phib, u, T, p, M):
    a = 3*np.exp(-2*(Ep(p,M) - u)/T)
    b = 1 + 3*(phi + phib*np.exp(-(Ep(p,M) - u)/T))*(np.exp(-(Ep(p,M) - u)/T)) + np.exp(-3*(Ep(p,M) - u)/T)
    return a/b

def dZplusM(phi, phib, u, T, p, M):
    a = -3*M*((phib + 2*phi*np.exp(-(Ep(p,M) + u)/T))*np.exp(-(Ep(p,M) + u)/T) + np.exp(-3*(Ep(p,M)+u)/T))
    b = Ep(p,M)*T*(1 + 3*(phib + phi*np.exp(-(Ep(p,M) + u)/T))*(np.exp(-(Ep(p,M) + u)/T)) + np.exp(-3*(Ep(p,M) + u)/T))
    return a/b

def dZminusM(phi, phib, u, T, p, M):
    a = -3*M*((phi + 2*phib*np.exp(-(Ep(p,M) - u)/T))*np.exp(-(Ep(p,M) - u)/T) + np.exp(-3*(Ep(p,M)-u)/T))
    b = Ep(p,M)*T*(1 + 3*(phi + phib*np.exp(-(Ep(p,M) - u)/T))*(np.exp(-(Ep(p,M) - u)/T)) + np.exp(-3*(Ep(p,M) - u)/T))
    return a/b

def dOmegaphi(phi,phib,u,T,M):
    a = dUphi(phi,phib,T)
    b = sp.integrate.quad(lambda k: (k**2)*(dZminusphi(phi, phib, u, T, k, M) + dZplusphi(phi, phib, u, T, k, M)), 0, np.inf)[0]
    return a - (1/pi2)*Nf*T*b

def dOmegaphib(phi,phib,u,T,M):
    a = dUphib(phi,phib,T)
    b = sp.integrate.quad(lambda k: (k**2)*(dZminusphib(phi, phib, u, T, k, M) + dZplusphib(phi, phib, u, T, k, M)), 0, np.inf)[0]
    return a - (1/pi2)*Nf*T*b

def dOmegaM(phi,phib,u,T,M):
    a = (M-m)/G
    b = sp.integrate.quad(lambda k: (k**2)*(dZminusM(phi, phib, u, T, k, M) + dZplusM(phi, phib, u, T, k, M)), 0, np.inf)[0]
    c = sp.integrate.quad(lambda k: (k**2)*dEp(k,M), 0, L)[0]
    return a - (1/pi2)*(Nf*T*b - 3*Nf*c)

def sistema(u,T,chuteinit):
    z = np.empty(3)
    phi = z[0]
    phib = z[1]
    M = z[2]
    def eqns(z):
        z[0] = dOmegaphi(phi,phib,u,T,M)
        z[1] = dOmegaphib(phi,phib,u,T,M)
        z[2] = dOmegaM(phi,phib,u,T,M)
        return z
    result = sp.optimize.fsolve(eqns, chuteinit)
    return result
