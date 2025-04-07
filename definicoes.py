import numpy as np
import scipy as sp
from constants import *
from scipy.optimize import fsolve


def Ep(p,M):
    return np.sqrt(p**2 + M**2)

def U(phi, phib, T):
    return (T**4)*(-0.5*(a0 + a1*(To/T) + a2*(To/T)**2 + a3*(To/T)**3)*phib*phi - (0.16667*b3)*(phi**3 + phib**3) + (0.25*b4)*(phib*phi)**2)

def zminus(phi, phib, u, T, p, M):
    return np.log(1 + 3*(phi + phib*np.exp((-1/T)*(Ep(p,M) - u)))*(np.exp((-1/T)*(Ep(p,M) - u))) + np.exp((-3/T)*(Ep(p,M) - u)))

def zplus(phi, phib, u, T, p, M):
    return np.log(1 + 3*(phib + phi*np.exp((-1/T)*(Ep(p,M) + u)))*(np.exp((-1/T)*(Ep(p,M) + u))) + np.exp((-3/T)*(Ep(p,M) + u)))

def potencial(phi, phib, u, T, M):
    a = sp.integrate.quad(lambda k: (k**2)*(zminus(phi, phib, u, T, k, M) + zplus(phi, phib, u, T, k, M)), 0, np.inf)
    b = Ep(L,M)*(0.0689736 + 0.081375*(M**2) - (0.125*(M**3)*np.arcsin(L/M))/(np.sqrt(1 + (L**2)/(M**2))))
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
    a = -3*M*(np.exp(-3*(Ep(p,M)+u)/T) + phi*np.exp(-2*(Ep(p,M)+u)/T) + np.exp(-(Ep(p,M)+u)/T)*(phib + phi*np.exp(-(Ep(p,M)+u)/T)))/(Ep(p,M)*T)
    b = 1 + 3*(phib + phi*np.exp(-(Ep(p,M) + u)/T))*(np.exp(-(Ep(p,M) + u)/T)) + np.exp(-3*(Ep(p,M) + u)/T)
    return a/b

def dZminusM(phi, phib, u, T, p, M):
    a = -3*M*(np.exp(-3*(Ep(p,M)-u)/T) + phib*np.exp(-2*(Ep(p,M)-u)/T) + np.exp(-(Ep(p,M)-u)/T)*(phi + phib*np.exp(-(Ep(p,M)-u)/T)))/(Ep(p,M)*T)
    b = 1 + 3*(phi + phib*np.exp(-(Ep(p,M) - u)/T))*(np.exp(-(Ep(p,M) - u)/T)) + np.exp(-3*(Ep(p,M) - u)/T)
    return a/b

def dOmegaphi(phi,phib,u,T,M):
    a = dUphi(phi,phib,T)
    b = sp.integrate.quad(lambda k: (k**2)*(dZminusphi(phi, phib, u, T, k, M) + dZplusphi(phi, phib, u, T, k, M)), 0, np.inf)[0]
    return a - (Nf*T*b)/pi2

def dOmegaphib(phi,phib,u,T,M):
    a = dUphib(phi,phib,T)
    b = sp.integrate.quad(lambda k: (k**2)*(dZminusphib(phi, phib, u, T, k, M) + dZplusphib(phi, phib, u, T, k, M)), 0, np.inf)[0]
    return a - (Nf*T*b)/pi2

def dOmegaM(phi,phib,u,T,M):
    a = (M-m)/(2*G)
    b = sp.integrate.quad(lambda k: (k**2)*(dZminusM(phi, phib, u, T, k, M) + dZplusM(phi, phib, u, T, k, M)), 0, np.inf)[0]
    c = sp.integrate.quad(lambda k: (k**2)*dEp(k,M), 0, L)[0]
    return a - (Nf*T*b + 3*Nf*c)/pi2

def d2OmegaM(phi,phib,u,T,M):               #
    a = phi + phib*np.exp(-Ep(u,M)/T)       #
    b = phi*np.exp(-Ep(u,M)/T) + phib       #
    return -3*(a*b)/(Ep(u,M)*T)             #
                                            # só para ter algo definido, não são essas funções de fato.
def d3OmegaM(phi,phib,u,T,M):               #       
    a = phi + phib*np.exp(-Ep(u,M)/T)       #
    b = phi*np.exp(-Ep(u,M)/T) + phib       #
    return 3*(a*b)/(Ep(u,M)*T)**2           #

def CEPeqns(phi,phib,u,T,M):
    a = dOmegaphi(phi,phib,u,T,M)
    b = dOmegaphib(phi,phib,u,T,M)
    c = dOmegaM(phi,phib,u,T,M)
    d = d2OmegaM(phi,phib,u,T,M)
    e = d3OmegaM(phi,phib,u,T,M)
    return [a,b,c,d,e]

