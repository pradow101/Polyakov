from definicoes import *
import scipy as sp
from sympy import *
from constants import *


phi = symbols('phi')
phib = symbols('phib')
M = symbols('M')
u = 0.34
T = symbols('T')


F0 = ((M-m)**2)/(2*G) - (1/(pi2))*(Nf*T)*(sp.integrate.quad(lambda k: (k**2)*(np.log(1 + 3*(phi + phib*np.exp((-1/T)*(np.sqrt(k**2 + M**2) - u)))*(np.exp((-1/T)*(np.sqrt(k**2 + M**2) - u))) + np.exp((-3/T)*(np.sqrt(k**2 + M**2) - u))) + np.log(1 + 3*(phib + phi*np.exp((-1/T)*(np.sqrt(k**2 + M**2) + u)))*(np.exp((-1/T)*(np.sqrt(k**2 + M**2) + u))) + np.exp((-3/T)*(np.sqrt(k**2 + M**2) + u))), 0, np.inf))) - (1/(pi2))*3*Nf*(sp.integrate.quad(lambda k: k**2 * np.sqrt(k**2 + M**2), 0, L)) + (T**4)*(-0.5*(a0 + a1*(To/T) + a2*(To/T)**2 + a3*(To/T)**3)*phib*phi - (0.16667*b3)*(phi**3 + phib**3) + (0.25*b4)*(phib*phi)**2)


np.diff(F0, phi)