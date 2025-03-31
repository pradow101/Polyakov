from definicoes import *
import matplotlib.pyplot as plt


def system_equations(variables, u, T):
    phi, phib, M = variables
    return [dOmegaphi(phi, phib, u, T, M), 
            dOmegaphib(phi, phib, u, T, M), 
            dOmegaM(phi, phib, u, T, M)]

def solve_system(u, T, initial_guess):
    try:
        solution = fsolve(system_equations, initial_guess, args=(u, T))
        return solution
    except Exception as e: #Has to have this try/except to avoid errors in the fsolve function
        print(f"Error during solution for T={T}: {e}")
        return None

def solverTrange(u, T_vals, chuteinit):
    results = []
    chuteatual = chuteinit
    for T in T_vals:
        solution = solve_system(u, T, chuteatual)
        phi, phib, M = solution
        chuteatual = solution
        results.append((T, phi, phib, M))
    return results
    
## Implementar certo as precisões para as soluções e arrumar fsolve para continuar tentando até achar solução correta.


T_vals = np.linspace(0.045, 0.08, 150)
chuteinit = [0.01, 0.01, 0.3]
u = 0.34
solutions = solverTrange(u, T_vals, chuteinit)

phi_vals = []
phib_vals = []
M_vals = []

for T, phi, phib, M in solutions:
    phi_vals.append(phi)
    phib_vals.append(phib)
    M_vals.append(M)

#plt.plot(T_vals, phi_vals,label='phi')
#plt.plot(T_vals, phib_vals,label='phib')
#plt.plot(T_vals, M_vals,label='M')
#plt.legend()


yplot = []
y2plot = []
for i in range(len(T_vals)):
    a = T_vals[i]
    b = phi_vals[i]
    c = phib_vals[i]
    d = M_vals[i]
    yi = dOmegaM(b,c,u,a,d)
    yib = dOmegaphi(b,c,u,a,d)
    yplot.append(yi)
    y2plot.append(yib)

plt.plot(T_vals, yplot,label='dOmegaM')
plt.plot(T_vals, y2plot,label='dOmegaphi')
plt.show()