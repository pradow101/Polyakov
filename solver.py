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
    results = {}
    chuteatual = chuteinit
    for T in T_vals:
        solucao = solve_system(u, T, chuteatual)
        results[T] = solucao
        chuteatual = solucao
    return results

T_vals = np.linspace(0.01, 0.3, 150)
chuteinit = [0.01, 0.01, 0.3]
u = 0
solutions = solverTrange(u, T_vals, chuteinit)

phi_vals = []
phib_vals = []
M_vals = []

for T, solution in solutions.items():
    phi, phib, M = solution
    phi_vals.append(phi)
    phib_vals.append(phib)
    M_vals.append(M)

for i in range(len(M_vals)):
    M_vals[i] = M_vals[i]/0.32552503427962565

plt.plot(T_vals, phi_vals,label='phi')
#plt.plot(T_vals, phib_vals,label='phib')
plt.plot(T_vals, M_vals,label='M')
plt.legend()
plt.show()