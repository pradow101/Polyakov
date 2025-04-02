from definicoes import *
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# def system_equations(variables, u, T):
#     phi, phib, M = variables
#     return [dOmegaphi(phi, phib, u, T, M), 
#             dOmegaphib(phi, phib, u, T, M), 
#             dOmegaM(phi, phib, u, T, M)]

# def solve_system(u, T, initial_guess):
#     try:
#         solution = fsolve(system_equations, initial_guess, args=(u, T))
#         return solution
#     except Exception as e: #Has to have this try/except to avoid errors in the fsolve function
#         print(f"Error during solution for T={T}: {e}")
#         return None

# def solverTrange(u, T_vals, chuteinit):
#     results = []
#     chuteatual = chuteinit
#     for T in T_vals:
#         solution = solve_system(u, T, chuteatual)
#         phi, phib, M = solution
#         chuteatual = solution
#         results.append((T, phi, phib, M))
#     return results
    
## Implementar certo as precisões para as soluções e arrumar fsolve para continuar tentando até achar solução correta.


# T_vals = np.linspace(0.045, 0.08, 150)
# chuteinit = [0.01, 0.01, 0.3]
# u = 0.34
# solutions = solverTrange(u, T_vals, chuteinit)

# phi_vals = []
# phib_vals = []
# M_vals = []

# for T, phi, phib, M in solutions:
#     phi_vals.append(phi)
#     phib_vals.append(phib)
#     M_vals.append(M)



def Minpotencial(u,T,chuteinit):
    vars[0] = chuteinit[0]
    vars[1] = chuteinit[1]
    vars[2] = chuteinit[2]
    try:
        solution = minimize(lambda vars: potencial(vars,u,T), chuteinit)
        phi, phib, M = solution.x
        return solution
    except Exception as e:
        print(f"Error during minimization for T={T}: {e}")
        return None
    
def solverTrangeMin(u, T_vals, chuteinit):
    results = []
    chuteatual = chuteinit
    for T in T_vals:
        solution = Minpotencial(u, T, chuteatual)
        phi, phib, M = solution
        chuteatual = solution
        results.append((T, phi, phib, M))
    return results

T_vals = np.linspace(0.01,0.3,30)
chuteinit = [0.01, 0.01, 0.3]
u = 0.34
phi_vals_min = []
phib_vals_min = []
M_vals_min = []
solutionsmin = solverTrangeMin(u, T_vals, chuteinit)
for T, phi, phib, M in solutionsmin:
    phi_vals_min.append(phi)
    phib_vals_min.append(phib)
    M_vals_min.append(M)


print(phib_vals_min)
chuteinit = [0.1, 0.1, 0.3] 
u = 0.34
T = 0.15
sol = Minpotencial(u, T, chuteinit)  ##VERY dependent on the initial guess
print(sol.x)
a = dOmegaM(sol.x[0],sol.x[1],u,T,sol.x[2])
print(a)
