from definicoes import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def system_equations(variables, u, T):
    phi, phib, M = variables
    return [dOmegaphi(phi, phib, u, T, M), 
            dOmegaphib(phi, phib, u, T, M), 
            dOmegaM(phi, phib, u, T, M)]

def solve_system(u, T, initial_guess):
    try:
        solution = fsolve(system_equations, initial_guess, args=(u, T))
        return solution
    except Exception as e:
        print(f"Error during solution for T={T}: {e}")
        return None

def solve_for_T_range(u, T_values, initial_guess):
    results = {}
    current_guess = initial_guess

    for T in T_values:
        solution = solve_system(u, T, current_guess)
        if solution is not None:
            results[T] = solution
            current_guess = solution  # Update the initial guess
        else:
            results[T] = None #store None when solution fails.
            current_guess = initial_guess #reset initial guess
    return results

# Example usage:
u_value = 0.2  # Fixed u value
T_range = np.linspace(0.01, 0.3, 100)  # Range of T values
initial_guess = [0.01, 0.01, 0.3]  # Initial guess

solutions = solve_for_T_range(u_value, T_range, initial_guess)

for T, solution in solutions.items():
    if solution is not None:
        phi, phib, M = solution
        print(f"T = {T}: phi = {phi}, phib = {phib}, M = {M}")
    else:
        print(f"T = {T}: Solution not found.")



# Plot the results
T_values = []
phi_values = []
phib_values = []
M_values = []

for T, solution in solutions.items():
    if solution is not None:
        phi, phib, M = solution
        T_values.append(T)
        phi_values.append(phi)
        phib_values.append(phib)
        M_values.append(M)

# Plotting
plt.figure(figsize=(12, 6))


plt.plot(T_values, phi_values)
plt.xlabel("T")
plt.ylabel("phi")
plt.title("phi vs T")


plt.plot(T_values, phib_values)
plt.xlabel("T")
plt.ylabel("phib")
plt.title("phib vs T")

plt.plot(T_values, M_values)
plt.xlabel("T")
plt.ylabel("M")
plt.title("M vs T")

plt.tight_layout()
plt.show()
