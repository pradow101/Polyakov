from solver import *
import matplotlib.pyplot as plt


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

# plt.plot(T_vals, np.log10(np.abs(yplot)),label='dOmegaM')
# plt.plot(T_vals, np.log10(np.abs(y2plot)),label='dOmegaphi')


# plt.plot(T_vals, phi_vals,label='phi')
# plt.plot(T_vals, phib_vals,label='phib')
plt.plot(T_vals, M_vals,label='M')
plt.scatter(T_vals, phi_vals,s=0.5)
plt.legend()
plt.show()
