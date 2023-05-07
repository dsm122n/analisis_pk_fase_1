import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Define the function to be optimized
def func(params, x, y):
    a, b = params
    return (z - (a * x + b * y))**2

# Generate some noisy data
xdata = np.linspace(0, 1, 10)
ydata = np.linspace(0, 1, 10)
X, Y = np.meshgrid(xdata, ydata)
Z = 3 * X + 2 * Y + 0.5 * np.random.randn(*X.shape)

# Define the initial guess for the parameters
init_params = [1, 1]

# Define a function to compute the sum of the squared residuals
def residuals(params, x, y, z):
    return np.sum(func(params, x, y, z))

# Use the Levenberg-Marquardt algorithm to optimize the parameters
params_history = [init_params]
for i in range(20):
    result = leastsq(residuals, params_history[-1], args=(X, Y, Z), ftol=1e-5)
    params_history.append(result[0])

# Create an animation that shows the optimization progress
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

def update_plot(params):
    ax.clear()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.plot_surface(X, Y, Z, cmap='viridis')
    ax.set_title(f'Residuals: {residuals(params, X, Y, Z):.3f}')
    x_fit = np.linspace(0, 1, 100)
    y_fit = np.linspace(0, 1, 100)
    X_fit, Y_fit = np.meshgrid(x_fit, y_fit)
    Z_fit = params[0] * X_fit + params[1] * Y_fit
    ax.plot_surface(X_fit, Y_fit, Z_fit, alpha=0.5, cmap='coolwarm')

ani = FuncAnimation(fig, update_plot, frames=params_history, repeat=False)

plt.show()



import libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

import data
asc = pd.read_csv("covariables_asc.csv")
only tac == tac1 patients
asc_filtered = asc[asc['tac'] == "tac1"]
asc_filtered = asc_filtered[asc_filtered['px'] == "p01"]
asc_filtered = asc_filtered.drop(['tac', 'px'], axis=1)
asc_filtered = asc_filtered.dropna(axis=1)

def modelo_ajuste_vc_tac1(V, Cl):
    """
    This function is a translation from r code to python.
    """
    k = Cl / V
    t = np.array([0, 15, 30, 60, 90, 120, 180])
    t_cambio = np.array([30, 90])
    R = 20.2362
    conc = t.copy()
    conc[t <= t_cambio[0]] = (R / Cl) * (1 - np.exp(-k * t[t <= t_cambio[0]] / 60.0))
    conc[(t <= t_cambio[1]) & (t > t_cambio[0])] = (R / Cl) * (
        1 - np.exp(-k * (t_cambio[0] / 60.0))
    ) * np.exp(-k * (t[(t <= t_cambio[1]) & (t > t_cambio[0])] - t_cambio[0]) / 60.0) + (
        (R / 6.0) / Cl
    ) * (1 - np.exp(-k * (t[(t <= t_cambio[1]) & (t > t_cambio[0])] - t_cambio[0]) / 60.0))
    conc[t > t_cambio[1]] = (R / Cl) * (
        1 - np.exp(-k * (t_cambio[0] / 60.0))
    ) * np.exp(
        -k * (t[t > t_cambio[1]] - t_cambio[0]) / 60.0
    ) + (
        (R / 6.0) / Cl
    ) * (
        1 - np.exp(-k * ((t_cambio[1] - t_cambio[0]) / 60.0))
    ) * np.exp(
        -k * (t[t > t_cambio[1]] - t_cambio[1]) / 60.0
    ) + (
        0.0 / Cl
    ) * (
        1 - np.exp(-k * (t[t > t_cambio[1]] - t_cambio[1]) / 60.0)
    )
    return conc

def modelo_ajuste_vc_tac2(V, Cl):
    