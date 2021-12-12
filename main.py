import pylab

import numpy as np
import matplotlib.pyplot as plt

Nx = 100
Ny = 100
Nt = 100
T = 0

x = np.linspace (0, 1, Nx)
y = np.linspace (0, 2, Ny)
t = np.linspace (0, T, Nt)

hx = x[1] - x[0]
hy = y[1] - y[0]
tau = t[1] - t[0]

gamma_x = 9 * tau / hx**2
gamma_y = 9 * tau / hy**2


u = np.zeros((Nx, Ny, 2 * Nt + 1))


for i in range(0, Nx):
    for j in range(0, Ny):
        u[i, j, 0] = np.sin(np.pi * x[i] / 2) * (4 - y[j]*y[j])


def F_1(i1, i2, j):
    return 0.5 * gamma_y * (u[i1, i2 - 1, j - 1] + u[i1, i2 + 1, j - 1]) + (1 - gamma_y) * u[i1, i2, j - 1] + 0.5 * tau * 2 * x[i1] * y[i2] * (tau * (j + 1) / 2)


def F_2(i1, i2, j):
    return 0.5 * gamma_x * (u[i1 - 1, i2, j - 1] + u[i1 + 1, i2, j - 1]) + (1 - gamma_x) * u[i1, i2, j - 1] + 0.5 * tau * 2 * x[i1] * y[i2] * (tau * (j + 1) / 2)


def progonka_x(i2, j):
    d = np.zeros(Nx)
    sigma = np.zeros(Nx)
    d[1] = 0
    sigma[1] = 0
    A = 0.5 * gamma_x
    B = 1 + gamma_x
    C = 0.5 * gamma_x
    for m in range(1, Nx - 1):
        Fm = (-1) * F_1(m, i2, j)
        d[m + 1] = C / (B - A * d[m])
        sigma[m + 1] = (Fm - A * sigma[m]) / (A * d[m] - B)
        u[Nx - 1, i2, j] = sigma[-1] / (1 - d[-1])
    for m in range(Nx - 1, 0, -1):
        u[m - 1, i2, j] = d[m] * u[m, i2, j] + sigma[m]


def progonka_y(i1, j):
    d = np.zeros(Ny)
    sigma = np.zeros(Ny)
    d[1] = 1
    sigma[1] = 0
    A = 0.5 * gamma_y
    B = 1 + gamma_y
    C = 0.5 * gamma_y
    for m in range(1, Ny - 1):
        Fm = (-1) * F_2(i1, m, j)
        d[m + 1] = C / (B - A * d[m])
        sigma[m + 1] = (Fm - A * sigma[m]) / (A * d[m] - B)
        u[i1, Ny - 1, j] = 0
    for m in range(Ny - 1, 0, -1):
        u[i1, m - 1, j] = d[m] * u[i1, m, j] + sigma[m]

for j in (range (1, 2 * Nt, 2)):
    for i2 in range(1, Ny - 1):
        progonka_x(i2, j)
    for i1 in range(1, Nx - 1):
        progonka_y(i1, j + 1)

fig = plt.figure(figsize= (8, 6))
plt.pcolormesh(y, x, u[:,:,-1], cmap="inferno")
plt.colorbar()
plt.ylabel('x')
plt.xlabel('y')
plt.title ('Значение функции U в момент времени t = ' + str(T))


# times=[0, 0.25, 0.5, 0.75, 1., 1.5]
# fig = plt.figure(figsize = (18, 9))
# for i in range(len(times)):
#     plt.subplot(2, 3, i+1)
#     r = int(times[i] / T * (2 * Nt))
#     plt.pcolormesh(y, x, u[:, :, r], cmap = 'inferno', vmin = 0, vmax = 0.08)
#     plt.colorbar()
#     plt.ylabel('х', fontsize = 8)
#     plt.xlabel('y', fontsize = 8)
#     plt.title('Численное решение U в момент времени t = ' + str(times[i]), fontsize = 8)
# plt.tight_layout()
