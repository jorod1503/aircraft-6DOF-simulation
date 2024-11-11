import matplotlib.pyplot as plt
import numpy as np

from numerical_integrators import numerical_integration_methods
from goverment_equations import flat_earth_eom

#------------------------------------------------------------------------
#Parte 1: Inicialización de la simulación
#------------------------------------------------------------------------

#Definir vehículo (esfera)
r_sphere = 0.
#Definir condiciones iniciales
u0_bf_mps = []
v0_bf_mps = []
w0_bf_mps = []
p0_bf_rps = []
q0_bf_rps = []
r0_bf_rps = []
phi0_rad = []
theta0_rad = []
psi0_rad = []
p10_n_m = []
p20_n_m = []
p30_n_m = []

#Asignar condiciones iniciales al arreglo
x0 = np.array([
    u0_bf_mps,
    v0_bf_mps,
    w0_bf_mps,
    p0_bf_rps,
    q0_bf_rps,
    r0_bf_rps,
    phi0_rad,
    theta0_rad,
    psi0_rad,
    p10_n_m,
    p20_n_m,
    p30_n_m
])

#Transponiendo el arreglo
x0 = x0.transpose(); nx0 = x0.size

#Definiendo  condiciones de tiempo
t0_s = 0
tf_s = 10
h_s = 0.01

#---------------------------------------------------------------------------------------
#Parte 2: Aproximación numérica de la solución de las ecuaciones que rigen la dinámica
#---------------------------------------------------------------------------------------

#Pre localizar la solución del arreglo
t_s = np.arange(t0_s, tf_s + h_s, h_s); nt_s = t_s.size
x = np.empty((nx0, nt_s), dtype = float)

#Asignación de la condición inicial, x0 a la solución del arreglo x
x[:,0] = x0; 

#Ejecución del Euler mejorado
t_s, x = numerical_integration_methods.forward_euler(flat_earth_eom.flat_earth_eom, t_s, x, h_s)

#Acciones de post procesamiento de datos

#---------------------------------------------------------------------------------------
#Parte 3: Graficar datos
#---------------------------------------------------------------------------------------

#Crear sub gráficos y definir arreglo de figuras
fig, (ax1, ax2) = plt.subplots(1,2, figsize = (10,6)) #1 fila 2 columnas

#Graficar la línea 1 en el primer subgráfico
ax1.plot(t_s, x[0,:], label = 'Line 1')
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('Line 1 Value')
ax1.set_tittle('Line 1')
ax1.grid(True)

#Graficar la línea 2 en el primer subgráfico
ax1.plot(t_s, x[:,1], label = 'Line 2')
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('Line 2 Value')
ax1.set_tittle('Line 2')
ax1.grid(True)

plt.show()
