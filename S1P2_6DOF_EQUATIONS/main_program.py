import matplotlib.pyplot as plt
import numpy as np
import math

from numerical_integrators import numerical_integration_methods
from goverment_equations import flat_earth_eom

#------------------------------------------------------------------------
#Parte 1: Inicialización de la simulación
#------------------------------------------------------------------------

#Definir vehículo (esfera)
r_sphere_m = 0.08
m_sphere_kg = 5
J_sphere_kgm2 = 0.4*m_sphere_kg*r_sphere_m**2

amod = {"m_kg": 1, \
        "Jxz_b_kgm2": 0, \
        "Jxx_b_kgm2": J_sphere_kgm2, \
        "Jyy_b_kgm2": J_sphere_kgm2, \
        "Jzz_b_kgm2": J_sphere_kgm2   } # Aircraft model properties 

#Definir las condiciones iniciales (de aeronave en rutina de trim)
u0_bf_mps   = 0
v0_bf_mps   = 0
w0_bf_mps   = 0
p0_bf_rps   = 0
q0_bf_rps   = 0
r0_bf_rps   = 0
phi0_rad    = 90*math.pi/180
theta0_rad  = 90*math.pi/180
psi0_rad    = 0
p10_n_m     = 0
p20_n_m     = 0
p30_n_m     = 0

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

#Transponiendo el arreglo x0 = x0.transpose();
nx0 = x0.size

#Definiendo  condiciones de tiempo
t0_s = 0
tf_s = 10
h_s = 0.005

#---------------------------------------------------------------------------------------
#Parte 2: Aproximación numérica de la solución de las ecuaciones que rigen la dinámica
#---------------------------------------------------------------------------------------

#Pre localizar la solución del arreglo
t_s = np.arange(t0_s, tf_s + h_s, h_s); nt_s = t_s.size
x = np.zeros((nx0, nt_s))#, dtype = float)

#Asignación de la condición inicial, x0 a la solución del arreglo x
x[:, 0] = x0

#Ejecución del Euler mejorado
t_s, x = numerical_integration_methods.forward_euler(flat_earth_eom.flat_earth_eom, t_s, x, h_s, amod)

#Acciones de post procesamiento de datos

#---------------------------------------------------------------------------------------
#Parte 3: Graficar datos
#---------------------------------------------------------------------------------------

#Crear sub gráficos y definir arreglo de figuras
fig, axes = plt.subplots(2,4, figsize = (10,6)) 
fig.set_facecolor('black')

#Velocidad en x u^b_CM/n
axes[0, 0].plot(t_s, x[0,:], color='yellow')
axes[0, 0].set_xlabel('Time [s]', color='white')
axes[0, 0].set_ylabel('u [m/s]', color='white')
axes[0, 0].grid(True)
axes[0, 0].set_facecolor('black')
axes[0, 0].tick_params(colors='white')

#Velocidad en y v^b_CM/n
axes[0, 1].plot(t_s, x[1,:], color='yellow')
axes[0, 1].set_xlabel('Time [s]', color='white')
axes[0, 1].set_ylabel('v [m/s]', color='white')
axes[0, 1].grid(True)
axes[0, 1].set_facecolor('black')
axes[0, 1].tick_params(colors='white')

#Velocidad en z w^b_CM/n
axes[0, 2].plot(t_s, x[2,:], color='yellow')
axes[0, 2].set_xlabel('Time [s]', color='white')
axes[0, 2].set_ylabel('w [m/s]', color='white')
axes[0, 2].grid(True)
axes[0, 2].set_facecolor('black')
axes[0, 2].tick_params(colors='white')

#ángulo de roll
axes[0, 3].plot(t_s, x[6,:], color='yellow')
axes[0, 3].set_xlabel('Time [s]', color='white')
axes[0, 3].set_ylabel('phi [rad]', color='white')
axes[0, 3].grid(True)
axes[0, 3].set_facecolor('black')
axes[0, 3].tick_params(colors='white')

#Rate de roll p^b_b/n
axes[1, 0].plot(t_s, x[3,:], color='yellow')
axes[1, 0].set_xlabel('Time [s]', color='white')
axes[1, 0].set_ylabel('p [r/s]', color='white')
axes[1, 0].grid(True)
axes[1, 0].set_facecolor('black')
axes[1, 0].tick_params(colors='white')

#Rate de pitch q^b_b/n
axes[1, 1].plot(t_s, x[4,:], color='yellow')
axes[1, 1].set_xlabel('Time [s]', color='white')
axes[1, 1].set_ylabel('q [r/s]', color='white')
axes[1, 1].grid(True)
axes[1, 1].set_facecolor('black')
axes[1, 1].tick_params(colors='white')

#Rate de yaw r^b_b/n
axes[1, 2].plot(t_s, x[5,:], color='yellow')
axes[1, 2].set_xlabel('Time [s]', color='white')
axes[1, 2].set_ylabel('r [r/s]', color='white')
axes[1, 2].grid(True)
axes[1, 2].set_facecolor('black')
axes[1, 2].tick_params(colors='white')

#Ángulo pitch theta
axes[1, 3].plot(t_s, x[7,:], color='yellow')
axes[1, 3].set_xlabel('Time [s]', color='white')
axes[1, 3].set_ylabel('theta [rad]', color='white')
axes[1, 3].grid(True)
axes[1, 3].set_facecolor('black')
axes[1, 3].tick_params(colors='white')

#No hay gráfica del ángulo de yaw

plt.tight_layout()
#plt.savefig('imagenes_guardadas/sphere_drop_test_1.png')
plt.show()
