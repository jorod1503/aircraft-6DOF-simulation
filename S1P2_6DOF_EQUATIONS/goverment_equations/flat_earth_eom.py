import math
import numpy as np

def flat_earth_eom(t, x, amod):
    '''
    Para nombrar las variables se utiliza la siguiente convención: <nombre de la variable>_<sistema de coordenadas si aplica>_<unidades>.
    Por ejemplo, en el pitch rate q en el marco de referencia fijo en el cuerpo bf, con unidades de radianes por segundo sería q_b_rps
    Argumentos:
        t - tiempo [s], escalar
        x - vector de estado en tempo t [diversas unidades], arreglo numpy
            x[0] = u_b_mps, velocidad en x
            x[1] = v_b_mps, velocidad en y
            x[2] = w_b_mps, velocidad en z
            x[3] = p_b_rps, velocidad angular roll
            x[4] = q_b_rps, velocidad angular pitch
            x[5] = r_b_rps, velocidad angular yaw
            x[6] = phi_rad, ángulo roll
            x[7] = theta_rad, ángulo pitch
            x[8] = psi_rad, ángulo yaw
            x[9] = p1_n_m, posición x
            x[10] = p2_n_m, posición y
            x[11] = p3_n_m, posición z
    amod = datos del modelo de la aeronave almacenado como diccionario, contiene múltiples parámetros
    
    Respuestas:
        dx - derivada respecto al tiempo para cada estado en x
    '''
    #Pre localización del lado izquierdo de la ecuación
    dx = np.empty((12,),dtype=float)

    #Asignación de variables de estado a nombres de variables
    u_b_mps = x[0]
    v_b_mps = x[1]
    w_b_mps = x[2]
    p_b_rps = x[3]
    q_b_rps = x[4]
    r_b_rps = x[5]
    phi_rad = x[6]
    theta_rad = x[7]
    psi_rad = x[8]
    p1_n_m = x[9]
    p2_n_m = x[10]
    p3_n_m = x[11]

    #Obtención de momentos de masa e inercia
    m_kg =  amod["m_kg"]
    Jxz_b_kgm2 = amod["Jxz_b_kgm2"]
    Jxx_b_kgm2 = amod["Jxx_b_kgm2"]
    Jyy_b_kgm2 = amod["Jyy_b_kgm2"]
    Jzz_b_kgm2 = amod["Jzz_b_kgm2"]

    #Cálculo de datos de aire (Mach, altitud, AoA, AoS)

    #Modelo atmosférico

    #Gravedad
    gz_n_mps2 = 9.81

    #Convirtiendo gravedad a sistema de coordenadas fijo en el cuerpo
    gx_b_mps2 = -math.sin(theta_rad) * gz_n_mps2
    gy_b_mps2 = math.sin(phi_rad) * math.cos(theta_rad) * gz_n_mps2
    gz_b_mps2 = math.cos(phi_rad) * math.cos(theta_rad) * gz_n_mps2

    #Fuerzas externas
    Fx_b_kgmps2 = 0
    Fy_b_kgmps2 = 0
    Fz_b_kgmps2 = 0

    #Momentos externos
    l_b_kgm2ps2 = 0
    m_b_kgm2ps2 = 0
    n_b_kgm2ps2 = 0

    #ecuación del denominador en roll y yaw
    Den = Jxx_b_kgm2 * Jzz_b_kgm2 - Jxz_b_kgm2**2

    #Velocidad en eje x, estado u_b_mps
    dx[0] = 1 / m_kg * Fx_b_kgmps2 + gx_b_mps2 - w_b_mps * q_b_rps + v_b_mps * r_b_rps

    #Velocidad en eje y, estado v_b_mps
    dx[1] = 1 / m_kg * Fy_b_kgmps2 + gy_b_mps2 - u_b_mps * r_b_rps + w_b_mps * p_b_rps

    #Velocidad en eje z, estado w_b_mps
    dx[2] = 1 / m_kg * Fz_b_kgmps2 + gz_b_mps2 - v_b_mps * p_b_rps + u_b_mps * q_b_rps

    #Ecuaciones de roll, estado p_b_rps
    dx[3] = (Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2 + Jzz_b_kgm2) * p_b_rps * q_b_rps - \
            (Jzz_b_kgm2 * (Jzz_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2**2) * q_b_rps * r_b_rps + \
            Jzz_b_kgm2 * l_b_kgm2ps2 + \
            Jxz_b_kgm2 * n_b_kgm2ps2)/Den
    
    #Ecuaciones de roll, estado p_b_rps
    dx[4] = ((Jzz_b_kgm2 - Jxx_b_kgm2) * p_b_rps * r_b_rps - \
            Jxz_b_kgm2 * (p_b_rps**2 - r_b_rps**2) + m_b_kgm2ps2)/Jyy_b_kgm2
    
    #Ecuaciones de roll, estado p_b_rps
    dx[5] = ((Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2 )+ Jxz_b_kgm2**2) * p_b_rps * q_b_rps + \
            Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2 + Jzz_b_kgm2**2) * q_b_rps * r_b_rps + \
            Jxz_b_kgm2 * l_b_kgm2ps2 + \
            Jxz_b_kgm2 * n_b_kgm2ps2)/Den
    
    #Ecuaciones cinemáticas de Euler
    dx[6] = p_b_rps + math.sin(phi_rad)*math.tan(theta_rad)*q_b_rps + \
                      math.cos(phi_rad)*math.tan(theta_rad)*r_b_rps 
    dx[7] = math.cos(phi_rad)*q_b_rps - \
            math.sin(phi_rad)*r_b_rps
    dx[8] = math.sin(phi_rad)/math.cos(theta_rad)*q_b_rps + \
            math.cos(phi_rad)/math.cos(theta_rad)*r_b_rps

    #Ecuaciones de posición (navegación)
    dx[9] = 0
    dx[10] = 0
    dx[11] = 0

    return dx

