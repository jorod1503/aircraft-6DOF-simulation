import numpy as np

def forward_euler(f, t_s, x, h_s,amod):

    '''
    Método de Euler mejorado:

    Argumentos de entrada:
        f: función del lado derecho de la ecuación diferencial (dx/dt = f(t,x))
        t_s: Vector de puntos al tiempo en cel cuan las soluciones serán aproximadas
        x: Aproximación numérica de datos  a la ecuación diferencial f
        h_s: Tamaño de muestra en segundos

    Respuestas:
        t_s: Vector de puntos en tiempo en el cual las soluciones numéricas fueron aproximadas
        x: Solución numérica aproximada de los datos de la ecuación diferencial f
    '''

    for i in range(1, len(t_s)):
        x[:,i] = x[:,i-1] + h_s * f(t_s[i-1], x[:,i-1],amod)
    
    return t_s, x