import numpy as np
import matplotlib.pyplot as plt

# Constantes del oscilador
G, W = 0.5, 2 * np.pi

def x_ex(t): return np.exp(-G * t) * np.cos(W * t)
def v_ex(t): return -np.exp(-G * t) * (G * np.cos(W * t) + W * np.sin(W * t))
def a_ex(t): return np.exp(-G * t) * ((G**2 - W**2) * np.cos(W * t) + 2 * G * W * np.sin(W * t))

def calcular_derivadas(dt):
    n_puntos = int(round(3.0 / dt)) + 1
    t = np.linspace(0.0, 3.0, n_puntos)
    x = x_ex(t)

    # 1. Primera derivada (velocidad) por rebanadas
    v_adel = (x[1:] - x[:-1]) / dt
    v_atra = (x[1:] - x[:-1]) / dt
    v_cent = (x[2:] - x[:-2]) / (2 * dt)
    v_cinco = (-x[4:] + 8*x[3:-1] - 8*x[1:-3] + x[:-4]) / (12 * dt)

    # 2. Segunda derivada (aceleración) central
    a_cent = (x[2:] - 2*x[1:-1] + x[:-2]) / (dt**2)

    # Errores máximos frente a valores exactos
    err_v_adel = np.max(np.abs(v_adel - v_ex(t[:-1])))
    err_v_atra = np.max(np.abs(v_atra - v_ex(t[1:])))
    err_v_cent = np.max(np.abs(v_cent - v_ex(t[1:-1])))
    err_v_cinco = np.max(np.abs(v_cinco - v_ex(t[2:-2])))
    err_a_cent = np.max(np.abs(a_cent - a_ex(t[1:-1])))

    return t, x, v_cent, a_cent, [err_v_adel, err_v_atra, err_v_cent, err_v_cinco, err_a_cent]

# Ejecución para dt = 0.01 y dt = 0.005
t_01, x_01, v_c01, a_c01, errs_01 = calcular_derivadas(0.01)
_, _, _, _, errs_005 = calcular_derivadas(0.005)

# Guardar derivadas_oscilador.dat con doble precisión completa (16 decimales)
datos_interiores = np.column_stack((t_01[1:-1], x_01[1:-1], v_c01, a_c01))
np.savetxt("derivadas_oscilador.dat", datos_interiores, fmt="%.16e")

# Guardar errores_derivadas.dat
tabla_errores = np.array([
    [0.01, *errs_01],
    [0.005, *errs_005]
])
np.savetxt("errores_derivadas.dat", tabla_errores, fmt="%.16e")

# Generar gráfica derivadas_oscilador.png
t_int = t_01[1:-1]
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

ax1.plot(t_int, v_c01, 'b-', label='v num (central)', alpha=0.8)
ax1.plot(t_int, v_ex(t_int), 'r--', label='v exacta', alpha=0.7)
ax1.set_ylabel('Velocidad [m/s]')
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3)

ax2.plot(t_int, a_c01, 'g-', label='a num (central)', alpha=0.8)
ax2.plot(t_int, a_ex(t_int), 'm--', label='a exacta', alpha=0.7)
ax2.set_xlabel('Tiempo [s]')
ax2.set_ylabel('Aceleración [m/s²]')
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3)

plt.suptitle(r'Oscilador amortiguado: Velocidad y Aceleración ($\Delta t = 0.01$ s)')
plt.tight_layout()
plt.savefig('derivadas_oscilador.png', dpi=150)
plt.close()
