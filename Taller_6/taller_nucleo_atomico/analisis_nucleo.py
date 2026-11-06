"""Taller 6 -- Módulo de análisis (escriba aquí sus cuatro funciones)."""
import numpy as np
import matplotlib.pyplot as plt

def leer_valle(archivo="taller_nucleo_atomico/valle.dat"):
    """Lee valle.dat y devuelve A, Z y E_B."""
    datos = np.loadtxt(archivo)
    return datos[:, 0].astype(int), datos[:, 1].astype(int), datos[:, 2]

def encontrar_maximo_enlace(A, Z, E_B):
    """Encuentra el núcleo más ligado (E_B más negativo) con un ciclo for."""
    min_E_B = E_B[0]
    best_A, best_Z = A[0], Z[0]
    for a_i, z_i, eb_i in zip(A, Z, E_B):
        if eb_i < min_E_B:
            min_E_B = eb_i
            best_A = a_i
            best_Z = z_i
    return best_A, best_Z, min_E_B

def guardar_maximo_python(best_A, best_Z, min_E_B, archivo="taller_nucleo_atomico/maximo_python.dat"):
    """Guarda maximo_python.dat con el formato pedido."""
    with open(archivo, "w") as f:
        f.write(f"{best_A:4d} {best_Z:4d} {min_E_B:16.8f}\n")

def graficar_curva_enlace(A, E_B, experimentales, archivo="taller_nucleo_atomico/curva_enlace.png"):
    """Grafica la curva del modelo y los puntos experimentales."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(A, -E_B, "-", color="#162A54", label="Modelo gota líquida")

    exp_A = [e[2] for e in experimentales]
    exp_BA = [e[3] for e in experimentales]
    ax.scatter(exp_A, exp_BA, color="#B02A20", label="Datos experimentales (AME2020)", zorder=5)

    ax.set_xlabel(r"$A$ (Número de nucleones)")
    ax.set_ylabel(r"$B/A = -E_B$ [MeV]")
    ax.set_title("Curva de energía de enlace por nucleón")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend()
    plt.tight_layout()
    plt.savefig(archivo, dpi=300)
    plt.show()
