"""Operaciones comunes que mantienen los ejemplos centrados en el modelo."""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy import signal


def escalon(t, inicio=0):
    """Escalón unitario, igual a uno a partir de ``inicio``."""
    return np.asarray(np.asarray(t) >= inicio, dtype=float)


def graficar(t, series, *, titulo="", xlabel="Tiempo (s)", ylabel="Amplitud"):
    """Graficar series con estilos sólido, raya-punto y punteado en ciclo."""
    fig, ax = plt.subplots()
    estilos = ("-", "-.", ":")
    for indice, (nombre, valores) in enumerate(series.items()):
        ax.plot(t, valores, label=nombre, linestyle=estilos[indice % len(estilos)])
    ax.set(title=titulo, xlabel=xlabel, ylabel=ylabel)
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    return fig, ax


def campo_fase(f, limites_x, limites_y, *, puntos=(25, 21), titulo="",
               xlabel="Estado 1", ylabel="Estado 2"):
    """Representar un campo autónomo de dos estados con flechas normalizadas.

    El color indica la magnitud de las derivadas expresadas en las unidades
    elegidas para cada eje; con unidades distintas no es una rapidez física.
    """
    x = np.linspace(*limites_x, puntos[0])
    y = np.linspace(*limites_y, puntos[1])
    X, Y = np.meshgrid(x, y)
    U, V = f(X, Y)
    U, V = np.asarray(U, dtype=float), np.asarray(V, dtype=float)
    if U.shape != X.shape or V.shape != X.shape:
        raise ValueError("El campo debe devolver dos matrices con la forma de la malla.")
    if not np.all(np.isfinite(U)) or not np.all(np.isfinite(V)):
        raise ValueError("El campo debe contener derivadas finitas.")
    magnitud = np.hypot(U, V)
    norma = np.where(magnitud == 0, 1, magnitud)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    flechas = ax.quiver(X, Y, U / norma, V / norma, magnitud,
                        cmap="viridis", pivot="mid")
    fig.colorbar(flechas, ax=ax, label="Magnitud de derivadas (escala elegida)")
    ax.set(title=titulo, xlabel=xlabel, ylabel=ylabel)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig, ax


def retrato_fase(trayectorias, *, titulo="", xlabel="Estado 1", ylabel="Estado 2"):
    """Representar trayectorias bidimensionales y marcar sus condiciones iniciales."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, trayectoria in enumerate(trayectorias):
        trayectoria = np.asarray(trayectoria, dtype=float)
        if trayectoria.ndim != 2 or trayectoria.shape[1] != 2 or not len(trayectoria):
            raise ValueError("Cada trayectoria debe tener forma (instantes, 2).")
        if not np.all(np.isfinite(trayectoria)):
            raise ValueError("Las trayectorias deben contener valores finitos.")
        ax.plot(trayectoria[:, 0], trayectoria[:, 1], label=f"Trayectoria {i + 1}")
        ax.plot(trayectoria[0, 0], trayectoria[0, 1], "k.", markersize=6)
    ax.set(title=titulo, xlabel=xlabel, ylabel=ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    return fig, ax


def convolucion(t, entrada, impulso):
    """Aproximar la convolución causal continua con muestras uniformes.

    Se supone que ambas señales comienzan en t=0 y valen cero antes de él.
    """
    t = np.asarray(t, dtype=float)
    entrada = np.asarray(entrada, dtype=float)
    impulso = np.asarray(impulso, dtype=float)
    if t.ndim != 1 or len(t) < 2 or not np.isclose(t[0], 0):
        raise ValueError("El tiempo debe ser un vector que comience en cero.")
    dt = np.diff(t)
    if not np.all(np.isfinite(t)) or dt[0] <= 0 or not np.allclose(dt, dt[0]):
        raise ValueError("El tiempo debe tener paso constante y positivo.")
    if entrada.shape != t.shape or impulso.shape != t.shape:
        raise ValueError("Las señales deben tener la misma forma que el tiempo.")
    if not np.all(np.isfinite(entrada)) or not np.all(np.isfinite(impulso)):
        raise ValueError("Las señales deben contener valores finitos.")
    return signal.convolve(entrada, impulso, mode="full")[: len(t)] * dt[0]


def integrar(f, tiempo, inicial, *, metodo="rk4"):
    """Integrar x'=f(t, x) en una malla creciente con Euler o RK4.

    Retorna una matriz con una fila por instante y una columna por estado.
    ``f`` recibe siempre el estado como vector, incluso en problemas escalares.
    """
    tiempo = np.asarray(tiempo, dtype=float)
    inicial = np.atleast_1d(np.asarray(inicial, dtype=float))
    if tiempo.ndim != 1 or len(tiempo) < 2 or not np.all(np.isfinite(tiempo)):
        raise ValueError("El tiempo debe ser un vector finito de al menos dos puntos.")
    if np.any(np.diff(tiempo) <= 0):
        raise ValueError("El tiempo debe ser estrictamente creciente.")
    if inicial.ndim != 1 or inicial.size == 0 or not np.all(np.isfinite(inicial)):
        raise ValueError("La condición inicial debe ser un vector finito no vacío.")
    if metodo not in ("euler", "rk4"):
        raise ValueError("El método debe ser 'euler' o 'rk4'.")

    def evaluar(t, x):
        resultado = np.asarray(f(t, x), dtype=float)
        if resultado.shape != inicial.shape or not np.all(np.isfinite(resultado)):
            raise ValueError("La derivada debe ser un vector finito con la forma del estado.")
        return resultado

    estados = np.empty((len(tiempo), len(inicial)))
    estados[0] = inicial
    for i, h in enumerate(np.diff(tiempo)):
        t, x = tiempo[i], estados[i]
        k1 = evaluar(t, x)
        if metodo == "euler":
            estados[i + 1] = x + h * k1
        else:
            k2 = evaluar(t + h / 2, x + h * k1 / 2)
            k3 = evaluar(t + h / 2, x + h * k2 / 2)
            k4 = evaluar(t + h, x + h * k3)
            estados[i + 1] = x + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if not np.all(np.isfinite(estados[i + 1])):
            raise ValueError("La integración produjo estados no finitos.")
    return estados


def routh(coeficientes):
    """Tabla de Routh para coeficientes reales numéricos en orden descendente.

    Los casos especiales de pivote cero o fila nula requieren tratamiento
    simbólico adicional y se indican explícitamente mediante ValueError.
    """
    coef = [sp.sympify(c) for c in coeficientes]
    if len(coef) < 2 or any(c.is_number is not True or c.is_real is not True
                            or c.is_finite is not True for c in coef):
        raise ValueError("Se requieren coeficientes reales numéricos finitos.")
    if coef[0] == 0:
        raise ValueError("El coeficiente principal debe ser distinto de cero.")
    grado = len(coef) - 1
    columnas = (len(coef) + 1) // 2
    tabla = sp.zeros(grado + 1, columnas)
    for j, valor in enumerate(coef[::2]):
        tabla[0, j] = valor
    for j, valor in enumerate(coef[1::2]):
        tabla[1, j] = valor
    for i in range(2, grado + 1):
        pivote = tabla[i - 1, 0]
        if pivote == 0:
            raise ValueError("Pivote cero o fila nula: se requiere el caso especial de Routh.")
        for j in range(columnas - 1):
            tabla[i, j] = sp.simplify(
                (pivote * tabla[i - 2, j + 1]
                 - tabla[i - 2, 0] * tabla[i - 1, j + 1]) / pivote
            )
    if tabla[grado, 0] == 0:
        raise ValueError("Fila nula: se requiere el caso especial de Routh.")
    return tabla
