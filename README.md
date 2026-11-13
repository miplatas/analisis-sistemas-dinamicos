# Análisis y Simulación de Sistemas Dinámicos

Material de clase, ejemplos resueltos y plantillas de simulación para la Unidad de Aprendizaje **Análisis y Simulación de Sistemas Dinámicos** (5° semestre, Licenciatura en Ingeniería en Electromovilidad - FIME, UANL).

## Estructura de Ramas

* **`main`**: Documentación general del curso, programa analítico y enunciados.
* **`python`**: Todos los ejemplos y simulaciones resueltos en Python (`numpy`, `scipy`, `matplotlib`, `control`, `sympy`).
* **`scilab`**: Todos los ejemplos y simulaciones resueltos en Scilab (`.sce` / `.sci`).

## Contenido por Fases

1. **Fase 1:** Fundamentos y clasificación de sistemas dinámicos.
2. **Fase 2:** Análisis entrada-salida de sistemas lineales (tiempo y frecuencia, SISO y MIMO).
3. **Fase 3:** Métodos numéricos para la simulación de sistemas dinámicos no lineales (Euler y RK4).
4. **Fase 4:** Análisis de estabilidad (Routh-Hurwitz, linealización y Lyapunov).
5. **PIA:** Proyecto Integrador de Aprendizaje.

## Inicio rápido en la rama Python

Los notebooks presentan conceptos, ejemplos ejecutables y ejercicios de cada
tema de las cuatro fases. La biblioteca compartida `sistemas_dinamicos` reúne
operaciones repetitivas para que los ejemplos se centren en el modelo; también
se utilizan directamente NumPy, SciPy y SymPy cuando ayudan a mostrarlo.

Desde la raíz del repositorio, con Python 3.10 o posterior:

```bash
git switch python
python -m venv .venv
source .venv/bin/activate  # En Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python -m jupyter lab
```

La instalación editable permite importar la biblioteca desde los notebooks
de cualquier carpeta. En Jupyter o VS Code, selecciona el intérprete/kernel
del entorno `.venv` y ejecuta las celdas en orden. Si utilizas un kernel distinto,
instala las dependencias y la biblioteca también en ese entorno.

### Cuadernos

* **Fase 1:** [Señales y sistemas](Fase_01_Fundamentos_y_Modelado/1_1_Senales_y_Sistemas/tema.ipynb),
  [clasificación](Fase_01_Fundamentos_y_Modelado/1_2_Clasificacion_Propiedades/tema.ipynb),
  [ecuaciones y Laplace](Fase_01_Fundamentos_y_Modelado/1_3_Ecuaciones_Diferenciales_Laplace/tema.ipynb),
  [transferencia y estados](Fase_01_Fundamentos_y_Modelado/1_4_Funcion_Transferencia_Espacio_Estados/tema.ipynb).
* **Fase 2:** [impulso y convolución](Fase_02_Analisis_Sistemas_Lineales/2_1_Respuesta_Impulso_Convolucion/tema.ipynb),
  [transitorio y error](Fase_02_Analisis_Sistemas_Lineales/2_2_Respuesta_Transitoria_Error_SISO/tema.ipynb),
  [estados MIMO](Fase_02_Analisis_Sistemas_Lineales/2_3_Solucion_Ecuacion_Estados_MIMO/tema.ipynb),
  [Bode](Fase_02_Analisis_Sistemas_Lineales/2_4_Respuesta_Frecuencia_Bode/tema.ipynb).
* **Fase 3:** [Euler](Fase_03_Metodos_Numericos_No_Lineales/3_1_Metodo_Euler/tema.ipynb),
  [Runge-Kutta 4](Fase_03_Metodos_Numericos_No_Lineales/3_2_Metodo_Runge_Kutta_4/tema.ipynb),
  [error y estabilidad numérica](Fase_03_Metodos_Numericos_No_Lineales/3_3_Errores_Paso_Estabilidad_Numerica/tema.ipynb).
* **Fase 4:** [polos y Routh-Hurwitz](Fase_04_Analisis_Estabilidad/4_1_Polos_Ceros_Routh_Hurwitz/tema.ipynb),
  [equilibrios y linealización](Fase_04_Analisis_Estabilidad/4_2_Puntos_Equilibrio_Linealizacion/tema.ipynb),
  [Lyapunov](Fase_04_Analisis_Estabilidad/4_3_Estabilidad_Lyapunov/tema.ipynb).

Las carpetas de evidencias y del PIA quedan pendientes de sus enunciados y rúbricas.

### Biblioteca y pruebas

La biblioteca ofrece `escalon`, `graficar`, `convolucion`, `integrar` y
`routh`. La convolución supone señales causales muestreadas uniformemente
desde cero; `integrar` admite Euler explícito y RK4 con una malla creciente.
`routh` admite coeficientes reales numéricos y avisa cuando hay pivotes cero
o filas nulas, que requieren el tratamiento especial del criterio.

Para comprobar la biblioteca:

```bash
python -m unittest discover -s tests -v
```
