# Análisis y Simulación de Sistemas Dinámicos

Material de clase, ejemplos resueltos y plantillas de simulación para la Unidad de Aprendizaje **Análisis y Simulación de Sistemas Dinámicos** (5° semestre, Licenciatura en Ingeniería en Electromovilidad - FIME, UANL).

El material utiliza Python y notebooks de Jupyter, con bibliotecas como
NumPy, SciPy, Matplotlib, SymPy y python-control. La rama principal es `main`.

## Contenido por Fases

1. **Fase 1:** Fundamentos y clasificación de sistemas dinámicos.
2. **Fase 2:** Análisis entrada-salida de sistemas lineales (tiempo y frecuencia, SISO y MIMO).
3. **Fase 3:** Métodos numéricos para la simulación de sistemas dinámicos no lineales (Euler y RK4).
4. **Fase 4:** Análisis de estabilidad (Routh-Hurwitz, linealización y Lyapunov).
5. **PIA:** Proyecto Integrador de Aprendizaje.

## Inicio rápido

Los notebooks presentan conceptos, ejemplos ejecutables y ejercicios de cada
tema de las cuatro fases. La biblioteca compartida `sistemas_dinamicos` reúne
operaciones repetitivas para que los ejemplos se centren en el modelo; también
se utilizan directamente NumPy, SciPy y SymPy cuando ayudan a mostrarlo.

### Instalación local para estudiantes

Necesitas Python 3.10 o posterior y Git. Descarga el repositorio y entra
en la carpeta del proyecto:

```bash
git clone -b main https://github.com/miplatas/analisis-sistemas-dinamicos.git
cd analisis-sistemas-dinamicos
```

Si ya descargaste el repositorio, entra en su carpeta y actualiza `main`:

```bash
git switch main
git pull --ff-only origin main
```

Cada estudiante debe crear su propio entorno virtual. La carpeta `.venv`
no se incluye en el repositorio porque depende de la computadora y del sistema
operativo:

```bash
python -m venv .venv
```

Activa el entorno con el comando correspondiente a tu terminal:

```bash
# Linux o macOS
source .venv/bin/activate
```

```powershell
# Windows: PowerShell
.\.venv\Scripts\Activate.ps1
```

```bat
REM Windows: símbolo del sistema (CMD)
.venv\Scripts\activate.bat
```

Con el entorno activado, instala las dependencias y la biblioteca del curso:

```bash
python -m pip install -r requirements.txt
```

La instalación inicial necesita conexión a internet; después puedes ejecutar
los notebooks sin conexión. Para usarlos, elige una de estas opciones:

* **VS Code:** instala las extensiones Python y Jupyter, abre la carpeta del
  proyecto y un notebook. Al ejecutar una celda, selecciona **Python
  Environments** y el intérprete de `.venv` como kernel.
* **JupyterLab:** desde la terminal con el entorno activado, ejecuta:

```bash
python -m jupyter lab
```

La instalación editable permite importar la biblioteca desde los notebooks
de cualquier carpeta. Ejecuta las celdas en orden. Si utilizas un kernel distinto,
instala las dependencias y la biblioteca también en ese entorno. En futuras
sesiones no necesitas repetir la instalación: basta con activar `.venv` y abrir
los notebooks; repítela si cambian las dependencias.

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
