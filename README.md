# CNLS-Galerkin-Caso-de-Estudio-Numerico
Repositorio con el código utilizado para resolver la ecuación de Schrödinger no lineal acoplada mediante el método de Galerkin, incluyendo implementaciones, simulaciones y resultados del estudio presentado en la investigación.

# Descripción
Este repositorio contiene la implementación computacional completa del método de elementos finitos de Galerkin para la resolución numérica de la Ecuación de Schrödinger No Lineal (ESNL) acoplada. El trabajo reproduce y valida los resultados del paper de Ismail (2008) sobre simulación de solitones y sus interacciones no lineales.
La implementación incluye el esquema numérico de Galerkin P1 con discretización temporal implícita mediante la regla del punto medio, junto con el método de Newton para resolver el sistema no lineal resultante en cada paso temporal.
Sistema de ecuaciones

## Sistema de ecuaciones

El código resuelve el siguiente sistema acoplado de ecuaciones diferenciales parciales:

$$i \frac{\partial \Psi_1}{\partial t} + \frac{1}{2} \frac{\partial^2 \Psi_1}{\partial x^2} + (|\Psi_1|^2 + e|\Psi_2|^2) \Psi_1 = 0, \quad - \infty < x < \infty,$$
$$i \frac{\partial \Psi_2}{\partial t} + \frac{1}{2} \frac{\partial^2 \Psi_2}{\partial x^2} + (e|\Psi_1|^2 + |\Psi_2|^2) \Psi_2 = 0, \quad - \infty < x < \infty, $$

donde `e` es el coeficiente de modulación de fase cruzada, sujeto a condiciones de frontera de Neumann homogéneas.

# Características principales
### Propiedades numéricas

- Esquema conservativo: Preserva cantidades físicas fundamentales (masa, momento lineal, energía)
- Estabilidad incondicional: Validado mediante análisis de von Neumann
- Convergencia de segundo orden: O(h^2 + k^2) en espacio y tiempo respectivamente
- Procesamiento paralelo: Implementación optimizada con parfor para estudios de convergencia

# Validación

- Comparación con soluciones analíticas exactas para solitones simétricos
- Verificación de cantidades conservadas mediante cuadratura numérica
- Análisis sistemático de convergencia espacial y temporal

# Instalación 
git clone https://github.com/usuario/cnls-galerkin.git
cd cnls-galerkin

No se requieren pasos de instalación adicionales. Los archivos están listos para ejecutarse en MATLAB.

## Estructura del repositorio
```

├── nls.m                          # Solver principal (Galerkin + Newton) resuelve la ecuación con condiciones iniciales
├── nls_invariants.m               # Cálculo de cantidades conservadas
├── convergencia_espacial.m        # Análisis de convergencia espacial
├── convergencia_temporal.m        # Análisis de convergencia temporal
├── ResultadosESNL.mlx             # Notebook con todos los escenearios, casos, pruebas de convengencia y demás resultados
└── solucionExactaConservacion.mlx # Verificación simbólica de la conservacion del moemnto lineal y la energia 

```
