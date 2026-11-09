# Gas lift en un pad · más producción, o la misma con menos gas

[![Pruebas de física](https://github.com/juanlaigle-hub/gas-lift-pad/actions/workflows/pruebas.yml/badge.svg)](https://github.com/juanlaigle-hub/gas-lift-pad/actions/workflows/pruebas.yml)

**[Abrir la demo](https://juanlaigle-hub.github.io/gas-lift-pad/)** · datos simulados · [English summary](#english-summary)

Modelo físico de 6 pozos con gas lift que comparten la misma compresión. La demo muestra tres cosas que importan en la operación de un pad:

| Pestaña | Qué muestra | Resultado con los datos simulados |
|---|---|---|
| **1 · Repartir el gas** | Explica el problema en 3 pasos, anima el paso de repartir el gas en partes iguales a repartirlo según la respuesta de cada pozo, y cuenta en palabras qué pozos ceden gas y cuáles lo reciben. Se elige el objetivo: más petróleo con el mismo gas, o el mismo petróleo con menos gas | Con 120 kSm³/d: **+0,6 % de petróleo con el mismo gas, o el mismo petróleo con 11 % menos gas**. Con los pozos al inicio del período y en todo el rango de capacidad de la demo (40 a 300 kSm³/d), la ganancia de producción va de 0,4 a 2,6 % y el ahorro de gas de 10 a 54 %; el ahorro es mayor cuando sobra compresión y el reparto uniforme inyecta de más. Con los pozos ya declinados (control "Días de operación"), los beneficios son menores |
| **2 · Ensayo con declinación** | Un ensayo de escalones en un pozo que declina: la parábola ajustada a los puntos confunde declinación con respuesta al gas | Error típico del óptimo (mediana de 60 casos): 3,3 kSm³/d con la parábola y entre 0 y 1 kSm³/d con el modelo físico y la declinación estimada de la historia previa |
| **3 · Diagnóstico** | Compara lo medido (caudal de líquido, presión de casing, corte de agua) con lo que predice el modelo y estima la causa más probable | Reconoce parafina, inyección por una válvula superior y frac hit con entrada de agua; con anomalías leves duda, y en operación normal da falsas alarmas en ~12 % de las mediciones |

El problema de la pestaña 2 es el que describe el trabajo técnico SPE-212577-MS (2023) sobre gas lift en Vaca Muerta: las curvas de optimización no fueron concluyentes porque los pozos tardan en estabilizarse y la declinación distorsiona los resultados, y la optimización debe hacerse a nivel de pad.

## Modelo

- **Reservorio:** IPR de Vogel; declinación hiperbólica de Arps de la productividad.
- **Columna de producción:** integrada en tramos de 50 m con Runge–Kutta de 2º orden; gravedad con la fracción de líquido del modelo de deriva de Zuber–Findlay (C₀ = 1,2, velocidad de deriva de Harmathy) y fricción de Darcy (Haaland). Por encima de la válvula operadora, gas de formación más gas inyectado.
- **Punto de operación:** cruce de la IPR con la columna, por bisección.
- **Presión de casing:** presión del tubing en la válvula más la caída en la válvula, menos el peso de la columna de gas en el anular.
- **Reparto óptimo:** programación dinámica sobre las curvas de los pozos (exacta en la grilla y válida aunque un pozo necesite un mínimo de gas para arrancar).
- **Diagnóstico:** cada anomalía se simula con el mismo modelo para obtener su huella en los residuos; para cada causa se estima la severidad que mejor explica lo medido y se calcula su probabilidad.

### Simplificaciones

- Un solo régimen de flujo en toda la columna; en un proyecto real se usan correlaciones por régimen o un modelo mecanístico calibrado.
- Todo el gas de formación se considera libre; z constante y temperatura lineal.
- Pozo vertical equivalente hasta el talón, sin la rama horizontal.
- Estado estacionario: no representa el cabeceo ni la descarga.
- Pozos, fluidos y declinación representativos, no de un yacimiento en particular.

## Pruebas

`tests/fisica.test.js` verifica en cada cambio: IPR, aumento de presión con la profundidad, forma de la curva de gas lift (sin gas no fluye, máximo interior, el exceso de gas resta), consistencia del punto de operación, presión de casing en rango, declinación, optimalidad del reparto (comparado con búsqueda exhaustiva), recuperación del óptimo en el ensayo con declinación y signos de las huellas del diagnóstico.

```bash
node --test tests/fisica.test.js
```

## Archivos

- `index.html`: la demo (abre sin servidor).
- `modelo.js`: núcleo físico, compartido por la página y las pruebas.
- `tests/`: pruebas automáticas.

## English summary

Physics-based demo of a gas-lifted pad (6 wells sharing one compressor): optimal allocation of injection gas (dynamic programming over well performance curves built from a Vogel IPR and a drift-flux tubing model), the effect of production decline on step-rate tests (a fitted parabola misplaces the optimum; the physical model with decline estimated from prior history recovers it), and model-residual diagnosis of paraffin, injection through an upper valve and frac hits. Simulated data.

---

[Juan Laiglecia](https://www.linkedin.com/in/juanlaiglecia) · [Ingeniería VR](https://www.ingenieriavr.com) · © Todos los derechos reservados.
