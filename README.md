# Proyecto MN — Sistema No Lineal de Tanques Acoplados
## Infraestructura de Experimentos en MATLAB

## Estructura de archivos

```
proyecto_tanques/
├── functions/
│   ├── campo_vectorial.m    % f(h) con clipping de raices negativas y saturacion por Hmax
│   ├── calc_jacobian.m      % Jacobiana analitica 2x2 de f(h)
│   ├── newton_multivar.m    % Newton-Raphson multivariable -> equilibrio h*
│   ├── solver_rk4.m         % Integrador explicito RK4 (soporta Qin(t) variable)
│   └── solver_implicit.m    % Euler Implicito + Newton-Raphson embebido (robusto)
├── main_in_mlx.mlx         % Orquesta los 4 experimentos de manera modular
├── main.m         % Orquesta los 4 experimentos (ver mas abajo)
├── data/                    % Se crea automaticamente: .mat con resultados de cada experimento
└── figures/                 % Se crea automaticamente: .png de cada figura (200 dpi)
```

## Como ejecutarlo

1. Abrir MATLAB y situar la carpeta `functions/` en el path (el script ya hace
   `addpath('functions')` automaticamente, basta con ejecutar `main.m`
   desde la carpeta raiz del proyecto).
2. Ejecutar `main.m` de principio a fin (Run, o seccion por seccion
   con Ctrl+Enter si se trabaja en el editor normal).
3. Para entregar el archivo como Live Script (`.mlx`), tal como exigen los
   lineamientos del curso: abrir main.m` en MATLAB y usar
   Save As → MATLAB Live Code File (*.mlx). Las secciones `%% ...` ya
   estan delimitadas para que el Live Editor las reconozca como celdas
   independientes con su salida (graficos y texto) intercalada.

## Resumen de los 4 experimentos (valores ya validados)

| Experimento | Metodo(s) | Resultado clave |
|---|---|---|
| 1. Trayectoria base | Newton-Raphson + RK4 | h* = [1.0000, 1.5625] m; RK4 reproduce h* con error < 1e-9 |
| 2. Escalon de Qin (t=500s) | RK4 con Qin(t) | Nuevo equilibrio h* = [2.5600, 4.0000] m, alcanzado por la simulacion |
| 3. Rigidez (c1=0.5, c2=0.001, dt=4s) | RK4 vs Euler Implicito | RK4 oscila en [-0.16, 0.00] m (inestable); Implicito converge a 0.0100 m con ~2-3 iteraciones Newton/paso |
| 4. Inundacion (Qin=0.15, Hmax=3m) | Euler Implicito | Tanque 1 desborda en t≈35.0 s, Tanque 2 en t≈83.5 s; volumen perdido ≈ 92.9 m³ (T1) y 24.5 m³ (T2) en 1500 s |

## Notas de diseño relevantes

- **Experimento 3**: con los parametros base del sistema (c1=0.05) los
  autovalores de la Jacobiana en el equilibrio son demasiado pequenos
  (|lambda| ~ 0.025 [1/s]) para desestabilizar a RK4 con dt=4s
  (dt·lambda ≈ 0.1, muy por debajo del limite de estabilidad ≈ 2.785).
  Por ello se incremento c1 a 0.5 (orificio de salida del Tanque 1 mas
  grande, dinamica rapida) manteniendo c2=0.001 (valvula del Tanque 2
  casi cerrada, dinamica lenta), maximizando la separacion de escalas
  temporales y logrando dt·lambda1 = -10, claramente fuera de la region
  de estabilidad de RK4. Esto se documenta dentro de main.m`.
- **Experimento 4**: la saturacion por desborde introduce una
  discontinuidad genuina (no un error numerico) en el campo vectorial
  exactamente en h=Hmax. El Newton-Raphson embebido en
  `solver_implicit.m` incorpora backtracking (line search) y deteccion
  de estancamiento para converger de forma robusta y precisa incluso en
  ese instante de contacto con la pared del tanque.
