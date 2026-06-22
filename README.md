# Dinámica Atmosférica sobre la Península Ibérica

**Autor:** Jesús Moral  


---

## Objetivo

Resolver numéricamente las **ecuaciones primitivas de la dinámica atmosférica**  sobre la
Península Ibérica, para entenderlas mejor y practicar, usando datos reales de Open-Meteo (modelo GFS) como condiciones
iniciales y de contorno.

---

## Física del modelo

### Coordenadas isobáricas y presión vs. altitud

El modelo trabaja en **coordenadas isobáricas (x, y, p, t)**: la vertical no es la
altitud geométrica sino la **presión** del aire.

```
Altitud   Presión típica   k (nivel solver)
─────────────────────────────────────────────
  ~100 m   1000 hPa   →   k = 0  (nivel más bajo del solver)
  ~750 m    925 hPa   →   k = 1
 ~1500 m    850 hPa   →   k = 2
 ~3000 m    700 hPa   →   k = 3
 ~5500 m    500 hPa   →   k = 4
 ~9200 m    300 hPa   →   k = 5  (tope del dominio)
```

> **Por qué no hay nivel de superficie:** La superficie terrestre no es una superficie
> isobárica — cada punto tiene una presión distinta según su elevación (desde ~800 hPa en
> alta montaña hasta ~1015 hPa en costa). Incluirla en un esquema de diferencias finitas
> isobáricas crearía gradientes artificiales de presión incompatibles con el modelo. El
> nivel 1000 hPa (~100 m de altitud) es el proxy de "capa baja" del dominio.

> **Regla fundamental:** cuanto más alto → menos presión. La columna de aire que tienes
> por encima pesa menos. Por eso 300 hPa es el nivel más ALTO y 1000 hPa el más BAJO.

La consecuencia directa en el código es que **k aumenta con la altitud y decrece con la
presión**: `LEVELS_hPa = [1000, 925, 850, 700, 500, 300]`.

### Velocidad vertical isobárica ω

En coordenadas isobáricas, la velocidad vertical es `ω = dp/dt` [Pa/s]:

| ω | Significado físico |
|---|-------------------|
| ω > 0 | Hundimiento: la parcela se mueve hacia mayor p → baja de altitud |
| ω < 0 | Ascenso: la parcela se mueve hacia menor p → sube de altitud |
| ω = 0 | Sin movimiento vertical |

Esto es el **opuesto** del convenio habitual de velocidad vertical w [m/s], donde w > 0
significa ascenso.

### Ecuaciones prognosticadas

```
∂u/∂t = −u ∂u/∂x − v ∂u/∂y − ω ∂u/∂p + fv − ∂Φ/∂x + ν∇²u
∂v/∂t = −u ∂v/∂x − v ∂v/∂y − ω ∂v/∂p − fu − ∂Φ/∂y + ν∇²v
∂T/∂t = −u ∂T/∂x − v ∂T/∂y − ω(∂T/∂p − Rd·T/(cp·p)) + ν∇²T
```

### Ecuaciones diagnósticas (post-proceso cada paso)

| Campo | Cálculo |
|-------|---------|
| ω (continuidad) | ∂ω/∂p = −(∂u/∂x + ∂v/∂y), integrada de p=300 hPa (ω=0) hacia abajo |
| Φ (hidrostático) | Φₖ = Φₖ₋₁ − Rd·T̄·ln(pₖ/pₖ₋₁), anclada en Φ₁₀₀₀ fijo |

---

## Dominio y malla

```
Longitud : −10° a +5°   paso 0.75°   →  21 columnas
Latitud  :  35.5° a 44°  paso 0.75°   →  12 filas
Total    :  21 × 12 = 252 nodos horizontales
Vertical :  7 niveles (superficie + 1000/925/850/700/500/300 hPa)
```

**Resolución horizontal:** ~83 km en latitud, ~60–67 km en longitud (depende de la latitud).  
**CFL:** Δx/U_max = 64 000 m / 50 m/s ≈ 1280 s → se usa Δt = 600 s (CFL ≈ 0.50).

---

## Estructura del proyecto

```
meteo/
├── crear_malla.py               # Genera la malla regular y data.db
├── datos_openmeteo.py           # Descarga datos GFS por nodo → data.db
├── condiciones_iniciales.py     # Construye campo 3D → condiciones_iniciales_fd.csv
├── solver_fd.py                 # Solver diferencias finitas upwind → resultados_fd.csv
├── plot_resultados_fd.py        # Visualiza resultados → PNGs por nivel
├── data.db                      # Base de datos (nodos + niveles isobáricos)
├── malla_rectangular.png        # Mapa con nodos tierra/mar
├── condiciones_iniciales_fd.csv # Campo inicial (1512 filas = 252 nodos × 6 niveles)
├── resultados_fd.csv            # Campos finales tras 1h (1512 filas)
└── resultado_fd_nivel_*hPa.png  # Mapas por nivel: ΔT, |V|, Φ, ω
```

---


## Base de datos — `data.db`

```sql
nodes (252 filas) — FIJA tras crear_malla.py
  id       INTEGER PK
  i        INTEGER    índice columna (lon, 0=oeste)
  j        INTEGER    índice fila    (lat, 0=sur)
  lon      REAL       Longitud [°E]
  lat      REAL       Latitud [°N]
  elev_m   REAL       Elevación [m s.n.m.]
  is_land  INTEGER    1=tierra, 0=mar

pressure_levels (1500 filas) — sobreescrita en cada ejecución de datos_openmeteo.py
  node_id   INTEGER FK → nodes
  timestamp TEXT
  level_hPa INTEGER    1000 / 925 / 850 / 700 / 500 / 300
  T_C       REAL       Temperatura [°C]
  geop_m    REAL       Altura geopotencial [m]
  u_ms      REAL       Viento zonal [m/s]
  v_ms      REAL       Viento meridional [m/s]
  omega_pas REAL       Velocidad vertical ω [Pa/s]
  dewp_C    REAL       Punto de rocío [°C]
  q         REAL       Humedad específica [kg/kg]
  PRIMARY KEY (node_id, timestamp, level_hPa)
```

---

## Flujo de ejecución

```
── Preparación (ejecutar una sola vez) ───────────────────────────────────────
1.  python3 crear_malla.py           → data.db  +  malla_rectangular.png

── Datos frescos (repetir cuando se quieran datos actualizados) ──────────────
2.  python3 datos_openmeteo.py       → data.db (surface + pressure_levels)
3.  python3 condiciones_iniciales.py → condiciones_iniciales_fd.csv

── Simulación ────────────────────────────────────────────────────────────────
4.  python3 solver_fd.py             → resultados_fd.csv
5.  python3 plot_resultados_fd.py    → resultado_fd_*.png
```

---

## Limitaciones del modelo

| Limitación | Motivo |
|-----------|--------|
| Sin radiación (Q=0) | Ecuación termodinámica seca adiabática |
| Sin nivel de superficie | La superficie no es isobárica — dominio empieza en 1000 hPa |
| Humedad pasiva | q se advecta pero no retroalimenta la dinámica |
| BCs laterales rígidas | Flujo entrante/saliente sin absorción de ondas |
| Sin actualización de p_sfc | La presión de superficie se mantiene constante |
| Dominio pequeño (15°×8.5°) | Efectos de contorno importantes en simulaciones >2h |

---

## Dependencias

```bash
pip install numpy matplotlib cartopy requests shapely
```

No requiere FEniCS ni ningún solver externo. Python 3.10+ recomendado.

---

## APIs utilizadas

| API | Endpoint | Auth | Uso |
|-----|----------|------|-----|
| Open-Meteo Forecast | `api.open-meteo.com/v1/forecast` | Sin clave | T, viento, ω, geopotencial por nodo |
| Open-Meteo Elevation | `api.open-meteo.com/v1/elevation` | Sin clave | Elevación del terreno en batch |
