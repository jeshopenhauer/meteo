# Dinámica Atmosférica sobre la Península Ibérica

**Autor:** Jesús Moral  
**Asignatura:** Meteorología Física

---

## Objetivo

Resolver numéricamente las **ecuaciones primitivas de la dinámica atmosférica** sobre la
Península Ibérica usando datos reales de Open-Meteo (modelo GFS) como condiciones
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

## Scripts

### `crear_malla.py`

Genera la malla rectangular regular y crea `data.db` desde cero.

| | Detalle |
|--|---------|
| **Input** | Ninguno (parámetros en el código) |
| **Output** | `data.db` (nodos), `malla_rectangular.png` |

- Máscara tierra/mar: Natural Earth 50m (shapely)
- Elevaciones: API Open-Meteo Elevation en batch de 100 nodos
- Resultado: 252 nodos, ~141 tierra, ~111 mar, elevación 0–1971 m

---

### `datos_openmeteo.py`

Descarga todos los campos meteorológicos para cada nodo desde Open-Meteo GFS (28 km).
Una petición HTTP por nodo, 3 reintentos en caso de timeout. Sobreescribe los datos
previos en cada ejecución.

**Fuente:** `https://api.open-meteo.com/v1/forecast`

| | Detalle |
|--|---------|
| **Input** | `data.db` → tabla `nodes` |
| **Output** | `data.db` → tablas `surface` y `pressure_levels` |

**Variables de superficie:**

| Campo | Descripción | Unidad |
|-------|-------------|--------|
| `T_2m` | Temperatura a 2 m | °C |
| `p_sfc` | Presión de superficie | hPa |
| `u_10m`, `v_10m` | Viento a 10 m (componentes) | m/s |
| `dewp_2m` | Punto de rocío a 2 m | °C |
| `q_2m` | Humedad específica | kg/kg |

**Variables isobáricas** (1000/925/850/700/500/300 hPa):

| Campo | Descripción | Unidad |
|-------|-------------|--------|
| `T_C` | Temperatura | °C |
| `geop_m` | Altura geopotencial | m |
| `u_ms`, `v_ms` | Viento (componentes) | m/s |
| `omega_pas` | Velocidad vertical ω | Pa/s |
| `dewp_C` | Punto de rocío | °C |
| `q` | Humedad específica | kg/kg |

> **Importante sobre niveles subterráneos:** Para nodos de montaña, el nivel 1000 hPa
> puede estar por debajo del terreno. Open-Meteo rellena estos niveles con valores
> extrapolados del modelo GFS (geopotencial coherente, ~100 m, NO la elevación del
> terreno). Esto garantiza gradientes de presión realistas y evita inestabilidades.

---

### `condiciones_iniciales.py`

Construye el campo inicial 3D (252 nodos × 7 niveles) leyendo `data.db`.

| | Detalle |
|--|---------|
| **Input** | `data.db` → tablas `nodes`, `surface`, `pressure_levels` |
| **Output** | `condiciones_iniciales_fd.csv` (1764 filas) |

**Algoritmo por nodo:**

1. `k=0` (superficie): T_2m, p_sfc, u/v 10m, q_2m de la tabla `surface`; Φ = g·elev_m
2. `k=1..6`: datos isobáricos de `pressure_levels`; **Φ = g·geop_m directamente** (no
   integración hidrostática desde superficie, para evitar artefactos orográficos)
3. Niveles subterráneos (p_nivel ≥ p_sfc): T, u, v, q de superficie; Φ del modelo GFS
4. Nodos sin datos: interpolación nearest-neighbor de los vecinos más próximos

**Formato de salida:**

| Columna | Descripción | Unidad |
|---------|-------------|--------|
| `node_id`, `i`, `j` | Índices del nodo | — |
| `k` | Índice de nivel (0=sfc, 1=1000hPa … 6=300hPa) | — |
| `lon`, `lat` | Coordenadas geográficas | ° |
| `level_hPa` | Nivel isobárico (vacío en k=0) | hPa |
| `T_K` | Temperatura | K |
| `p_Pa` | Presión | Pa |
| `u_ms`, `v_ms` | Viento zonal y meridional | m/s |
| `omega_pas` | Velocidad vertical isobárica | Pa/s |
| `Phi` | Geopotencial | m²/s² |
| `q` | Humedad específica | kg/kg |

---

### `solver_fd.py`

Resuelve las ecuaciones primitivas mediante **diferencias finitas con esquema upwind**
y **Euler explícito** en tiempo.

| | Detalle |
|--|---------|
| **Input** | `condiciones_iniciales_fd.csv` |
| **Output** | `resultados_fd.csv` (1764 filas, misma estructura) |

**Parámetros numéricos:**

| Parámetro | Valor | Justificación |
|-----------|-------|---------------|
| Δt | 600 s | CFL ≈ 0.50 con U_max=50 m/s, Δx=64 km |
| T_END | 3600 s | 1 hora de simulación (6 pasos) |
| ν (difusión) | 10⁵ m²/s | Estabilización numérica (τ_difusión ≈ 11 h) |

**Condiciones de contorno:**

| Frontera | Variable | Tipo |
|----------|----------|------|
| 4 paredes laterales | u, v, T | Dirichlet (valores CI fijos) |
| Tope 300 hPa | u, v | Dirichlet (CI fijo) |
| Tope 300 hPa | T | Libre (T puede evolucionar) |
| Tope 300 hPa | ω | = 0 (tapa rígida) |
| 1000 hPa | ω | Libre (calculado por continuidad) |
| Superficie k=0 | Todos | BC fija (no evoluciona en el solver) |

**Diagnósticos cada paso temporal:**

1. **ω**: integra `∂ω/∂p = −∇·V` desde k=5 (300 hPa, ω=0) hacia k=0 (1000 hPa)
2. **Φ**: integra hidrostático `Φₖ = Φₖ₋₁ − Rd·T̄·ln(pₖ/pₖ₋₁)` desde k=0 (fijo) hacia k=5

---

### `plot_resultados_fd.py`

Genera mapas en proyección Lambert Conformal para cada nivel.

```
python3 plot_resultados_fd.py           # todos los niveles + superficie
python3 plot_resultados_fd.py 850 500   # sólo 850 y 500 hPa
python3 plot_resultados_fd.py sfc       # sólo superficie
```

| | Detalle |
|--|---------|
| **Input** | `condiciones_iniciales_fd.csv`, `resultados_fd.csv` |
| **Output** | `resultado_fd_superficie.png`, `resultado_fd_nivel_{1000,925,850,700,500,300}hPa.png` |

**Subplots por figura isobárica:**

| Posición | Variable | Colormap |
|----------|----------|----------|
| [0,0] | ΔT anomalía temperatura (°C) | RdBu_r (divergente) |
| [0,1] | Velocidad \|V\| (m/s) + quiver | YlOrRd |
| [1,0] | Geopotencial Φ/g (m) | viridis |
| [1,1] | Velocidad vertical ω (Pa/s) | BrBG_r (divergente) |

**Subplots figura superficie:**

| Posición | Variable |
|----------|----------|
| [0,0] | ΔT_2m (siempre 0 — BC fija) |
| [0,1] | T_2m inicial con flechas viento 10m |
| [1,0] | Presión de superficie p_sfc (hPa) |
| [1,1] | Elevación del terreno (m) |

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
