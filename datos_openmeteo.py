"""
datos_openmeteo.py
Descarga datos meteorológicos de Open-Meteo (GFS, 28 km) para cada nodo
de la malla rectangular y los almacena en data.db.

Variables:
  Superficie (k=0): temperature_2m, surface_pressure, windspeed_10m,
                    winddirection_10m, dewpoint_2m
  Niveles isobáricos: 1000 / 925 / 850 / 700 / 500 / 300 hPa
                      → T, geopotential_height, dewpoint,
                        windspeed, winddirection, vertical_velocity

Ejecución:
  python datos_openmeteo.py          ← hora actual (índice 0 del forecast)
  python datos_openmeteo.py 3        ← hora +3 del forecast
"""

import math
import sqlite3
import datetime
import sys
import requests
from pathlib import Path

DB_FILE       = Path("data.db")
FORECAST_URL  = "https://api.open-meteo.com/v1/forecast"
NIVELES_ISO   = [1000, 925, 850, 700, 500, 300]

Rd = 287.058
g  = 9.80665


# ---------------------------------------------------------------------------
# Termodinámica
# ---------------------------------------------------------------------------

def calc_q(Td_C: float, p_hPa: float) -> float | None:
    if Td_C is None or p_hPa is None:
        return None
    try:
        es = 6.112 * math.exp((17.67 * Td_C) / (Td_C + 243.5))
        return 0.622 * es / (p_hPa - 0.378 * es)
    except Exception:
        return None


def calc_uv(wspd_kmh: float, wdir_deg: float) -> tuple[float, float]:
    if wspd_kmh is None or wdir_deg is None:
        return 0.0, 0.0
    ms  = wspd_kmh / 3.6
    rad = math.radians(wdir_deg)
    return -ms * math.sin(rad), -ms * math.cos(rad)


# ---------------------------------------------------------------------------
# Open-Meteo — una petición por nodo
# ---------------------------------------------------------------------------

VARS_ISO = ",".join(
    v for grupo in [
        [f"temperature_{n}hPa"         for n in NIVELES_ISO],
        [f"geopotential_height_{n}hPa" for n in NIVELES_ISO],
        [f"dewpoint_{n}hPa"            for n in NIVELES_ISO],
        [f"windspeed_{n}hPa"           for n in NIVELES_ISO],
        [f"winddirection_{n}hPa"       for n in NIVELES_ISO],
        [f"vertical_velocity_{n}hPa"   for n in NIVELES_ISO],
    ]
    for v in grupo
)
def obtener_datos_nodo(lat: float, lon: float) -> dict:
    r = requests.get(FORECAST_URL, params={
        "latitude":      lat,
        "longitude":     lon,
        "hourly":        VARS_ISO,
        "forecast_days": 1,
        "timezone":      "UTC",
    }, timeout=20)
    r.raise_for_status()
    return r.json()["hourly"]


# ---------------------------------------------------------------------------
# Base de datos
# ---------------------------------------------------------------------------

def leer_nodos(con: sqlite3.Connection) -> list[dict]:
    rows = con.execute(
        "SELECT id, i, j, lon, lat, elev_m, is_land FROM nodes ORDER BY id"
    ).fetchall()
    return [{"id": r[0], "i": r[1], "j": r[2],
             "lon": r[3], "lat": r[4], "elev_m": r[5],
             "is_land": r[6]} for r in rows]


def limpiar_datos(con: sqlite3.Connection) -> None:
    """Vacía la tabla de datos isobáricos (no la malla)."""
    con.execute("DELETE FROM pressure_levels")
    con.commit()


def guardar_niveles(con: sqlite3.Connection, node_id: int,
                    ts: str, h: dict, hora: int) -> None:
    for n in NIVELES_ISO:
        T     = h[f"temperature_{n}hPa"][hora]
        geop  = h[f"geopotential_height_{n}hPa"][hora]
        Td    = h[f"dewpoint_{n}hPa"][hora]
        wspd  = h[f"windspeed_{n}hPa"][hora]
        wdir  = h[f"winddirection_{n}hPa"][hora]
        omega = h[f"vertical_velocity_{n}hPa"][hora]
        u, v  = calc_uv(wspd, wdir)
        q     = calc_q(Td, n)
        con.execute(
            "INSERT OR REPLACE INTO pressure_levels "
            "(node_id,timestamp,level_hPa,T_C,geop_m,u_ms,v_ms,"
            "omega_pas,dewp_C,q) "
            "VALUES (?,?,?,?,?,?,?,?,?,?)",
            (node_id, ts, n, T, geop, u, v, omega, Td, q)
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(hora: int = 0):
    con = sqlite3.connect(DB_FILE)
    nodos = leer_nodos(con)
    if not nodos:
        print("ERROR: data.db sin nodos. Ejecuta primero crear_malla.py")
        con.close()
        return

    ts = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    limpiar_datos(con)   # sobrescribe: borra todo antes de insertar

    N = len(nodos)
    n_tierra = sum(n["is_land"] for n in nodos)
    print("=" * 70)
    print(f"DESCARGA OPEN-METEO (GFS)  —  {N} nodos  ({n_tierra} tierra)")
    print(f"Timestamp: {ts}   Hora forecast: +{hora}h")
    print("=" * 70)
    print(f"\n{'#':>4}  {'ID':>5}  {'lon':>7}  {'lat':>6}  "
          f"{'T_1000(°C)':>11}  {'Z_1000(m)':>10}  {'u(m/s)':>7}")
    print("-" * 60)

    ok = err = 0
    for idx, nodo in enumerate(nodos, 1):
        for intento in range(10):          # hasta 3 reintentos por timeout
            try:
                h = obtener_datos_nodo(nodo["lat"], nodo["lon"])
                guardar_niveles(con, nodo["id"], ts, h, hora)
                con.commit()

                T    = h["temperature_1000hPa"][hora]
                geop = h["geopotential_height_1000hPa"][hora]
                u, _ = calc_uv(h["windspeed_1000hPa"][hora],
                               h["winddirection_1000hPa"][hora])
                tipo = "T" if nodo["is_land"] else "M"
                print(f"{idx:>4}  {nodo['id']:>5}{tipo}  "
                      f"{nodo['lon']:>7.2f}  {nodo['lat']:>6.2f}  "
                      f"{T:>11.1f}  {geop:>10.0f}  {u:>7.2f}")
                ok += 1
                break
            except Exception as e:
                if intento == 2:
                    print(f"{idx:>4}  {nodo['id']:>5}  "
                          f"{nodo['lon']:>7.2f}  {nodo['lat']:>6.2f}  "
                          f"ERROR tras 3 intentos: {e}")
                    err += 1

    con.close()
    print("\n" + "=" * 70)
    print(f"✓ Nodos OK: {ok}   Errores: {err}")
    print(f"✓ Timestamp guardado: {ts}")
    print("\nSiguiente paso:")
    print("  python condiciones_iniciales.py")
    print("=" * 70)


if __name__ == "__main__":
    hora = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    main(hora)
