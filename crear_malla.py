"""
crear_malla.py
Genera la malla rectangular regular sobre la Península Ibérica.

Grid: lon [-10°, 5°] × lat [35.5°, 44°]  paso 0.75°
  → 21 × 12 = 252 nodos totales
  → ~90–100 nodos en tierra (máscara Natural Earth 50m)

Salida: data.db  →  tablas nodes, surface, pressure_levels
"""

import sqlite3
import numpy as np
import requests
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from pathlib import Path

# ---------------------------------------------------------------------------
# Parámetros del grid
# ---------------------------------------------------------------------------
LON_MIN, LON_MAX, DLON = -10.0,  5.0, 0.5
LAT_MIN, LAT_MAX, DLAT =  35.5, 44.0, 0.5

LONS = np.arange(LON_MIN, LON_MAX + DLON / 2, DLON)
LATS = np.arange(LAT_MIN, LAT_MAX + DLAT / 2, DLAT)
NI, NJ = len(LONS), len(LATS)          # columnas (lon), filas (lat)

DB_FILE = Path("data.db")
ELEVATION_URL = "https://api.open-meteo.com/v1/elevation"


# ---------------------------------------------------------------------------
# Máscara tierra/mar (Natural Earth 50m)
# ---------------------------------------------------------------------------

def construir_mascara_tierra() -> sgeom.base.BaseGeometry:
    shpfile = shpreader.natural_earth(
        resolution="50m", category="physical", name="land"
    )
    land_geom = unary_union([rec.geometry
                             for rec in shpreader.Reader(shpfile).records()])
    return land_geom


def es_tierra(lon: float, lat: float, land: sgeom.base.BaseGeometry) -> bool:
    return land.contains(sgeom.Point(lon, lat))


# ---------------------------------------------------------------------------
# Elevaciones en batch (Open-Meteo Elevation API)
# ---------------------------------------------------------------------------

def obtener_elevaciones(lats_list: list, lons_list: list,
                        batch: int = 100) -> list:
    import time
    elevaciones = []
    n_batches = (len(lats_list) + batch - 1) // batch
    for idx, start in enumerate(range(0, len(lats_list), batch)):
        lat_b = lats_list[start:start + batch]
        lon_b = lons_list[start:start + batch]
        for intento in range(5):
            try:
                r = requests.get(ELEVATION_URL, params={
                    "latitude":  ",".join(f"{v:.4f}" for v in lat_b),
                    "longitude": ",".join(f"{v:.4f}" for v in lon_b),
                }, timeout=60)
                r.raise_for_status()
                elevaciones.extend(r.json()["elevation"])
                print(f"  batch {idx+1}/{n_batches} OK", flush=True)
                time.sleep(100)   # pausa entre batches para evitar rate-limit
                break
            except Exception as e:
                espera = 10 * 2**intento
                print(f"  batch {idx+1} intento {intento+1} error ({e}), "
                      f"reintentando en {espera}s…", flush=True)
                time.sleep(espera)
        else:
            raise RuntimeError(f"batch {idx+1} falló tras 5 intentos")
    return elevaciones


# ---------------------------------------------------------------------------
# Base de datos
# ---------------------------------------------------------------------------

def init_db(con: sqlite3.Connection) -> None:
    con.executescript("""
        DROP TABLE IF EXISTS pressure_levels;
        DROP TABLE IF EXISTS nodes;

        CREATE TABLE nodes (
            id       INTEGER PRIMARY KEY,
            i        INTEGER NOT NULL,     -- índice columna (lon, 0=oeste)
            j        INTEGER NOT NULL,     -- índice fila    (lat, 0=sur)
            lon      REAL    NOT NULL,
            lat      REAL    NOT NULL,
            elev_m   REAL    DEFAULT 0.0,
            is_land  INTEGER DEFAULT 1     -- 1=tierra, 0=mar
        );

        CREATE TABLE pressure_levels (
            node_id   INTEGER NOT NULL REFERENCES nodes(id),
            timestamp TEXT    NOT NULL,
            level_hPa INTEGER NOT NULL,   -- 1000, 925, 850, 700, 500, 300
            T_C       REAL,               -- temperatura [°C]
            geop_m    REAL,               -- altura geopotencial [m]
            u_ms      REAL,               -- viento zonal [m/s]
            v_ms      REAL,               -- viento meridional [m/s]
            omega_pas REAL,               -- velocidad vertical isobárica [Pa/s]
            dewp_C    REAL,               -- punto de rocío [°C]
            q         REAL,               -- humedad específica [kg/kg]
            PRIMARY KEY (node_id, timestamp, level_hPa)
        );
    """)
    con.commit()


# ---------------------------------------------------------------------------
# Visualización del grid
# ---------------------------------------------------------------------------

def plot_grid(nodos: list) -> None:
    PLATE = ccrs.PlateCarree()
    LAMBERT = ccrs.LambertConformal(
        central_longitude=-2.0, central_latitude=40.0,
        standard_parallels=(35, 45))

    fig, ax = plt.subplots(figsize=(14, 9),
                           subplot_kw={"projection": LAMBERT})
    ax.add_feature(cfeature.OCEAN.with_scale("50m"),
                   facecolor="#E8F4F8", zorder=0)
    ax.add_feature(cfeature.LAND.with_scale("50m"),
                   facecolor="#F5F5DC", zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"),
                   linewidth=0.6, zorder=10)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"),
                   linewidth=0.5, linestyle="--", zorder=10)
    ax.set_extent([-10.5, 5.5, 35.0, 44.5], crs=PLATE)

    lons_all = [n["lon"] for n in nodos]
    lats_all = [n["lat"] for n in nodos]

    ax.scatter(lons_all, lats_all, s=14, c="crimson", transform=PLATE, zorder=20)

    ax.set_title(
        f"Malla regular {DLON}° × {DLAT}°  —  {NI}×{NJ} = {NI*NJ} nodos",
        fontsize=13, fontweight="bold"
    )

    plt.tight_layout()
    fig.savefig("malla_rectangular.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.show()
    print("✓ malla_rectangular.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print(f"MALLA RECTANGULAR  {NI}×{NJ} = {NI*NJ} nodos")
    print(f"  lon: {LON_MIN}° a {LON_MAX}°  paso {DLON}°")
    print(f"  lat: {LAT_MIN}° a {LAT_MAX}°  paso {DLAT}°")
    print("=" * 65)

    # Máscara tierra
    print("\nCargando máscara tierra/mar…")
    land = construir_mascara_tierra()

    # Construir lista de nodos
    nodos = []
    node_id = 0
    all_lats, all_lons = [], []
    for j, lat in enumerate(LATS):
        for i, lon in enumerate(LONS):
            all_lats.append(lat)
            all_lons.append(lon)
            nodos.append({"id": node_id, "i": i, "j": j,
                          "lon": lon, "lat": lat,
                          "elev_m": 0.0,
                          "is_land": int(es_tierra(lon, lat, land))})
            node_id += 1

    n_tierra = sum(n["is_land"] for n in nodos)
    print(f"Nodos en tierra: {n_tierra}  /  mar: {len(nodos)-n_tierra}")

    # Elevaciones (Open-Meteo batch)
    print("\nDescargando elevaciones…")
    elevs = obtener_elevaciones(all_lats, all_lons)
    for n, e in zip(nodos, elevs):
        n["elev_m"] = max(0.0, float(e))      # 0 para nodos marinos

    # Guardar en DB
    print(f"\nCreando {DB_FILE}…")
    con = sqlite3.connect(DB_FILE)
    init_db(con)
    con.executemany(
        "INSERT INTO nodes (id,i,j,lon,lat,elev_m,is_land) VALUES (?,?,?,?,?,?,?)",
        [(n["id"], n["i"], n["j"], n["lon"], n["lat"],
          n["elev_m"], n["is_land"]) for n in nodos]
    )
    con.commit()
    con.close()

    print(f"✓ {DB_FILE}  →  {len(nodos)} nodos guardados")
    print(f"  Rango elevación tierra: "
          f"{min(n['elev_m'] for n in nodos if n['is_land']):.0f}–"
          f"{max(n['elev_m'] for n in nodos if n['is_land']):.0f} m")

    # Mapa
    print("\nGenerando mapa del grid…")
    plot_grid(nodos)

    print("\n" + "=" * 65)
    print("Siguiente paso:")
    print("  python datos_openmeteo.py   ← descarga campos meteorológicos")
    print("=" * 65)


if __name__ == "__main__":
    main()
