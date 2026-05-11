"""
eumet_descargar.py
Descarga el producto MSG SEVIRI más reciente y genera imágenes PNG
recortadas sobre la Península Ibérica en varios canales.

Uso:
    python eumet_descargar.py                  # último producto disponible
    python eumet_descargar.py --horas 6        # buscar en últimas 6h
    python eumet_descargar.py --idx 2          # segundo producto más reciente

Salida: PNG por canal en el directorio de trabajo
"""

import os
import argparse
import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")   # backend sin ventana — evita crash tkinter

from dotenv import load_dotenv
import eumdac
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

load_dotenv()
KEY    = os.getenv("EUMETSAT_KEY")
SECRET = os.getenv("EUMETSAT_SECRET")

COLLECTION_ID = "EO:EUM:DAT:MSG:HRSEVIRI"
DOWNLOAD_DIR  = Path("datos_sat")

# Extensión del dominio (Península Ibérica + margen)
LON_MIN, LON_MAX = -11.0,  6.0
LAT_MIN, LAT_MAX =  34.5, 45.0

PLATE   = ccrs.PlateCarree()
LAMBERT = ccrs.LambertConformal(
    central_longitude=-2.0, central_latitude=40.0,
    standard_parallels=(35, 45))

# Canales SEVIRI a visualizar y sus propiedades
CANALES = {
    "VIS006":  dict(cmap="gray",       label="Reflectancia VIS 0.6 µm [%]",  pct=(2, 98)),
    "IR_108":  dict(cmap="RdYlBu_r",   label="TB IR 10.8 µm [K]",            pct=(2, 98)),
    "WV_062":  dict(cmap="Blues_r",    label="TB Vapor de Agua 6.2 µm [K]",  pct=(2, 98)),
    "IR_039":  dict(cmap="hot_r",      label="TB IR 3.9 µm [K]",             pct=(5, 95)),
}


# ── Autenticación ─────────────────────────────────────────────────────────────

def conectar():
    token = eumdac.AccessToken((KEY, SECRET))
    print(f"Token: válido hasta {token.expiration}")
    return eumdac.DataStore(token)


# ── Búsqueda ──────────────────────────────────────────────────────────────────

def buscar_productos(datastore, horas=3):
    col   = datastore.get_collection(COLLECTION_ID)
    t_fin = datetime.datetime.now(datetime.timezone.utc)
    t_ini = t_fin - datetime.timedelta(hours=horas)
    prods = list(col.search(dtstart=t_ini, dtend=t_fin))
    prods.sort(key=lambda p: str(p), reverse=True)
    print(f"Productos encontrados: {len(prods)}")
    for i, p in enumerate(prods[:5]):
        try:
            print(f"  [{i}] {p}  ({p.sensing_start:%H:%M UTC})")
        except Exception:
            print(f"  [{i}] {p}")
    return prods


# ── Descarga ──────────────────────────────────────────────────────────────────

def descargar(producto, directorio: Path) -> list[Path]:
    directorio.mkdir(parents=True, exist_ok=True)
    archivos = []
    entries  = list(producto.entries)
    for entry in entries:
        destino = directorio / entry
        if destino.exists():
            print(f"  Caché: {entry}")
            archivos.append(destino)
            continue
        print(f"  Descargando {entry}…", end=" ", flush=True)
        with producto.open(entry=entry) as src, open(destino, "wb") as dst:
            bloque = src.read(1024 * 1024)
            total  = 0
            while bloque:
                dst.write(bloque)
                total  += len(bloque)
                bloque  = src.read(1024 * 1024)
        print(f"{total // (1024*1024):.0f} MB")
        archivos.append(destino)
    return archivos


# ── Lectura con satpy ─────────────────────────────────────────────────────────

def leer_satpy(archivos: list[Path]):
    from satpy import Scene

    nat_files = [str(f) for f in archivos if f.suffix == ".nat"]
    if not nat_files:
        raise SystemExit("No se encontró archivo .nat")

    print(f"\nLeyendo con satpy: {nat_files[0]}")
    scn = Scene(reader="seviri_l1b_native", filenames=nat_files)

    disponibles = scn.available_dataset_names()
    canales_leer = [c for c in CANALES if c in disponibles]
    print(f"Canales disponibles: {len(disponibles)}  |  A leer: {canales_leer}")

    scn.load(canales_leer)
    return scn, canales_leer


# ── Recorte al dominio ────────────────────────────────────────────────────────

def recortar(scn):
    from pyresample import create_area_def

    area = create_area_def(
        "iberia",
        {"proj": "lcc", "lon_0": -2.0, "lat_0": 40.0,
         "lat_1": 35.0, "lat_2": 45.0, "ellps": "WGS84"},
        area_extent=[-900_000, -600_000, 800_000, 600_000],
        resolution=3000,   # 3 km por pixel
    )
    return scn.resample(area)


# ── Visualización ─────────────────────────────────────────────────────────────

def plot_canal(scn_crop, canal, timestamp_str):
    cfg = CANALES[canal]
    try:
        data = scn_crop[canal].values.astype(float)
    except KeyError:
        print(f"  Canal {canal} no disponible tras recorte.")
        return

    data[data < -1e10] = np.nan

    fig, ax = plt.subplots(figsize=(10, 7),
                           subplot_kw={"projection": LAMBERT})
    ax.add_feature(cfeature.OCEAN.with_scale("50m"),
                   facecolor="#D0E8F5", zorder=0)
    ax.add_feature(cfeature.LAND.with_scale("50m"),
                   facecolor="#F5F0E8", zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"),
                   linewidth=0.8, zorder=10)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"),
                   linewidth=0.5, linestyle="--", zorder=10)
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=PLATE)

    vmin = float(np.nanpercentile(data, cfg["pct"][0]))
    vmax = float(np.nanpercentile(data, cfg["pct"][1]))

    # Coordenadas lon/lat del área recortada
    lons, lats = scn_crop[canal].attrs.get("area").get_lonlats()

    img = ax.pcolormesh(lons, lats, data,
                        cmap=cfg["cmap"],
                        norm=mcolors.Normalize(vmin=vmin, vmax=vmax),
                        transform=PLATE, shading="auto", zorder=5)

    plt.colorbar(img, ax=ax, orientation="horizontal",
                 pad=0.04, shrink=0.85, label=cfg["label"])
    ax.set_title(f"MSG SEVIRI — {canal}  |  {timestamp_str}",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()

    fname = f"seviri_{canal}_{timestamp_str.replace(' ','_').replace(':','')}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  ✓ {fname}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--horas", type=int, default=3)
    parser.add_argument("--idx",   type=int, default=0,
                        help="Índice del producto (0=más reciente)")
    args = parser.parse_args()

    if not KEY or not SECRET:
        raise SystemExit("ERROR: añade EUMETSAT_KEY y EUMETSAT_SECRET en .env")

    ds   = conectar()
    prods = buscar_productos(ds, args.horas)
    if not prods:
        raise SystemExit("Sin productos disponibles")

    prod = prods[args.idx]
    try:
        ts = prod.sensing_start.strftime("%Y-%m-%d %H:%M UTC")
    except Exception:
        ts = str(prod)
    print(f"\nDescargando producto: {prod}  ({ts})")

    archivos = descargar(prod, DOWNLOAD_DIR)

    print("\nProcesando con satpy…")
    scn, canales = leer_satpy(archivos)

    print("Recortando al dominio ibérico…")
    scn_crop = recortar(scn)

    print("Generando figuras…")
    for canal in canales:
        plot_canal(scn_crop, canal, ts)

    print("\nHecho.")


if __name__ == "__main__":
    main()
