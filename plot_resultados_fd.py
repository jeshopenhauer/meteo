"""
plot_resultados_fd.py
Visualiza los resultados del solver FD sobre la malla rectangular.
Genera un PNG por nivel isobárico con 4 subplots.
"""

import sys
import csv
import json
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path

CI_CSV  = Path("condiciones_iniciales_fd.csv")
RES_CSV = Path("resultados_fd.csv")

NI, NJ = 61, 35
LON_MIN, DLON = -10.0, 0.25
LAT_MIN, DLAT =  35.5, 0.25
LONS = np.array([LON_MIN + i * DLON for i in range(NI)])
LATS = np.array([LAT_MIN + j * DLAT for j in range(NJ)])
LEVELS_hPa = [1000, 925, 850, 700, 500, 300]

PLATE   = ccrs.PlateCarree()
LAMBERT = ccrs.LambertConformal(
    central_longitude=-2.0, central_latitude=40.0,
    standard_parallels=(35, 45))


def leer_csv(path, NK=6):
    """Devuelve arrays (NI,NJ,NK) para los 6 niveles isobáricos."""
    T   = np.full((NI, NJ, NK), np.nan)
    u   = np.full((NI, NJ, NK), np.nan)
    v   = np.full((NI, NJ, NK), np.nan)
    om  = np.full((NI, NJ, NK), np.nan)
    Phi = np.full((NI, NJ, NK), np.nan)
    with open(path) as f:
        for row in csv.DictReader(f):
            k_ci = int(row["k"])
            ii, jj = int(row["i"]), int(row["j"])
            kk = k_ci - 1
            if kk < 0 or kk >= NK:
                continue
            T  [ii, jj, kk] = float(row["T_K"])
            u  [ii, jj, kk] = float(row["u_ms"])
            v  [ii, jj, kk] = float(row["v_ms"])
            om [ii, jj, kk] = float(row["omega_pas"])
            Phi[ii, jj, kk] = float(row["Phi"])
    return T, u, v, om, Phi


def add_map_features(ax):
    ax.add_feature(cfeature.OCEAN.with_scale("50m"),
                   facecolor="#E8F4F8", zorder=0)
    ax.add_feature(cfeature.LAND.with_scale("50m"),
                   facecolor="#F5F5DC", zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"),
                   linewidth=0.7, zorder=10)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"),
                   linewidth=0.5, linestyle="--", zorder=10)
    ax.set_extent([-10.5, 5.5, 35.0, 44.5], crs=PLATE)


def pcm(ax, data, cmap, norm, title):
    """pcolormesh on regular lon/lat grid."""
    LON2D, LAT2D = np.meshgrid(LONS, LATS, indexing="ij")
    img = ax.pcolormesh(LON2D, LAT2D, data,
                        cmap=cmap, norm=norm,
                        transform=PLATE, zorder=5,
                        shading="auto")
    ax.set_title(title, fontsize=10, fontweight="bold")
    return img


def quiver_overlay(ax, u2d, v2d, stride=2):
    """Flechas de viento cada 'stride' nodos."""
    LON2D, LAT2D = np.meshgrid(LONS[::stride], LATS[::stride], indexing="ij")
    ax.quiver(LON2D, LAT2D, u2d[::stride, ::stride],
              v2d[::stride, ::stride],
              transform=PLATE, zorder=15,
              scale=500, width=0.003, color="k", alpha=0.7)


def plot_level(k, T0_k, Tf_k, u_k, v_k, om_k, Phi_k, level_hPa, tag="",
               ts_ci="", ts_res=""):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10),
                              subplot_kw={"projection": LAMBERT})
    subtitle = f"CI: {ts_ci}  →  +1h: {ts_res}" if ts_ci else ""
    fig.suptitle(
        f"Resultados solver FD — Nivel {level_hPa} hPa{tag}\n{subtitle}",
        fontsize=13, fontweight="bold", y=0.99
    )

    # --- ΔT anomalía temperatura
    dT = Tf_k - T0_k
    flat = dT[~np.isnan(dT)].ravel()
    # Escala basada en percentil 90 → resuelve variaciones de ~0.5°C sin saturar
    vmax = max(float(np.percentile(np.abs(flat), 90)), 0.5)
    norm_dT = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    ax = axes[0, 0]
    add_map_features(ax)
    dt_label = f"+1h  ({ts_res})" if ts_res else "+1h"
    img = pcm(ax, dT, "RdBu_r", norm_dT,
              f"ΔT {dt_label}  [{float(np.nanmin(dT)):.1f} / {float(np.nanmax(dT)):.1f}] °C")
    plt.colorbar(img, ax=ax, orientation="horizontal", pad=0.04,
                 label="ΔT [°C]", shrink=0.85, extend="both")

    # --- Velocidad del viento + flechas
    spd = np.sqrt(u_k**2 + v_k**2)
    norm_spd = mcolors.Normalize(0, float(np.nanmax(spd)))
    ax = axes[0, 1]
    add_map_features(ax)
    img = pcm(ax, spd, "YlOrRd", norm_spd,
              f"|V| (m/s)  [max {float(np.nanmax(spd)):.1f}]")
    quiver_overlay(ax, u_k, v_k)
    plt.colorbar(img, ax=ax, orientation="horizontal", pad=0.04,
                 label="|V| [m/s]", shrink=0.85)

    # --- Geopotencial Φ/g (metros)
    geop_m = Phi_k / 9.80665
    norm_phi = mcolors.Normalize(float(np.nanmin(geop_m)),
                                  float(np.nanmax(geop_m)))
    ax = axes[1, 0]
    add_map_features(ax)
    img = pcm(ax, geop_m, "viridis", norm_phi,
              f"Φ/g (m)  [{float(np.nanmin(geop_m)):.0f} / {float(np.nanmax(geop_m)):.0f}]")
    plt.colorbar(img, ax=ax, orientation="horizontal", pad=0.04,
                 label="Geopotencial [m]", shrink=0.85)

    # --- ω velocidad vertical isobárica
    om_abs = float(max(abs(float(np.nanmin(om_k))),
                       abs(float(np.nanmax(om_k))), 0.01))
    norm_om = mcolors.TwoSlopeNorm(vmin=-om_abs, vcenter=0, vmax=om_abs)
    ax = axes[1, 1]
    add_map_features(ax)
    img = pcm(ax, om_k, "BrBG_r", norm_om,
              f"ω (Pa/s)  [{float(np.nanmin(om_k)):.3f} / {float(np.nanmax(om_k)):.3f}]")
    plt.colorbar(img, ax=ax, orientation="horizontal", pad=0.04,
                 label="ω [Pa/s]  (+ hundimiento)", shrink=0.85)

    plt.tight_layout()
    fname = f"resultado_fd_nivel_{level_hPa}hPa.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  ✓ {fname}")


def main():
    if not CI_CSV.exists():
        print(f"ERROR: {CI_CSV} no encontrado"); return
    if not RES_CSV.exists():
        print(f"ERROR: {RES_CSV} no encontrado"); return

    # Carga metadatos de tiempo (generados por condiciones_iniciales.py)
    ts_ci = ts_res = ""
    meta_path = Path("ci_meta.json")
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        raw_ci = meta.get("timestamp_ci", "")
        t_end  = int(meta.get("T_END_s", 3600))
        if raw_ci:
            dt_ci  = datetime.datetime.fromisoformat(raw_ci.replace("Z", "+00:00"))
            dt_res = dt_ci + datetime.timedelta(seconds=t_end)
            ts_ci  = dt_ci.strftime("%Y-%m-%d %H:%M UTC")
            ts_res = dt_res.strftime("%Y-%m-%d %H:%M UTC")

    print("Leyendo campos…")
    T0, u0, v0, om0, Phi0 = leer_csv(CI_CSV)
    Tf, uf, vf, omf, Phif = leer_csv(RES_CSV)

    levels = list(sys.argv[1:]) if len(sys.argv) > 1 else [str(l) for l in LEVELS_hPa]
    level_set = set(levels)

    print("Generando figuras…")
    for kk, lhPa in enumerate(LEVELS_hPa):
        if str(lhPa) not in level_set:
            continue
        plot_level(
            kk,
            T0[:, :, kk] - 273.15,
            Tf[:, :, kk] - 273.15,
            uf[:, :, kk], vf[:, :, kk],
            omf[:, :, kk],
            Phif[:, :, kk],
            lhPa,
            ts_ci=ts_ci, ts_res=ts_res
        )

    # Resumen global
    dT_all = Tf - T0
    print(f"\nResumen global:")
    print(f"  ΔT isobárico: {float(np.nanmin(dT_all)):.2f} – {float(np.nanmax(dT_all)):.2f} K")
    print(f"  ω_max: {float(np.nanmax(np.abs(omf))):.4f} Pa/s")
    print(f"  |V|  : {float(np.nanmax(np.sqrt(uf**2+vf**2))):.1f} m/s")


if __name__ == "__main__":
    main()
