"""
condiciones_iniciales.py
Construye el campo 3D inicial sobre la malla rectangular (21×12×6 niveles isobáricos)
a partir de data.db (Open-Meteo).

Niveles isobáricos (k=1..6 en el CSV, kk=0..5 en el solver):
  k=1 : 1000 hPa
  k=2 :  925 hPa
  k=3 :  850 hPa
  k=4 :  700 hPa
  k=5 :  500 hPa
  k=6 :  300 hPa

Salida: condiciones_iniciales_fd.csv
  node_id, i, j, k, lon, lat, level_hPa, T_K, p_Pa, u_ms, v_ms, omega_pas, Phi, q
"""

import sqlite3
import math
import csv
import json
from pathlib import Path

DB_FILE = Path("data.db")
OUT_CSV = Path("condiciones_iniciales_fd.csv")

g   = 9.80665
NIVELES = [1000, 925, 850, 700, 500, 300]   # k=1..6


# ---------------------------------------------------------------------------
# Carga de datos
# ---------------------------------------------------------------------------

def cargar_nodos(con: sqlite3.Connection) -> dict:
    rows = con.execute(
        "SELECT id, i, j, lon, lat, elev_m, is_land FROM nodes ORDER BY id"
    ).fetchall()
    return {r[0]: {"id": r[0], "i": r[1], "j": r[2],
                   "lon": r[3], "lat": r[4],
                   "elev_m": r[5], "is_land": r[6]} for r in rows}


def cargar_niveles(con: sqlite3.Connection) -> dict:
    """Devuelve dict (node_id, level_hPa) → datos verticales."""
    ts = con.execute("SELECT MAX(timestamp) FROM pressure_levels").fetchone()[0]
    rows = con.execute(
        "SELECT node_id, level_hPa, T_C, geop_m, u_ms, v_ms, omega_pas, dewp_C, q "
        "FROM pressure_levels WHERE timestamp=?", (ts,)
    ).fetchall()
    return {(r[0], r[1]): {"T_C": r[2], "geop_m": r[3],
                            "u_ms": r[4], "v_ms": r[5],
                            "omega_pas": r[6], "dewp_C": r[7],
                            "q": r[8],
                            "Phi": g * r[3] if r[3] is not None else None
                            } for r in rows}


# ---------------------------------------------------------------------------
# Interpolación espacial para nodos sin datos
# ---------------------------------------------------------------------------

def interpolar_nivel(node_id: int, nivel: int, nodos: dict,
                     niveles: dict) -> dict | None:
    """Nearest-neighbor para nodos sin datos de un nivel isobárico."""
    nodo = nodos[node_id]
    best_dist = float("inf")
    best = None
    for (nid, lev), datos in niveles.items():
        if lev != nivel:
            continue
        n = nodos[nid]
        d = (n["lon"] - nodo["lon"])**2 + (n["lat"] - nodo["lat"])**2
        if d < best_dist:
            best_dist = d
            best = datos
    return best


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    con = sqlite3.connect(DB_FILE)
    nodos   = cargar_nodos(con)
    niveles = cargar_niveles(con)
    con.close()

    print(f"Nodos: {len(nodos)}  |  Niveles: {len(niveles)}")

    rows_out = []

    for node_id in sorted(nodos.keys()):
        nodo = nodos[node_id]

        # k=1..6 (niveles isobáricos)
        # Φ = g × geop_m de Open-Meteo directamente (evita artefactos orográficos)
        for k, nivel in enumerate(NIVELES, 1):
            p_nivel_Pa = nivel * 100.0

            if (node_id, nivel) in niveles:
                lv = niveles[(node_id, nivel)]
            else:
                lv = interpolar_nivel(node_id, nivel, nodos, niveles)
                if lv is None:
                    lv = {"T_C": 0.0, "u_ms": 0.0, "v_ms": 0.0,
                          "omega_pas": 0.0, "q": 0.0, "Phi": 0.0}
                    print(f"  [WARN] nodo {node_id} nivel {nivel} sin datos")

            T_K = (lv["T_C"] + 273.15) if lv["T_C"] is not None else 273.15
            u   = lv["u_ms"]      if lv["u_ms"]      is not None else 0.0
            v   = lv["v_ms"]      if lv["v_ms"]      is not None else 0.0
            om  = lv["omega_pas"] if lv["omega_pas"]  is not None else 0.0
            q   = lv["q"]         if lv["q"]          is not None else 0.0
            phi = lv.get("Phi") or 0.0

            rows_out.append({
                "node_id":    node_id,
                "i":          nodo["i"],
                "j":          nodo["j"],
                "k":          k,
                "lon":        nodo["lon"],
                "lat":        nodo["lat"],
                "level_hPa":  nivel,
                "T_K":        round(T_K, 4),
                "p_Pa":       p_nivel_Pa,
                "u_ms":       round(u, 4),
                "v_ms":       round(v, 4),
                "omega_pas":  round(om, 6),
                "Phi":        round(phi, 2),
                "q":          round(q, 8),
            })

    # Escritura CSV
    campos = ["node_id", "i", "j", "k", "lon", "lat", "level_hPa",
              "T_K", "p_Pa", "u_ms", "v_ms", "omega_pas", "Phi", "q"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(rows_out)

    # Guarda metadatos para el plot (timestamp CI + duración simulación)
    con2 = sqlite3.connect(DB_FILE)
    ts_ci = con2.execute("SELECT MAX(timestamp) FROM pressure_levels").fetchone()[0]
    con2.close()
    meta = {"timestamp_ci": ts_ci, "T_END_s": 3600}
    with open("ci_meta.json", "w") as mf:
        json.dump(meta, mf)

    n_niveles = len(NIVELES)
    print(f"\n✓ {OUT_CSV}  →  {len(rows_out)} filas  ({len(nodos)} nodos × {n_niveles} niveles)")
    T_vals   = [r["T_K"]  for r in rows_out]
    Phi_vals = [r["Phi"]  for r in rows_out]
    spd      = [math.hypot(r["u_ms"], r["v_ms"]) for r in rows_out]
    print(f"  T  rango: {min(T_vals):.1f} – {max(T_vals):.1f} K")
    print(f"  Φ  rango: {min(Phi_vals):.0f} – {max(Phi_vals):.0f} m²/s²")
    print(f"  |V| máx: {max(spd):.1f} m/s")
    print("\nSiguiente paso:")
    print("  python solver_fd.py")


if __name__ == "__main__":
    main()
