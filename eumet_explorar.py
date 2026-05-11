"""
eumet_explorar.py
Busca los productos MTG FCI más recientes sobre la Península Ibérica.
Muestra metadata de los productos disponibles y sus canales.

Uso:
    python eumet_explorar.py
    python eumet_explorar.py --horas 6      # últimas 6 horas (default: 3)
    python eumet_explorar.py --col EO:EUM:DAT:MTG:FCI-1C
"""

import os
import argparse
import datetime
import json
from dotenv import load_dotenv
import eumdac

load_dotenv()
KEY    = os.getenv("EUMETSAT_KEY")
SECRET = os.getenv("EUMETSAT_SECRET")

# Colecciones MTG/MSG de interés
COLECCIONES = {
    "MTG-FCI":  "EO:EUM:DAT:MTG:FCI-1C",
    "MSG-SEVIRI-HRV": "EO:EUM:DAT:MSG:HRSEVIRI",
    "MSG-SEVIRI-RSS": "EO:EUM:DAT:MSG:RSS",
}

# Bounding box Península Ibérica
LON_MIN, LON_MAX = -10.5,  5.5
LAT_MIN, LAT_MAX =  35.0, 44.5


def conectar():
    token = eumdac.AccessToken((KEY, SECRET))
    print(f"Token válido hasta: {token.expiration}")
    return eumdac.DataStore(token)


def explorar_coleccion(datastore, collection_id, horas=3):
    print(f"\n{'='*60}")
    print(f"Colección: {collection_id}")

    try:
        col = datastore.get_collection(collection_id)
        print(f"Título: {col.title}")
    except Exception as e:
        print(f"  ERROR accediendo a la colección: {e}")
        return

    t_fin   = datetime.datetime.now(datetime.timezone.utc)
    t_ini   = t_fin - datetime.timedelta(hours=horas)

    print(f"Buscando productos entre {t_ini:%Y-%m-%d %H:%M} y {t_fin:%Y-%m-%d %H:%M} UTC")
    print(f"Área: lon [{LON_MIN},{LON_MAX}] lat [{LAT_MIN},{LAT_MAX}]")

    try:
        productos = col.search(
            dtstart=t_ini,
            dtend=t_fin,
        )
        lista = list(productos)
    except Exception as e:
        print(f"  ERROR en búsqueda: {e}")
        return

    if not lista:
        print("  Sin productos en este intervalo.")
        return

    print(f"\nProductos encontrados: {len(lista)}")
    for p in lista[:10]:
        print(f"\n  ID  : {p}")
        try:
            print(f"  Hora: {p.sensing_start}  →  {p.sensing_end}")
        except Exception:
            pass
        try:
            entries = list(p.entries)
            print(f"  Ficheros ({len(entries)}):")
            for e in entries[:8]:
                print(f"    {e}")
        except Exception:
            pass

    if len(lista) > 10:
        print(f"\n  … y {len(lista)-10} productos más.")

    return lista


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--horas", type=int, default=3)
    parser.add_argument("--col",   type=str, default=None,
                        help="ID de colección (default: prueba todas)")
    args = parser.parse_args()

    if not KEY or not SECRET:
        raise SystemExit("ERROR: añade EUMETSAT_KEY y EUMETSAT_SECRET en .env")

    ds = conectar()

    if args.col:
        explorar_coleccion(ds, args.col, args.horas)
    else:
        for nombre, cid in COLECCIONES.items():
            print(f"\nProbando {nombre}…")
            explorar_coleccion(ds, cid, args.horas)


if __name__ == "__main__":
    main()
