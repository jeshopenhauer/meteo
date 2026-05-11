"""
eumet_auth.py
Verifica las credenciales EUMETSAT y muestra las colecciones disponibles.

Uso:
    python eumet_auth.py

Requiere .env con:
    EUMETSAT_KEY=...
    EUMETSAT_SECRET=...
"""

import os
from dotenv import load_dotenv
import eumdac

load_dotenv()

KEY    = os.getenv("EUMETSAT_KEY")
SECRET = os.getenv("EUMETSAT_SECRET")

if not KEY or not SECRET:
    raise SystemExit(
        "ERROR: crea el archivo .env con EUMETSAT_KEY y EUMETSAT_SECRET"
    )

print("Autenticando…")
token = eumdac.AccessToken((KEY, SECRET))
print(f"  Token válido hasta: {token.expiration}")

datastore = eumdac.DataStore(token)

print("\nColecciones disponibles (primeras 30):")
for i, col in enumerate(datastore.collections):
    print(f"  {col}")
    if i >= 29:
        print("  …(más colecciones disponibles)")
        break
