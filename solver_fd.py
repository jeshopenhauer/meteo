"""
solver_fd.py — FEniCS forma variacional + Euler implícito (Newton)
Ejecutar con:
  PKG_CONFIG_PATH=/home/almo/miniconda3/envs/fenics_env/lib/pkgconfig \
  /home/almo/miniconda3/envs/fenics_env/bin/python3 solver_fd.py

Malla 3D : 61×35×6 = 12810 nodos  (i_lon, j_lat, p_hPa)
Espacio  : MixedElement CG1 × 3  (u, v, T)
Tiempo   : Euler implícito Δt=300s, Newton con MUMPS
PGF      : forma fuerte directa ∫(∂Φ/∂x)·φ·J dx — sin integral de borde
ω        : diagnóstico explícito cada paso — integra ∂ω/∂p = −∇·V desde tope
Φ        : diagnóstico hidrostático explícito cada paso
BCs      : ninguna — fronteras abiertas; unicidad garantizada por (u−u⁰)/Δt

Escalado variacional (malla índice → físico, Jacobiano J = DX·DY·dp_c):
  tiempo/Coriolis   : × J
  advección-i       : × DY·dp_c        (DX cancela con 1/DX del gradiente)
  advección-j       : × DX·dp_c        (DY cancela con 1/DY del gradiente)
  advección-p       : × DX·DY          (dp_c cancela con 1/dp_c de ∂/∂p)
  adiabático κ·T/p  : × J              (dp_c NO cancela; κ = Rd/Cp ≈ 0.286)
"""

import sys, csv, numpy as np
from pathlib import Path

try:
    from fenics import *
    set_log_level(LogLevel.ERROR)
except ImportError:
    sys.exit("Ejecutar con: "
             "/home/almo/miniconda3/envs/fenics_env/bin/python3 solver_fd.py")

# ── Constantes ────────────────────────────────────────────────────────────────
Rd  = 287.058
Cp  = 1004.0
g   = 9.80665
Om  = 7.2921e-5
RE  = 6_371_000.0

# ── Grid ──────────────────────────────────────────────────────────────────────
LON_MIN, DLON = -10.0, 0.25
LAT_MIN, DLAT =  35.5, 0.25
NI, NJ = 61, 35
LONS = np.array([LON_MIN + i*DLON for i in range(NI)])
LATS = np.array([LAT_MIN + j*DLAT for j in range(NJ)])

LEVELS_hPa = np.array([1000., 925., 850., 700., 500., 300.])
LEVELS_Pa  = LEVELS_hPa * 100.
NK = len(LEVELS_hPa)

DX_arr = RE * np.cos(np.radians(LATS)) * np.radians(DLON)   # (NJ,) m
DY     = RE * np.radians(DLAT)                               # m
F_arr  = 2.*Om * np.sin(np.radians(LATS))                    # (NJ,) s⁻¹
DP_FAC = 100.0    # Pa por hPa

DT    = 60.0
T_END = 3600.0
NSTEP = int(T_END / DT)

CI_CSV  = Path("condiciones_iniciales_fd.csv")
OUT_CSV = Path("resultados_fd.csv")


# ── Malla 3D: (i, j, p_hPa) ──────────────────────────────────────────────────
# Vértice (i,j,k) → índice (j*NI + i)*NK + k

def build_mesh():
    horiz_tris = []
    for j in range(NJ - 1):
        for i in range(NI - 1):
            v00 = j*NI + i;      v10 = j*NI + (i+1)
            v01 = (j+1)*NI + i;  v11 = (j+1)*NI + (i+1)
            horiz_tris.append([v00, v10, v11])
            horiz_tris.append([v00, v11, v01])
    horiz_tris = np.array(horiz_tris)
    n_cells = len(horiz_tris) * (NK - 1) * 3
    n_vert  = NI * NJ * NK

    mesh = Mesh()
    ed   = MeshEditor()
    ed.open(mesh, "tetrahedron", 3, 3)
    ed.init_vertices(n_vert)
    ed.init_cells(n_cells)

    for j in range(NJ):
        for i in range(NI):
            nid = j*NI + i
            for k, p in enumerate(LEVELS_hPa):
                ed.add_vertex(int(nid*NK + k), Point(float(i), float(j), float(p)))

    cidx = 0
    for tri in horiz_tris:
        a, b, c = tri
        for k in range(NK - 1):
            v = [int(x*NK + k) for x in (a,b,c)] + [int(x*NK + k+1) for x in (a,b,c)]
            ed.add_cell(int(cidx),   [v[0], v[1], v[2], v[5]])
            ed.add_cell(int(cidx+1), [v[0], v[1], v[5], v[4]])
            ed.add_cell(int(cidx+2), [v[0], v[4], v[5], v[3]])
            cidx += 3

    ed.close()
    return mesh


print("Construyendo malla 3D…", flush=True)
mesh = build_mesh()

# Espacio mixto: (u, v, T)
elem = MixedElement([FiniteElement("CG", mesh.ufl_cell(), 1)] * 3)
W    = FunctionSpace(mesh, elem)

V_sc = FunctionSpace(mesh, "CG", 1)
v2d  = vertex_to_dof_map(V_sc)
def vert(i, j, k): return int((j*NI + i)*NK + k)


# ── Expresiones dependientes de posición ─────────────────────────────────────

class DxE(UserExpression):
    def eval(self, v, x):
        j = min(max(int(round(x[1])), 0), NJ-1); v[0] = float(DX_arr[j])
    def value_shape(self): return ()

class FcE(UserExpression):
    def eval(self, v, x):
        j = min(max(int(round(x[1])), 0), NJ-1); v[0] = float(F_arr[j])
    def value_shape(self): return ()

dx_e   = DxE(degree=0)
fc_e   = FcE(degree=0)
DY_c   = Constant(DY)
dp_c   = Constant(DP_FAC)
dt_c   = Constant(DT)
RdCp_c = Constant(Rd / Cp)     # κ = Rd/Cp ≈ 0.2857


# ── BCs ───────────────────────────────────────────────────────────────────────

class LateralBC(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], 0.) or near(x[0], float(NI-1))
                                or near(x[1], 0.) or near(x[1], float(NJ-1)))

class TopBC(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 300.)

lat_bc = LateralBC()
top_bc = TopBC()


# ── E/S ───────────────────────────────────────────────────────────────────────

def leer_ci():
    u=np.zeros((NI,NJ,NK)); v=np.zeros((NI,NJ,NK))
    T=np.full((NI,NJ,NK),273.15); Phi=np.zeros((NI,NJ,NK))
    q=np.zeros((NI,NJ,NK)); nid={}
    with open(CI_CSV) as f:
        for row in csv.DictReader(f):
            kk=int(row["k"])-1; ii,jj=int(row["i"]),int(row["j"])
            if kk<0 or kk>=NK: continue
            u[ii,jj,kk]=float(row["u_ms"]); v[ii,jj,kk]=float(row["v_ms"])
            T[ii,jj,kk]=float(row["T_K"]);  Phi[ii,jj,kk]=float(row["Phi"])
            q[ii,jj,kk]=float(row["q"]);    nid[(ii,jj)]=int(row["node_id"])
    return u,v,T,Phi,q,nid


def arr_to_fn(arr3d, fn_scalar):
    vals = np.zeros(V_sc.dim())
    for k in range(NK):
        for j in range(NJ):
            for i in range(NI):
                vals[v2d[vert(i,j,k)]] = arr3d[i,j,k]
    fn_scalar.vector().set_local(vals)
    fn_scalar.vector().apply("insert")


def fn_to_arr(fn_scalar):
    vals = fn_scalar.vector().get_local()
    out  = np.zeros((NI,NJ,NK))
    for k in range(NK):
        for j in range(NJ):
            for i in range(NI):
                out[i,j,k] = vals[v2d[vert(i,j,k)]]
    return out


# ── Diagnósticos explícitos ───────────────────────────────────────────────────

def diag_omega(u_arr, v_arr):
    """ω desde continuidad: ∂ω/∂p = −(∂u/∂x + ∂v/∂y), integrado desde tope."""
    om = np.zeros((NI, NJ, NK))
    for k in range(NK-2, -1, -1):
        # Divergencia a nivel medio k+0.5 (promedio de k y k+1)
        um = 0.5*(u_arr[:,:,k] + u_arr[:,:,k+1])
        vm = 0.5*(v_arr[:,:,k] + v_arr[:,:,k+1])

        # ∂u/∂x (central en interior, unilateral en bordes)
        du = np.zeros((NI, NJ))
        du[1:-1,:] = (um[2:,:] - um[:-2,:]) / 2.0
        du[0,:]    = um[1,:] - um[0,:]
        du[-1,:]   = um[-1,:] - um[-2,:]
        div_u = du / DX_arr[np.newaxis, :]    # DX_arr[j] broadcast sobre i

        # ∂v/∂y (central en interior, unilateral en bordes)
        dv = np.zeros((NI, NJ))
        dv[:,1:-1] = (vm[:,2:] - vm[:,:-2]) / 2.0
        dv[:,0]    = vm[:,1] - vm[:,0]
        dv[:,-1]   = vm[:,-1] - vm[:,-2]
        div_v = dv / DY

        dp = LEVELS_Pa[k] - LEVELS_Pa[k+1]       # > 0 siempre
        om[:,:,k] = om[:,:,k+1] - dp * (div_u + div_v)
    return om


def diag_phi(T_arr, Phi_bot):
    """Φ hidrostático: Φₖ = Φₖ₋₁ − Rd·T̄·ln(pₖ/pₖ₋₁), desde Φ(1000 hPa) fijo."""
    Phi = np.zeros((NI,NJ,NK)); Phi[:,:,0] = Phi_bot
    for k in range(1, NK):
        Tm = 0.5*(T_arr[:,:,k-1]+T_arr[:,:,k])
        Phi[:,:,k] = Phi[:,:,k-1] - Rd*Tm*np.log(LEVELS_Pa[k]/LEVELS_Pa[k-1])
    return Phi


# ── Forma variacional (3 incógnitas: u, v, T) ─────────────────────────────────

def build_form(sol, U_n, Om_fn, Phi_fn):
    """
    Euler implícito para (u, v, T). Sin condiciones Dirichlet laterales.
    ω y Φ son externos — diagnósticos del paso previo.

    PGF en forma FUERTE (sin IBP):
      ∫(∂Φ/∂x)·φ dΩ = ∫Phi_fn.dx(0)·φ·DY·dp dΩ
      → No genera integral de borde. No requiere Dirichlet en contornos laterales.
      → La unicidad está garantizada por el término de masa (u−u⁰)/Δt.

    Término adiabático: −ω·(Rd/Cp)·T/p  → escala con J (no con DX·DY)
    """
    u, v, T   = split(sol)
    u0,v0,T0  = split(U_n)
    phi_u, phi_v, phi_T = TestFunctions(W)

    p_Pa = SpatialCoordinate(mesh)[2] * dp_c   # p_hPa × 100 → Pa
    J    = dx_e * DY_c * dp_c                  # DX · DY · 100

    # ── Momento-x ─────────────────────────────────────────────────────────────
    F_u = ((u-u0)/dt_c * phi_u * J
           + u*u.dx(0) * phi_u * DY_c * dp_c
           + v*u.dx(1) * phi_u * dx_e * dp_c
           + Om_fn*u.dx(2) * phi_u * dx_e * DY_c
           - fc_e*v    * phi_u * J
           + Phi_fn.dx(0) * phi_u * DY_c * dp_c    # PGF-x forma fuerte
          ) * dx

    # ── Momento-y ─────────────────────────────────────────────────────────────
    F_v = ((v-v0)/dt_c * phi_v * J
           + u*v.dx(0) * phi_v * DY_c * dp_c
           + v*v.dx(1) * phi_v * dx_e * dp_c
           + Om_fn*v.dx(2) * phi_v * dx_e * DY_c
           + fc_e*u    * phi_v * J
           + Phi_fn.dx(1) * phi_v * dx_e * dp_c    # PGF-y forma fuerte
          ) * dx

    # ── Termodinámica ─────────────────────────────────────────────────────────
    F_T = ((T-T0)/dt_c * phi_T * J
           + u*T.dx(0) * phi_T * DY_c * dp_c
           + v*T.dx(1) * phi_T * dx_e * dp_c
           + Om_fn*T.dx(2) * phi_T * dx_e * DY_c       # ω·∂T/∂p_hPa · DX·DY
           - Om_fn*RdCp_c*T/p_Pa * phi_T * J           # −ω·κT/p · J
          ) * dx

    return F_u + F_v + F_T


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("="*65)
    print("SOLVER FEniCS — forma variacional (IBP) + Euler implícito")
    print(f"  incógnitas: (u,v,T)   ω,Φ diagnósticos explícitos")
    print(f"  Grid: {NI}×{NJ}×{NK}   DOF: {W.dim()}   Δt={DT:.0f}s")
    print("="*65)

    u0a,v0a,T0a,Phi0a,q,nid = leer_ci()
    Phi_bot = Phi0a[:,:,0].copy()

    u_ci  = Function(V_sc); arr_to_fn(u0a,  u_ci)
    v_ci  = Function(V_sc); arr_to_fn(v0a,  v_ci)
    T_ci  = Function(V_sc); arr_to_fn(T0a,  T_ci)

    Om_n  = Function(V_sc)      # ω del paso anterior
    Phi_n = Function(V_sc)      # Φ del paso anterior

    # ω inicial: consistente con nuestra fórmula FD aplicada a los CI
    # (evita el salto en paso 1 respecto al ω de Open-Meteo)
    om0a = diag_omega(u0a, v0a)
    arr_to_fn(om0a,   Om_n)
    arr_to_fn(Phi0a,  Phi_n)

    U_n = Function(W); sol = Function(W)
    for fn_sc, sub in ((u_ci,0),(v_ci,1),(T_ci,2)):
        assign(U_n.sub(sub), fn_sc)
    sol.assign(U_n)

    # Sin Dirichlet laterales: la PGF en forma fuerte elimina la necesidad.
    # El término (u-u0)/Δt garantiza unicidad sin BC adicionales.
    bcs = []

    F    = build_form(sol, U_n, Om_n, Phi_n)
    J_ac = derivative(F, sol)

    solver_prm = {
        "nonlinear_solver": "newton",
        "newton_solver": {
            "linear_solver":           "mumps",
            "maximum_iterations":      30,
            "relaxation_parameter":    0.7,
            "relative_tolerance":      1e-6,
            "absolute_tolerance":      1e-8,
            "error_on_nonconvergence": False,
        }
    }

    print("\nCompilando forma FEniCS (primera vez ~30 s)…", flush=True)
    print(f"\n{'paso':>5}  {'t(h)':>6}  {'T_min':>8}  {'T_max':>8}  "
          f"{'ΔT_min':>7}  {'|V|_max':>8}  {'ω_max':>8}")
    print("-"*65)

    T_arr = T0a.copy(); u_arr = u0a.copy()
    v_arr = v0a.copy(); w_arr = om0a.copy()

    for step in range(1, NSTEP+1):
        solve(F == 0, sol, bcs, J=J_ac, solver_parameters=solver_prm)

        parts = sol.split(deepcopy=True)
        u_arr = fn_to_arr(parts[0])
        v_arr = fn_to_arr(parts[1])
        T_arr = fn_to_arr(parts[2])

        # Diagnósticos explícitos para el siguiente paso
        w_arr   = diag_omega(u_arr, v_arr)
        Phi_new = diag_phi(T_arr, Phi_bot)

        arr_to_fn(w_arr,   Om_n)
        arr_to_fn(Phi_new, Phi_n)

        U_n.assign(sol)

        if np.isnan(T_arr).any() or T_arr.min() < 150 or T_arr.max() > 400:
            print(f"\n[ABORT] Inestabilidad en paso {step}"); break

        dT = T_arr - T0a
        print(f"{step:>5}  {step*DT/3600:>6.2f}  {T_arr.min():>8.2f}  "
              f"{T_arr.max():>8.2f}  {dT.min():>7.2f}  "
              f"{np.sqrt(u_arr**2+v_arr**2).max():>8.2f}  "
              f"{np.abs(w_arr).max():>8.5f}", flush=True)

    Phi_final = diag_phi(T_arr, Phi_bot)
    campos = ["node_id","i","j","k","lon","lat","level_hPa",
              "T_K","p_Pa","u_ms","v_ms","omega_pas","Phi","q"]
    with open(OUT_CSV,"w",newline="") as f:
        wr = csv.DictWriter(f, fieldnames=campos); wr.writeheader()
        for ii in range(NI):
            for jj in range(NJ):
                for kk in range(NK):
                    wr.writerow({
                        "node_id": nid.get((ii,jj), ii*NJ+jj),
                        "i":ii, "j":jj, "k":kk+1,
                        "lon":    round(LONS[ii],4),
                        "lat":    round(LATS[jj],4),
                        "level_hPa": int(LEVELS_hPa[kk]),
                        "T_K":       round(float(T_arr[ii,jj,kk]),4),
                        "p_Pa":      float(LEVELS_Pa[kk]),
                        "u_ms":      round(float(u_arr[ii,jj,kk]),4),
                        "v_ms":      round(float(v_arr[ii,jj,kk]),4),
                        "omega_pas": round(float(w_arr[ii,jj,kk]),6),
                        "Phi":       round(float(Phi_final[ii,jj,kk]),2),
                        "q":         round(float(q[ii,jj,kk]),8),
                    })

    dT_all = T_arr - T0a
    print(f"\n✓ {OUT_CSV}  {NI*NJ*NK} filas")
    print(f"  ΔT : {dT_all.min():.3f} – {dT_all.max():.3f} K")
    print(f"  |V|: {np.sqrt(u_arr**2+v_arr**2).max():.2f} m/s")
    print(f"  ω_max: {np.abs(w_arr).max():.4f} Pa/s")
    print("  python3 plot_resultados_fd.py")


if __name__ == "__main__":
    main()
