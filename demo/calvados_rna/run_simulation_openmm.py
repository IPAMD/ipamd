# polyU30 OpenMM 验证脚本 — 参数全部写死，对应 demo/calvados_rna/run_simulation.py
import math
import os
import sys
import time

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit

# ========== 体系（同 run_simulation.py）==========
N_RES = 30
BOX = 30.0
CENTER = 15.0
T_K = 293.0
DT_PS = 0.01
TOTAL_NS = 2.0
SNAP = 200
RG_FRAME_START = 51
RG_FRAME_END = 200
STEPS = int(TOTAL_NS * 1000 / DT_PS)       # 200000
REPORT = STEPS // SNAP                     # 1000，共 200 帧
RG_EVERY = REPORT
GPU = 3
SEED = 99
OUT = "openmm_simulation_20_mmol_L"

# ========== 力场（Calvados_rna + Calvados2 debye, T=293K, I=20mM）==========
KAPPA = 0.463644
DH_SHIFT = 0.039390084
Q_BB = -1.315298
EPS = 0.8368

os.makedirs(OUT, exist_ok=True)

# --- 构型：直线 backbone，侧链同 calvados_rna_from_sequence ---
bb = np.array([[CENTER, CENTER, CENTER + 0.59 * i] for i in range(N_RES)])
vecs = [- (bb[1] - bb[0])] + [bb[i + 1] - bb[i] for i in range(N_RES - 1)] + [-(bb[-1] - bb[-2])]
ref = np.array([1.0, 1.0, 1.0])
rs = np.zeros((N_RES, 3))
for i in range(N_RES):
    v1 = vecs[i] / np.linalg.norm(vecs[i])
    v2 = vecs[i + 1] / np.linalg.norm(vecs[i + 1])
    if np.dot(v1, v2) > 0:
        v1 = -v1
    if np.linalg.norm(np.cross(v1, v2)) > 1e-6:
        t = (v1 + v2) / np.linalg.norm(v1 + v2)
    else:
        t = np.cross(v1, np.cross(v1, ref))
        t /= np.linalg.norm(t)
    rs[i] = bb[i] + 0.45 * t

topo = app.Topology()
chain = topo.addChain()
residue = topo.addResidue("polyU", chain)
elem = app.Element.getBySymbol("C")
pos, masses = [], []
for i in range(N_RES):
    topo.addAtom("BB", elem, residue)
    pos.append(bb[i].tolist())
    masses.append(211.107 if i == 0 else 195.108 if i == N_RES - 1 else 194.1)
    topo.addAtom("RS", elem, residue)
    pos.append(rs[i].tolist())
    masses.append(111.1)
positions = unit.Quantity(pos, unit.nanometer)
masses = np.array(masses)

# --- System ---
system = mm.System()
system.setDefaultPeriodicBoxVectors(
    mm.Vec3(BOX, 0, 0), mm.Vec3(0, BOX, 0), mm.Vec3(0, 0, BOX)
)
for m in masses:
    system.addParticle(m * unit.amu)

bonds, angles = [], []
bf = mm.HarmonicBondForce()
for i in range(N_RES):
    b, r = 2 * i, 2 * i + 1
    bf.addBond(b, r, 0.54 * unit.nanometer, 2200.0)
    bonds.append((b, r))
    if i < N_RES - 1:
        b2 = 2 * (i + 1)
        bf.addBond(b, b2, 0.59 * unit.nanometer, 1400.0)
        bonds.append((b, b2))
system.addForce(bf)

af = mm.HarmonicAngleForce()
for i in range(1, N_RES - 1):
    a, b, c = 2 * (i - 1), 2 * i, 2 * (i + 1)
    af.addAngle(a, b, c, math.pi, 4.2)
    angles.append((a, b, c))
system.addForce(af)

excl = set()
for a, b in bonds:
    excl.add((min(a, b), max(a, b)))
for i, j, k in angles:
    for a, b in ((i, j), (i, k), (j, k)):
        excl.add((min(a, b), max(a, b)))

ah = mm.CustomNonbondedForce(
    "select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6-shift),"
    "4*eps*((s/r)^12-(s/r)^6-l*shift)+eps*(1-l)); "
    "s=0.5*(s1+s2); l=0.5*(l1+l2); "
    "shift=(0.5*(s1+s2)/rc)^12-(0.5*(s1+s2)/rc)^6"
)
ah.addGlobalParameter("eps", EPS)
ah.addGlobalParameter("rc", 2.0)
ah.addPerParticleParameter("s")
ah.addPerParticleParameter("l")
ah.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
ah.setCutoffDistance(2.0 * unit.nanometer)
for i in range(2 * N_RES):
    ah.addParticle([0.6954, 0.0] if i % 2 == 0 else [0.5959, 1.18])
for a, b in excl:
    ah.addExclusion(a, b)
system.addForce(ah)

dh = mm.CustomNonbondedForce("q*(exp(-kappa*r)/r-shift); q=q1*q2")
dh.addGlobalParameter("kappa", KAPPA)
dh.addGlobalParameter("shift", DH_SHIFT)
dh.addPerParticleParameter("q")
dh.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
dh.setCutoffDistance(4.0 * unit.nanometer)
for i in range(2 * N_RES):
    dh.addParticle([Q_BB if i % 2 == 0 else 0.0])
for a, b in excl:
    dh.addExclusion(a, b)
system.addForce(dh)

st = mm.CustomBondForce(
    "select(step(r-2^(1/6)*s),n*4*eps*l*((s/r)^12-(s/r)^6-shift),"
    "n*4*eps*((s/r)^12-(s/r)^6-l*shift)+n*eps*(1-l)); "
    "shift=(s/rc_nb)^12-(s/rc_nb)^6"
)
st.addGlobalParameter("eps", EPS)
st.addGlobalParameter("rc_nb", 2.0)
st.addPerBondParameter("s")
st.addPerBondParameter("l")
st.addPerBondParameter("n")
st.setUsesPeriodicBoundaryConditions(True)
for i in range(N_RES - 1):
    st.addBond(2 * i + 1, 2 * i + 3, [0.4, 1.18, 15.0])
system.addForce(st)

# --- 模拟 ---
integrator = mm.LangevinMiddleIntegrator(
    T_K * unit.kelvin, 0.001 / unit.picosecond, DT_PS * unit.picosecond
)
integrator.setRandomNumberSeed(SEED)
# Change platform to 'Reference' for CPU execution as 'CUDA' was not found.
platform = mm.Platform.getPlatformByName("Reference")
sim = app.Simulation(
    topo, system, integrator, platform,
    # The 'properties' argument is specific to certain platforms (like CUDA).
    # It's removed when using 'Reference' to avoid errors.
    # {"DeviceIndex": str(GPU), "CudaPrecision": "mixed"},
)
sim.context.setPositions(positions)

print(f"Using platform: {platform.getName()} | kappa={KAPPA} | q_BB={Q_BB}")
print("Minimize...")
sim.minimizeEnergy(maxIterations=10_000)
with open(os.path.join(OUT, "em.pdb"), "w", encoding="utf-8") as f:
    app.PDBFile.writeFile(topo, sim.context.getState(getPositions=True).getPositions(), f)

sim.context.setVelocitiesToTemperature(T_K * unit.kelvin, SEED)
dcd = os.path.join(OUT, "trajectory.dcd")
sim.reporters.append(app.DCDReporter(dcd, REPORT))
sim.reporters.append(
    app.StateDataReporter(
        os.path.join(OUT, "simulation.log"), REPORT,
        step=True, potentialEnergy=True, temperature=True, speed=True,
    )
)
sim.reporters.append(
    app.StateDataReporter(sys.stdout, REPORT, step=True, potentialEnergy=True, temperature=True)
)

print(f"MD {STEPS} steps ({STEPS * DT_PS / 1000} ns)...")
rg_path = os.path.join(OUT, "rg.csv")
with open(rg_path, "w", encoding="utf-8") as rg_f:
    rg_f.write("time_ns,step,Rg_nm\n")

def calc_rg(sim):
    p = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
    p = np.asarray(p, dtype=float)
    if p.ndim == 3:
        p = p[0]
    com = np.average(p, axis=0, weights=masses)
    return math.sqrt(np.average(np.sum((p - com) ** 2, axis=1), weights=masses))

done, t0, last_rg = 0, time.time(), float("nan")
while done < STEPS:
    chunk = min(RG_EVERY, STEPS - done)
    sim.step(chunk)
    done += chunk

    last_rg = calc_rg(sim)
    time_ns = done * DT_PS / 1000
    with open(rg_path, "a", encoding="utf-8") as rg_f:
        rg_f.write(f"{time_ns:.4f},{done},{last_rg:.6f}\n")

    el = time.time() - t0
    pct = 100 * done / STEPS
    spd = done / el if el else 0
    eta = (STEPS - done) / spd / 60 if spd else 0
    n = int(40 * done / STEPS)
    sys.stdout.write(
        f"\r[{('█' * n).ljust(40, '░')}] {pct:5.1f}%  {done}/{STEPS}  "
        f"{spd:,.0f} step/s  ETA {eta:.1f} min  Rg {last_rg:.4f} nm  "
    )
    sys.stdout.flush()
print()

# --- 后处理：帧 51–200 平均 Rg（同 ipamd target_frame='51-200'）---
# OpenMM 8 的 DCDFile 仅用于写入；Rg 已在模拟中每 REPORT 步记录，直接读 rg.csv
rgs = []
with open(rg_path, encoding="utf-8") as f:
    next(f)  # header
    for line in f:
        if line.startswith("mean"):
            continue
        _, step, rg = line.strip().split(",")
        frame = int(step) // REPORT
        if RG_FRAME_START <= frame <= RG_FRAME_END:
            rgs.append(float(rg))
mean_rg = float(np.mean(rgs))
print(f"Mean Rg (frames {RG_FRAME_START}-{RG_FRAME_END}): {mean_rg:.4f} nm")
with open(rg_path, "a", encoding="utf-8") as f:
    f.write(f"mean_{RG_FRAME_START}_{RG_FRAME_END},{mean_rg}\n")
print(f"rg.data = {mean_rg}")
