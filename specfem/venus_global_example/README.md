# venus_specfem

Generate SPECFEM2D run folders for a global axisymmetric Venus model with
a fluid atmosphere and a liquid core. The standard experiment is a pair
of runs on the same mesh, one with Venus air and one without, whose
difference shows how the atmosphere changes what a surface seismometer
records. Air receivers at 1–60 km show what reaches balloon altitudes.

The tool only writes files. It never runs the mesher or the solver and
never submits jobs: you copy the run folders to your cluster and submit
them yourself.

## Default experiment

| Setting | Value |
|---|---|
| Shortest period | 20 s (Ricker f0 = 0.02 Hz) |
| Mesh | ~1.7 M elements, `multiplication_factor = 50` |
| Air | 0–100 km with realistic sound speed; absorbing sponge from 70 to 100 km |
| Simulated time | 2,000 s |
| Cost per run | ~90–275 core-hours, roughly 0.5–1.6 h on 168 cores |

## Files

| File | What it is |
|---|---|
| `venus_config.py` | The only file you normally edit |
| `make_venus_run.py` | The generator (Python standard library only, Python ≥ 3.7) |
| `compare_runs.py` | Overlays seismograms from two runs (needs numpy and matplotlib) |
| `fortran/define_external_model.f90` | Model reader you install once in SPECFEM2D |
| `data/venus_interior.txt` | Solid planet and core, region-tagged, SI units |
| `data/venus_air_profile.txt` | Approximate VIRA atmosphere (replace it if you have a published table) |
| `templates/create_mesh_Venus.F90` | Stand-alone mesher, adapted from the SPECFEM2D global AK135F example with Venus layering |
| `templates/Par_file`, `templates/SOURCE` | Copied from SPECFEM2D's `global_Earth_ak135f/DATA`; replace them with your SPECFEM2D version's copies if the solver reports missing parameters |

## Steps

### 1. Install the model reader in SPECFEM2D (once)

The model comes from a text file in each run folder
(`DATA/venus_model.txt`) instead of arrays compiled into the solver, so
changing the model never requires editing Fortran or recompiling.

```bash
cd /path/to/specfem2d
cp src/specfem2D/define_external_model.f90 \
   src/specfem2D/define_external_model.f90.orig
cp /path/to/venus_specfem/fortran/define_external_model.f90 src/specfem2D/
make xspecfem2D
```

The routine only interpolates between table rows of an element's own
region, so values are never blended across the core–mantle boundary or
the ground–air interface.

### 2. Edit `venus_config.py`

Most defaults are already set for the standard experiment. Set the paths
and the cluster block.

```python
SPECFEM_BIN = "/path/to/specfem2d/bin"
MIN_PERIOD_S = 20.0                     # sets mesh size and f0
CASES = {"air": 1.0, "noair": 1.0e-3}   # air-density factor per run
SPONGE = True                           # 70-100 km, Q down to 3
SOURCE_TYPE = "force"                   # vertical force on the axis
SOURCE_DEPTH_KM = 10.0
SIM_DURATION_S = 2000.0
GROUND_STATION_DEG = [5, 10, 20, 30, 45, 60, 90]
AIR_STATION_ALT_KM = [1, 10, 30, 50, 60]
NPROC = 168
SLURM = dict(nodes=5, ntasks=168, partition="standard", time="12:00:00")
```

### 3. Generate the run folders

```bash
python make_venus_run.py venus_config.py
# writes runs/venus_T20/air/ and runs/venus_T20/noair/, each with:
#   DATA/Par_file  DATA/SOURCE  DATA/STATIONS  DATA/venus_model.txt
#   MESH/Nummaterial_velocity_file_NO_MUD
#   create_mesh_venus.F90   run_summary.txt
#   run_venus_cluster.sh    run_venus_local.sh
```

Read `run_summary.txt` first. It lists the element size against the size
needed, the time step, the number of steps and the estimated cost. The
script stops with a message if a mesher setting would fail the mesher's
own checks.

### 4. Run on the cluster or a laptop

Both scripts compile and run the stand-alone mesher, then run
`xmeshfem2D` and `xspecfem2D`. Run them from inside the run folder. Each
one first sets `NPROC` in `DATA/Par_file` to its own process count, so
the same folder works on either machine.

On the cluster (edit the `#SBATCH` lines if needed):

```bash
cd runs/venus_T20/air   && sbatch run_venus_cluster.sh
cd ../noair             && sbatch run_venus_cluster.sh
```

On a laptop or workstation, use long periods (`MIN_PERIOD_S` of 60–100 s)
and set `LOCAL_NPROC`, plus `LOCAL_SPECFEM_BIN` if SPECFEM2D lives in a
different place there:

```bash
cd runs/venus_T60/air   && ./run_venus_local.sh > run.log 2>&1
cd ../noair             && ./run_venus_local.sh > run.log 2>&1
```

`run_summary.txt` gives the expected run time for both.

### 5. Check the first few minutes of the solver log

These lines confirm the setup took effect before you spend the full
allocation.

```text
Venus table rows read:          ...        # new routine is in the binary
rho min/max (kg/m3)  :   7.4E-05  ...      # SI units, air reaches the top
Qkappa min           :   3.00              # sponge is present
Total number of acoustic elements = ...    # air + core are fluid
viscoacoustic attenuation:  T
*** Maximum suggested time step = ...      # must be larger than DT
*** Max CFL stability condition ... = ...  # must be below ~0.5
```

### 6. Compare the two runs

```bash
python compare_runs.py runs/venus_T20 G010 G030 G060 A010H010 A010H050
# writes runs/venus_T20/figures/<station>_<component>.png
```

Each figure overlays the two runs and plots their difference, for
vertical and radial velocity and for pressure.

## How the model and mesh are set up

The mesher is adapted from SPECFEM2D's global AK135F axisymmetric example
(`EXAMPLES/applications/axisymmetric_examples/global_Earth_ak135f`). Its
layers are shifted to represent Venus: the top layer is the atmosphere,
then the upper mantle, the lower mantle, and a liquid core. Each choice
below says what the generator does and why.

### 1. Realistic sound speed in the air

Venus air carries sound at ~410 m/s at the ground and ~210 m/s at 90 km.
The sound speed controls how strongly the ground couples to the air and
where energy is trapped.

```text
# DATA/venus_model.txt  (region radius rho vp vs qkappa qmu)
4 6051900.000 6.633060e+01 406.4900 0.0000 9999.0000 9999.0000
4 6101900.000 1.618290e+00 289.7800 0.0000 9999.0000 9999.0000
4 6141900.000 1.000200e-03 210.7800 0.0000    4.5000 9999.0000
```

The values come from `data/venus_air_profile.txt`, an approximate VIRA
profile.

### 2. Absorbing sponge at the top of the air

An untreated top edge of a fluid acts as a rigid lid and reflects
everything. Combined with density falling ~1,000,000× with height, that
makes the air ring. A layer whose Q falls to 3 absorbs waves before they
reach the top. SPECFEM2D's PML only damps along Cartesian x/z directions,
so it can't follow a curved planetary boundary, and Stacey conditions are
disabled for AXISYM.

```text
# Par_file
ATTENUATION_VISCOACOUSTIC = .true.     # attenuation in fluids (the sponge)
PML_BOUNDARY_CONDITIONS   = .false.
# venus_model.txt: 1/Q ramps linearly from 0 at 70 km to 1/3 at 100 km
```

### 3. Air resolved at 20 s, without wasting elements

Air elements must be no wider than ~3.3 km for 20 s waves at 207 m/s.
Shifting the Earth layers also squeezes the original 56×m radial
elements into Venus's 738 km upper mantle. Those elements end up 10×
flatter than wide, and they set the time step. Reducing the interior
radial counts makes runs ~10× cheaper.

```fortran
integer, parameter :: multiplication_factor = 50
integer, parameter :: nspec_rad_670_surf           = 126   ! air (was 15*m)
integer, parameter :: nspec_rad_CMB_670            = 400   ! was 56*m
integer, parameter :: nspec_rad_doubling_OC_to_CMB = 400   ! was 24*m
integer, parameter :: nspec_rad_ICB_to_doubling_OC = 400   ! was 40*m
```

### 4. Mesher radii match the model table

If the mesher's core boundary sits even 1 m off the table's, elements on
one side sample the other side's properties. The generator reads every
boundary radius from the model table.

```fortran
double precision, parameter :: RICB = 3352753.000d0   ! core-mantle boundary
double precision, parameter :: RCMB = 5313568.000d0   ! mid-mantle boundary
```

### 5. Solver settings

Newmark time stepping is used; LDDRK costs several times more per
simulated second for no benefit here. The time step is estimated from
the mesh and calibrated against the solver's own stability check.

```text
time_stepping_scheme = 1        # Newmark
DT                   = 0.012    # ~80% of the estimated maximum (0.015)
NSTEP                = 166700   # 2000 s
```

### 6. Source and stations

A shallow source puts energy into surface waves, where the air coupling
happens. On the symmetry axis only a vertical force or an explosion is
valid. Ground stations sit 20 m below the surface, so they are clearly in
rock instead of on the rock–air boundary. Air stations sit at fixed
altitudes below the sponge.

```text
# SOURCE
zs = 6041900.0     anglesource = 0.     f0 = 0.02
# STATIONS (G = ground, A = air: A010H050 = 10 degrees out, 50 km up)
G030      VE  3025940.00  5241081.82  0.0  0.0
A010H050  VE  1059583.82  6009198.43  0.0  0.0
```

### 7. A clean "no atmosphere" reference

Removing the air layer would require a different mesh. Making the air
1000× lighter keeps the mesh identical, and the pressure the air exerts
on the ground scales with density, so the only difference between the
runs is the atmosphere.

```python
CASES = {"air": 1.0, "noair": 1.0e-3}
```

## Reading the results

- Treat late arrivals in the air with care even with the sponge.
  `run_summary.txt` gives the vertical travel time to the top of the air
  (~360 s). Any reflection that survives the sponge would arrive about
  twice that long after each arrival.
- Pressure records (`PRE.semp`) at air stations, divided by the vertical
  ground velocity below them, give the ground-to-air coupling. At the
  Venus surface, ρc ≈ 2.7×10⁴ Pa·s/m.
- Amplitude grows with height as the air thins. Compare the 1, 10, 30, 50
  and 60 km stations at the same distance.
- The model is purely acoustic in the air: no gravity, wind or viscosity.
  That is a reasonable approximation at these periods, well above Venus's
  ~2 mHz acoustic cutoff.

## Troubleshooting

| Symptom | Fix |
|---|---|
| Log says "external model: ak135-f" or has no "Venus table" line | The binary wasn't rebuilt with the new routine. Redo step 1. |
| "some parameters are missing in your Par_file" | Your SPECFEM2D version differs from the Par_file template. Copy your version's example Par_file to `templates/Par_file` and regenerate. |
| The run blows up (NaN, or an exploding "Max norm") | Set `DT_OVERRIDE` to about 80% of "Maximum suggested time step" and regenerate. |
| The mesher runs out of memory | Set `MIN_PERIOD_S = 30` (m ≈ 34). This also makes runs ~3× cheaper. |
| The job hits its wall-time limit | Use a longer partition, or shorten `SIM_DURATION_S`. |

## Re-extracting a hard-coded model

To convert a `define_external_model.f90` that stores the model as
hard-coded `*_ak135(i) = ...` arrays into `data/venus_interior.txt`:

```bash
python make_venus_run.py --extract-table path/to/define_external_model.f90
```
