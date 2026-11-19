"""
Settings for one Venus global AXISYM experiment.

Edit this file, then run:  python make_venus_run.py venus_config.py
Every case in CASES becomes its own run folder under OUTPUT_DIR/EXPERIMENT.
Nothing is run or submitted. Each run folder gets two scripts, run from
inside that folder:
    run_venus_cluster.sh  ->  sbatch run_venus_cluster.sh
    run_venus_local.sh    ->  ./run_venus_local.sh > run.log 2>&1
"""

# --- where things live ------------------------------------------------------
SPECFEM_BIN = "/home/bhchow/REPOS/specfem2d/bin"   # rebuilt with the new
                                                     # define_external_model
FC = "gfortran"            # compiler for the stand-alone mesher
OUTPUT_DIR = "runs"
EXPERIMENT = "venus_T60"   # folder name for this set of runs

# --- resolution -------------------------------------------------------------
MIN_PERIOD_S = 60.0        # shortest period resolved; Ricker f0 = 1/(2.5*T)
POINTS_PER_WAVELENGTH = 5.0
MULTIPLICATION_FACTOR = None   # None = choose automatically from the period
INTERIOR_RADIAL_FACTOR = 8     # radial element count factor in the mantle
                               # (was 56/24/40; 8 keeps elements well shaped)

# --- atmosphere -------------------------------------------------------------
AIR_TOP_KM = 100.0         # top of the mesh above the ground
SPONGE = True              # absorbing Q sponge at the top of the air
SPONGE_BOTTOM_KM = 70.0    # sponge from here to AIR_TOP_KM
SPONGE_Q_MIN = 3.0         # Q at the very top (1/Q ramps up linearly)

# case name -> factor on air density. 1e-3 makes the air too light to
# push on the ground, which acts as the "no atmosphere" reference run
# on exactly the same mesh.
CASES = {"air": 1.0, "noair": 1.0e-3}

# --- source -----------------------------------------------------------------
SOURCE_TYPE = "force"      # "force" (vertical) or "explosion"
SOURCE_DEPTH_KM = 10.0     # below the ground, on the symmetry axis
SIM_DURATION_S = 2000.0    # simulated time

# DT is estimated conservatively. After the first run, read the stability
# lines in OUTPUT_FILES/output_solver.txt and set DT_OVERRIDE if needed.
DT_OVERRIDE = None

# --- stations (angles in degrees from the source) ---------------------------
SEISMOTYPE = "2,4"         # 2 = velocity, 4 = pressure
GROUND_STATION_DEG = [5, 10, 20, 30, 45, 60, 90]
AIR_STATION_DEG = [10, 30, 60]
AIR_STATION_ALT_KM = [1, 10, 30, 50, 60]   # keep below SPONGE_BOTTOM_KM

# --- cluster ----------------------------------------------------------------
NPROC = 168
SLURM = dict(
    job_name="venus",
    nodes=5,
    ntasks=168,
    partition="t2standard",     # pick a partition with enough wall time
    time="12:00:00",
    output="log.txt",
)

# --- laptop / workstation (run_venus_local.sh) -------------------------------
LOCAL_NPROC = 4            # MPI processes on the laptop
LOCAL_SPECFEM_BIN = None   # laptop SPECFEM2D bin dir; None = SPECFEM_BIN
LOCAL_FC = None            # laptop Fortran compiler; None = FC
LOCAL_MPIRUN = "mpirun"    # "" only if SPECFEM2D was built without MPI
                           # (then LOCAL_NPROC must be 1)
# Laptop runs: use long periods, e.g. MIN_PERIOD_S = 60 to 100. The cost
# line in run_summary.txt gives the expected laptop run time.
