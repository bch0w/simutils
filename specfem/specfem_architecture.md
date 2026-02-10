# SPECFEM3D Architecture Overview
2026-01-30 Written entirely by GPT-5 through VSCode after asking it to introspect SPECFEM code to help explain code architecture to a new developer. Not reviewed.

## Overview
SPECFEM3D is a spectral-element method (SEM) code for simulating wave propagation in 3D heterogeneous media. It solves the elastic, acoustic, poroelastic, or coupled wave equations on unstructured hexahedral meshes, primarily for seismic applications but also applicable to other wave phenomena.

## Code Execution Flowchart

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              SPECFEM3D Workflow                                 │
└─────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│  xmeshfem3D     │    │ xdecompose_mesh │    │xgenerate_databases│    │   xspecfem3D   │
│  (Mesh Gen)     │    │  (Partition)    │    │   (Setup DB)      │    │   (Solver)     │
└─────────────────┘    └─────────────────┘    └─────────────────┘    └─────────────────┘
         │                       │                       │                       │
         │                       │                       │                       │
         ▼                       ▼                       ▼                       ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   MESH-default/ │    │ DATABASES_MPI/ │    │ DATABASES_MPI/ │    │ OUTPUT_FILES/   │
│   (Mesh files)  │    │ (Partitioned   │    │ (Precomputed   │    │ (Results)       │
│                 │    │  mesh files)   │    │  databases)    │    │                 │
└─────────────────┘    └─────────────────┘    └─────────────────┘    └─────────────────┘
         ▲                       ▲                       ▲                       │
         │                       │                       │                       │
         └───────────────────────┼───────────────────────┼───────────────────────┘
                                 │                       │
                                 ▼                       ▼
                    ┌─────────────────┐    ┌─────────────────┐
                    │   DATA/Par_file │    │   DATA/Par_file │
                    │ (Configuration) │    │ (Configuration) │
                    └─────────────────┘    └─────────────────┘

┌─────────────────────────────────────────────────────────────────────────────────┐
│                          Internal Solver Flow                                  │
└─────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────┐
│   xspecfem3D    │
│   (Main)        │
└─────────────────┘
         │
         ▼
┌─────────────────┐     ┌─────────────────┐
│  initialize_    │     │  read_mesh_     │
│  simulation()   │     │  databases()    │
└─────────────────┘     └─────────────────┘
         │                       │
         └──────────┬────────────┘
                    ▼
         ┌─────────────────┐
         │ setup_GLL_      │
         │ points()        │
         └─────────────────┘
                    ▼
         ┌─────────────────┐
         │ setup_sources_  │
         │ receivers()     │
         └─────────────────┘
                    ▼
         ┌─────────────────┐
         │ prepare_        │
         │ timerun()       │
         └─────────────────┘
                    ▼
         ┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
         │  iterate_time() │────▶│ compute_forces_ │────▶│ compute_add_    │
         │  (Time Loop)    │     │ viscoelastic()  │     │ sources()       │
         └─────────────────┘     └─────────────────┘     └─────────────────┘
                    ▲                       │                       │
                    │                       ▼                       ▼
                    │            ┌─────────────────┐     ┌─────────────────┐
                    │            │ compute_stacey_ │     │ write_          │
                    │            │ viscoelastic()  │     │ seismograms()   │
                    │            └─────────────────┘     └─────────────────┘
                    │
                    ▼
         ┌─────────────────┐
         │ finalize_       │
         │ simulation()    │
         └─────────────────┘

┌─────────────────────────────────────────────────────────────────────────────────┐
│                          Module Dependencies                                  │
└─────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│  specfem_par    │◄────┤ specfem_par_    │◄────┤ specfem_par_    │
│  (Main params)  │     │ elastic         │     │ acoustic        │
└─────────────────┘     └─────────────────┘     └─────────────────┘
         ▲                       ▲                       ▲
         │                       │                       │
         ▼                       ▼                       ▼
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│ shared_par      │     │ pml_par         │     │ coupling        │
│ (Shared utils)  │     │ (PML)           │     │ modules         │
└─────────────────┘     └─────────────────┘     └─────────────────┘
         ▲
         │
         ▼
┌─────────────────┐
│  constants      │
│  (Core consts)  │
└─────────────────┘

┌─────────────────────────────────────────────────────────────────────────────────┐
│                            Data Flow                                           │
└─────────────────────────────────────────────────────────────────────────────────┘

Par_file ─────────────────────────────────────────────────────────────────────────►
    │
    ├─► xmeshfem3D ──► Mesh files ──► xdecompose_mesh ──► Partitioned mesh ──►
    │                                                                                │
    └─► xgenerate_databases ──► Databases ──► xspecfem3D ──► Wavefields ──► Output
                                                                       │
                                                                       ▼
                                                            Seismograms/Snapshots
```

## Key Components and Workflow

### 1. Mesh Generation (`src/meshfem3D/`)
- **Purpose**: Creates or reads 3D hexahedral meshes
- **Main executable**: `xmeshfem3D`
- **Key features**:
  - Supports internal mesher for simple geometries
  - Reads external meshes (e.g., from CUBIT/GEOCUBIT)
  - Handles topography, oceans, and complex boundaries
  - Outputs mesh files to `MESH-default/` or similar

### 2. Mesh Decomposition (`src/decompose_mesh/`)
- **Purpose**: Partitions the mesh for parallel computation
- **Main executable**: `xdecompose_mesh`
- **Key features**:
  - Uses partitioning libraries (SCOTCH, METIS, PATOH)
  - Distributes elements across MPI processes
  - Optimizes load balancing and communication

### 3. Database Generation (`src/generate_databases/`)
- **Purpose**: Precomputes simulation data structures
- **Main executable**: `xgenerate_databases`
- **Key features**:
  - Sets up Gauss-Lobatto-Legendre (GLL) integration points
  - Computes mass matrices and stiffness matrices
  - Assigns material properties from models (1D, 3D, tomography)
  - Prepares absorbing boundaries and coupling surfaces
  - Outputs databases to `OUTPUT_FILES/DATABASES_MPI/`

### 4. Wave Solver (`src/specfem3D/`)
- **Purpose**: Time-steps the wave equation
- **Main executable**: `xspecfem3D`
- **Key features**:
  - Explicit time integration (Newmark or LDDRK schemes)
  - Supports multiple physics: elastic, acoustic, poroelastic, coupled
  - Handles attenuation, anisotropy, gravity
  - Computes seismograms and snapshots
  - Adjoint simulations for tomography

## Common Structures and Patterns

### Modules and Parameters
- **`constants`**: Fundamental constants (GLL points NGLLX=5, I/O units, precision)
- **`shared_input_parameters`**: Simulation parameters from `DATA/Par_file`
- **`specfem_par`**: Main solver parameters and arrays
- **`shared_par`**: Shared utilities and I/O routines

### Spectral Element Method Basics
- **Elements**: Hexahedral (8 or 27 nodes)
- **Integration**: Gauss-Lobatto-Legendre quadrature (typically 5x5x5 points per element)
- **Basis functions**: Lagrange polynomials at GLL points
- **Global assembly**: MPI-based parallel assembly of sparse matrices

### Parallelism
- **MPI**: Domain decomposition across processes
- **OpenMP**: Thread-level parallelism within processes (optional)
- **GPU**: CUDA/HIP acceleration for compute-intensive kernels

### Physics Modules
- **Elastic**: Displacement formulation, stress-strain
- **Acoustic**: Pressure potential formulation
- **Poroelastic**: Solid-fluid coupling
- **Coupling**: Interfaces between different media types

### I/O and Data Formats
- **Database files**: Binary MPI I/O for mesh and precomputed data
- **Seismograms**: ASCII, binary, SU, ASDF formats
- **Snapshots**: VTK, HDF5 for visualization
- **Advanced I/O**: ADIOS2 for high-performance I/O

## Interconnections

1. **Mesh → Decomposition**: Mesh files are partitioned for parallel processing
2. **Decomposition → Databases**: Partition info guides database setup
3. **Databases → Solver**: Precomputed arrays enable efficient time-stepping
4. **Shared Modules**: Common utilities (GLL, MPI, I/O) used across all components
5. **Parameter Flow**: `Par_file` configures all stages consistently

## Key Algorithms

- **Time Integration**: Explicit schemes with CFL stability constraints
- **Boundary Conditions**: Absorbing (PML/CPML), free surface, coupling
- **Source Modeling**: Moment tensors, force sources, external time functions
- **Receiver Processing**: Interpolation to station locations
- **Adjoint Methods**: Reverse time migration for sensitivity kernels

## Deep Code Insights for Developers

### Core Data Structures

#### Mesh Representation (`specfem_par`)
```fortran
integer :: NSPEC_AB, NGLOB_AB  ! Number of spectral elements and global points
integer, dimension(:,:,:,:), allocatable :: ibool  ! Maps local (i,j,k,ispec) to global iglob
real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,ystore,zstore  ! Global coordinates
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
  xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore  ! Derivatives for mapping
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: jacobianstore  ! Jacobian determinants
```

The `ibool` array is crucial - it maps the 4D local element indexing (i,j,k,ispec) to 1D global node indexing, enabling efficient sparse matrix operations.

#### Material Properties
```fortran
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kappastore,mustore  ! Lame parameters
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore  ! Density
```

These are stored per GLL point, allowing for heterogeneous media.

#### Wavefields
```fortran
real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ,veloc,accel  ! Displacement, velocity, acceleration
```

Wavefields are stored globally, with MPI communication handling boundaries.

### Main Time Loop Structure (`iterate_time.F90`)

```fortran
do it = it_begin,it_end
  ! Update wavefields (Newmark or LDDRK)
  call update_displ_Newmark()
  
  ! Compute forces for each physics domain
  if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_forward_calling()
  if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_calling()
  if (POROELASTIC_SIMULATION) call compute_forces_poroelastic_calling()
  
  ! Add sources
  call compute_add_sources_viscoelastic(accel)
  
  ! Apply boundary conditions
  call compute_stacey_viscoelastic()
  
  ! Write seismograms
  call write_seismograms()
enddo
```

The loop alternates between updating wavefields and computing right-hand-side forces.

### Spectral Element Computation (`compute_forces_viscoelastic.F90`)

The core SEM computation follows this pattern for each element:

```fortran
do ispec = 1, NSPEC_AB
  ! Loop over GLL points
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        ! Compute strain from displacement gradients
        duxdxl = 0._CUSTOM_REAL
        duxdyl = 0._CUSTOM_REAL
        duxdzl = 0._CUSTOM_REAL
        ! ... (compute derivatives using hprime matrices)
        
        ! Compute stress from strain and material properties
        sigma_xx = lambda * (duxdxl + duydyl + duzdzl) + 2._CUSTOM_REAL * mu * duxdxl
        ! ... (similar for other components)
        
        ! Compute force contribution
        ! ... (integrate against test functions)
      enddo
    enddo
  enddo
enddo
```

This implements the weak form: ∫ σ : ∇v dΩ = boundary terms + source terms.

### GLL Basis and Integration

```fortran
! GLL points and weights (from constants.h.in)
double precision, dimension(NGLLX) :: xigll,wxgll

! Derivative matrices
real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx  ! dL_i/dξ at GLL points
```

The `hprime` matrices enable efficient computation of derivatives without explicit differentiation.

### MPI Communication

```fortran
! Assembly buffers
real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh
integer, dimension(:), allocatable :: request_send_vector_ext_mesh

! Synchronization
call assemble_MPI_vector()
```

MPI handles element boundaries between processes, ensuring continuity.

### Source Implementation (`compute_add_sources_viscoelastic.F90`)

```fortran
! For each source
do isource = 1, nsources_local
  ! Get source time function
  stf = get_stf_viscoelastic(it, isource)
  
  ! Add to acceleration at source location
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        iglob = ibool(i,j,k,ispec_selected_source(isource))
        accel(:,iglob) = accel(:,iglob) + stf * sourcearrays(:,i,j,k,isource)
      enddo
    enddo
  enddo
enddo
```

Sources are distributed across GLL points using precomputed `sourcearrays`.

### Boundary Conditions

- **Absorbing boundaries**: Implemented via Stacey conditions or PML
- **Free surface**: Zero traction (σ·n = 0)
- **Coupling**: Continuity of displacement/traction at material interfaces

### Adjoint Simulations

For tomography, adjoint sources are added similarly to forward sources, but with reversed time:

```fortran
! Read adjoint sources
call read_adjoint_sources_ASDF()

! Add to adjoint acceleration
accel_adj(:,iglob) = accel_adj(:,iglob) + stf_adj * adj_source(:,irec)
```

### GPU Acceleration

GPU kernels parallelize over elements and GLL points:

```fortran
! CUDA kernel launch
call compute_forces_viscoelastic_cuda(Mesh_pointer, iphase, backward_simulation)
```

### Key Modification Points for Developers

1. **Adding new physics**: Extend `compute_forces_*` routines and add corresponding modules
2. **New boundary conditions**: Modify `compute_stacey_*` or add new boundary routines
3. **Custom sources**: Extend `compute_add_sources_*` or `get_stf_*`
4. **Material models**: Modify `get_model` in generate_databases
5. **I/O formats**: Add new writers in `write_output_*`
6. **Time schemes**: Extend `update_displ_*` routines

### Performance Considerations

- **Memory layout**: Arrays are optimized for cache locality (GLL points contiguous)
- **Vectorization**: FORCE_VECTORIZATION flag enables SIMD optimizations
- **Load balancing**: Mesh partitioning affects parallel efficiency
- **I/O bottlenecks**: ADIOS/HDF5 used for large-scale output

This architecture enables scalable simulations on thousands of cores while maintaining accuracy through high-order spectral elements.

## Key Components and Workflow

### 1. Mesh Generation (`src/meshfem3D/`)
- **Purpose**: Creates or reads 3D hexahedral meshes
- **Main executable**: `xmeshfem3D`
- **Key features**:
  - Supports internal mesher for simple geometries
  - Reads external meshes (e.g., from CUBIT/GEOCUBIT)
  - Handles topography, oceans, and complex boundaries
  - Outputs mesh files to `MESH-default/` or similar

### 2. Mesh Decomposition (`src/decompose_mesh/`)
- **Purpose**: Partitions the mesh for parallel computation
- **Main executable**: `xdecompose_mesh`
- **Key features**:
  - Uses partitioning libraries (SCOTCH, METIS, PATOH)
  - Distributes elements across MPI processes
  - Optimizes load balancing and communication

### 3. Database Generation (`src/generate_databases/`)
- **Purpose**: Precomputes simulation data structures
- **Main executable**: `xgenerate_databases`
- **Key features**:
  - Sets up Gauss-Lobatto-Legendre (GLL) integration points
  - Computes mass matrices and stiffness matrices
  - Assigns material properties from models (1D, 3D, tomography)
  - Prepares absorbing boundaries and coupling surfaces
  - Outputs databases to `OUTPUT_FILES/DATABASES_MPI/`

### 4. Wave Solver (`src/specfem3D/`)
- **Purpose**: Time-steps the wave equation
- **Main executable**: `xspecfem3D`
- **Key features**:
  - Explicit time integration (Newmark or LDDRK schemes)
  - Supports multiple physics: elastic, acoustic, poroelastic, coupled
  - Handles attenuation, anisotropy, gravity
  - Computes seismograms and snapshots
  - Adjoint simulations for tomography

## Common Structures and Patterns

### Modules and Parameters
- **`constants`**: Fundamental constants (GLL points NGLLX=5, I/O units, precision)
- **`shared_input_parameters`**: Simulation parameters from `DATA/Par_file`
- **`specfem_par`**: Main solver parameters and arrays
- **`shared_par`**: Shared utilities and I/O routines

### Spectral Element Method Basics
- **Elements**: Hexahedral (8 or 27 nodes)
- **Integration**: Gauss-Lobatto-Legendre quadrature (typically 5x5x5 points per element)
- **Basis functions**: Lagrange polynomials at GLL points
- **Global assembly**: MPI-based parallel assembly of sparse matrices

### Parallelism
- **MPI**: Domain decomposition across processes
- **OpenMP**: Thread-level parallelism within processes (optional)
- **GPU**: CUDA/HIP acceleration for compute-intensive kernels

### Physics Modules
- **Elastic**: Displacement formulation, stress-strain
- **Acoustic**: Pressure potential formulation
- **Poroelastic**: Solid-fluid coupling
- **Coupling**: Interfaces between different media types

### I/O and Data Formats
- **Database files**: Binary MPI I/O for mesh and precomputed data
- **Seismograms**: ASCII, binary, SU, ASDF formats
- **Snapshots**: VTK, HDF5 for visualization
- **Advanced I/O**: ADIOS2 for high-performance I/O

## Interconnections

1. **Mesh → Decomposition**: Mesh files are partitioned for parallel processing
2. **Decomposition → Databases**: Partition info guides database setup
3. **Databases → Solver**: Precomputed arrays enable efficient time-stepping
4. **Shared Modules**: Common utilities (GLL, MPI, I/O) used across all components
5. **Parameter Flow**: `Par_file` configures all stages consistently

## Key Algorithms

- **Time Integration**: Explicit schemes with CFL stability constraints
- **Boundary Conditions**: Absorbing (PML/CPML), free surface, coupling
- **Source Modeling**: Moment tensors, force sources, external time functions
- **Receiver Processing**: Interpolation to station locations
- **Adjoint Methods**: Reverse time migration for sensitivity kernels

## Deep Code Insights for Developers

### Core Data Structures

#### Mesh Representation (`specfem_par`)
```fortran
integer :: NSPEC_AB, NGLOB_AB  ! Number of spectral elements and global points
integer, dimension(:,:,:,:), allocatable :: ibool  ! Maps local (i,j,k,ispec) to global iglob
real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,ystore,zstore  ! Global coordinates
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
  xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore  ! Derivatives for mapping
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: jacobianstore  ! Jacobian determinants
```

The `ibool` array is crucial - it maps the 4D local element indexing (i,j,k,ispec) to 1D global node indexing, enabling efficient sparse matrix operations.

#### Material Properties
```fortran
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kappastore,mustore  ! Lame parameters
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore  ! Density
```

These are stored per GLL point, allowing for heterogeneous media.

#### Wavefields
```fortran
real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ,veloc,accel  ! Displacement, velocity, acceleration
```

Wavefields are stored globally, with MPI communication handling boundaries.

### Main Time Loop Structure (`iterate_time.F90`)

```fortran
do it = it_begin,it_end
  ! Update wavefields (Newmark or LDDRK)
  call update_displ_Newmark()
  
  ! Compute forces for each physics domain
  if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_forward_calling()
  if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_calling()
  if (POROELASTIC_SIMULATION) call compute_forces_poroelastic_calling()
  
  ! Add sources
  call compute_add_sources_viscoelastic(accel)
  
  ! Apply boundary conditions
  call compute_stacey_viscoelastic()
  
  ! Write seismograms
  call write_seismograms()
enddo
```

The loop alternates between updating wavefields and computing right-hand-side forces.

### Spectral Element Computation (`compute_forces_viscoelastic.F90`)

The core SEM computation follows this pattern for each element:

```fortran
do ispec = 1, NSPEC_AB
  ! Loop over GLL points
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        ! Compute strain from displacement gradients
        duxdxl = 0._CUSTOM_REAL
        duxdyl = 0._CUSTOM_REAL
        duxdzl = 0._CUSTOM_REAL
        ! ... (compute derivatives using hprime matrices)
        
        ! Compute stress from strain and material properties
        sigma_xx = lambda * (duxdxl + duydyl + duzdzl) + 2._CUSTOM_REAL * mu * duxdxl
        ! ... (similar for other components)
        
        ! Compute force contribution
        ! ... (integrate against test functions)
      enddo
    enddo
  enddo
enddo
```

This implements the weak form: ∫ σ : ∇v dΩ = boundary terms + source terms.

### GLL Basis and Integration

```fortran
! GLL points and weights (from constants.h.in)
double precision, dimension(NGLLX) :: xigll,wxgll

! Derivative matrices
real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx  ! dL_i/dξ at GLL points
```

The `hprime` matrices enable efficient computation of derivatives without explicit differentiation.

### MPI Communication

```fortran
! Assembly buffers
real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh
integer, dimension(:), allocatable :: request_send_vector_ext_mesh

! Synchronization
call assemble_MPI_vector()
```

MPI handles element boundaries between processes, ensuring continuity.

### Source Implementation (`compute_add_sources_viscoelastic.F90`)

```fortran
! For each source
do isource = 1, nsources_local
  ! Get source time function
  stf = get_stf_viscoelastic(it, isource)
  
  ! Add to acceleration at source location
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        iglob = ibool(i,j,k,ispec_selected_source(isource))
        accel(:,iglob) = accel(:,iglob) + stf * sourcearrays(:,i,j,k,isource)
      enddo
    enddo
  enddo
enddo
```

Sources are distributed across GLL points using precomputed `sourcearrays`.

### Boundary Conditions

- **Absorbing boundaries**: Implemented via Stacey conditions or PML
- **Free surface**: Zero traction (σ·n = 0)
- **Coupling**: Continuity of displacement/traction at material interfaces

### Adjoint Simulations

For tomography, adjoint sources are added similarly to forward sources, but with reversed time:

```fortran
! Read adjoint sources
call read_adjoint_sources_ASDF()

! Add to adjoint acceleration
accel_adj(:,iglob) = accel_adj(:,iglob) + stf_adj * adj_source(:,irec)
```

### GPU Acceleration

GPU kernels parallelize over elements and GLL points:

```fortran
! CUDA kernel launch
call compute_forces_viscoelastic_cuda(Mesh_pointer, iphase, backward_simulation)
```

### Key Modification Points for Developers

1. **Adding new physics**: Extend `compute_forces_*` routines and add corresponding modules
2. **New boundary conditions**: Modify `compute_stacey_*` or add new boundary routines
3. **Custom sources**: Extend `compute_add_sources_*` or `get_stf_*`
4. **Material models**: Modify `get_model` in generate_databases
5. **I/O formats**: Add new writers in `write_output_*`
6. **Time schemes**: Extend `update_displ_*` routines

### Performance Considerations

- **Memory layout**: Arrays are optimized for cache locality (GLL points contiguous)
- **Vectorization**: FORCE_VECTORIZATION flag enables SIMD optimizations
- **Load balancing**: Mesh partitioning affects parallel efficiency
- **I/O bottlenecks**: ADIOS/HDF5 used for large-scale output

This architecture enables scalable simulations on thousands of cores while maintaining accuracy through high-order spectral elements.
