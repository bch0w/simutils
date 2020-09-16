! =============================================================================
! 11.09.2020
!
! This is an edited version of the sem_model_slice program from SPECFEM3D.
! Provided in a modified state to Bryant Chow by Kai Wang. Bryant has edited for
! code formatting, intended for use with a Python-based wrapper Semslicer.
! As such, there are Python formatting curly-brackets included which will need
! to be set by Semslicer.
!
! NOTE: This version is considerably different from my current checked out
!       version of Specfem (75e1785c912b996ac438deaa5b88137ac7705b0a), but was
!       used by Kai Wang as recently as 2018 (?) so I've opted to use it.
!
! github.com/bch0w/simutils/model_slice
!
! =============================================================================
program sem_model_slice

  implicit none

  include 'mpif.h'
  ! Usually found in specfem3d/setup/constants.h
  include '{set_path_to_constants}/constants.h'
  ! Usually found in specfem3d/setup/precision.h
  include '{set_path_to_precision}/precision.h'
  !include 'values_from_mesher.h'

  ! Hard code the max number of GLL points expected
  integer, parameter :: NMAXPTS = 10000000

  integer :: ier, sizeprocs, myrank, ios, i, j, k, ispec, iglob, ipt, npts
  character(len=500) :: xyz_infile, topo_dir, model_dir, data_name, &
          outfile, local_data_file, prname
  real(kind=CUSTOM_REAL), dimension(NMAXPTS) :: x, y, z, v, vmin, vall, &
          distmin, dist, distall
  integer, dimension(NMAXPTS) :: ispec_min, ix_min, iy_min, iz_min
  real(kind=CUSTOM_REAL), dimension(2, NMAXPTS) :: in, out

  !!!
  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(:,:,:,:), allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore, ystore, zstore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: vstore
  !!!

  ! Two implementations for topography:
  ! true  --> replace "air" points with NaN (vertical cross sections)
  ! false --> take the closest value to the "air" points (horizontal x-section)
  logical, parameter :: TOPOGRAPHY = .{set_topography}.

  character(len=100) :: topo_file
  integer, dimension(NX, NY) :: itopo
  double precision :: elevation

  ! MPI initialization
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, sizeprocs, ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)

  ! Input arguments
  call getarg(1, xyz_infile)
  call getarg(2, topo_dir)
  call getarg(3, model_dir)
  call getarg(4, data_name)
  call getarg(5, outfile)

  if(myrank==0)then
    write(*,*) trim(xyz_infile)
    write(*,*) trim(topo_dir)
    write(*,*) trim(model_dir)
    write(*,*) trim(data_name)
    write(*,*) trim(outfile)
  endif

  ! Read points to be interpolated
  open(11, file=xyz_infile, iostat=ios)
  i = 0
  do while (1 == 1)
    i = i + 1
    read(11, *, iostat = ios) x(i), y(i), z(i)
    if (ios /= 0) exit
  enddo
  close(11)
  npts = i - 1
  if (myrank == 0) then
    write(*, *) 'Total number of points = ', npts
    if (npts > NMAXPTS .or. npts <= 0) call &
            exit_mpi(myrank, 'max npts exceeded')
  endif

  ! Kai added read mesh to allocate grid points
  write(prname, '(a,i6.6,a)') trim(topo_dir)//'/proc',myrank,'_'
  open(unit=27, file=trim(prname)//'external_mesh.bin', status='old', &
          action='read', form='unformatted', iostat=ier)
  read(27) NSPEC_AB
  read(27) NGLOB_AB

  ! Allocates mesh arrays to [xyz]store
  allocate(ibool(NGLLX, NGLLY, NGLLZ, NSPEC_AB), stat=ier)
  if (ier /= 0) stop 'Error allocating array ibool'
  allocate(xstore(NGLOB_AB), ystore(NGLOB_AB), zstore(NGLOB_AB), stat=ier)
  if (ier /= 0) stop 'Error allocating array [xyz]store'
  allocate(vstore(NGLLX, NGLLY, NGLLZ, NSPEC_AB), stat=ier)
  if (ier /= 0) stop 'Error allocating array vstore'

  ! read ibool file and global point arrays
  read(27) ibool
  read(27) xstore
  read(27) ystore
  read(27) zstore
  close(27)

  !write(*,*)'myrank,NGLOB_AB,NSPEC_AB=',myrank,NGLOB_AB,NSPEC_AB
  !!!!

  ! Read in model data as vstore
  write(prname,'(a,i6.6,a)') trim(model_dir)//'/proc', myrank,'_'
  local_data_file = trim(prname) // trim(data_name)//'.bin'
  open(unit = 27, file=local_data_file, status='old', iostat = ios, &
          form ='unformatted')
  if (ios /= 0) call exit_mpi(myrank, &
          'Error reading model data file '//trim(local_data_file))
  read(27) vstore
  close(27)

  ! Search for local minimum-distance point
  distmin(1:npts) = HUGEVAL

  do ispec = 1, NSPEC_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i, j, k, ispec)
          dist(1: npts) = &
                  dsqrt((x(1:npts) - dble(xstore(iglob))) ** 2 &
                  +(y(1:npts)-dble(ystore(iglob))) ** 2 &
                  +(z(1:npts)-dble(zstore(iglob))) ** 2)
          do ipt=1,npts
            if(dist(ipt) < distmin(ipt)) then
              distmin(ipt) = dist(ipt)
              ispec_min(ipt) = ispec
              ix_min(ipt) = i
              iy_min(ipt) = j
              iz_min(ipt) = k
              vmin(ipt) = vstore(i, j, k, ispec)
            endif
          enddo
        enddo
      enddo
    enddo
    ! end of loop on all the elements in current slice
  enddo

  ! Frees up memory
  deallocate(ibool)
  deallocate(xstore,ystore,zstore,vstore)

  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if (myrank == 0) print *, 'Done looping over global points ...'

  ! choose the minimum value

  !  write(prname,'(a,i6.6,a)') 'OUTPUT_FILES/in_',myrank,'.txt'
  !  open(33,file=prname)
  !  do i = 1, npts
  !    write(33,*) i, myrank, distmin(i), ispec_min(i), ix_min(i), iy_min(i), iz_min(i),vmin(i)
  !  enddo
  !  close(33)

  do i = 1, npts
    in(1, i) = distmin(i)
    in(2, i) = myrank    ! myrank is coerced to a double
  enddo
  call MPI_REDUCE(in,out,npts,CUSTOM_MPI_2REAL,MPI_MINLOC,0,MPI_COMM_WORLD,ier)

  !  if (myrank == 0)  then
  !   open(33,file='OUTPUT_FILES/out.txt')
  !   do i = 1, npts
  !     write(33,*) i, out(1,i), out(2,i)
  !  enddo
  !   close(33)
  !  endif

  call MPI_BCAST(out, 2*npts, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ier)

  v(1:npts) = 0
  dist(1:npts) = 0.

  do i = 1, npts
    if (myrank == nint(out(2,i))) then
      v(i) = vmin(i)
!      if (GLL_INTERPOLATION) call xeg_search(x(i),y(i),z(i),&
!                 ispec_min(i),ix_min(i),iy_min(i),iz_min(i),v(i))
      dist(i) = distmin(i)
    endif
  enddo

  call MPI_REDUCE(v, vall, npts, CUSTOM_MPI_TYPE, MPI_SUM, 0, &
          MPI_COMM_WORLD, ier)
  call MPI_REDUCE(dist, distall, npts, CUSTOM_MPI_TYPE, MPI_SUM, 0, &
          MPI_COMM_WORLD, ier)

  ! Write out the final .XYZ file containing location and model value
  if (myrank == 0) then
    if (TOPOGRAPHY) then
      topo_file='/ibrixfs1/home/lqy/cmt3d/test_dir/'//trim(TOPO_FILE_SOCAL)
      call read_topo(itopo, NX , NY, topo_file)
    endif

    print *, 'Writing output file ...'
    open(12, file=outfile, status='unknown')
    do i = 1, npts
      ! Check to see if values above topography are being queried
      if (TOPOGRAPHY) then
         call topo_value(itopo, dble(x(i)), dble(y(i)), elevation)
         if (elevation < z(i)) vall(i)= -1000.00
      endif
      write(12,*) x(i), y(i), z(i), vall(i), distall(i)
    enddo
    close(12)
  endif

  call MPI_FINALIZE(ier)

end program sem_model_slice

!------------------------------------------------

subroutine read_topo(itopo_bathy_basin, NX, NY, topo_file)
  ! ===========================================================================
  ! Read topography file which should be provided as a 1 column file in
  ! the variabel 'topo_file'
  ! ===========================================================================
  implicit none
  include "constants.h"
  integer NX, NY

  ! Use integer array to store topography values
  integer itopo(NX, NY)
  character(len=100) topo_file

  integer ix,iy

  ! Initiate topography as 0 to start
  itopo(:,:) = 0

  open(unit=13, file=topo_file, status='old', action='read')
  do iy = 1, NY
    do ix = 1, NX
      read(13, *) itopo(ix, iy)
    enddo
  enddo
  close(13)

end subroutine read_topo

!------------------------------------------------

subroutine topo_value(itopo, x, y, elevation)
  ! ===========================================================================
  ! Query topographic height based on X and Y coordinate values
  ! ===========================================================================
  implicit none

  integer, dimension(NX, NY) :: itopo
  double precision :: x, y
  double precision elevation

  double precision :: long, lat
  integer :: icornerlong, icornerlat
  double precision :: long_corner, lat_corner, ratio_xi, ratio_eta
  integer, parameter :: UTM_PROJECTION_ZONE  = {set_utm_proj}
  logical, parameter :: SUPPRESS_UTM_PROJECTION  = .{set_suppress_utm_proj}.

  ! To avoid requring constants.h file lat/lon values
  ! Need to set the origin corner of the topography file
  double precision, parameter :: TOPO_ORIG_LON = {set_topo_orig_lon}
  double precision, parameter :: TOPO_ORIG_LAT = {set_topo_orig_lat}

  ! Call utm_geo to (maybe) convert lon lat to UTM
  call utm_geo(long, lat, x, y, UTM_PROJECTION_ZONE, IUTM2LONGLAT, &
          SUPPRESS_UTM_PROJECTION)

  ! Get coordinate of corner in bathy/topo model
  icornerlong = &
          int((long - ORIG_LONG_TOPO_SOCAL) / &DEGREES_PER_CELL_TOPO_SOCAL) + 1
  icornerlat = &
          int((lat - ORIG_LAT_TOPO_SOCAL) / DEGREES_PER_CELL_TOPO_SOCAL) + 1

  ! Avoid edge effects and extend with identical point if outside model
  if(icornerlong < 1) icornerlong = 1
  if(icornerlong > NX_TOPO_SOCAL-1) icornerlong = NX_TOPO_SOCAL-1
  if(icornerlat < 1) icornerlat = 1
  if(icornerlat > NY_TOPO_SOCAL-1) icornerlat = NY_TOPO_SOCAL-1

  ! Compute coordinates of corner
  long_corner = ORIG_LONG_TOPO_SOCAL + &
          (icornerlong - 1) * DEGREES_PER_CELL_TOPO_SOCAL
  lat_corner = ORIG_LAT_TOPO_SOCAL + &
          (icornerlat - 1) * DEGREES_PER_CELL_TOPO_SOCAL

  ! Compute ratio for interpolation
  ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO_SOCAL
  ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO_SOCAL

  ! Avoid edge effects
  if(ratio_xi < 0.) ratio_xi = 0.
  if(ratio_xi > 1.) ratio_xi = 1.
  if(ratio_eta < 0.) ratio_eta = 0.
  if(ratio_eta > 1.) ratio_eta = 1.

  ! interpolate elevation at current point
  elevation = &
        itopo(icornerlong, icornerlat) * (1. - ratio_xi) * (1. - ratio_eta) + &
        itopo(icornerlong + 1, icornerlat) * ratio_xi * (1. - ratio_eta) + &
        itopo(icornerlong + 1, icornerlat + 1) * ratio_xi * ratio_eta + &
        itopo(icornerlong, icornerlat + 1) * (1. - ratio_xi) * ratio_eta


end subroutine topo_value


