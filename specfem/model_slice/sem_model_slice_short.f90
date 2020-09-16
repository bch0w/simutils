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
! REQUIRES:
!   exit_mpi.f90
!
! RUBRIC:
!   sem_model_slice xyz_infile model_dir data_name outfile
!   where:
!   xyz_infile: path to the input for a regular grid, which is a 3 column double
!               precision file with columns representing x, y, z respectively.
!
!
! =============================================================================
program sem_model_slice

  implicit none

  include 'mpif.h'
  ! Usually found in specfem3d/setup/constants.h
  include '/home/chowbr/project/bchow/specfem/specfem3d/setup/constants.h'
  ! Usually found in specfem3d/setup/precision.h
  include '/home/chowbr/project/bchow/specfem/specfem3d/setup/precision.h'
  !include 'values_from_mesher.h'

  ! Hard code the max number of GLL points expected
  integer, parameter :: NMAXPTS = 10000000

  integer :: ier, sizeprocs, myrank, ios, i, j, k, ispec, iglob, ipt, npts
  character(len=500) :: xyz_infile, model_dir, data_name, &
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
  !!

  integer, dimension(NX, NY) :: itopo
  double precision :: elevation

  ! MPI initialization
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, sizeprocs, ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)

  ! Input arguments
  call getarg(1, xyz_infile)
  call getarg(2, model_dir)
  call getarg(3, data_name)
  call getarg(4, outfile)

  if(myrank==0)then
    write(*,*) trim(xyz_infile)
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
  write(prname, '(a,i6.6,a)') trim(model_dir)//'/proc',myrank,'_'
  open(unit=27, file=trim(prname)//'external_mesh.bin', status='old', &
          action='read', form='unformatted', iostat=ier)
  read(27) NSPEC_AB
  read(27) NGLOB_AB

  ! Allocates mesh arrays to [xyz]store for unstructured grid
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

  ! Read in chosen model data as variable vstore
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
    print *, 'Writing output file ...'
    open(12, file=outfile, status='unknown')
    do i = 1, npts
      write(12,*) x(i), y(i), z(i), vall(i), distall(i)
    enddo
    close(12)
  endif

  call MPI_FINALIZE(ier)

end program sem_model_slice

!------------------------------------------------------------------------------

subroutine exit_MPI(myrank, error_msg)
  ! ============================================================================
  ! A subroutine that ends the MPI process gracefully and exits MPI
  ! ============================================================================
  implicit none
  include 'mpif.h'
  include '/home/chowbr/project/bchow/specfem/specfem3d/setup/constants.h'
  include '/home/chowbr/project/bchow/specfem/specfem3d/setup/precision.h'

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  integer myrank
  character(len=*) error_msg

  integer ier
  character(len=80) outputname
  character(len=150) OUTPUT_FILES

  ! write error message to screen
  write(*, *) error_msg(1:len(error_msg))
  write(*, *) 'Error detected, aborting MPI... proc ', myrank

  ! write error message to file
  write(outputname, "('/error_message',i6.6,'.txt')") myrank
  open(unit=IERROR, file='OUTPUT_FILES/'//outputname, status='unknown')
  write(IERROR, *) error_msg(1:len(error_msg))
  write(IERROR, *) 'Error detected, aborting MPI... proc ', myrank
  close(IERROR)

  ! close output file
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  ! stop all the MPI processes, and exit
  ! on some machines, MPI_FINALIZE needs to be called before MPI_ABORT
  call MPI_FINALIZE(ier)
  call MPI_ABORT(ier)
  stop 'error, program ended in exit_MPI'

end subroutine exit_MPI

