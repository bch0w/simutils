!========================================================================
!
! Venus 1-D model reader for SPECFEM2D (MODEL = external)
!
! Drop-in replacement for src/specfem2D/define_external_model.f90.
! Install once, rebuild xspecfem2D, and from then on change the model by
! editing the text file DATA/venus_model.txt (written by
! make_venus_run.py). No more Fortran edits for model changes.
!
! File format (lines starting with '#' are comments):
!   first data line : number of rows N
!   next N lines    : region  radius[m]  rho[kg/m3]  vp[m/s]  vs[m/s]
!                     qkappa  qmu
!
! Rows must be grouped by region and sorted by increasing radius inside
! each region. "region" is the material flag the mesher writes for each
! element (MESH/Material_*). A GLL point is only ever interpolated
! between rows of its own element's region, so values are never mixed
! across a discontinuity (core-mantle boundary, ground-air interface).
! Radii outside a region's range are clamped to its end values.
!
!========================================================================

  subroutine define_external_model(coord,material_element,ibool, &
                                   rho,vp,vs,QKappa_attenuation, &
                                   Qmu_attenuation, &
                                   c11,c12,c13,c15,c23,c25,c33,c35,c55)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN,IIN, &
                       IN_DATA_FILES,MAX_STRING_LEN

  use specfem_par, only: myrank,nspec,nglob

  implicit none

  double precision, dimension(NDIM,nglob), intent(in) :: coord
  integer, dimension(nspec), intent(in) :: material_element
  integer, dimension(NGLLX,NGLLZ,nspec), intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: &
    rho,vp,vs
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: &
    QKappa_attenuation,Qmu_attenuation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: &
    c11,c15,c13,c33,c35,c55,c12,c23,c25

  ! largest region flag allowed in the table
  integer, parameter :: MAX_REGION = 20

  character(len=MAX_STRING_LEN) :: filename, line
  integer :: nrow, irow, ier, ireg, i, j, ispec, iglob, k, klo, khi
  integer, dimension(MAX_REGION) :: ifirst, ilast
  integer, dimension(:), allocatable :: reg
  double precision, dimension(:), allocatable :: rad, trho, tvp, tvs, &
                                                 tqk, tqm
  double precision :: r, frac

  filename = trim(IN_DATA_FILES)//'venus_model.txt'

  if (myrank == 0) then
    write(IMAIN,*) '  external model: Venus 1-D table ',trim(filename)
    call flush_IMAIN()
  endif

  ! ---- read the table -------------------------------------------------
  open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) call stop_the_code('Error opening DATA/venus_model.txt')

  ! skip comments, then read the number of rows
  nrow = 0
  do
    read(IIN,'(a)',iostat=ier) line
    if (ier /= 0) call stop_the_code('venus_model.txt: missing row count')
    line = adjustl(line)
    if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
    read(line,*,iostat=ier) nrow
    if (ier /= 0 .or. nrow < 2) &
      call stop_the_code('venus_model.txt: bad row count')
    exit
  enddo

  allocate(reg(nrow),rad(nrow),trho(nrow),tvp(nrow),tvs(nrow), &
           tqk(nrow),tqm(nrow))

  irow = 0
  do while (irow < nrow)
    read(IIN,'(a)',iostat=ier) line
    if (ier /= 0) call stop_the_code('venus_model.txt: too few rows')
    line = adjustl(line)
    if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
    irow = irow + 1
    read(line,*,iostat=ier) reg(irow),rad(irow),trho(irow),tvp(irow), &
                            tvs(irow),tqk(irow),tqm(irow)
    if (ier /= 0) call stop_the_code('venus_model.txt: unreadable row')
    if (reg(irow) < 1 .or. reg(irow) > MAX_REGION) &
      call stop_the_code('venus_model.txt: region flag out of range')
  enddo
  close(IIN)

  ! first/last row of each region (rows of one region must be contiguous)
  ifirst(:) = 0
  ilast(:) = 0
  do irow = 1,nrow
    ireg = reg(irow)
    if (ifirst(ireg) == 0) then
      ifirst(ireg) = irow
    else if (ilast(ireg) /= irow-1) then
      call stop_the_code('venus_model.txt: region rows not contiguous')
    endif
    ilast(ireg) = irow
  enddo

  ! ---- assign values at every GLL point ------------------------------
  c11(:,:,:) = 0._CUSTOM_REAL
  c13(:,:,:) = 0._CUSTOM_REAL
  c15(:,:,:) = 0._CUSTOM_REAL
  c33(:,:,:) = 0._CUSTOM_REAL
  c35(:,:,:) = 0._CUSTOM_REAL
  c55(:,:,:) = 0._CUSTOM_REAL
  c12(:,:,:) = 0._CUSTOM_REAL
  c23(:,:,:) = 0._CUSTOM_REAL
  c25(:,:,:) = 0._CUSTOM_REAL

  do ispec = 1,nspec
    ireg = material_element(ispec)
    if (ireg < 1 .or. ireg > MAX_REGION) &
      call stop_the_code('mesh region flag not found in venus_model.txt')
    if (ifirst(ireg) == 0) &
      call stop_the_code('mesh region flag not found in venus_model.txt')
    klo = ifirst(ireg)
    khi = ilast(ireg)

    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        r = sqrt(coord(1,iglob)**2 + coord(2,iglob)**2)

        if (r <= rad(klo) .or. klo == khi) then
          k = klo + 1
          frac = 0.d0
        else if (r >= rad(khi)) then
          k = khi
          frac = 1.d0
        else
          k = klo + 1
          do while (rad(k) < r)
            k = k + 1
          enddo
          if (rad(k) > rad(k-1)) then
            frac = (r - rad(k-1)) / (rad(k) - rad(k-1))
          else
            frac = 1.d0
          endif
        endif
        if (klo == khi) k = klo

        if (k == klo) then
          rho(i,j,ispec) = real(trho(k),kind=CUSTOM_REAL)
          vp(i,j,ispec) = real(tvp(k),kind=CUSTOM_REAL)
          vs(i,j,ispec) = real(tvs(k),kind=CUSTOM_REAL)
          QKappa_attenuation(i,j,ispec) = real(tqk(k),kind=CUSTOM_REAL)
          Qmu_attenuation(i,j,ispec) = real(tqm(k),kind=CUSTOM_REAL)
        else
          rho(i,j,ispec) = real(trho(k-1) + frac*(trho(k)-trho(k-1)), &
                                kind=CUSTOM_REAL)
          vp(i,j,ispec) = real(tvp(k-1) + frac*(tvp(k)-tvp(k-1)), &
                               kind=CUSTOM_REAL)
          vs(i,j,ispec) = real(tvs(k-1) + frac*(tvs(k)-tvs(k-1)), &
                               kind=CUSTOM_REAL)
          QKappa_attenuation(i,j,ispec) = &
            real(tqk(k-1) + frac*(tqk(k)-tqk(k-1)),kind=CUSTOM_REAL)
          Qmu_attenuation(i,j,ispec) = &
            real(tqm(k-1) + frac*(tqm(k)-tqm(k-1)),kind=CUSTOM_REAL)
        endif
      enddo
    enddo
  enddo

  if (myrank == 0) then
    write(IMAIN,*) '  Venus table rows read: ',nrow
    write(IMAIN,*) '  rho min/max (kg/m3)  : ',minval(rho),maxval(rho)
    write(IMAIN,*) '  vp  min/max (m/s)    : ',minval(vp),maxval(vp)
    write(IMAIN,*) '  Qkappa min           : ',minval(QKappa_attenuation)
    call flush_IMAIN()
  endif

  deallocate(reg,rad,trho,tvp,tvs,tqk,tqm)

  end subroutine define_external_model
