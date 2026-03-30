! Module exposing:
!   joint_bin_sum       - 2D weighted joint histogram of two variables over time,
!                         weighted by a time-varying field (nlon x nlat x ntime)
!   joint_bin_sumarea   - same as joint_bin_sum but weighted by a static cell-area
!                         field (nlon x nlat), independent of time
!   joint_bin_count     - same as joint_bin_sumarea but with unit weights
!                         (equivalent to counting samples per 2D bin)
!   bin_sum             - 1D weighted histogram of var1 over time,
!                         with var2 accumulated as weight in var1 bins
!   bin_sumarea         - same as bin_sum but weighted by static cell area
!                         (nlon x nlat), independent of time
!   bin_count           - same as bin_sumarea but with unit weights
!                         (equivalent to counting samples per var1 bin)
!   compute_cell_areas  - spherical grid cell areas from lon/lat cell-centre arrays
! bin_index is a private helper and is not accessible outside the module.
module hist
  implicit none
  private
  public :: joint_bin_sum
  public :: joint_bin_sumarea
  public :: joint_bin_count
  public :: bin_sum
  public :: bin_sumarea
  public :: bin_count
  public :: compute_cell_areas

contains

  ! joint_bin_sum
  !
  ! Computes a 2D joint histogram of var1 and var2, weighted by var3.
  ! All three input arrays share the same (nlon, nlat, ntime) shape.
  ! Bin edges for each variable are passed by the caller; the number of bins
  ! is derived as size(var1_edges)-1 and size(var2_edges)-1 respectively.
  ! Values outside the edge range are clipped to the nearest edge before
  ! binning so that no data point is discarded.
  !
  ! Arguments:
  !   var1       - first  variable to bin,  shape (nlon, nlat, ntime)
  !   var2       - second variable to bin,  shape (nlon, nlat, ntime)
  !   var3       - weight field to accumulate, shape (nlon, nlat, ntime)
  !   var1_edges - monotonically increasing bin edges for var1 (nbin_var1+1)
  !   var2_edges - monotonically increasing bin edges for var2 (nbin_var2+1)
  !   varbin     - output 2D histogram, shape (nbin_var1, nbin_var2), allocated here
  subroutine joint_bin_sum(var1, var2, var3, var1_edges, var2_edges, varbin)
    implicit none

    ! Inputs
    real(8), intent(in) :: var1(:,:,:)       ! first  variable  (nlon x nlat x ntime)
    real(8), intent(in) :: var2(:,:,:)       ! second variable  (nlon x nlat x ntime)
    real(8), intent(in) :: var3(:,:,:)       ! weight field     (nlon x nlat x ntime)
    real(8), intent(in) :: var1_edges(:)     ! bin edges for var1 (nbin_var1+1 values)
    real(8), intent(in) :: var2_edges(:)     ! bin edges for var2 (nbin_var2+1 values)

    ! Output
    real(8), allocatable, intent(out) :: varbin(:,:)  ! accumulated histogram (nbin_var1 x nbin_var2)

    ! Local variables
    integer :: i, ix, iy
    integer :: nlon, nlat, ntime
    integer :: nbin_var1, nbin_var2
    integer :: i1, i2          ! bin indices for current grid point
    real(8) :: v1, v2          ! clipped values of var1 / var2

    ! Derive grid dimensions from var1
    nlon  = size(var1, 1)
    nlat  = size(var1, 2)
    ntime = size(var1, 3)

    ! Guard: all three fields must be co-located on the same grid
    if (size(var2, 1) /= nlon .or. size(var2, 2) /= nlat .or. size(var2, 3) /= ntime) then
      error stop 'joint_bin_sum: var2 shape must match var1 shape'
    end if

    if (size(var3, 1) /= nlon .or. size(var3, 2) /= nlat .or. size(var3, 3) /= ntime) then
      error stop 'joint_bin_sum: var3 shape must match var1 shape'
    end if

    if (size(var1_edges) < 2 .or. size(var2_edges) < 2) then
      error stop 'joint_bin_sum: each edges array must have at least 2 values'
    end if

    ! Number of bins = number of edges minus one
    nbin_var1 = size(var1_edges) - 1
    nbin_var2 = size(var2_edges) - 1

    ! Allocate and zero-initialise the output histogram
    allocate(varbin(nbin_var1, nbin_var2))
    varbin = 0.0d0

    ! Loop over time steps and accumulate var3 into the appropriate bin
    do i = 1, ntime
      do iy = 1, nlat
        do ix = 1, nlon
          ! Clip values to the valid edge range before binning
          v1 = max(var1_edges(1), min(var1_edges(size(var1_edges)), var1(ix,iy,i)))
          v2 = max(var2_edges(1), min(var2_edges(size(var2_edges)), var2(ix,iy,i)))

          ! Find the 1-based bin index for each variable
          i1 = bin_index(v1, var1_edges)
          i2 = bin_index(v2, var2_edges)

          ! Accumulate weight only for valid bin indices
          if (i1 >= 1 .and. i1 <= nbin_var1 .and. i2 >= 1 .and. i2 <= nbin_var2) then
            varbin(i1, i2) = varbin(i1, i2) + var3(ix, iy, i)
          end if
        end do
      end do
    end do

  end subroutine joint_bin_sum

  ! joint_bin_sumarea
  !
  ! Identical to joint_bin_sum, except the weight field (cell_area) is
  ! time-independent with shape (nlon, nlat).  The same cell area is therefore
  ! applied at every time step.
  !
  ! Arguments:
  !   var1       - first  variable to bin,  shape (nlon, nlat, ntime)
  !   var2       - second variable to bin,  shape (nlon, nlat, ntime)
  !   cell_area  - static weight field,     shape (nlon, nlat)
  !   var1_edges - monotonically increasing bin edges for var1 (nbin_var1+1)
  !   var2_edges - monotonically increasing bin edges for var2 (nbin_var2+1)
  !   varbin     - output 2D histogram, shape (nbin_var1, nbin_var2), allocated here
  subroutine joint_bin_sumarea(var1, var2, cell_area, var1_edges, var2_edges, varbin)
    implicit none

    ! Inputs
    real(8), intent(in) :: var1(:,:,:)       ! first  variable  (nlon x nlat x ntime)
    real(8), intent(in) :: var2(:,:,:)       ! second variable  (nlon x nlat x ntime)
    real(8), intent(in) :: cell_area(:,:)    ! static weight field (nlon x nlat)
    real(8), intent(in) :: var1_edges(:)     ! bin edges for var1 (nbin_var1+1 values)
    real(8), intent(in) :: var2_edges(:)     ! bin edges for var2 (nbin_var2+1 values)

    ! Output
    real(8), allocatable, intent(out) :: varbin(:,:)  ! accumulated histogram (nbin_var1 x nbin_var2)

    ! Local variables
    integer :: i, ix, iy
    integer :: nlon, nlat, ntime
    integer :: nbin_var1, nbin_var2
    integer :: i1, i2          ! bin indices for current grid point
    real(8) :: v1, v2          ! clipped values of var1 / var2

    ! Derive grid dimensions from var1
    nlon  = size(var1, 1)
    nlat  = size(var1, 2)
    ntime = size(var1, 3)

    ! Guard: var2 must match var1 in all three dimensions
    if (size(var2, 1) /= nlon .or. size(var2, 2) /= nlat .or. size(var2, 3) /= ntime) then
      error stop 'joint_bin_sumarea: var2 shape must match var1 shape'
    end if

    ! Guard: cell_area must match the horizontal grid of var1
    if (size(cell_area, 1) /= nlon .or. size(cell_area, 2) /= nlat) then
      error stop 'joint_bin_sumarea: cell_area horizontal shape must match var1'
    end if

    if (size(var1_edges) < 2 .or. size(var2_edges) < 2) then
      error stop 'joint_bin_sumarea: each edges array must have at least 2 values'
    end if

    ! Number of bins = number of edges minus one
    nbin_var1 = size(var1_edges) - 1
    nbin_var2 = size(var2_edges) - 1

    ! Allocate and zero-initialise the output histogram
    allocate(varbin(nbin_var1, nbin_var2))
    varbin = 0.0d0

    ! Loop over time steps and accumulate the static cell_area into the bin
    ! determined by (var1, var2) at each grid point and time step
    do i = 1, ntime
      do iy = 1, nlat
        do ix = 1, nlon
          ! Clip values to the valid edge range before binning
          v1 = max(var1_edges(1), min(var1_edges(size(var1_edges)), var1(ix,iy,i)))
          v2 = max(var2_edges(1), min(var2_edges(size(var2_edges)), var2(ix,iy,i)))

          ! Find the 1-based bin index for each variable
          i1 = bin_index(v1, var1_edges)
          i2 = bin_index(v2, var2_edges)

          ! Accumulate the static cell area only for valid bin indices
          if (i1 >= 1 .and. i1 <= nbin_var1 .and. i2 >= 1 .and. i2 <= nbin_var2) then
            varbin(i1, i2) = varbin(i1, i2) + cell_area(ix, iy)
          end if
        end do
      end do
    end do

  end subroutine joint_bin_sumarea

  ! joint_bin_count
  !
  ! Identical to joint_bin_sumarea, except every sample has weight 1.0.
  ! This subroutine therefore returns a simple 2D count histogram over
  ! (var1, var2), accumulated across all grid points and time steps.
  !
  ! Arguments:
  !   var1       - first  variable to bin,  shape (nlon, nlat, ntime)
  !   var2       - second variable to bin,  shape (nlon, nlat, ntime)
  !   var1_edges - monotonically increasing bin edges for var1 (nbin_var1+1)
  !   var2_edges - monotonically increasing bin edges for var2 (nbin_var2+1)
  !   varbin     - output 2D histogram, shape (nbin_var1, nbin_var2), allocated here
  subroutine joint_bin_count(var1, var2, var1_edges, var2_edges, varbin)
    implicit none

    ! Inputs
    real(8), intent(in) :: var1(:,:,:)       ! first  variable  (nlon x nlat x ntime)
    real(8), intent(in) :: var2(:,:,:)       ! second variable  (nlon x nlat x ntime)
    real(8), intent(in) :: var1_edges(:)     ! bin edges for var1 (nbin_var1+1 values)
    real(8), intent(in) :: var2_edges(:)     ! bin edges for var2 (nbin_var2+1 values)

    ! Output
    real(8), allocatable, intent(out) :: varbin(:,:)  ! count histogram (nbin_var1 x nbin_var2)

    ! Local variables
    integer :: i, ix, iy
    integer :: nlon, nlat, ntime
    integer :: nbin_var1, nbin_var2
    integer :: i1, i2          ! bin indices for current grid point
    real(8) :: v1, v2          ! clipped values of var1 / var2

    ! Derive grid dimensions from var1
    nlon  = size(var1, 1)
    nlat  = size(var1, 2)
    ntime = size(var1, 3)

    ! Guard: var2 must match var1 in all three dimensions
    if (size(var2, 1) /= nlon .or. size(var2, 2) /= nlat .or. size(var2, 3) /= ntime) then
      error stop 'joint_bin_count: var2 shape must match var1 shape'
    end if

    if (size(var1_edges) < 2 .or. size(var2_edges) < 2) then
      error stop 'joint_bin_count: each edges array must have at least 2 values'
    end if

    ! Number of bins = number of edges minus one
    nbin_var1 = size(var1_edges) - 1
    nbin_var2 = size(var2_edges) - 1

    ! Allocate and zero-initialise the output histogram
    allocate(varbin(nbin_var1, nbin_var2))
    varbin = 0.0d0

    ! Loop over time steps and increment the selected bin by one
    do i = 1, ntime
      do iy = 1, nlat
        do ix = 1, nlon
          ! Clip values to the valid edge range before binning
          v1 = max(var1_edges(1), min(var1_edges(size(var1_edges)), var1(ix,iy,i)))
          v2 = max(var2_edges(1), min(var2_edges(size(var2_edges)), var2(ix,iy,i)))

          ! Find the 1-based bin index for each variable
          i1 = bin_index(v1, var1_edges)
          i2 = bin_index(v2, var2_edges)

          ! Increment count only for valid bin indices
          if (i1 >= 1 .and. i1 <= nbin_var1 .and. i2 >= 1 .and. i2 <= nbin_var2) then
            varbin(i1, i2) = varbin(i1, i2) + 1.0d0
          end if
        end do
      end do
    end do

  end subroutine joint_bin_count

  ! bin_sum
  !
  ! Similar to joint_bin_sum, but bins only var1 and accumulates var2 as the
  ! weight in each var1 bin.  The output is therefore a 1D weighted histogram.
  !
  ! Arguments:
  !   var1       - variable to bin, shape (nlon, nlat, ntime)
  !   var2       - weight field to accumulate, shape (nlon, nlat, ntime)
  !   var1_edges - monotonically increasing bin edges for var1 (nbin_var1+1)
  !   varbin     - output 1D histogram, length nbin_var1, allocated here
  subroutine bin_sum(var1, var2, var1_edges, varbin)
    implicit none

    ! Inputs
    real(8), intent(in) :: var1(:,:,:)       ! variable to bin (nlon x nlat x ntime)
    real(8), intent(in) :: var2(:,:,:)       ! weight field    (nlon x nlat x ntime)
    real(8), intent(in) :: var1_edges(:)     ! bin edges for var1 (nbin_var1+1 values)

    ! Output
    real(8), allocatable, intent(out) :: varbin(:)    ! weighted histogram (nbin_var1)

    ! Local variables
    integer :: i, ix, iy
    integer :: nlon, nlat, ntime
    integer :: nbin_var1
    integer :: i1
    real(8) :: v1

    ! Derive grid dimensions from var1
    nlon  = size(var1, 1)
    nlat  = size(var1, 2)
    ntime = size(var1, 3)

    ! Guard: var2 must match var1 in all three dimensions
    if (size(var2, 1) /= nlon .or. size(var2, 2) /= nlat .or. size(var2, 3) /= ntime) then
      error stop 'bin_sum: var2 shape must match var1 shape'
    end if

    if (size(var1_edges) < 2) then
      error stop 'bin_sum: var1_edges must have at least 2 values'
    end if

    ! Number of bins = number of edges minus one
    nbin_var1 = size(var1_edges) - 1

    ! Allocate and zero-initialise the output histogram
    allocate(varbin(nbin_var1))
    varbin = 0.0d0

    ! Loop over time steps and accumulate var2 into var1 bins
    do i = 1, ntime
      do iy = 1, nlat
        do ix = 1, nlon
          ! Clip var1 to the valid edge range before binning
          v1 = max(var1_edges(1), min(var1_edges(size(var1_edges)), var1(ix,iy,i)))

          ! Find the 1-based bin index for var1
          i1 = bin_index(v1, var1_edges)

          ! Accumulate weight only for valid bin indices
          if (i1 >= 1 .and. i1 <= nbin_var1) then
            varbin(i1) = varbin(i1) + var2(ix, iy, i)
          end if
        end do
      end do
    end do

  end subroutine bin_sum

  ! bin_sumarea
  !
  ! Similar to joint_bin_sumarea, but bins only var1 and accumulates a static
  ! cell-area weight for each (ix, iy) across all time steps.
  !
  ! Arguments:
  !   var1       - variable to bin, shape (nlon, nlat, ntime)
  !   cell_area  - static weight field, shape (nlon, nlat)
  !   var1_edges - monotonically increasing bin edges for var1 (nbin_var1+1)
  !   varbin     - output 1D histogram, length nbin_var1, allocated here
  subroutine bin_sumarea(var1, cell_area, var1_edges, varbin)
    implicit none

    ! Inputs
    real(8), intent(in) :: var1(:,:,:)       ! variable to bin (nlon x nlat x ntime)
    real(8), intent(in) :: cell_area(:,:)    ! static weight   (nlon x nlat)
    real(8), intent(in) :: var1_edges(:)     ! bin edges for var1 (nbin_var1+1 values)

    ! Output
    real(8), allocatable, intent(out) :: varbin(:)    ! weighted histogram (nbin_var1)

    ! Local variables
    integer :: i, ix, iy
    integer :: nlon, nlat, ntime
    integer :: nbin_var1
    integer :: i1
    real(8) :: v1

    ! Derive grid dimensions from var1
    nlon  = size(var1, 1)
    nlat  = size(var1, 2)
    ntime = size(var1, 3)

    ! Guard: cell_area must match horizontal dimensions of var1
    if (size(cell_area, 1) /= nlon .or. size(cell_area, 2) /= nlat) then
      error stop 'bin_sumarea: cell_area horizontal shape must match var1'
    end if

    if (size(var1_edges) < 2) then
      error stop 'bin_sumarea: var1_edges must have at least 2 values'
    end if

    ! Number of bins = number of edges minus one
    nbin_var1 = size(var1_edges) - 1

    ! Allocate and zero-initialise the output histogram
    allocate(varbin(nbin_var1))
    varbin = 0.0d0

    ! Loop over time steps and accumulate static cell_area into var1 bins
    do i = 1, ntime
      do iy = 1, nlat
        do ix = 1, nlon
          ! Clip var1 to the valid edge range before binning
          v1 = max(var1_edges(1), min(var1_edges(size(var1_edges)), var1(ix,iy,i)))

          ! Find the 1-based bin index for var1
          i1 = bin_index(v1, var1_edges)

          ! Accumulate weight only for valid bin indices
          if (i1 >= 1 .and. i1 <= nbin_var1) then
            varbin(i1) = varbin(i1) + cell_area(ix, iy)
          end if
        end do
      end do
    end do

  end subroutine bin_sumarea

  ! bin_count
  !
  ! Similar to joint_bin_count, but bins only var1 and increments by one for
  ! each sample.  The result is a simple count histogram over var1.
  !
  ! Arguments:
  !   var1       - variable to bin, shape (nlon, nlat, ntime)
  !   var1_edges - monotonically increasing bin edges for var1 (nbin_var1+1)
  !   varbin     - output 1D histogram, length nbin_var1, allocated here
  subroutine bin_count(var1, var1_edges, varbin)
    implicit none

    ! Inputs
    real(8), intent(in) :: var1(:,:,:)       ! variable to bin (nlon x nlat x ntime)
    real(8), intent(in) :: var1_edges(:)     ! bin edges for var1 (nbin_var1+1 values)

    ! Output
    real(8), allocatable, intent(out) :: varbin(:)    ! count histogram (nbin_var1)

    ! Local variables
    integer :: i, ix, iy
    integer :: nlon, nlat, ntime
    integer :: nbin_var1
    integer :: i1
    real(8) :: v1

    ! Derive grid dimensions from var1
    nlon  = size(var1, 1)
    nlat  = size(var1, 2)
    ntime = size(var1, 3)

    if (size(var1_edges) < 2) then
      error stop 'bin_count: var1_edges must have at least 2 values'
    end if

    ! Number of bins = number of edges minus one
    nbin_var1 = size(var1_edges) - 1

    ! Allocate and zero-initialise the output histogram
    allocate(varbin(nbin_var1))
    varbin = 0.0d0

    ! Loop over time steps and increment selected var1 bin by one
    do i = 1, ntime
      do iy = 1, nlat
        do ix = 1, nlon
          ! Clip var1 to the valid edge range before binning
          v1 = max(var1_edges(1), min(var1_edges(size(var1_edges)), var1(ix,iy,i)))

          ! Find the 1-based bin index for var1
          i1 = bin_index(v1, var1_edges)

          ! Increment count only for valid bin indices
          if (i1 >= 1 .and. i1 <= nbin_var1) then
            varbin(i1) = varbin(i1) + 1.0d0
          end if
        end do
      end do
    end do

  end subroutine bin_count

  ! bin_index
  !
  ! Returns the 1-based bin index of x within the sorted edges array using
  ! binary search. Returns 0 if x is outside [edges(1), edges(end)].
  ! The last bin is closed on both sides: x == edges(n+1) maps to bin n.
  integer function bin_index(x, edges)
    implicit none
    real(8), intent(in) :: x        ! value to locate
    real(8), intent(in) :: edges(:) ! monotonically increasing bin edges
    integer :: lo, hi, mid, n

    n = size(edges) - 1  ! number of bins

    ! Out-of-range check
    if (x < edges(1) .or. x > edges(n+1)) then
      bin_index = 0
      return
    end if

    ! Special case: x falls exactly on the last edge -> last bin
    if (x == edges(n+1)) then
      bin_index = n
      return
    end if

    ! Binary search for the interval [edges(mid), edges(mid+1)) containing x
    lo = 1
    hi = n
    do while (lo <= hi)
      mid = (lo + hi) / 2
      if (x >= edges(mid) .and. x < edges(mid+1)) then
        bin_index = mid
        return
      else if (x < edges(mid)) then
        hi = mid - 1
      else
        lo = mid + 1
      end if
    end do

    bin_index = 0  ! should not be reached for in-range values
  end function bin_index

  ! compute_cell_areas
  !
  ! Computes the area of every grid cell on a perfect sphere of radius R_earth,
  ! given 1D arrays of cell-centre longitudes and latitudes in degrees.
  !
  ! Cell boundaries are defined as midpoints between adjacent cell centres.
  ! The outermost boundaries are obtained by symmetric extrapolation from the
  ! two nearest centres.  Latitude edges are clamped to [-90, +90] degrees to
  ! avoid exceeding the poles.
  !
  ! The area of a spherical cell is:
  !   A(i,j) = R^2 * |Δλ(i)| * |sin(φ_north(j)) - sin(φ_south(j))|
  ! where Δλ is the longitudinal width and φ are latitude boundaries, all in
  ! radians.  This expression is exact for a sphere.
  !
  ! Arguments:
  !   lon       - 1D array of cell-centre longitudes [degrees], length nlon (>= 2)
  !   lat       - 1D array of cell-centre latitudes  [degrees], length nlat (>= 2)
  !   cell_area - output cell areas [m^2], shape (nlon, nlat), allocated here
  subroutine compute_cell_areas(lon, lat, cell_area)
    implicit none

    ! Inputs
    real(8), intent(in) :: lon(:)   ! cell-centre longitudes [degrees]
    real(8), intent(in) :: lat(:)   ! cell-centre latitudes  [degrees]

    ! Output
    real(8), allocatable, intent(out) :: cell_area(:,:)  ! cell areas [m^2], (nlon x nlat)

    ! Physical and mathematical constants
    real(8), parameter :: R_earth = 6.371d6        ! mean radius of the Earth [m]
    real(8), parameter :: pi      = acos(-1.0d0)   ! pi
    real(8), parameter :: deg2rad = pi / 180.0d0   ! degrees-to-radians conversion factor

    ! Local variables
    integer :: i, j, nlon, nlat
    real(8), allocatable :: lon_edges(:)  ! longitude cell boundaries [radians], length nlon+1
    real(8), allocatable :: lat_edges(:)  ! latitude  cell boundaries [radians], length nlat+1
    real(8) :: dlon              ! longitudinal width of the current cell [radians]
    real(8) :: sin_south, sin_north  ! sine of the southern / northern latitude boundaries

    nlon = size(lon)
    nlat = size(lat)

    if (nlon < 2 .or. nlat < 2) then
      error stop 'compute_cell_areas: lon and lat must each contain at least 2 values'
    end if

    ! --- Longitude edges [radians] ---
    ! Interior edges: midpoint of two adjacent cell centres.
    ! Boundary edges: symmetric extrapolation from the two outermost centres.
    allocate(lon_edges(nlon+1))
    lon_edges(1)      = deg2rad * (lon(1)    - 0.5d0 * (lon(2)    - lon(1)))
    do i = 2, nlon
      lon_edges(i)    = deg2rad * 0.5d0 * (lon(i-1) + lon(i))
    end do
    lon_edges(nlon+1) = deg2rad * (lon(nlon) + 0.5d0 * (lon(nlon) - lon(nlon-1)))

    ! --- Latitude edges [radians] ---
    ! Same midpoint / extrapolation convention; edges are clamped to ±pi/2
    ! to prevent the grid from exceeding the poles.
    allocate(lat_edges(nlat+1))
    lat_edges(1)      = deg2rad * (lat(1)    - 0.5d0 * (lat(2)    - lat(1)))
    lat_edges(1)      = max(-0.5d0 * pi, lat_edges(1))   ! clamp to south pole
    do j = 2, nlat
      lat_edges(j)    = deg2rad * 0.5d0 * (lat(j-1) + lat(j))
    end do
    lat_edges(nlat+1) = deg2rad * (lat(nlat) + 0.5d0 * (lat(nlat) - lat(nlat-1)))
    lat_edges(nlat+1) = min( 0.5d0 * pi, lat_edges(nlat+1))  ! clamp to north pole

    ! --- Area computation ---
    ! A(i,j) = R^2 * |delta_lon| * |sin(lat_north) - sin(lat_south)|
    allocate(cell_area(nlon, nlat))
    do j = 1, nlat
      sin_south = sin(lat_edges(j))
      sin_north = sin(lat_edges(j+1))
      do i = 1, nlon
        dlon = abs(lon_edges(i+1) - lon_edges(i))  ! always positive by construction
        cell_area(i, j) = R_earth**2 * dlon * abs(sin_north - sin_south)
      end do
    end do

  end subroutine compute_cell_areas

end module hist
