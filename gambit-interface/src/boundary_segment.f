module Boundary_Segment

  use Constants, only: RNP
  use Basic_IO, only: new_unit
  use Spectral_Element_Vertex
  use Spectral_Element_Edge
  use Boundary_Curve
  use Boundary_Condition
  implicit none
  private

  type, public :: BoundaryEdge
     integer :: e = 0  ! associated spectral element edge
  end type BoundaryEdge

  type, public :: BoundarySegment
    character(len=25) :: name = ''   ! boundary name (setup by import)
    integer           :: curve = 0   ! associated boundary curve
    integer           :: bcond = 0   ! associated boundary condition
    type(BoundaryEdge), pointer :: edge(:) => null()
  end type BoundarySegment

  public :: ConstructBoundarySegments
  real(RNP), parameter :: def_epsilon = 1e-5_RNP

! Conventions
! 1. The list edge(:)%e is sorted, that means edge(i)%e is the neighbor
!    of edge(i+1)%e
! 2. Because the element direction is defined in an anticlockwise direction
!    and the edge vertices are defined by the associated element with the lower
!    element id the curve parameter "t" grows continuously. That means:
!    edge(i)%geo%t(2) > edge(i)%geo%t(1)
!    edge(bseg%edge(i)%e)%geo%t > edge(bseg%edge(i-1)%e)%geo%t.
! 3. The list bseg%edge(:)%e is defined in an anticlockwise direction too.

contains

!===============================================================================

subroutine ConstructBoundarySegments(vertex,edge,bseg,curve,bcond)
  type(SpectralElementVertex), intent(inout)	:: vertex(:)
  type(SpectralElementEdge), intent(inout)	:: edge(:)
  type(BoundarySegment), intent(inout)		:: bseg(:)
  type(BoundaryCurve), pointer			:: curve(:)
  type(BoundaryCondition), pointer		:: bcond(:)
  integer, allocatable				:: bList(:)
  integer					:: i, j

  write(*,'(A)',advance='NO') '  Construct boundaries from edges ...'

! init boundary curves =======================================================

  if (size(bseg) < size(ControlBoundaryCurve)) then
    STOP 'Error: More than one curve per boundary is not allowed yet.'
  else
    allocate( curve(size(bseg)) )
    curve%name = bseg%name
  end if

  ! assign control data to boundary curves
  ! has to be modified if number of curves > number of boundary segments
  do i = 1, size(curve)
    bseg(i)%curve = i
    do j = 1, size(ControlBoundaryCurve)
      if (curve(i)%name  == ControlBoundaryCurve(j)%name) &
      curve(i) = ControlBoundaryCurve(j)
    end do
  end do

  ! delete boundary curve control data
  deallocate( ControlBoundaryCurve )

  ! setup boundary curve link by boundary name and user definitions
  ! has to be modified if number of curves > number of boundary segments
  do i = 1, size(bseg)
    bseg(i)%curve = i
  end do

! setup edges ================================================================

  ! count associated boundary edges
  allocate( bList(1:size(bseg)) )
  bList = 0
  do i = 1, size(edge)
    if (edge(i)%b(1) /= 0) bList(edge(i)%b(1)) = bList(edge(i)%b(1)) + 1
  end do

  ! assign edges to boundary segments
  do i = 1, size(bseg)
    allocate( bseg(i)%edge(bList(i)) )
  end do
  bList = 0
  do i = 1, size(edge)
    if (edge(i)%b(1) /= 0) then
      bseg(edge(i)%b(1))%edge(bList(edge(i)%b(1))+1)%e = i
      bList(edge(i)%b(1)) = bList(edge(i)%b(1)) + 1
    end if
  end do

  ! sort edge list if necessary
  do i = 1, size(bseg)
    do j = 1, size(bseg(i)%edge) - 1
      if ( SharedNode(edge(bseg(i)%edge(j)%e),         &
                      edge(bseg(i)%edge(j+1)%e)) == 0) &
        call SortEdgeList(edge,bseg(i))
    end do
  end do

  ! setup edge list direction
  do i = 1, size(bseg)
    call SetupEdgeListDirection(edge,bseg(i))
  end do

  ! check edge - boundary correlation
  do i = 1, size(bseg)
    do j = 1, size(bseg(i)%edge)
      edge(bseg(i)%edge(j)%e)%b(2) = j
      if (edge(bseg(i)%edge(j)%e)%b(1) /= i) &
         STOP 'Error: Edge%b(1) isnt equal boundary segment.'
    end do
  end do

  ! setup boundary segment curve parameter
  do i = 1, size(bseg)
    call SetupBSCurveParameter(vertex,edge,bseg(i))
  end do

  ! setup boundary segment "curved" flag
  do i = 1, size(bseg)
    do j = 1, size(bseg(i)%edge)
      if (curve(i)%geom > 1) edge(bseg(i)%edge(j)%e)%geo%curved = .true.
    end do
  end do

  write(*,'(A)') '  done'

  write(*,'(A)') repeat('*',45)
  write(*,'(A)') '  Setup boundary geometry'
  write(*,'(A,/)') repeat('*',45)

  ! setup boundary segment control points
  ! has to be modified if number of curves > number of boundary segments
  do i = 1, size(bseg)
    write(*,'(3A)',advance='NO') '  ', bseg(i)%Name(1:10), ': '
    call SetupBSControlPoints(vertex,edge,bseg(i),curve(i))
    call SetupBSSpline(curve(i))
  end do

  ! move boundary edge vertices onto the corresponding spline
  ! has to be modified if number of curves > number of boundary segments
  do i = 1, size(bseg)
    call SetupBSVertices(vertex,edge,bseg(i),curve(i))
  end do

  write(*,'(/,A)') repeat('*',45)
  write(*,'(A)',advance='YES') '  done'
  write(*,'(A)') repeat('*',45)

! init boundary conditions ===================================================

  write(*,'(A)') '  Setup boundary conditions'
  write(*,'(A,/)') repeat('*',45)

  if (size(bseg) < size(ControlBoundaryCondition)) then
    STOP 'Error: More than one bc per boundary is not allowed yet.'
  else
    allocate( bcond(size(bseg)) )
    bcond%name = bseg%name
  end if

  ! assign control data to boundary conditions
  ! has to be modified if number of conditions > number of boundary segments
  do i = 1, size(bcond)
    bseg(i)%bcond = i
    do j = 1, size(ControlBoundaryCondition)
      if (bcond(i)%name  == ControlBoundaryCondition(j)%name) &
      bcond(i) = ControlBoundaryCondition(j)
    end do
  end do

  ! delete boundary condition control data
  deallocate( ControlBoundaryCondition )

  ! setup boundary condition link by boundary name and user definitions
  ! has to be modified if number of conditions > number of boundary segments
  do i = 1, size(bseg)
    bseg(i)%bcond = i
  end do

  do i = 1, size(bseg)
!    write(*,'(3A)',advance='NO') '  ', bseg(i)%Name(1:10), ': '

    write(*, '(4A)',advance='YES') 'boundary segment: ', bcond(i)%name(1:25),&
          &'      key: ', bcond(i)%key

!    select case(bcond(i)%type)
!    case(0)
!      write(*,'(A)',advance='YES') 'no treatment'
!    case(1)
!      write(*,'(A)',advance='YES') 'reflective wall'
!    case default
!      STOP 'Error: This type of boundary condition is not implemented yet.'
!    end select

  end do

  write(*,'(/,A)') repeat('*',45)
  write(*,'(A)',advance='YES') '  done'
  write(*,'(A)') repeat('*',45)

end subroutine ConstructBoundarySegments

!===============================================================================

subroutine SetupBSCurveParameter(vertex,edge,bseg)
  type(SpectralElementVertex), intent(in)	::  vertex(:)
  type(SpectralElementEdge), intent(inout)	::  edge(:)
  type(BoundarySegment), intent(inout)		::  bseg
  real(RNP)					::  t(0:size(bseg%Edge))
  real(RNP)					::  ts, x(2,2)
  integer					::  i, K

  K = size(bseg%Edge)
  t = 0
  do i = 1, K
    x(1,:) = vertex(edge(bseg%edge(i)%e)%v(1))%x
    x(2,:) = vertex(edge(bseg%edge(i)%e)%v(2))%x
    ts     = sqrt( (x(2,1)-x(1,1))**2 + (x(2,2)-x(1,2))**2 )
    t(i)   = t(i-1) + ts
  end do
  t = t / t(K)

  if (K == 1) then
    edge(bseg%edge(1)%e)%geo%t(1) = 0
    edge(bseg%edge(1)%e)%geo%t(2) = 1
    return
  end if

  forall(i=1:K)
    edge(bseg%edge(i)%e)%geo%t(1) = t(i-1)
    edge(bseg%edge(i)%e)%geo%t(2) = t(i)
  end forall

end subroutine SetupBSCurveParameter

!===============================================================================

subroutine SetupBSControlPoints(vertex,edge,bseg,curve)
  type(SpectralElementVertex), intent(in) :: vertex(:)
  type(SpectralElementEdge), intent(in)	  :: edge(:)
  type(BoundarySegment), intent(in)	  :: bseg
  type(BoundaryCurve), intent(inout)	  :: curve
  integer				  :: i, n, unit, ios, K, e
  real(RNP)				  :: x(2), d(2), eps = def_epsilon
  real(RNP), allocatable		  :: xc(:,:)
  character(len=25)			  :: file

  n = 0
  if (LEN_TRIM(curve%file) > 0) n = 1

  if (size(bseg%edge) == 1) then
    allocate(curve%xc(0:1,2))
    e = bseg%edge(1)%e
    curve%xc(0,:) = vertex(edge(e)%v(1))%x
    curve%xc(1,:) = vertex(edge(e)%v(2))%x
    write(*,'(A)',advance='YES') 'use gridpoints'
    return
  end if

  select case(n)
  case(0) 	! use grid points as control points
  write(*,'(A)',advance='YES') 'use gridpoints'
    K = size(bseg%edge)
    allocate(curve%xc(0:K,2))
    curve%xc(0,:) = vertex(edge(bseg%edge(1)%e)%v(1))%x
    do i = 1, K
      curve%xc(i,:) = vertex(edge(bseg%edge(i)%e)%v(2))%x
    end do

  case(1)	! use control points from file
    file = curve%file
    write(*,'(2A)',advance='YES') 'use file: ', file
    K = - 1
    unit = new_unit()
    open(unit, IOSTAT=ios, file=file, action='READ', status='OLD')
    if (ios /= 0) STOP 'An error occurred while opening the file.'
    do
      read(unit, *, iostat=ios)
      if (ios /= 0) exit
      K = K + 1
    end do
    rewind(unit)
    allocate( curve%xc(0:K,2), xc(0:K,2) )
    do i = 0, K
      read(unit, *) xc(i,1), xc(i,2)
    end do
    close(unit)

    ! setup control points direction
    x    = vertex(edge(bseg%edge(1)%e)%v(1))%x
    d(1) = sqrt((xc(0,1)-x(1))**2 + (xc(0,2)-x(2))**2)
    x    = vertex(edge(bseg%edge(size(bseg%edge))%e)%v(2))%x
    d(2) = sqrt((xc(0,1)-x(1))**2 + (xc(0,2)-x(2))**2)
    if (d(1) <= eps) forall(i=0:K) curve%xc(i,:) = xc(i,:)
    if (d(2) <= eps) forall(i=K:0:-1) curve%xc(i,:) = xc(K-i,:)
    if (ALL(d(:) <= eps)) then
      write(*,'(2A)') repeat(' ',14), '(closed curve)'
      forall(i=0:K) curve%xc(i,:) = xc(i,:)
    else
      if (d(1) <= eps) forall(i=0:K) curve%xc(i,:) = xc(i,:)
      if (d(2) <= eps) forall(i=K:0:-1) curve%xc(i,:) = xc(K-i,:)
    end if
    if (minval(d) > eps) &
      STOP 'Error: Points from file doesnt fit to gridpoints.'
  end select

end subroutine SetupBSControlPoints

!===============================================================================

subroutine SortEdgeList(edge,bseg)
  type(SpectralElementEdge), intent(in)   :: edge(:)
  type(BoundarySegment), intent(inout)    :: bseg
  integer				  :: i, j, v
  integer, allocatable			  :: e(:), n(:)
  logical				  :: inner

  allocate(e(size(bseg%edge)),n(size(bseg%edge)))

  ! search the vertex which is not owned by two edges and use the corresponding
  ! edge as start edge
  e = bseg%edge(:)%e
  n = 0
  do i = 1, size(e)
    v = edge(e(i))%v(1)
    inner = .false.
    do j = 1, size(e)
      if (i /= j) then
        if (SharedNode(edge(e(i)),edge(e(j))) == v) then
          inner = .true.
          exit
        end if
      end if
    end do
    if (inner .EQV. .false.) then
      n(1) = e(i)
      e(i) = 0
      exit
    end if
  end do
  if (n(1) == 0) then ! closed curve
    n(1) = e(1)
    e(1) = 0
  end if

  ! sort the edges by shared nodes
  do i = 1, size(e) - 1
    do j = 1, size(e)
      if (e(j) /= 0 .AND. n(i) /= e(j)) then
        if (SharedNode(edge(n(i)),edge(e(j))) /= 0) then
          n(i+1) = e(j)
          e(j) = 0
          exit
        end if
      end if
    end do
    if (n(i+1) == 0) STOP 'Error: Boundary edges are not contiguous.'
  end do
  bseg%edge(:)%e = n
  write(*,'(/,/,3A,/)') '  ', bseg%name(1:10), ': Boundary edges resorted'

end subroutine SortEdgeList

!===============================================================================

subroutine SetupEdgeListDirection(edge,bseg)
  type(SpectralElementEdge), intent(in)   :: edge(:)
  type(BoundarySegment), intent(inout)    :: bseg
  integer				  :: i, e(2), v
  integer, allocatable			  :: elist(:)
  if (size(bseg%edge) < 2) return
  e(1) = bseg%edge(1)%e
  e(2) = bseg%edge(2)%e
  v    = edge(e(1))%v(1)
  if (SharedNode(edge(e(1)),edge(e(2))) == v) then
    allocate(elist(size(bseg%edge)))
    eList = bseg%edge(:)%e
    forall(i=1:size(elist)) bseg%edge(i)%e = elist(size(elist)-i+1)
    deallocate(elist)
  end if
end subroutine SetupEdgeListDirection

!===============================================================================

subroutine SetupBSVertices(vertex,edge,bseg,curve)
  type(SpectralElementVertex), intent(inout)	:: vertex(:)
  type(SpectralElementEdge), intent(inout)	:: edge(:)
  type(BoundarySegment), intent(inout)		:: bseg
  type(BoundaryCurve), intent(inout)		:: curve
  integer					:: i, v(2)
  real(RNP)					:: x(2), xi
  do i = 1, size(bseg%Edge)
    v  = edge(bseg%Edge(i)%e)%v
    xi = -1
    call BSPoint(curve,edge(bseg%edge(i)%e)%geo%t,x,xi)
    vertex(v(1))%x = x
    xi  = 1
    call BSPoint(curve,edge(bseg%edge(i)%e)%geo%t,x,xi)
    vertex(v(2))%x = x
  end do
end subroutine SetupBSVertices

!===============================================================================

end module Boundary_Segment
