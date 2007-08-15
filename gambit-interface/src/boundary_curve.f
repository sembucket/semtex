module Boundary_Curve

  use Constants
  implicit none
  private

  type, public :: StrLine
    real(RNP), pointer :: tc(:) => null()
  end type StrLine

  type, public :: BSpline
    integer ::  P = 3  ! polynomial degree
    integer, pointer :: tc(:) => null()
  end type BSpline

  type, public :: Cubic
    integer            :: n = 1 ! 0: isoparametric, 1: curve length based
    real(RNP), pointer :: tc(:) => null()
    real(RNP), pointer :: xcd(:,:) => null()
  end type Cubic

! Änderung 15.12.2006-----
  type, public :: CircArc
    real(RNP), pointer :: tc(:) => null()
  end type CircArc
!-------------------------

  type, public :: BoundaryCurve
    character(len=25)      :: name = ''   ! boundary curve name
    character(len=25)      :: file = ''   ! file which contains control points
    integer                :: geom    = 1       ! boundary geometry
! Änderung 20.12.2006------
    real(RNP)              :: x0 = 0.0          ! Center coordinates x-value
    real(RNP)              :: y0 = 0.0          ! Center coordinates y-value
    real(RNP), allocatable :: r0(:)             ! Arc radius for control points(0:Nmax)
!--------------------------
    real(RNP), pointer	   :: xc(:,:) => null() ! control points (0:Nmax)
    type(StrLine), pointer :: sLine   => null()
    type(BSpline), pointer :: bSpline => null()
    type(Cubic), pointer   :: Cubic   => null()
! Änderung 15.12.2006------
    type(CircArc), pointer :: cArc    => null()
!--------------------------
  end type BoundaryCurve

  type(BoundaryCurve), public, save, allocatable :: ControlBoundaryCurve(:)

  public :: GetBndryCurveControlData
  public :: BSPoint
  public :: SetupBSSpline

! 1. Geometry type
!    1 = Straight
!    2 = Bernstein polynomial approximation
!    3 = B-spline approximation
!    4 = Cubic spline interpolation
!    5 = Arc segments
! 2. BoundaryEdge(i) and BoundaryEdge(i+1) are neighbors

contains

!===============================================================================

subroutine GetBndryCurveControlData(unit)
  integer, intent(in) :: unit
  integer	      :: i, N, geom, ios
  character(len=25)   :: name, file

! Änderung 20.12.2006--------
  real(RNP)           :: x0, y0

  namelist/boundary_curve/ name, file, geom, x0, y0
!  namelist/boundary_curve/ name, file, geom
!----------------------------

  ! count number of user boundary conditions
  N = 0
  rewind(unit)
  do
    name = ''
    read(unit,nml=boundary_curve,IOSTAT=ios)
    if (ios /= 0) exit
    if (LEN_TRIM(name) > 0) N = N + 1
  end do
  allocate( ControlBoundaryCurve(N) )
  rewind(unit)

  do i = 1, N
    name = ControlBoundaryCurve(i)%name
    file = ControlBoundaryCurve(i)%file
    geom = ControlBoundaryCurve(i)%geom
! Änderung 20.12.2006------------------
    x0 = ControlBoundaryCurve(i)%x0
    y0 = ControlBoundaryCurve(i)%y0
!--------------------------------------
    do
      name = ''
      read(unit,nml=boundary_curve)
      if (LEN_TRIM(name) > 0) exit
    end do
    ControlBoundaryCurve(i)%name = name
    ControlBoundaryCurve(i)%file = file
    ControlBoundaryCurve(i)%geom = geom
! Änderung 20.12.2006------------------
    ControlBoundaryCurve(i)%x0 = x0
    ControlBoundaryCurve(i)%y0 = y0
!--------------------------------------
  end do

end subroutine GetBndryCurveControlData

!===============================================================================

subroutine BSPoint(curve, te, x, xi)
  type(BoundaryCurve), intent(in) :: curve
  real(RNP), intent(out)	  :: x(2)
  real(RNP), intent(in)		  :: xi, te(2) ! -1 <= xi <= 1
  real(RNP), pointer		  :: tc(:), xc(:,:), xcd(:,:)
  real(RNP)		:: t	! gloabal t value
  integer, pointer	:: ti(:)
  integer		:: P

  xc => curve%xc
  t  =  te(1) + ( te(2) - te(1) ) * (xi+one) * half
  if (te(1) > te(2)) &
  t  =  te(2) + ( te(1) - te(2) ) * (xi+one) * half

  select case(curve%geom)
  case(1)
    tc => curve%sLine%tc
    call StraightLine(tc, xc(:,1), xc(:,2), t, x)
  case(2)
    call BernsteinSpline(xc(:,1), xc(:,2), t, x)
  case(3)
    P = curve%BSpline%P
    ti => curve%BSpline%tc
    call B_Spline(ti, xc(:,1), xc(:,2), P, t, x)
  case(4)
    tc   => curve%Cubic%tc
    xcd  => curve%Cubic%xcd
    call CubicSplineVector(tc, xc, xcd, t, x)
! Änderung 20.12.2006--------

  case(5)
! Eigentlich nicht nötig, da diese subroutine nur von 
! ExportBoundarySplineToGnuplot aufgerufen wird
    tc => curve%cArc%tc
    call StraightLine(tc, xc(:,1), xc(:,2), t, x)

!----------------------------
  case default
    STOP 'Error: Boundary geometry is not defined.'
  end select

end subroutine BSPoint

!===============================================================================

subroutine SetupBSSpline(curve)
  type(BoundaryCurve), intent(inout) :: curve
  integer, pointer		     :: ti(:)
  real(RNP), pointer		     :: tc(:), xc(:,:), xcd(:,:)
  real(RNP)			     :: x0, y0
  integer			     :: K, P, n

  K = size(curve%xc(:,1)) - 1
  select case(curve%Geom)
  case(1)
    write(*,'(A14,A)',advance='YES') '', 'use straight interpolation'
    allocate( curve%sLine )
    allocate( curve%sLine%tc(0:K) )
    tc => curve%sLine%tc
    xc => curve%xc
    call InitStraightLineNodes(tc, xc(:,1), xc(:,2))
  case(2)
    write(*,'(A14,A)',advance='YES') '', 'use bernstein approximation'
  case(3)
    write(*,'(A14,A)',advance='YES') '', 'use b-spline approximation'
    allocate( curve%bSpline )
    P = curve%bSpline%P
    allocate ( curve%BSpline%tc(0:K+P+1) )
    ti => curve%BSpline%tc
    call InitB_SplineNodes(ti, K, P)
  case(4)
    write(*,'(A14,A)',advance='YES') '', 'use cubic spline interpolation'
    allocate( curve%Cubic )
    allocate( curve%Cubic%tc(0:K) )
    allocate( curve%Cubic%xcd(0:K,2) )
    tc  => curve%Cubic%tc
    xc  => curve%xc
    xcd => curve%Cubic%xcd
    n   = curve%Cubic%n
    call InitCubicSplineNodes(tc, xc(:,1), xc(:,2), n)
    call InitCubicSpline(tc, xc(:,1), xcd(:,1))
    call InitCubicSpline(tc, xc(:,2), xcd(:,2))

! Änderung 20.12.2006 ------------
  case(5)
    write(*,'(A14,A)',advance='YES') '', 'use circular arc segments'
    allocate( curve%cArc )
    allocate( curve%cArc%tc(0:K) )
    tc => curve%cArc%tc
    xc => curve%xc
    call InitStraightLineNodes(tc, xc(:,1), xc(:,2))
    allocate( curve%r0(0:K) )
    x0 = curve%x0
    y0 = curve%y0
    call CircularArc(xc(:,1), xc(:,2), curve%r0(0:K), x0, y0)
!---------------------------------

  case default
    STOP 'Error: Boundary geometry is not defined.'
  end select

end subroutine SetupBSSpline

!===============================================================================

subroutine InitCubicSplineNodes(tc, xc, yc, n)
  real(RNP), intent(inout)	:: tc(0:)
  real(RNP), intent(in)		:: xc(0:), yc(0:)
  real(RNP)			:: Kreal
  integer, intent(in)		:: n
  integer			:: i, K

  K = size(tc) - 1
  Kreal = K
  select case(n)
  case(0) 			! equidistant
    forall(i=0:K) tc(i) = i / Kreal
  case(1:)			! arc length
    tc(0) = 0
    do i = 1, K
      tc(i) = sqrt((xc(i)-xc(i-1))**2 + (yc(i)-yc(i-1))**2) + tc(i-1)
    end do
    tc = tc / tc(K)
  end select

end subroutine InitCubicSplineNodes

!===============================================================================

subroutine InitCubicSpline(tc, qc, qd)
  real(RNP), intent(in)		:: tc(0:), qc(0:)
  real(RNP), intent(inout)	:: qd(0:)
  integer			:: i, K
  real(RNP)			:: p, qn, sig, un, qp1, qpn
  real(RNP), allocatable	:: u(:)

  K = size(qc) - 1
  qp1 = 1.e30_RNP
  qpn = 1.e30_RNP
  allocate( u(0:K) )
  if (qp1 > 0.99e30_RNP) then
    qd(0) = zero
    u(0)  = zero
  else
    qd(0) = -half
    u(0)  = (3._RNP/(tc(1)-tc(0))) * ((qc(1)-qc(0)) / (tc(1)-tc(0))-qp1)
  endif
  do i = 1, K - 1
    sig   = (tc(i)-tc(i-1)) / (tc(i+1)-tc(i-1))
    p     = sig*qd(i-1) + 2._RNP
    qd(i) = (sig-one) / p
    u(i)  = (6._RNP*((qc(i+1)-qc(i)) / (tc(i+1)-tc(i))-(qc(i)-qc(i-1)) &
            / (tc(i)-tc(i-1)))       / (tc(i+1)-tc(i-1))-sig*u(i-1)) / p
  enddo
  if (qpn > 0.99e30_RNP) then
    qn = zero
    un = zero
  else
    qn = half
    un = (3._RNP/(tc(K)-tc(K-1))) * (qpn-(qc(K)-qc(K-1)) / (tc(K)-tc(K-1)))
  endif
  qd(K) = (un-qn*u(K-1)) / (qn*qd(K-1)+one)

  do i = K - 1, 0, -1
    qd(i) = qd(i) * qd(i+1) + u(i)
  enddo
  deallocate(u)

end subroutine InitCubicSpline

!===============================================================================

subroutine CubicSplineVector(ts, xc, xcd, t, x)
  real(RNP), intent(in)		:: ts(0:), xc(0:,:), xcd(0:,:), t
  real(RNP), intent(inout)	:: x(:)

  x(1) = CubicSpline(ts, xc(0:,1), xcd(0:,1), t)
  x(2) = CubicSpline(ts, xc(0:,2), xcd(0:,2), t)

end subroutine CubicSplineVector

!===============================================================================

function CubicSpline(tc, qc, qcd, t) result(q)
  real(RNP), intent(in)		:: tc(0:), qc(0:), qcd(0:), t
  integer			:: K, ka, khi, klo
  real(RNP)			:: a, b, h, q

  K = size(tc) - 1
  klo = 0
  khi = K
  do while(khi - klo > 1 )
    ka = (khi + klo) / 2
    if(tc(ka) > t) then
      khi = ka
    else
      klo = ka
    end if
  end do
  h = tc(khi) - tc(klo)
  if (h == 0) STOP 'Error: Bad tc input in splint.'
  a = (tc(khi) - t) / h
  b = (t - tc(klo)) / h
  q = a*qc(klo)+b*qc(khi) + ((a**3-a)*qcd(klo)+(b**3-b)*qcd(khi))*(h**2)/6._RNP

end function CubicSpline

!===============================================================================

subroutine InitB_SplineNodes(ti, K, P)
  integer, intent(in)		:: K, P
  integer, intent(inout)	:: ti(0:)
  integer			:: i

  do i = 0, K + P + 1
    if (i <= P)                  ti(i) = ( 0         )
    if (P + 1 <= i .AND. i <= K) ti(i) = ( i - P     )
    if (i > K)                   ti(i) = ( K - P + 1 )
  end do

end subroutine InitB_SplineNodes

!===============================================================================

subroutine B_Spline(ti, xc, yc, P, t, x)
  integer, intent(in)			:: P, ti(0:)
  real(RNP), intent(in)			:: xc(0:), yc(0:), t
  real(RNP), intent(inout)		:: x(2)
  integer				:: K, i, j, d
  real(RNP), allocatable		:: b(:,:)
  real(RNP)				:: ta, KP

  K = size(ti) - P - 2
  allocate( b(0:K+P,0:P) )
  KP = K - P + 1
  ta = t * KP
  b = 0
  do i = 0, K + P
    b(i,0) = 0
    if ( ti(i) <= ta .AND. ta < ti(i + 1) ) b(i,0) = 1
  end do

  do j = 1, P
    do i = 0, K + P - j
      b(i,j) = 0
      d = ti(i+j) - ti(i)
      if (d /= 0) b(i,j) = b(i,j) + (ta-ti(i)) * b(i,j-1) / d
      d = ti(i+j+1) - ti(i+1)
      if (d /= 0) b(i,j) = b(i,j) + (ti(i+j+1)-ta) * b(i+1,j-1) / d
    end do
  end do
  if (t == 1) b(K,P) = 1 ! the one and only error in this code, dont know why
  x = 0
  do i = 0, K
    x(1) = x(1) + xc(i) * b(i,P)
    x(2) = x(2) + yc(i) * b(i,P)
  end do
  deallocate(b)

end subroutine B_Spline

!===============================================================================

subroutine InitStraightLineNodes(tc, xc, yc)
  real(RNP), intent(inout)	:: tc(0:)
  real(RNP), intent(in)		:: xc(0:), yc(0:)
  integer			:: i, K

  K = size(tc)-1
  tc(0) = 0
  do i = 1, K
    tc(i) = sqrt((xc(i)-xc(i-1))**2 + (yc(i)-yc(i-1))**2) + tc(i-1)
  end do
  tc = tc / tc(K)

end subroutine InitStraightLineNodes

!===============================================================================

subroutine StraightLine(tc, xc, yc, t, x)
  real(RNP), intent(in)			:: tc(0:), xc(0:), yc(0:), t
  real(RNP), intent(inout)		:: x(2)
  integer				:: i

  x = 0
  do i = 1, size(tc)-1
    if ( tc(i-1) <= t .AND. t <= tc(i) ) then
      x(1) = xc(i-1) + (xc(i)-xc(i-1)) * (t-tc(i-1)) / (tc(i)-tc(i-1))
      x(2) = yc(i-1) + (yc(i)-yc(i-1)) * (t-tc(i-1)) / (tc(i)-tc(i-1))
    end if
  end do

end subroutine StraightLine

!===============================================================================

subroutine BernsteinSpline(xc, yc, t, x)
  real(RNP), intent(in)		:: xc(0:), yc(0:), t
  real(RNP), intent(out)	:: x(2)
  integer			:: K, i, j
  real(RNP), allocatable	:: xb(:,:), yb(:,:)

  K = size(yc) - 1
  allocate ( xb(0:K,0:K), yb(0:K,0:K) )
  xb(:,0) = xc
  yb(:,0) = yc
  do j = 1, K
    do i = 0, K - j
      xb(i,j) = (1 - t) * xb(i,j-1) + t * xb(i+1,j-1)
      yb(i,j) = (1 - t) * yb(i,j-1) + t * yb(i+1,j-1)
    end do
  end do
  x(1) = xb(0,K)
  x(2) = yb(0,K)
  deallocate(xb,yb)

end subroutine BernsteinSpline

!===============================================================================

! Änderung 20.12.2006

subroutine CircularArc(xc, yc, r0, x0, y0)
  real(RNP), intent(in) 	:: xc(0:), yc(0:), x0, y0
  real(RNP), intent(out)	:: r0(0:)
  integer			:: i, K

  K = size(yc)-1
  do i = 0, K
    r0(i) = sqrt( (xc(i)-x0)**2 + (yc(i)-y0)**2 )
  end do

end subroutine CircularArc

!===============================================================================

end module Boundary_Curve
