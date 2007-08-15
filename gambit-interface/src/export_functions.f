module Export_Functions

  use Constants
  use Basic_IO, only: new_unit
  use Spectral_Element_Vertex
  use Spectral_Element_Edge
  use Spectral_Element
  use Boundary_Segment
  use Boundary_Curve

! Änderung 07.09.2006
  use Boundary_Condition

  implicit none
  private

  type, public :: ExportControlData
    character(len=80) :: casename = 'default'
    integer           :: Nexport  = 0 ! number of time steps to export variables
    logical	      :: on_start = .false. ! write before time integration ?
    logical	      :: on_end   = .false. ! write after time integration ?
  end type ExportControlData

  type(ExportControlData), public, save :: ExportData

  public :: GetExportControlData
  public :: ExportGridToTecplot
  public :: ExportBoundarySplineToGnuplot
  public :: ExportBoundaryVerticesToGnuplot
  public :: ExportBoundaryControlPointsToGnuplot

  public :: ExportGridToSemtex !Änderung 04.09.2006

  character(len=20), parameter :: fmt1 = '(1ES13.5)'
  character(len=20), parameter :: fmt2 = '(2ES13.5)'
  character(len=20), save :: fmtT = ''

contains

!===============================================================================

subroutine GetExportControlData(unit)
  integer, intent(in) :: unit
  character(len=80)   :: casename
  integer	      :: ios
  namelist/exportcase/ casename
  casename = ExportData%casename
  rewind(unit)
  read(unit,nml=exportcase,IOSTAT=ios)
  ExportData%casename = casename
end subroutine GetExportControlData

!===============================================================================

subroutine ExportGridToTecplot(elmt,vertex)
  type(SpectralElement), intent(in)	  :: elmt(:)
  type(SpectralElementVertex), intent(in) :: vertex(:)
  character(len=25)			  :: file
  integer				  :: unit, i, j, ios

  file = exportdata%casename
  write(file(len_trim(file)+1:),'(A)') '_grid.plt'

  write(*,'(A)')  '  Writing the gridfile:'
  write(*,'(2A)') '  ... ', file
  write(*,'(A)')  repeat('*',45)

  unit = new_unit()
  open(unit, IOSTAT=ios, file=file, action='WRITE', status='REPLACE')
  if (ios /= 0) STOP 'An error occurred while opening the Tecplotfile.'

  write(unit,'(A)') 'TITLE = "2D-SEM-Results"'
  write(unit,'(A)') 'VARIABLES  = X, Y'
  write(unit,'(A,I0,A,I0,A)') 'ZONE T="elmt-sorted",N=',4*size(elmt),',E=' &
			,size(elmt),',ET=QUADRILATERAL,F=FEPOINT'

  do i = 1, size(elmt)
    do j = 1, 4
      write(unit,fmt2) vertex(elmt(i)%v(j))%x
    end do
  end do
  do i = 1, size(elmt)
    write(unit,'(4(I0,1x))') 4*(i-1)+1,4*(i-1)+2,4*(i-1)+3,4*(i-1)+4
  end do

  write(unit,'(/,A)') 'TITLE = "2D-SEM-Grid-orig"'
  write(unit,'(A)') 'VARIABLES  = X, Y'
  write(unit,'(A,I0,A,I0,A)') 'ZONE T="global-sorted", N=',size(vertex),',E=' &
			,size(elmt),',ET=QUADRILATERAL, F=FEPOINT'

  do i = 1, size(vertex)
    write(unit,fmt2) vertex(i)%x
  end do
  do i = 1, size(elmt)
    write(unit,'(4(I0,1x))') elmt(i)%v
  end do
  close(unit)

end subroutine ExportGridToTecplot

!===============================================================================

subroutine ExportBoundarySplineToGnuplot(edge,bseg,curve)
  type(SpectralElementEdge), intent(in)	:: edge(:)
  type(BoundarySegment), intent(in)	:: bseg
  type(BoundaryCurve), intent(in)	:: curve
  character(len=80)			:: file
  integer				:: i, j, unit, ios, K
  real(RNP)				:: t, Kreal, x(2)

  K     = 50
  Kreal = K
  file = exportdata%casename
  write(file(len_trim(file)+1:),'(3A)') &
  '_bndry-spline_', bseg%name(1:len_trim(bseg%name)), '.dat'

  write(*,'(2A)')   '  Writing the boundary spline:'
  write(*,'(2A)') '  ... ', file
  write(*,'(A)')  repeat('*',45)
  unit = new_unit()
  open(unit, IOSTAT=ios, file=file, action='WRITE', status='REPLACE')
  if (ios /= 0) STOP 'An error occurred while opening the file.'
  do i = 1, size(bseg%Edge)
    do j = -K, K
      t = j / Kreal
      call BSPoint(curve,edge(bseg%edge(i)%e)%geo%t,x,t)
      write(unit,fmt2) x
    end do
  end do
  close(unit)
  write(*,'(A,/,A)') '  done', repeat('*',45)
end subroutine ExportBoundarySplineToGnuplot

!===============================================================================

subroutine ExportBoundaryVerticesToGnuplot(vertex, edge, bseg)
  type(SpectralElementVertex), intent(in) :: vertex(:)
  type(SpectralElementEdge), intent(in)	  :: edge(:)
  type(BoundarySegment), intent(in)	  :: bseg
  character(len=80)			  :: file
  integer				  :: i, unit, ios

  file = exportdata%casename
  write(file(len_trim(file)+1:),'(3A)') &
  '_bndry-vertices_', bseg%name(1:len_trim(bseg%name)), '.dat'

  write(*,'(2A)')   '  Writing the boundary vertices:'
  write(*,'(2A)') '  ... ', file
  write(*,'(A)')  repeat('*',45)
  unit = new_unit()
  open(unit, IOSTAT=ios, file=file, action='WRITE', status='REPLACE')
  if (ios /= 0) STOP 'An error occurred while opening the file.'
  write(unit,fmt2) vertex(edge(bseg%edge(1)%e)%v(1))%x
  do i = 1, size(bseg%Edge)
    write(unit,fmt2) vertex(edge(bseg%Edge(i)%e)%v(2))%x
  end do
  close(unit)
  write(*,'(A,/,A)') '  done', repeat('*',45)

end subroutine ExportBoundaryVerticesToGnuplot

!===============================================================================

subroutine ExportBoundaryControlPointsToGnuplot(bseg,curve)
  type(BoundarySegment), intent(in) :: bseg
  type(BoundaryCurve), intent(in)   :: curve
  character(len=80)		    :: file = ''
  integer			    :: i, unit, ios

  file = exportdata%casename
  write(file(len_trim(file)+1:),'(3A)') &
  '_bndry-ctrl-points_', bseg%name(1:len_trim(bseg%name)), '.dat'

  write(*,'(2A)') '  Writing the boundary control points:'
  write(*,'(2A)') '  ... ', file
  write(*,'(A)')  repeat('*',45)
  unit = new_unit()
  open(unit, IOSTAT=ios, file=file, action='WRITE', status='REPLACE')
  if (ios /= 0) STOP 'An error occurred while opening the file.'
  do i = 0, size(curve%xc(:,1)) - 1
    write(unit,fmt2) curve%xc(i,:)
  end do
  close(unit)
  write(*,'(A,/,A)') '  done', repeat('*',45)

end subroutine ExportBoundaryControlPointsToGnuplot

!===============================================================================

subroutine ExportGridToSemtex(elmt,vertex,edge,bseg,curve,bcond) !Änderung 04.09.2006

  type(SpectralElement), intent(in)	  :: elmt(:)
  type(SpectralElementVertex), intent(in) :: vertex(:)
  type(SpectralElementEdge), intent(in)	  :: edge(:)

  type(BoundarySegment), pointer		:: bseg(:)
  type(BoundaryCurve), pointer			:: curve(:)
  type(BoundaryCondition), pointer		:: bcond(:)

  character(len=25)			  :: file, curvetype
  integer				  :: unit, i, ii, j, jj, k, ios, element, p1, p2, n_edges, i2, kk, mm
  character, allocatable                  :: temp(:)
  integer,   allocatable                  :: temp_no(:)
  logical                                 :: flag, flag2
  real(RNP)                               :: radius, xc, yc, distance

  file = exportdata%casename
  write(file(len_trim(file)+1:),'(A)')


  write(*,'(A)')  '  Writing the gridfile:'
  write(*,'(2A)') '  ... ', file
  write(*,'(A)')  repeat('*',45)

  unit = new_unit()
  open(unit, IOSTAT=ios, file=file, action='WRITE', status='REPLACE')
  if (ios /= 0) STOP 'An error occurred while opening the SemTex file.'

write(unit,'(A)')  repeat('#',80)
write(unit,'(2A)') '# session file: ', trim(file)
write(unit,'(A)')  repeat('#',80)
write(unit,'(A)') ''

write(unit,'(A)') '<USER>'
write(unit,'(A)') '# Initial conditions go here'
write(unit,'(A)') '</USER>'
write(unit,'(A)') ''

write(unit,'(A)') '<FIELDS>'
write(unit,'(A)') '# Specification of field variables'
write(unit,'(A)') '</FIELDS>'
write(unit,'(A)') ''

write(unit,'(A)') '<TOKENS>'
write(unit,'(A)') '# Solver parameters go here'
write(unit,'(A)') '</TOKENS>'
write(unit,'(A)') ''

! Determine the number of different boundary condition keys ( =groups)
allocate(temp(1:size(bcond)))

! All boundaries marked with "p" have to be detected

! if (bcond(1)%key == 'p') then k = 0

!temp(1) = bcond(1)%key
!k = 1
!flag = .true.

!do i = 1, size(bcond)
!  do j = 1, k
!    if (bcond(i)%key == temp(j)) then
! "flag" is true when the given key has previously been found
!      flag = .true.
!      exit
!    else
!      flag = .false.
!    end if
!  end do
!  if (flag .eqv. .false.) then
!    k = k + 1
!    temp(k) = bcond(i)%key
!  end if
!end do

!if (bcond(1)%key =='p')  then
!  k = 0
!else
!  temp(1) = bcond(1)%key
!  k = 1
!end

k = 0

do i = 1, size(bcond)
! Periodic boundaries have no group, so they are not counted
  if (bcond(i)%key == 'p') then
    cycle
  else
    if (k == 0) then
      temp(1) = bcond(i)%key
      k = 1
      flag = .true.
    else
      do j = 1, k
        if (bcond(i)%key == temp(j)) then
          ! "flag" is true when key has previously been found
          flag = .true.
          exit
        else
          flag = .false.
        end if
      end do
    end if
    if (flag .eqv. .false.) then
      k = k + 1
      temp(k) =  bcond(i)%key
    end if
  end if
end do

if (k > 0) then

write(unit,'(A,I0,A)') '<GROUPS NUMBER=', k,'>'

do i = 1, k
      write(unit,'(I6,A6,A5,A)')  i, temp(i), '', 'group_name'
end do

write(unit,'(A)') '</GROUPS>'
write(unit,'(A)') ''

write(unit,'(A,I0,A)') '<BCS NUMBER=', k, '>'

do i = 1, k
      write(unit,'(I6,A6,A5,A)')  i, temp(i), '', '3'
      write(unit,'(A)') '                 <D>     u = 0      </D>'
      write(unit,'(A)') '                 <D>     v = 0      </D>'
      write(unit,'(A)') '                 <H>     p = 0      </H>'
end do

deallocate(temp)

write(unit,'(A)') '</BCS>'
write(unit,'(A)') ''

end if

write(unit,'(A,I0,A)') '<NODES NUMBER=',size(vertex),'>'

do i = 1, size(vertex)
    write(unit,'(I6,3ES13.5)') i , vertex(i)%x , 0.
end do

write(unit,'(A)') '</NODES>'
write(unit,'(A)') ''

write(unit,'(A,I0,A)') '<ELEMENTS NUMBER=',size(elmt),'>'

do i = 1, size(elmt)
  write(unit,'(I6,A,4I6,A)') i,'   <Q>', elmt(i)%v, '   </Q>' 
end do

write(unit,'(A)') '</ELEMENTS>'
write(unit,'(A)') ''

! Determine the number of edges that lie on a boundary: edge(i)%l(2) must be zero (second element)
!k = 0
!do i = 1, size(edge)
!  if (edge(i)%l(2) == 0 ) then
!  k = k + 1
!  end if
!end do

allocate(temp_no(1:size(bseg)))
temp_no = 0
kk = 0

! Detrmine the number of surfaces
k = 0
do i = 1, size(bcond)
  ! Check if boundary is periodic or not
  if (bcond(i)%key == 'p') then
    ! Boundary is periodic, check if it is recorded in "temp_no"
    flag = .false.
    do mm = 1, size(temp_no)
      if (i == temp_no(mm) ) then
        ! Boundary has been recorded before
        flag = .true.
        exit
      end if
    end do
    ! If boundary has not been recorded before ...
    if (flag .eqv. .false.) then
      ! Look for the complementary boundary of "i"
      do ii = 1, size(bcond)
        if (bcond(i)%link == bcond(ii)%name) then
          ! Corresponding boundary "ii" found: record it in temp_no
          kk = kk + 1
          temp_no(kk) = ii
          exit
        end if
      end do
      ! Add edges to total
      k = k + size(bseg(i)%edge)
    end if
  else
    ! Boundary is not periodic, just add all the edges to total
    k = k + size(bseg(i)%edge)
  end if
end do

! diagnostic

!write(unit,'(7A8)') '  bseg','  edge','   key','  name','  link','vert_1','vert_2'
!do i = 1, size(bseg)
!  do j = 1, size(bseg(i)%edge)
!    write(unit,'(I7,A,I7,A,A7,A,A7,A,A7,A,I7,A,I7,A,I7,A,I7)') i,' ', bseg(i)%edge(j)%e,' '&
!    &, trim(bcond(i)%key),' ', trim(bcond(i)%name),' ', trim(bcond(i)%link) &
!    &,' ', edge( bseg(i)%edge(j)%e )%v(1),' ', edge( bseg(i)%edge(j)%e )%v(2) &
!    &,' ', edge( bseg(i)%edge(j)%e )%b(1),' ', edge( bseg(i)%edge(j)%e )%b(2)
!  end do
!end do
!write(unit,'(A)') ''

!write(unit,*) temp_no
!write(unit,*) ''

write(unit,'(A,I0,A)') '<SURFACES NUMBER=', k,'>'

k = 1
do i = 1, size(edge)
  ! Check if edge(i) is a boundary edge: edge(i)%l(2) must be 0
  if (edge(i)%l(2) == 0 ) then
    ! Determine edge number (local numbering 1..4)
    j = 1
    do
      if (elmt( edge(i)%l(1) )%v(j) == edge(i)%v(1)) then 
        exit
      else
        j = j + 1
      end if
    end do
    ! If edge is part of a periodic boundary, the corresponding edge must be found
    if (bcond(edge(i)%b(1))%key == 'p') then

! Boundary is periodic, check if it is recorded in "temp_no"
    flag2 = .false.
    do mm = 1, size(temp_no)
      if (edge(i)%b(1) == temp_no(mm) ) then
        ! Boundary has been recorded before
        flag2 = .true.
        exit
      end if
    end do

if (flag2 .eqv. .false.) then
! New boundary, add it to file

      ! Looking for corresponding boundary
      p1 = edge(i)%b(1) !<- Number of first boundary condition/segment
      flag = .false.
      do ii = 1, size(bcond)
        if (bcond(edge(i)%b(1))%link == bcond(ii)%name) then
          if (bcond(ii)%key == 'p') flag = .true.
          p2 = ii !<- Number of second boundary condition/segment
          exit
        end if
      end do
      if (flag .eqv. .false.) then
         write(*,'(A,25A,A)') 'Error: Link not set correctly for periodic boundary "',&
         & trim(bcond(p1)%name), '"'
         write(*,'(A,25A,A)') '... or linked boundary "', trim(bcond(ii)%name),'" is not set to periodic.'
         stop
      else
      ! "p2" is the number of the corresponding periodic boundary
      ! It is assumed that both boundaried have the same number of edges,
      ! but the corresponding edges are numbered in opposite order
      n_edges =  size(bseg(p1)%edge)
      ! edge(i)%b(2)
      ! Corresponding "local" edge number on p2 is (-ii + n_edges + 1)
      i2 = -edge(i)%b(2) + n_edges + 1
      ! "bseg(p2)%edge(i2)%e" is global edge number on segment p2 for "local" edge number i2
      ! Determine edge number (local numbering 1..4 in element), "jj" is edge in local number
      jj = 1
      do
        if (elmt( edge( bseg(p2)%edge(i2)%e )%l(1) )%v(jj) == edge( bseg(p2)%edge(i2)%e )%v(1)) then 
          exit
        else
          jj = jj + 1
        end if
      end do
      ! Write entry for periodic edge
      write(unit,'(3I6,A,2I6,A)') k, edge(i)%l(1), j, '   <P>   ',&
      & edge( bseg(p2)%edge(i2)%e )%l(1), jj, '   </P>'
      k  = k + 1
      end if

end if

    else
      ! Entry for nonperiodic edge
      ! Assign boundary condition key to edge(i)
      write(unit,'(3I6,3A)') k, edge(i)%l(1), j, '   <B>        ', bcond( edge(i)%b(1) )%key, '      </B>'
      k = k + 1
    end if
  end if
end do

deallocate(temp_no)

write(unit,'(A)') '</SURFACES>'

! Determine the number of edges that lie on curved boundaries
k = 0
do i = 1, size(curve)
  if ( curve(i)%geom /= 1) then
  k = k + size(bseg(i)%edge)
  end if
end do

write(unit,'(A)') ''

if (k > 0) then

write(unit,'(A,I0,A)') '<CURVES NUMBER=', k,'>'

k = 0

! Go through all boundary segments == go through all specified curves
! Index of boundary segment == index of boundary curve ?
!        write(unit,'(3A6)') '#    k', '  elmt', '  edge' 
do i = 1, size(bseg)
! Determine what kind of curve (spline, arc, straight, etc.) the given boundary segment has
select case (curve(i)%geom)
  case(1)
  ! Straight.. do nothing
  case(4)
  ! Cubic splines
  curvetype = 'SPLINE'

  ! Go through all edges on given boundary segment
  do ii = 1, size(bseg(i)%edge)
  k = k + 1

  ! Determine edge number (local numbering 1..4)
    j = 1
    do
      if (elmt( edge(bseg(i)%edge(ii)%e)%l(1) )%v(j) == edge(bseg(i)%edge(ii)%e)%v(1)) then 
        exit
      else
        j = j + 1
      end if
    end do
! If no external file has been specified for the spline, it has to be 
! generated from element nodes. This is done afterwards in 'converter.f'
! by calling 'ExportBoundarySplineToGnuplot'.
    if (curve(i)%file == '') then 
        write(unit,'(3I6,10A)') k, edge(bseg(i)%edge(ii)%e)%l(1), j ,&
        '      <',trim(curvetype),'>   ',trim(exportdata%casename),'_bndry-spline_', &
        trim(bseg(i)%name),'.dat','   </',trim(curvetype),'>'
    else
        write(unit,'(3I6,7A)') k, edge(bseg(i)%edge(ii)%e)%l(1), j ,&
        '      <',trim(curvetype),'>   ',trim(curve(i)%file),'   </',trim(curvetype),'>'
    end if
  end do

   case(5)
   ! Circular Arc
   curvetype = 'ARC'

  ! Go through all edges on given boundary segment
  do ii = 1, size(bseg(i)%edge)
  k = k + 1

  ! Determine edge number (local numbering 1..4)
    j = 1
    do
      if (elmt( edge(bseg(i)%edge(ii)%e)%l(1) )%v(j) == edge(bseg(i)%edge(ii)%e)%v(1)) then 
        exit
      else
        j = j + 1
      end if
    end do
! Radius for edge is average of the radii for each of the two vertices belonging to that edge
! This assumes that both radii have approximately the same value
        radius = 0.5_RNP * (curve(i)%r0(ii-1) + curve(i)%r0(ii))
! Determine centroid of element attached to edge
        element = edge(bseg(i)%edge(ii)%e)%l(1)
        xc = 0.0_RNP
        yc = 0.0_RNP
        do jj = 1, 4
          xc = xc + vertex(elmt(element)%v(jj))%x(1)
          yc = yc + vertex(elmt(element)%v(jj))%x(2)
        end do
        xc = xc / 4.0_RNP
        yc = yc / 4.0_RNP
! Determine distance of element centroid from center of circle
        distance = sqrt((xc-curve(i)%x0)**2 + (yc-curve(i)%y0)**2)
! Set sign of radius, positive if convex, negative if concave
        if (distance > radius) radius = -radius
        write(unit,'(3I6,3A,ES13.5,3A)') k, edge(bseg(i)%edge(ii)%e)%l(1), j,&
        '      <',trim(curvetype),'>   ', radius,'   </',trim(curvetype),'>'

!    Ignore this
!        edge(bseg(i)%edge(ii)%e)%v(1),'   </',trim(curvetype),'>',&
!        curve(i)%r0(edge(bseg(i)%edge(ii)%e)%v(1),'   </',trim(curvetype),'>',&
!        vertex(edge(bseg(i)%edge(ii)%e)%v(1))%x,'   </',trim(curvetype),'>'
  end do

  case default
  ! Unknown kind of curve, treat boundary curve as straight line
  end select
end do

write(unit,'(A)') '</CURVES>'

end if
  close(unit)

end subroutine ExportGridToSemtex

!===============================================================================


end module Export_Functions
