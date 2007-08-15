program Converter

  use Constants
  use Spectral_Element_Vertex
  use Spectral_Element_Edge
  use Spectral_Element
  use Boundary_Segment
  use Boundary_Curve
  use Boundary_Condition
  use Data_Structure
  use Export_Functions
  use Control_Data
  use Import_Grid

  implicit none

  type(SpectralElementVertex), pointer	:: vertex(:)
  type(SpectralElementEdge), pointer	:: edge(:)
  type(SpectralElement), pointer	:: element(:)
  type(BoundarySegment), pointer	:: bseg(:)
  type(BoundaryCurve), pointer	        :: curve(:)
  type(BoundaryCondition), pointer      :: bcond(:)

! Ändegung  08.09.2006
! character(len=25), parameter :: casefile = 'testcase.inp'
  character(len=80)                     :: casefile
  integer                               :: i

!===============================================================================
! Import
!===============================================================================

! Ändegung  08.09.2006
  write(*,'(A)') 'specify input file (full name, including extension)'
  read(*,'(A)') casefile

  call GetControlData(casefile) ! controldata from inputfile
  call ImportGrid(vertex,edge,element,bseg) ! mesh import
  do i = 1, size(edge) ! initialisation
    allocate(edge(i)%geo)
  end do
  call ConstructBoundarySegments(vertex,edge,bseg,curve,bcond) ! boundary cnstr.
  call CheckDataStructure(vertex,edge,element,bseg) ! final check

!===============================================================================
! Export
!===============================================================================

!  do i = 1, size(bseg)
!    call ExportBoundaryVerticesToGnuplot(vertex, edge, bseg(i))
!    call ExportBoundarySplineToGnuplot(edge, bseg(i),curve(i))
!    call ExportBoundaryControlPointsToGnuplot(bseg(i),curve(i))
!  end do
!  call ExportGridToTecplot(element,vertex)

! Änderung 12.09.2006
  do i = 1, size(bseg)
  if ((curve(i)%geom == 4) .and. ( curve(i)%file == '')) then 
    call ExportBoundarySplineToGnuplot(edge, bseg(i),curve(i))
  end if
  end do

  call ExportGridToSemtex(element,vertex,edge,bseg,curve,bcond)

!===============================================================================

end program Converter
