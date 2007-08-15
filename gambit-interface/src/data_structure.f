module Data_Structure

  use Constants
  use Basic_IO, only: new_unit
  use Spectral_Element_Vertex
  use Spectral_Element_Edge
  use Spectral_Element
  use Boundary_Segment

  implicit none
  private

  public :: CheckDataStructure
  public :: CheckElementEdge
  public :: CheckElementEdges
  public :: CheckElementOrientations
  public :: CheckBoundaryEdges

contains

!===============================================================================

subroutine CheckDataStructure(vertex,edge,elmt,bseg)
  type(SpectralElementVertex), intent(in) :: vertex(:)
  type(SpectralElementEdge), intent(in)	  :: edge(:)
  type(SpectralElement), intent(in)	  :: elmt(:)
  type(BoundarySegment), intent(in)       :: bseg(:)

  write(*,'(A)',advance='NO')   '  Checking the data structure ...     '

  call CheckElementEdges(edge,elmt)
  call CheckElementOrientations(elmt,vertex)
  call CheckBoundaryEdges(edge,bseg)

  write(*,'(A)') 'done'
  write(*,'(A)') repeat('*',45)

end subroutine CheckDataStructure

!===============================================================================

subroutine CheckElementEdges(edge,elmt)
  type(SpectralElementEdge), intent(in)	:: edge(:)
  type(SpectralElement), intent(in)	:: elmt(:)
  integer				:: i
  do i = 1, size(elmt)
    call CheckElementEdge(edge,elmt(i),i)
  end do
end subroutine CheckElementEdges

!===============================================================================

subroutine CheckElementEdge(edge,elmt,elmtID)
  type(SpectralElementEdge), intent(in)	:: edge(:)
  type(SpectralElement), intent(in)	:: elmt
  integer, intent(in)			:: elmtID
  integer				:: i, j

! check connectivity
  do i = 1, 4
    if (elmt%e(i) == 0) STOP 'Error: Some edges are not defined.'
    j = i + 1
    if (j == 5) j = 1
    if (edge(elmt%e(i))%v(1) /= elmt%v(i) .AND. &
        edge(elmt%e(i))%v(2) /= elmt%v(i)) &
      STOP 'Error: Some edges are not connected.'

    if (edge(elmt%e(i))%v(1) /= elmt%v(j) .AND. &
        edge(elmt%e(i))%v(2) /= elmt%v(j)) &
      STOP 'Error: Some edges are not connected.'
  end do

  ! check element id's
  do i = 1, 4
    if (edge(elmt%e(i))%l(1) /= elmtID .AND. &
        edge(elmt%e(i))%l(2) /= elmtID )     &
      STOP 'Error: Edge-elmt correlation doesnt fit.'
    if (edge(elmt%e(i))%b(1) /= 0      .AND. &
        edge(elmt%e(i))%l(1) /= elmtID .AND. &
        edge(elmt%e(i))%l(2) /= 0 )          &
      STOP 'Error: BoundaryEdge-elmt correlation doesnt fit.'
  end do

end subroutine CheckElementEdge

!===============================================================================

subroutine CheckElementOrientations(elmt,vertex)
  type(SpectralElement), intent(in)	        :: elmt(:)
  type(SpectralElementVertex), intent(in)	:: vertex(:)
  integer					:: i

  do i = 1, size(elmt)
    if (ElementOrientation(vertex(elmt(i)%v)).EQV..FALSE.) &
      STOP 'Error: Element orientation doesnt fit.'
  end do

end subroutine CheckElementOrientations

!===============================================================================

subroutine CheckBoundaryEdges(edge,bseg)
  type(SpectralElementEdge), intent(in)	::  edge(:)
  type(BoundarySegment), intent(in)    	::  bseg(:)
  integer				::  i, j, v(2,2), e(2)
  logical				::  connected = .FALSE.

  ! check boundary edge neighborhood
  do i = 1, size(bseg)
    do j = 1, size(bseg(i)%Edge) - 1
      connected = .FALSE.
      e(1)   = bseg(i)%Edge(j)%e
      e(2)   = bseg(i)%Edge(j+1)%e
      v(1,1) = edge(e(1))%v(1)
      v(1,2) = edge(e(1))%v(2)
      v(2,1) = edge(e(2))%v(1)
      v(2,2) = edge(e(2))%v(2)
      if ( (v(1,1) == v(2,1)) .OR. &
           (v(1,1) == v(2,2)) .OR. &
           (v(1,2) == v(2,1)) .OR. &
           (v(1,2) == v(2,2)) ) connected = .TRUE.
      if (connected.EQV..FALSE.) STOP 'Error: Boundary edges are not connected.'
    end do
  end do

  ! check boundary edge definition
  do i = 1, size(edge)
    if ( edge(i)%b(1) /= 0 ) then
      if( bseg(edge(i)%b(1))%edge(edge(i)%b(2))%e /= i ) &
        STOP 'Error: Boundary edge definition doesnt fit.'
    end if
  end do

end subroutine CheckBoundaryEdges

!===============================================================================

end module Data_Structure
