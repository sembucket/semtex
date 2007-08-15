module Spectral_Element
  use Constants
  use Spectral_Element_Edge
  use Spectral_Element_Vertex
  implicit none
  private

  type, public :: SpectralElement
     integer ::  v(4) = 0  ! vertex IDs
     integer ::  e(4) = 0  ! edge IDs
  end type SpectralElement

  public :: FindUndefinedEdge
  public :: ElementOrientation
  public :: FlipElementOrientation

contains

!===============================================================================

integer function FindUndefinedEdge(element) result(N)
  type(SpectralElement), intent(in) ::  element
  integer ::  i
  N = 0
  do i = 1, 4
    if (element%e(i) == 0) then
       N = i
       exit
    end if
  end do
end function FindUndefinedEdge

!===============================================================================

logical function ElementOrientation(vertex) result(t1)
  type(SpectralElementVertex), intent(in) ::  vertex(4)
  real(RNP) ::  a(2), b(2)
  logical ::  t2

! .TRUE.  :: anticlockwise orientation
! .FALSE. :: clockwise orientation

  t1 = .TRUE.
  t2 = .TRUE.

  a = vertex(2)%x - vertex(1)%x
  b = vertex(3)%x - vertex(1)%x
  if (a(1)*b(2) - a(2)*b(1) < 0) t1 = .FALSE.

  a = vertex(3)%x - vertex(1)%x
  b = vertex(4)%x - vertex(1)%x
  if (a(1)*b(2) - a(2)*b(1) < 0) t2 = .FALSE.

  if (t1.NEQV.t2) STOP 'Error: Different face orientations.'

end function ElementOrientation

!===============================================================================

subroutine FlipElementOrientation(element)
  type(SpectralElement), intent(inout) ::  element
  element%v = element%v((/ 2,1,4,3 /))
  element%e = element%e((/ 1,4,3,2 /))
end subroutine FlipElementOrientation

!===============================================================================

end module Spectral_Element
