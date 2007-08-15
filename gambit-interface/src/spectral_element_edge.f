module Spectral_Element_Edge

  use Constants
  use Spectral_Element_Vertex
  use Boundary_Curve
  implicit none
  private

  type, public :: SpectralElementEdgeData
     logical   :: curved = .false. ! whether edge is curved or straight
     real(RNP) :: t(2)   = 0       ! curve parameter
  end type SpectralElementEdgeData

  type(SpectralElementEdgeData), public, save :: defaultEdgeData

  type, public :: SpectralElementEdge
     integer :: b(2) = 0  ! boundary ID (0 if none)
     integer :: v(2) = 0  ! vertex IDs
     integer :: l(2) = 0  ! element IDs
     type(SpectralElementEdgeData), pointer :: geo => null()
  end type SpectralElementEdge

  public :: SharedNode

contains

!===============================================================================

integer function SharedNode(edge1, edge2) result(v)
  type(SpectralElementEdge), intent(in) ::  edge1, edge2

  v = 0
  if (edge1%v(1) == edge2%v(1)) v = edge1%v(1)
  if (edge1%v(1) == edge2%v(2)) v = edge1%v(1)
  if (edge1%v(2) == edge2%v(1)) v = edge1%v(2)
  if (edge1%v(2) == edge2%v(2)) v = edge1%v(2)

end function SharedNode

!===============================================================================

end module Spectral_Element_Edge
