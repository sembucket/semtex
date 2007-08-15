module Spectral_Element_Vertex

  use Constants
  implicit none
  private

  type, public :: SpectralElementVertex
     real(RNP) ::  x(2) ! vertex coordinates
  end type SpectralElementVertex

end module Spectral_Element_Vertex
