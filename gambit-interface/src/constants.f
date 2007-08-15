module Constants                                                              !#

! kind parameters                                                             !#
  integer, parameter ::  RLP = selected_real_kind( 5) ! real low,             !#
  integer, parameter ::  RNP = selected_real_kind( 8) ! normal, and           !#
  integer, parameter ::  RHP = selected_real_kind( 8) ! high precision        !#

! numeric constants                                                           !#
  real(RHP), parameter ::  pi = 3.1415926535897932384626433832795029_RHP      !#
  real(RNP), parameter ::  one = 1._RNP, zero = 0._RNP, half = 0.5_RNP        !#

end module Constants
