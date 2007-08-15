module Boundary_Condition

  use Constants
  implicit none
  private

  type, public :: BoundaryCondition
    character(len=25)      :: name = '' ! boundary condition name

! Änderung 08.01.2007
    integer                :: type = 0  ! boundary condition type
    character              :: key = '' ! SEMTEX boundary group key
    character(len=25)      :: link = '' ! Associated bc name, only for periodic boundaries

  end type BoundaryCondition

  type(BoundaryCondition), public, save, allocatable :: ControlBoundaryCondition(:)

  public :: GetBndryConditionControlData

contains

!===============================================================================

subroutine GetBndryConditionControlData(unit)
  integer, intent(in) :: unit

! Änderung 08.01.2007
  integer	      :: i, N, type, ios
  character(len=25)   :: name, link
  character           :: key

  namelist/boundary_condition/ name, type, key, link

  ! count number of user boundary conditions
  N = 0
  rewind(unit)
  do
    name = ''
    read(unit,nml=boundary_condition,IOSTAT=ios)
    if (ios /= 0) exit
    if (LEN_TRIM(name) > 0) N = N + 1
  end do
  allocate( ControlBoundaryCondition(N) )
  rewind(unit)

  do i = 1, N
    name = ControlBoundaryCondition(i)%name
!   Änderung 08.01.2007
    type = ControlBoundaryCondition(i)%type
    key  = ControlBoundaryCondition(i)%key
    link = ControlBoundaryCondition(i)%link
    do
      name = ''
      read(unit,nml=boundary_condition)
      if (LEN_TRIM(name) > 0) exit
    end do
    ControlBoundaryCondition(i)%name = name
!   Änderung 08.01.2007
    ControlBoundaryCondition(i)%type = type
    ControlBoundaryCondition(i)%key  = key
    ControlBoundaryCondition(i)%link = link
  end do

end subroutine GetBndryConditionControlData

!===============================================================================

end module Boundary_Condition
