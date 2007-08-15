module Control_Data

  use Constants
  use Basic_IO, only: new_unit
  use Spectral_Element_Vertex
  use Spectral_Element_Edge
  use Spectral_Element
  use Boundary_Curve
  use Boundary_Condition
  use Import_Grid
  use Export_Functions
  implicit none
  private

  public :: GetControlData

contains

!===============================================================================

subroutine GetControlData(casefile)
  character(len=*), intent(in)  :: casefile
  integer			:: unit, ios

  write(*,'(A)') repeat('*',45)
  write(*,'(2A)')   '  Using the casefile: ' , casefile
  write(*,'(A)') repeat('*',45)
  unit = new_unit()
  open(unit, IOSTAT=ios, file=casefile, action='READ', status='OLD')
  if (ios /= 0) STOP 'Error: Could not open the casefile.'
  call GetImportControlData(unit)
  call GetBndryCurveControlData(unit)
  call GetBndryConditionControlData(unit)
  call GetExportControlData(unit)
  close(unit)

end subroutine GetControlData

!===============================================================================

end module  Control_Data
