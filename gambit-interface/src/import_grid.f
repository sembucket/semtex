module Import_Grid

  use Constants
  use Basic_IO
  use Data_Structure
  use Spectral_Element_Vertex
  use Spectral_Element_Edge
  use Boundary_Segment
  use Spectral_Element
  implicit none
  private

  type, public :: ImportControlData
    character(len=80)	::  File = ''	 ! file name
    integer		::  Type = 0	 ! 0 = Fluent/UNS
    real(RNP)		::  Scale(2) = 1 ! scaling vector
  end type ImportControlData

  type(ImportControlData), save :: ImportData

! Constants for gambit grid import
  character(len=25)	::  def_BoundaryName = 'Default_Boundary'
  character(len=5)	::  sep = '() #/'	! separators
  integer		::  BZId = 3		! boundary zone ID
  integer		::  IZId = 2		! interior zone ID
  integer		::  MZN  = 1000		! maximum zone numbers
  integer		::  def_BoundaryType = 1 ! default boundary type

  public :: GetImportControlData
  public :: ImportGrid

contains

!===============================================================================

subroutine GetImportControlData(unit)
  integer, intent(in)		 :: unit
  integer			 :: type, ios
  character(len=80)		 :: file
  real(RNP)			 :: xScale, yScale

  namelist/grid/ file, type, xScale, yScale

  file   = ImportData%file
  type   = ImportData%type
  xScale = ImportData%scale(1)
  yScale = ImportData%scale(2)

  rewind(unit)
  read(unit,nml=grid,IOSTAT=ios)

  ImportData%file     = file
  ImportData%type     = type
  ImportData%scale(1) = xScale
  ImportData%scale(2) = yScale

end subroutine GetImportControlData

!===============================================================================

subroutine ImportGrid(vertex,edge,elmt,bseg)
  type(SpectralElementVertex), pointer :: vertex(:)
  type(SpectralElementEdge), pointer   :: edge(:)
  type(SpectralElement), pointer       :: elmt(:)
  type(BoundarySegment), pointer       :: bseg(:)

  select case(ImportData%Type)
  case(0) 	! Fluent/UNS - Format
    call ImportFluentUNS(vertex,edge,elmt,bseg)
  case(1:)	! not implemented
    STOP 'Error: This file format is not supported.'
  end select

end subroutine ImportGrid

!===============================================================================

subroutine ImportFluentUNS(vertex,edge,elmt,bseg)
  type(SpectralElementVertex), pointer :: vertex(:)
  type(SpectralElementEdge), pointer   :: edge(:)
  type(SpectralElement), pointer       :: elmt(:)
  type(BoundarySegment), pointer       :: bseg(:)

  call ReadinMeshFile(vertex,edge,elmt,bseg)
  call ConstructElements(vertex,edge,elmt)

end subroutine ImportFluentUNS

!===============================================================================

subroutine ReadinMeshFile(vertex,edge,elmt,bseg)
  type(SpectralElementVertex), pointer :: vertex(:)
  type(SpectralElementEdge), pointer   :: edge(:)
  type(SpectralElement), pointer       :: elmt(:)
  type(BoundarySegment), pointer       :: bseg(:)
  character(len=80)		       :: file
  character(len=25)		       :: tmp(5), hexa = '(Z25)', fmt
  integer			       :: unit, ios, i, pos
  integer			       :: GZN, BZN, ztype, GZL(MZN)
  integer			       :: Nstart, Nend, N, data(5)
  real(RNP)			       :: vScale(2) !Scale Factor

! setup basic information
  vScale = ImportData%Scale
  file	 = ImportData%File

  write(*,'(A)')   '  Reading in the mesh:'
  write(*,'(2A)') '  ... ', file
  write(*,'(A,/)') repeat('*',45)

! open infile
  unit = new_unit()
  open(unit, IOSTAT=ios, file=file, action='READ', status='OLD')
  if (ios /= 0) STOP 'An error occurred while reading the gridfile.'

! read dimension
  call search_key(unit, '(2', set=' ')
  read(unit,'(Z1)') N
  if (N /= 2) STOP 'Error: Inputfile contains no 2d-data.'

! read vertex information
  call search_key(unit, '(10', set='(', iostat=ios)
  if (ios /= 0) STOP 'Error: no vertices exists.'
  read(unit,*) tmp(1:4)
  read(tmp(3),hexa) N
  write(*,'(A,I0)') '  Number of vertices:   ', N

! init vertex data
  allocate(vertex(1:N))

  call search_key(unit, '(10', set='(', iostat=ios)
  if (ios /= 0) STOP 'Error: irregular vertice section.'
  read(unit,*)

! read vertex coordinates
  do i = 1, N
    vertex(i)%x(1) = 0
    read(unit,*) vertex(i)%x(1), vertex(i)%x(2)
  end do

! scale vertex coordinates
  do i = 1, N
    vertex(i)%x(1) = vertex(i)%x(1) * vScale(1)
    vertex(i)%x(2) = vertex(i)%x(2) * vScale(2)
  end do

! read edge data
  BZN = 0 ! number of boundary zones
  GZL = 0 ! list of gambit zones

  call search_key(unit, '(13', set='(')
  if (ios /= 0) STOP 'Error: no edges exists.'

  read(unit,*) tmp(1:3)
  read(tmp(3),hexa) N
  write(*,'(A,I0)') '  Number of edges:      ', N

! init edge data
  allocate( edge(1:N) )

! switch through unknown number of edge zones
  do
    call search_key(unit, '(13', set='(', iostat=ios)
    if (ios /= 0) exit

    read(unit,*) tmp(1:4)
    read(tmp(1),hexa) GZN    ! Gambit zone number
    read(tmp(2),hexa) Nstart
    read(tmp(3),hexa) Nend
    read(tmp(4),hexa) ztype  ! zone type

    if (ztype == BZId) then
      BZN = BZN + 1
      if ( BZN <= MZN ) then; GZL(GZN) = BZN; else
        STOP 'Error: Too many zone numbers.'
      end if
    end if

    ! read connection
    do i = Nstart, Nend
      read(unit,*) tmp
      read(tmp,hexa) data
      if (data(1) /= 2) STOP 'Error: Irregular edge description.'
      edge(i)%b(1) = BZN
      edge(i)%b(2) = 0
      edge(i)%v(1) = data(2)
      edge(i)%v(2) = data(3)
      if (minval(data(4:5)) == 0) then
        edge(i)%l(1) = maxval(data(4:5))
        edge(i)%l(2) = 0
      else
        if (data(4) > data(5)) then
          edge(i)%l(1) = data(5)
          edge(i)%l(2) = data(4)
        else
          edge(i)%l(1) = data(4)
          edge(i)%l(2) = data(5)
        end if
      end if

      ! sort element list
      if (edge(i)%l(1) == 0) then
        if (edge(i)%l(2) == 0) &
        STOP 'Error: No elements are associated to the edge.'
        edge(i)%l(1) = edge(i)%l(2)
        edge(i)%l(2) = 0
      end if

      ! set interior id
      if (ztype == IZId) edge(i)%b(1) = 0

      ! check zone id
      if (ztype /= BZId .AND. ztype /= IZId) Stop 'Error: Undefined edge zone.'
    end do

    ! exit if NE elements read in
    if (Nend == N) exit
  end do

! read element information
  call search_key(unit, '(12', set='(', iostat=ios)
  if (ios /= 0) STOP 'Error: no elements exists.'

  read(unit,*) tmp(1:3)
  read(tmp(3),hexa) N

  write(*,'(A,I0)') '  Number of elements:   ', N

! init element data
  allocate(elmt(1:N))

  call search_key(unit, '(12', set='(', iostat=ios)
  if (ios /= 0) STOP 'Error: irregular element section.'

  read(unit,*) tmp
  pos = scan(tmp(5),sep) - 1

  read(tmp(1:4),     hexa) data(1:4)
  read(tmp(5)(1:pos),hexa) data(5)

  if (data(5) /= 3) STOP 'Error: no quadrangular elements defined.'

! init boundary segments
  allocate(bseg(1:BZN))
  bseg%Name = def_BoundaryName

! read zone information
  i = 0
  do
    call search_key(unit, '(45', set='(', iostat=ios)
    if (ios /= 0) exit

    read(unit,*) tmp(1:3)
    pos = scan(tmp(3),sep) - 1

    write(fmt,'(A,I0,A)') '(I', LEN_TRIM(tmp(1)), ')'
    read(tmp(1),fmt) i
    if ( GZL(i) /= 0 ) then
      bseg(GZL(i))%Name = tmp(3)(1:pos)
    end if
  end do
  write(*,'(A,I0)') '  Number of boundaries: ', BZN

  close(unit)

  write(*,'(/,A)') repeat('*',45)
  write(*,'(A)') '  done'
  write(*,'(A)') repeat('*',45)

end subroutine ReadinMeshFile

!===============================================================================

subroutine ConstructElements(vertex,edge,elmt)
  type(SpectralElementVertex), intent(in)  :: vertex(:)
  type(SpectralElementEdge), intent(inout) :: edge(:)
  type(SpectralElement), intent(inout)	   :: elmt(:)
  integer				   :: i, elmtID

  write(*,'(A)',advance='NO') '  Construct elements from edges ...   '

! find all edges owned by the element
  do i = 1, size(edge)
    elmtID = edge(i)%l(1)
    if (elmtID /= 0) then
      if ( FindUndefinedEdge(elmt(elmtID)) == 0 ) &
        STOP 'Error: All edges are already defined.'
    elmt(elmtID)%e(FindUndefinedEdge(elmt(elmtID))) = i
    end if
    elmtID = edge(i)%l(2)
    if (elmtID /= 0) then
      if ( FindUndefinedEdge(elmt(elmtID)) == 0 ) &
        STOP 'Error: All edges are already defined.'
      elmt(elmtID)%e(FindUndefinedEdge(elmt(elmtID))) = i
    end if
  end do

! sort edges in an anticlockwise or clockwise direction
  do i = 1, size(elmt)
    call SortElementEdge(edge(elmt(i)%e),elmt(i))
    call CheckElementEdge(edge,elmt(i),i)
  end do

! ensure element orientation
  do i = 1, size(elmt)
    if (ElementOrientation(vertex(elmt(i)%v)).EQV..FALSE.) &
      call FlipElementOrientation(elmt(i))
  end do

! check element orientation
  do i = 1, size(elmt)
    if (ElementOrientation(vertex(elmt(i)%v)).EQV..FALSE.) &
      STOP 'Error: Element orientation doesnt fit.'
  end do

! ensure edge normal direction
  do i = 1, size(elmt)
    call SetEdgeOrientation(edge,elmt(i),i)
  end do

  write(*,'(A)') 'done'
  write(*,'(A)') repeat('*',45)

end subroutine ConstructElements

!===============================================================================

subroutine SetEdgeOrientation(edge,elmt,elmtID)
  type(SpectralElementEdge), intent(inout) :: edge(:)
  type(SpectralElement), intent(in)	   :: elmt
  integer, intent(in)			   :: elmtID
  integer				   :: i, v(2)

  do i = 1, 4
    if (minval(edge(elmt%e(i))%l) == elmtID .OR. &
        minval(edge(elmt%e(i))%l) == 0) then
      if (edge(elmt%e(i))%v(1) /= elmt%v(i)) then
        v = edge(elmt%e(i))%v
        edge(elmt%e(i))%v(1) = v(2)
        edge(elmt%e(i))%v(2) = v(1)
      end if
    end if
  end do

end subroutine SetEdgeOrientation

!===============================================================================

subroutine SortElementEdge(edge,elmt)
  type(SpectralElementEdge), intent(in) :: edge(4)
  type(SpectralElement), intent(inout)  :: elmt
  integer				:: i, j
  integer				:: eo(4), en(4), eg(4), v(4)

! set up connectivity
  eo    = (/1, 2, 3, 4/)
  en    = 0
  en(1) = 1
  eo(1) = 0
  v(1:2) = edge(1)%v(1:2)

  do i = 2, 4
    do j = 2, 4
      if (eo(j) /= 0) then
        if (v(i) == edge(j)%v(1) .OR. v(i) == edge(j)%v(2)) then
          en(i) = j
          eo(j) = 0
          if (i /= 4) then
            v(i+1) =  edge(j)%v(1)
            if (v(i) == v(i+1)) v(i+1) = edge(j)%v(2)
          end if
          exit
        end if
      end if
    end do
  end do

  eg = elmt%e(en)
  elmt%e = eg
  elmt%v = v

end subroutine SortElementEdge

!===============================================================================

end module Import_Grid
