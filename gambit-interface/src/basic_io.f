module Basic_IO

  implicit none
  private
!
!--- control ------------------------------------------------------------------

! use [non]advancing io for "get" from default unit
  character(len=*), parameter ::  adv = 'NO' ! 'YES' or 'NO'

!
!--- public variables ---------------------------------------------------------

  integer, public ::  end_of_file = 0
  integer, public ::  end_of_record = 0
!
!--- public procedures --------------------------------------------------------

  public ::  new_unit
  ! get the number of an existing, not connected unit.

  public ::  search_key

  public ::  get
  interface get
    module procedure  &
        get_d_i,      & !  get default integer from default unit
        get_u_i,      & !  get default integer from specified unit
        get_d_rs,     & !  get single precision real from default unit
        get_u_rs,     & !  get single precision real from specified unit
        get_d_rd,     & !  get double precision real from default unit
        get_u_rd,     & !  get double precision real from specified unit
        get_d_l,      & !  get logical from default unit
        get_u_l         !  get logical real from specified unit
  end interface

!
!--- private module variables -------------------------------------------------

  integer, parameter ::  single = selected_real_kind(6)
  integer, parameter ::  double = selected_real_kind(13)
  character(len = max(digits(1), digits(1d0))) ::  buffer
  logical ::  initialized = .false.

contains

!--- module procedures --------------------------------------------------------

subroutine Initialize

  integer ::  unit, i

  if (initialized) return
  initialized = .true.
!
! end_of_file
!
  unit = 99
  open(unit, status='scratch', position='append')
  read(unit, *, iostat=end_of_file)
  rewind(unit)

!
! end_of_record
!
!  Änderung 13.03.2007
  write(unit,*)
  backspace(unit)
  read(unit, '(I2)', advance='NO', iostat=end_of_record) i
  close(unit)

end subroutine Initialize

!------------------------------------------------------------------------------

function new_unit()  result(unit)
  integer ::  unit
  logical ::  exists, opened
  integer ::  ios
  if (.not. initialized) call Initialize
  unit = 0
  do
      unit = unit + 1
      inquire( unit, exist=exists, opened=opened, iostat=ios )
      if (exists .and. .not. opened .and. ios == 0) return
  end do
  unit = -1
end function new_unit

!------------------------------------------------------------------------------

subroutine search_key( unit, key, set, iostat )

  integer, intent(in) ::  unit
  ! Unit number of the search file, which must be opened for
  ! formatted sequential read.

  character(len=*), intent(in) ::  key
  ! Keyword expression to search for.

  character(len=*), intent(in),  optional ::  set
  ! Set of separator characters. Ignored if unit is opened for
  ! unformatted read.

  integer, intent(out), optional ::  iostat

  logical   ::  open, ready
  integer   ::  ios
  character(len=20) ::  acc, fm, rd

  if (.not. initialized) call Initialize

  inquire(unit, opened=open, access=acc, form=fm, read=rd)

  ready = open .and. index(acc,'SEQUENTIAL') > 0 .and. index(rd, 'YES') > 0
  if (.not. ready) then
      write(*,*) '*** search_key failed!'
      write(*,*) 'unit ', unit,' not ready for sequential read:'
      write(*,*) 'opened = ', open
      write(*,*) 'access = "'//trim(acc)//'"  ("SEQUENTIAL" required)'
      write(*,*) 'read   = "'//trim(rd) //'"  ("YES" required)'
      stop
  end if

  if (index(fm,'UNFORMATTED') > 0) then
      call search_keyU( unit, adjustl(trim(key)), ios=ios )
  else if (present(set)) then
      if (len(set) == 0) then
          write(*,*) '*** search_key failed: set present but empty!'
          stop
      end if
      call search_keyF( unit, adjustl(trim(key)), set, ios )
  else
      call search_keyF( unit, adjustl(trim(key)), ios=ios )
  end if

  if (present(iostat)) then
      iostat = ios
  else if (ios /= 0) then
      write(*,*) '*** search_key: search for "'//key//'" failed:'
      if (ios == end_of_file) then
          write(*,*) 'Reached end of file.'
      else
          write(*,*) 'unit   = ', unit
          write(*,*) 'iostat = ', ios
          stop
      end if
  end if

end subroutine search_key

!..............................................................................

subroutine search_keyU( unit, key, ios )
  integer,          intent(in)  ::  unit
  character(len=*), intent(in)  ::  key
  integer,          intent(out) ::  ios
  character(len=len(key)) ::  b
  do
      read(unit, iostat=ios) b
      if (ios == end_of_file) exit
      if (b == key) exit
  end do
end subroutine search_keyU

!..............................................................................

subroutine search_keyF( unit, key, set, ios )

  integer,          intent(in)           ::  unit
  character(len=*), intent(in)           ::  key
  character(len=*), intent(in), optional ::  set
  integer,          intent(out)          ::  ios

  integer   ::  i, l
  character ::  c
  character(len=len(key)) ::  b

  l = len(key)
  search_records: do
      b = ' '
      i = 0
      search_key: do
          read(unit,fmt='(A)',advance='NO',iostat=ios) c
          if (ios == 0) then
              if (i == 0) then
                  if (c == ' ') cycle search_key
                  b(1:1) = c
              else
                  b = b(1:i)//c
              end if
              i = i + 1
              if (i >= l) exit search_key
          else if (ios == end_of_record) then
              cycle search_records
          else
              exit search_records
          end if
      end do search_key
      if (b /= key) cycle search_records
      if (present(set)) then
          search_separator: do
              read(unit,fmt='(A)',advance='NO',iostat=ios) c
              if (ios == 0) then
                  if (index(set,c) > 0) exit search_records !!! success !!!
              else if (ios == end_of_record) then
                  cycle search_records
              else
                  exit search_records
              end if
          end do search_separator
      else
          read(unit,'(A)',advance='YES')
          exit search_records !!! success !!!
      end if
  end do search_records

end subroutine search_keyF

!------------------------------------------------------------------------------

subroutine get_d_i( var, prompt, default, fmt )

  integer,          intent(inout)        ::  var
  character(len=*), intent(in)           ::  prompt
  integer,          intent(in), optional ::  default
  character(len=*), intent(in), optional ::  fmt

  integer ::  ios

  do
      write(*,'(A)',advance='NO') prompt
      if (present(default)) then
          write(*,'(A)',advance='NO') '['
          if (present(fmt)) then
              write(*,fmt,advance='NO') default
          else
              write(buffer,*) default
              write(*,'(A)',advance='NO')  trim(adjustl(buffer))
          end if
          write(*,'(A)',advance='NO') ']  '
      end if
      write(*,'(A)',advance=adv) ' '
      read(*,'(A)') buffer
      buffer = adjustl(buffer)
      if (len_trim(buffer) == 0 .and. present(default)) then
          var = default
          return
      end if
      read(buffer, *, iostat=ios) var
      if (ios <= 0)  exit
  end do

end subroutine get_d_i

!------------------------------------------------------------------------------

subroutine get_u_i( unit, var, iostat )

  integer, intent(in)  ::  unit
  integer, intent(out) ::  var
  integer, intent(out), optional ::  iostat

  integer ::  ios

  read(unit,'(A)') buffer
  buffer = adjustl(buffer)
  read(buffer, *, iostat=ios) var
  if (present(iostat)) then
    iostat = ios; return
  else if (ios /= 0) then
    write(*,*) ' Error No.', ios, &
               ' during GET of integer from UNIT ', unit
    stop
  end if

end subroutine get_u_i

!------------------------------------------------------------------------------

subroutine get_d_rd( var, prompt, default, fmt )

  real(double),     intent(out)          ::  var
  character(len=*), intent(in)           ::  prompt
  real(double),     intent(in), optional ::  default
  character(len=*), intent(in), optional ::  fmt

  integer ::  ios, pos

  do
      write(*,'(A)',advance='NO') prompt
      if (present(default)) then
          write(*,'(A)',advance='NO') '['
          if (present(fmt)) then
              write(*,fmt,advance='NO') default
          else
              write(buffer,*) default
              write(*,'(A)',advance='NO')  trim(adjustl(buffer))
          end if
          write(*,'(A)',advance='NO') ']  '
      end if
      write(*,'(A)',advance=adv) ' '
      read(*,'(A)') buffer
      buffer = adjustl(buffer)
      if (len_trim(buffer) == 0 .and. present(default)) then
          var = default
          return
      end if
      !
      ! add decimal point if necessary
      !
      if (scan(buffer,'.') == 0)  then
         buffer = adjustl(buffer);  pos = scan(buffer,'deDE')
         if ( pos == 0 )  then
             buffer = trim(buffer)//'.'
         else
             buffer = buffer(1:pos-1) // '.' // buffer(pos:)
         end if
      end if
      read(buffer, *, iostat=ios) var
      if (ios <= 0)  exit
  end do

end subroutine get_d_rd

!------------------------------------------------------------------------------

subroutine get_u_rd( unit, var, iostat )

  integer,      intent(in)  ::  unit
  real(double), intent(out) ::  var
  integer,      intent(out), optional ::  iostat

  integer ::  ios, pos

  read(unit,'(A)') buffer
  buffer = adjustl(buffer)
  !
  ! add decimal point if necessary
  !
  if (scan(buffer,'.') == 0)  then
     buffer = adjustl(buffer);  pos = scan(buffer,'deDE')
     if ( pos == 0 )  then
         buffer = trim(buffer)//'.'
     else
         buffer = buffer(1:pos-1) // '.' // buffer(pos:)
     end if
  end if
  read(buffer, *, iostat=ios) var
  if (present(iostat)) then
    iostat = ios; return
  else if (ios /= 0) then
    write(*,*) ' Error No.', ios, &
               ' during GET of double precision from UNIT ', unit
    stop
  end if

end subroutine get_u_rd

!------------------------------------------------------------------------------

subroutine get_d_rs( var, prompt, default, fmt )

  real(single),     intent(out)          ::  var
  character(len=*), intent(in)           ::  prompt
  real(single),     intent(in), optional ::  default
  character(len=*), intent(in), optional ::  fmt

  real(double) :: dvar, ddefault

  if (present(default))  then
    ddefault = default
    if (present(fmt))  then
       call get_d_rd( dvar, prompt, default=ddefault, fmt=fmt )
    else
       call get_d_rd( dvar, prompt, default=ddefault )
    end if
  else
    call get_d_rd( dvar, prompt )
  end if
  var = dvar

end subroutine get_d_rs

!------------------------------------------------------------------------------

subroutine get_u_rs( unit, var, iostat )

  integer,      intent(in)  ::  unit
  real(single), intent(out) ::  var
  integer,      intent(out), optional ::  iostat

  integer ::  ios
  real(double) ::  dvar

  call get_u_rd( unit, dvar, iostat=ios )
  if (ios == 0)  var = dvar
  if (present(iostat)) then
    iostat = ios; return
  else if (ios /= 0) then
    write(*,*) ' Error No.', ios, &
               ' during GET of real from UNIT ', unit
    stop
  end if

end subroutine get_u_rs

!------------------------------------------------------------------------------

subroutine get_d_l( var, prompt, default, fmt )

  logical,          intent(inout)        ::  var
  character(len=*), intent(in)           ::  prompt
  logical,          intent(in), optional ::  default
  character(len=*), intent(in), optional ::  fmt

  integer ::  ios

  do
      write(*,'(A)',advance='NO') prompt
      if (present(default)) then
          write(*,'(A)',advance='NO') '['
          if (present(fmt)) then
              write(*,fmt,advance='NO') default
          else
              write(*,'(L1)',advance='NO') default
          end if
          write(*,'(A)',advance='NO') ']  '
      end if
      write(*,'(A)',advance=adv) ' '
      read(*,'(A)') buffer
      buffer = adjustl(buffer)
      if (len_trim(buffer) == 0 .and. present(default)) then
          var = default
          return
      end if
      read(buffer, '(L2)', iostat=ios) var
      if (ios <= 0)  exit
  end do

end subroutine get_d_l

!------------------------------------------------------------------------------

subroutine get_u_l( unit, var, iostat )

  integer, intent(in)  ::  unit
  logical, intent(out) ::  var
  integer, intent(out), optional ::  iostat

  integer ::  ios

  read(unit,'(A)') buffer
  buffer = adjustl(buffer)
  read(buffer, '(L2)', iostat=ios) var
  if (present(iostat)) then
    iostat = ios; return
  else if (ios /= 0) then
    write(*,*) ' Error No.', ios, &
               ' during GET of logical from UNIT ', unit
    stop
  end if

end subroutine get_u_l

!------------------------------------------------------------------------------

end module Basic_IO
