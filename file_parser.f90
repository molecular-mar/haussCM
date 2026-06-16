module file_parser
  use omp_lib
  use iso_fortran_env, only: real64
  implicit none
contains
    subroutine parse_xyz(filename, natoms, symbols, coords)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: natoms
      character(len=:), allocatable :: comment
      character(len=2), dimension(:), allocatable, intent(out) :: symbols
      real(real64), dimension(:,:), allocatable, intent(out) :: coords

      integer :: file_unit = 10
      integer :: i, ios
      character(len=256) :: line

      open(unit=file_unit, file=filename, status='OLD', action='READ', iostat=ios)
      if (ios /= 0) then
       print *, "Error opening file: ", TRIM(filename)
       stop
      end if
      
      read(file_unit, *, iostat=ios) natoms
      if (ios /= 0) return
      
      read(file_unit, '(A)', iostat=ios) line
      if (ios /= 0) return
      comment = trim(adjustl(line))

      allocate(symbols(natoms))
      allocate(coords(natoms, 3))

      do i = 1, natoms
        read(file_unit, *, iostat=ios) symbols(i), coords(i,1), coords(i,2), coords(i,3)
        if (ios /= 0) exit
      end do

    close(file_unit)
    end subroutine parse_xyz
end module file_parser
