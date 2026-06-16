program hcm_main
  use hcm_core
  use file_parser
  use iso_fortran_env, only: real64
  implicit none
  integer :: NN, i
  real(real64), allocatable :: R(:,:)
  character(len=2), dimension(:), allocatable :: symbols
  character(len=200) :: cmd
  character(len=30) :: infile, outfile
  real(real64) :: DF = 0.5, DT = 0.5, DP = 0.5, DDT = 0
  real(real64) :: cidx
  logical :: file_exists
  integer :: NT = 0
  integer :: ios

  call get_command_argument(1, infile)
! call get_command_argument(2, outfile)

  inquire(file='inchi.dat', exist=file_exists)

  if (file_exists) then
    open(unit=11, file='inchi.dat', status='old', action='read')
    read(11,*)
    read(11,*) DF, DT, DP, DDT, NT
    close(11)
  end if

  call parse_xyz(infile, NN, symbols, R)

  print *, 'Starting chiral index computation...'

  cidx = compute_chiral_index(NN, R, DF, DT, DP, DDT, NT)

!  open(unit=20, file=outfile, status='replace')
!  write(20,*) 'Cindex=', cidx
!  close(20)

  print *, 'Cindex=', cidx

end program hcm_main
