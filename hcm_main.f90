program hcm_main
  use hcm_core
  use file_parser
  use iso_fortran_env, only: real64
  implicit none
  integer :: NN, i
  integer :: num_args
  real(real64), allocatable :: R(:,:)
  character(len=2), dimension(:), allocatable :: symbols
  character(len=200) :: cmd, arg
  character(len=256) :: infile, outfile
  !character(len=256) :: home_dir
  !character(len=256) :: arg_input
  !character(len=3) :: print_xyz
  real(real64) :: DF = 0.5, DT = 0.5, DP = 0.5, DDT = 0
  real(real64) :: cidx
  logical :: file_exists, print_xyz = .false.
  integer :: NT = 0
  integer :: ios

! Command line arguments handling

  num_args = command_argument_count() 

  i = 1
    do while (i < num_args)
        call get_command_argument(i, arg)

        select case (trim(arg))
        
        case ("-xyz")
          print_xyz = .true.
          i = i + 1

        case("-h", "--help")
          call print_usage()
          stop

        case default
            if (arg(1:1) == "-") then
              print *, "Unknown option: ", trim(arg)
              call print_usage()
              stop
            else 
              print *, "Unknown argument: ", trim(arg), " (only one file at a time is allowed)"
              call print_usage()
              stop
            end if

        end select
    end do

  call get_command_argument(num_args, infile)

!  arg_input = adjustl(arg_input)

!  print *, infile

!  if (infile == '~') then
!    call get_environment_variable("HOME", home_dir)
!        ! Combine the system home directory with everything after the '~'
!    infile = trim(home_dir) // trim(infile(2:))
!  else
!    infile = arg_input
!  end if

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

  cidx = compute_chiral_index(NN, R, DF, DT, DP, DDT, NT, print_xyz, symbols)

!  open(unit=20, file=outfile, status='replace')
!  write(20,*) 'Cindex=', cidx
!  close(20)

  print *, 'Cindex=', cidx

end program hcm_main

subroutine print_usage()
        print *, "Usage: ./haussCM [-xyz] INPUT_FILE"
        print *, "Options:"
        print *, "  -h, --help      Show this help message"
        print *, "  -xyz            Print xyz best fit coordinates for the enantiomers."
        print *, "                  Centered reference: 'centered_ref.xyz'."
        print *, "                  Best fit obtained enantiomer: 'enan_fit.xyz'."
        print *, "To modify the step sizes for the rotational exploration," 
        print *, "create a file named 'inchi.dat' in the current directory"
        print *, "with the following format:"
        print *, "DF DT DP DDT NT"
        print *, "DF, DT, DP -> Euler angle steps in degrees (default: 0.5 0.5 0.5)"
        print *, "DDT, NT -> Not used currently. Translation parameters (default: 0 0)."
end subroutine print_usage