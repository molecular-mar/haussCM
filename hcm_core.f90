module hcm_core
  use omp_lib
  use iso_fortran_env, only: real64
  implicit none
contains
  ! Helper: compute normalized Hausdorff (sabmax) between `Rcopy` and rotated `R1`
  subroutine best_coords(NN, R1, ii, jj, kk, DF, DT, DP, SYM)
    use file_parser
    implicit none
    integer, intent(in) :: NN, ii, jj, kk
    real(real64), intent(in) :: R1(NN,3)
    real(real64), allocatable :: Rbest(:,:)
    real(real64) :: alpha, beta, gamma
    real(real64) :: sabmax
    real(real64) :: pi, OGRA, DELTAX, DELTAY, DELTAZ
    integer :: i,j,k,ia
    real(real64) :: RM(3,3)
    real(real64), intent(in) :: DF, DT, DP
    real(real64), allocatable :: R2(:,:), D1(:), D1MIN(:)
    real(real64) :: AAMAX, BBMAX, abmax

    character(len=30) :: filetest
    character(len=2), intent(in) :: SYM(:)
    !character(len=2), allocatable :: symbols(:) 

    allocate(Rbest(NN,3), D1(NN), D1MIN(NN))
  
    filetest = 'centered_ref.xyz'
    !allocate(symbols(NN))
    !symbols = [ ('6', ia = 1, NN) ]
    call write_xyz(filetest, NN, SYM, R1)

    pi = 4.0_real64*atan(1.0_real64)
    OGRA = (2.0_real64*pi)/360.0_real64
    DELTAX = OGRA*DF
    DELTAY = OGRA*DT
    DELTAZ = OGRA*DP

    gamma = (ii - 1) * DELTAZ
    beta = (jj - 1) * DELTAY
    alpha = (kk - 1) * DELTAX

    RM(1,1) = cos(gamma)*cos(alpha) - cos(beta)*sin(alpha)*sin(gamma)
    RM(1,2) = cos(gamma)*sin(alpha) + cos(beta)*cos(alpha)*sin(gamma)
    RM(1,3) = sin(gamma)*sin(beta)
    RM(2,1) = -sin(gamma)*cos(alpha) - cos(beta)*sin(alpha)*cos(gamma)
    RM(2,2) = -sin(gamma)*sin(alpha) + cos(beta)*cos(alpha)*cos(gamma)
    RM(2,3) = cos(gamma)*sin(beta)
    RM(3,1) = sin(beta)*sin(alpha)
    RM(3,2) = -sin(beta)*cos(alpha)
    RM(3,3) = cos(beta)

    do i=1,NN
      do j=1,3
        Rbest(i,j) = 0.0_real64
      end do
      do j=1,3
        do k=1,3
          Rbest(i,j) = Rbest(i,j) + RM(j,k) * R1(i,k)
        end do
      end do
    end do

    filetest = 'enan_fit.xyz'
    call write_xyz(filetest, NN, SYM, Rbest)

  end subroutine

  function rotation_sabmax(NN, Rcopy, R1, alpha, beta, gamma, AN) result(sabmax)
    implicit none
    integer, intent(in) :: NN
    real(real64), intent(in) :: Rcopy(NN,3), R1(NN,3)
    real(real64), intent(in) :: alpha, beta, gamma, AN
    real(real64) :: sabmax
    integer :: i,j,k
    real(real64) :: RM(3,3)
    real(real64), allocatable :: R2(:,:), D1(:), D1MIN(:)
    real(real64) :: AAMAX, BBMAX, abmax

    allocate(R2(NN,3), D1(NN), D1MIN(NN))

    RM(1,1) = cos(gamma)*cos(alpha) - cos(beta)*sin(alpha)*sin(gamma)
    RM(1,2) = cos(gamma)*sin(alpha) + cos(beta)*cos(alpha)*sin(gamma)
    RM(1,3) = sin(gamma)*sin(beta)
    RM(2,1) = -sin(gamma)*cos(alpha) - cos(beta)*sin(alpha)*cos(gamma)
    RM(2,2) = -sin(gamma)*sin(alpha) + cos(beta)*cos(alpha)*cos(gamma)
    RM(2,3) = cos(gamma)*sin(beta)
    RM(3,1) = sin(beta)*sin(alpha)
    RM(3,2) = -sin(beta)*cos(alpha)
    RM(3,3) = cos(beta)

    do i=1,NN
      do j=1,3
        R2(i,j) = 0.0_real64
      end do
      do j=1,3
        do k=1,3
          R2(i,j) = R2(i,j) + RM(j,k) * R1(i,k)
        end do
      end do
    end do

    ! S -> R2
    do i=1,NN
      do j=1,NN
        D1(j) = sqrt( (Rcopy(i,1)-R2(j,1))**2 + (Rcopy(i,2)-R2(j,2))**2 + (Rcopy(i,3)-R2(j,3))**2 )
      end do
      D1MIN(i) = minval(D1)
    end do
    AAMAX = maxval(D1MIN)

    ! R2 -> S
    do j=1,NN
      do i=1,NN
        D1(i) = sqrt( (Rcopy(i,1)-R2(j,1))**2 + (Rcopy(i,2)-R2(j,2))**2 + (Rcopy(i,3)-R2(j,3))**2 )
      end do
      D1MIN(j) = minval(D1)
    end do
    BBMAX = maxval(D1MIN)

    abmax = max(AAMAX, BBMAX)
    sabmax = abmax / AN

    deallocate(R2, D1, D1MIN)
  end function rotation_sabmax

  function compute_chiral_index(NN, R, DF, DT, DP, DDT, NT, XYZ_P, SYM) result(best)
    use file_parser
    implicit none
    integer, intent(in) :: NN, NT
    real(real64), intent(in) :: R(NN,3)
    real(real64), intent(in) :: DF, DT, DP, DDT
    logical, intent(in) :: XYZ_P
    real(real64) :: best
    character(len=2), dimension(:), allocatable, intent(in) :: SYM

    integer :: i,j,k,ii,jj,kk,local_ii,local_jj,local_kk,BEST_II,BEST_JJ,BEST_KK,ia
    integer :: NA,NB,NC
    integer :: nmd, nm, nme, ij, ik, il
    real(real64) :: pi, OGRA, DELTAX, DELTAY, DELTAZ
    real(real64) :: dw, w1, w2, w3
    real(real64) :: ACC, ANTEMAX, local_antemax
    real(real64) :: X, Y, Z  ! cumulative Euler angles (persistent across iterations)
    real(real64), allocatable :: RA(:,:), R1(:,:), R2(:,:)
    real(real64), allocatable :: D1(:), D1MIN(:)
    real(real64), allocatable :: DA1(:), DA2(:)
    real(real64) :: AN, AUX
    real(real64) :: alpha,beta,gamma
    real(real64) :: RM(3,3)
    real(real64) :: AAMAX, BBMAX, AMAX, BMAX
    real(real64) :: abmax, sabmax
    real(real64) :: sum1, sum2, sum3, asum1, asum2, asum3
    real(real64), allocatable :: Rcopy(:,:)  ! local copy of R for centering

    character(len=30) :: filetest
    !character(len=2), allocatable :: symbols(:) 

    allocate(RA(NN,3), R1(NN,3), R2(NN,3), D1(NN), D1MIN(NN), DA1(NN), DA2(NN), Rcopy(NN,3))

    ! Copy R to local array
    do i = 1, NN
      do j = 1, 3
        Rcopy(i,j) = R(i,j)
      end do
    end do
    
    !filetest = 'test_coords1.xyz'
    !allocate(symbols(NN))
    !symbols = [ ('6', ia = 1, NN) ]
    !call write_xyz(filetest, NN, symbols, Rcopy)


    ! Center coordinates relative to center of mass
    sum1 = 0.0_real64
    sum2 = 0.0_real64
    sum3 = 0.0_real64
    do i = 1, NN
      sum1 = sum1 + Rcopy(i,1)
      sum2 = sum2 + Rcopy(i,2)
      sum3 = sum3 + Rcopy(i,3)
    end do
    asum1 = sum1 / dble(NN)
    asum2 = sum2 / dble(NN)
    asum3 = sum3 / dble(NN)
    do i = 1, NN
      Rcopy(i,1) = Rcopy(i,1) - asum1
      Rcopy(i,2) = Rcopy(i,2) - asum2
      Rcopy(i,3) = Rcopy(i,3) - asum3
    end do

    !filetest = 'test_coords2.xyz'
    !call write_xyz(filetest, NN, symbols, Rcopy)

    ! Mirror image (using centered coordinates)
    do i=1,NN
      do j=1,3
        RA(i,j) = - Rcopy(i,j)
      end do
    end do
    
    !filetest = 'test_coords3_enan.xyz'
    !call write_xyz(filetest, NN, symbols, RA)

    ! Compute diameter AN (max pairwise distance) using centered coordinates
    do i=1,NN
      do j=1,NN
        DA1(j) = sqrt( (Rcopy(i,1)-Rcopy(j,1))**2 + (Rcopy(i,2)-Rcopy(j,2))**2 + (Rcopy(i,3)-Rcopy(j,3))**2 )
      end do
      DA2(i) = maxval(DA1)
    end do
    AN = maxval(DA2)
    if (AN == 0.0_real64) AN = 1.0_real64
    print *, "Computed diameter = ", AN

    pi = 4.0_real64*atan(1.0_real64)
    OGRA = (2.0_real64*pi)/360.0_real64
    DELTAX = OGRA*DF
    DELTAY = OGRA*DT
    DELTAZ = OGRA*DP

    NA = max(1, int(360.0_real64/DF + 0.5_real64))
    NB = max(1, int(180.0_real64/DT + 0.5_real64))
    NC = max(1, int(360.0_real64/DP + 0.5_real64))

    ! Initialize cumulative angle variables
    !X = -DELTAX
    !Y = -DELTAY
    !Z = -DELTAZ

    ! translation grid setup 
    if (DDT == 0.0_real64 .or. NT == 0) then
      nmd = 0
    else
      nmd = 2*NT
    end if
    dw = DDT
    nm = NT
    w1 = - (NT+1) * dw
    w2 = w1
    w3 = w1
    nme = nmd + 1

    ACC = 1.0_real64

    do ij = 1, nme
      w1 = w1 + dw
      ! set first coordinate of R1
      do i = 1, NN
        R1(i,1) = RA(i,1) + w1
      end do
      w2 = - (nm+1) * dw
      do ik = 1, nme
        w2 = w2 + dw
        do i = 1, NN
          R1(i,2) = RA(i,2) + w2
        end do
        w3 = - (nm+1) * dw
        do il = 1, nme
          w3 = w3 + dw
          do i = 1, NN
            R1(i,3) = RA(i,3) + w3
          end do

          ! per-translation best (ANTEMAX)
          ANTEMAX = 1.0_real64

          !$OMP PARALLEL DEFAULT(SHARED) & 
          !$OMP PRIVATE(ii,jj,kk,X,Y,Z,alpha,beta,gamma) &
          !$OMP PRIVATE(local_antemax, local_ii, local_jj, local_kk)

          local_antemax = HUGE(local_antemax)
          local_ii = -1
          local_jj = -1
          local_kk = -1

          !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
          do ii = 1, NA
            do jj = 1, NB
              do kk = 1, NC
                Z = (ii - 1) * DELTAZ
                Y = (jj - 1) * DELTAY
                X = (kk - 1) * DELTAX
                sabmax = rotation_sabmax(NN, Rcopy, R1, X, Y, Z, AN)
                if (sabmax < local_antemax) then 
                    local_antemax = sabmax
                    local_ii = ii
                    local_jj = jj
                    local_kk = kk
                end if
              end do
            end do
          end do
          !$OMP END DO NOWAIT

          !$OMP CRITICAL
          if (local_antemax < ANTEMAX) then
            ANTEMAX = local_antemax
            BEST_II = local_ii
            BEST_JJ = local_jj
            BEST_KK = local_kk
          end if
          !$OMP END CRITICAL

          !$OMP END PARALLEL

          ! Track minimum across translations
          if (ANTEMAX < ACC) ACC = ANTEMAX

        end do
      end do
    end do

    if (XYZ_P) call best_coords(NN, RA, BEST_II, BEST_JJ, BEST_KK, DF, DT, DP, SYM)

    best = ACC

    deallocate(RA, R1, R2, D1, D1MIN, DA1, DA2, Rcopy)
  end function compute_chiral_index

end module hcm_core
