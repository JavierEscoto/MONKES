module API_Example_Adjoint_DKE_BTD_Solution_Legendre

  use Adjoint_DKE_BTD_Solution_Legendre
  
  implicit none
  
  private
  
  public :: Adjoint_Monoenergetic_Database_Input
  
  
  contains  
  
  
  
  ! *** Computes a monoenergetic coefficients database reading the input file
  ! "monkes_input.parameters"
  subroutine Adjoint_Monoenergetic_Database_Input
     integer, parameter :: N_nu_max = 500, N_E_r_max = 500
     integer, parameter :: M_theta_max = 500, M_zeta_max = 500, M_xi_max = 500
     real :: nu(N_nu_max), E_r(N_E_r_max)
     integer :: N_theta(M_theta_max), N_zeta(M_zeta_max) , N_xi(M_xi_max) 
     integer :: N_nu, N_E_r, M_theta, M_zeta, M_xi, ierr 
     namelist /parameters/N_theta, N_zeta, N_xi, nu, E_r       
     
     ! *** Read input parameters from "monkes.parameters" file
      N_theta = -1 ; N_zeta = -1 ; N_xi = -1 ; nu = -1d14 ; E_r = -1d14 
     open(1, file= "monkes_input.parameters", status="old") 
     read(1, nml=parameters, iostat=ierr)     
     close(1)  
     
     ! Count number of collisionalities and radial electric field to be
     ! includede in the database
     M_theta = count(N_theta > 0)  
     M_zeta = count(N_zeta > 0)  
     M_xi = count(N_xi > 1)  
     N_nu = count(nu > 0) ; N_E_r = count(E_r /= -1d14) 
     
     where( mod(N_theta(1:M_theta),2) == 0 )  N_theta(1:M_theta) = N_theta(1:M_theta) + 1
     where( mod(N_zeta(1:M_zeta),2) == 0 )    N_zeta(1:M_zeta) = N_zeta(1:M_zeta) + 1
     
     call Adjoint_Monoenergetic_Database_Scan( nu(1:N_nu), E_r(1:N_E_r),       &
                                               N_theta(1:M_theta),             &
                                               N_zeta(1:M_zeta),               &
                                               N_xi(1:M_xi) )
     
  end subroutine
  
  ! ********************************************************************
  ! Monoenergetic_Database_Scan( nu, E_r, N_theta, N_zeta, N_xi, DD, DD_33_Sp )
  !   Computes the monoenergetic coefficients D(i,j) for the set of specified
  !   collisionalities, radial electric field, resolution in theta, zeta, xi
  ! 
  ! - INPUTS:
  !     nu(N_nu): 
  !       Real vector of collision frequency divided by speed in [m^-1]
  !     E_r(N_E_r):  
  !       Real vector of radial electric field divided by speed in [kV s /m^2]
  !     N_theta(M_theta): 
  !       Integer vector (>0) of number of discretization points in the 
  !       Boozer poloidal angle theta. Must be an odd number for numerical stability. 
  !     N_zeta(M_zeta): 
  !       Integer vector (>0) of number of discretization points in the 
  !       Boozer toroidal angle zeta. Must be an odd number for numerical stability. 
  !     N_xi(M_xi):  
  !       Integer vector (>0) of number of Legendre modes for the pitch-angle cosine xi. 
  !
  ! - OPTIONAL OUTPUTS: (WILL BE NEEDED FOR MOMENTUM CONSERVATION)
  !     DD(3,3, N_nu, N_E_r, M_theta, M_zeta, M_xi) : 
  !       Real array of the monoenergetic geometric coefficients \hat{D}_{ij}. The
  !       normalization used for the monoenergetic coefficients is explained
  !       in [Escoto, NF (2023?)] (Provisional)
  !     DD_33_Sp(N_nu, N_E_r, M_theta, M_zeta, M_xi) : 
  !       Real array of the Sptizer conductivity monoenergetic geometric coefficient.
  !
  ! - OUTPUTS GIVEN IN OUTPUT FILE "monkes_Monoenergetic_Database.dat":
  !     D(3,3, N_nu, N_E_r, M_theta, M_zeta, M_xi) : 
  !       Real array of the monoenergetic geometric coefficients \hat{D}_{ij}. The
  !       normalization used for the monoenergetic coefficients is explained
  !       in [Escoto, NF (2023?)] (Provisional)
  !     D_33_Sp(N_nu, N_E_r, M_theta, M_zeta, M_xi) : 
  !       Real array of the Sptizer conductivity monoenergetic geometric coefficient.
  ! ********************************************************************
  subroutine Adjoint_Monoenergetic_Database_Scan( nu, E_r, N_theta, N_zeta, N_xi, DD, DD_33_Sp  )
     real, intent(in) :: nu(:), E_r(:) 
     integer, intent(in) :: N_theta(:), N_zeta(:), N_xi(:) 
     real, optional, intent(out) :: DD( 3, 3, size(nu), size(E_r), size(N_theta), size(N_zeta), size(N_xi) )
     real, optional, intent(out) :: DD_33_Sp( size(nu), size(E_r), size(N_theta), size(N_zeta), size(N_xi) )
     
     real :: D( 3, 3, size(nu), size(E_r), size(N_theta), size(N_zeta), size(N_xi) )
     real :: D_33_Sp( size(nu), size(E_r), size(N_theta), size(N_zeta), size(N_xi) )
  
     integer :: i, j, ii, jj, kk, N_nu, N_E_r, M_theta, M_zeta, M_xi    
     character(len=500) :: file_path
     real :: t_clock, rate, t0, t1
     integer :: c0, c1, c_rate 
     
     ! Initialize compiler internal clock for computing wall-clock time
     call system_clock(count_rate=c_rate) ; rate = real(c_rate)
     
     ! Vectors sizes
     N_nu = size(nu) ; N_E_r = size(E_r)
     M_theta = size(N_theta) ; M_zeta = size(N_zeta) ; M_xi = size(N_xi)
     
     ! Location for output
     file_path = "monkes_Adjoint_Monoenergetic_Database.dat"     
     ! Open output file and write header
     open(21, file=trim(file_path))
     write(21,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_33_Spitzer ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     do j = 1, N_E_r ! Loop electric field value
        do i = 1, N_nu ! Loop collisionality value   
           do kk = 1, M_xi ! Loop number of Legendre modes
              do ii = 1, M_theta ! Loop number of theta points
                 do jj = 1, M_zeta ! Loop number of zeta points                  
                 
                    call system_clock(c0) ; call cpu_time(t0) 
                    call Solve_BTD_Adjoint_DKE_Legendre( N_theta(ii),          &
                                                 N_zeta(jj),           &
                                                 N_xi(kk),             &
                                                 nu(i), E_r(j),        &
                                                 D(:,:,i,j,ii,jj,kk),  &
                                                 D_33_Sp(i,j,ii,jj,kk) )                 
                    call system_clock(c1) ; call cpu_time(t1)
                    
                    ! Wall-clock time in seconds
                    t_clock = ( c1 - c0 ) / rate   
                                  
                    ! Writing results on "monkes_Monoenergetic_Database.plt"
                    write(21,'(9999e)') nu(i), E_r(j), &
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_33_Sp(i,j,ii,jj,kk), &
                                        t_clock, t1-t0 
                 end do  
              end do  
           end do  
        end do  
     end do
     close(21)
     
     ! Monoenergetic coefficients as output
     if ( present(DD) ) DD = D
     if ( present(DD_33_Sp) ) DD_33_Sp = D_33_Sp
  end subroutine 
   

end module
