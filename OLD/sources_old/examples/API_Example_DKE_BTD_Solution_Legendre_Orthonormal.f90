module API_Example_DKE_BTD_Solution_Legendre_Orthonormal

  use DKE_BTD_Solution_Legendre_Orthonormal
  
  implicit none
  
  private
  
  public :: Monoenergetic_Database_Example_1_Orthonormal
  public :: DKE_zeta_Convergence_Example_1_Orthonormal
  public :: Monoenergetic_Database_Input_Orthonormal
  
  
  
  contains  
   
 
  subroutine Monoenergetic_Database_Example_1_Orthonormal
    integer, parameter :: N_theta = 23, N_zeta = 55, N_xi = 140
    integer, parameter :: N_nu = 36, N_E_rho = 10
    real, parameter :: nu_min = 1d-5 , nu_max = 3e2!300 
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma(3,3,N_nu,N_E_rho), Gamma_33_Spitzer(N_nu,N_E_rho) 
    real :: nu_factor, E_rho_factor, D11_factor, D31_factor
    character(len=500) :: folder_path, file_path
    character(len=4) nu_char, E_rho_char, N_theta_char, N_zeta_char, N_xi_char
    integer :: i, j
     
    real :: t0, t1, rate 
    integer :: c0, c1, c_rate   
      
    ! *** Vector of collisionality
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu  
    
    ! *** Vector of radial electric field E_rho / v where rho = r/a.
    E_rho = [ 0d0, 1d-5, 3d-5, 1d-4, 3d-4, 1d-3, 3d-3, 1d-2, 3d-2, 1d-1 ]
    
    N_theta_char = int_2_str(N_theta)
    N_zeta_char = int_2_str(N_zeta)
    N_xi_char = int_2_str(N_xi)
    
    folder_path = "results/Monoenergetic_Database_Example_Orthonormal_"&
                // trim(N_theta_char) //"_" &
                // trim(N_zeta_char) //"_"// trim(N_xi_char)
    
    file_path = trim(folder_path) // "/Monoenergetic"
    call execute_command_line( "mkdir "//trim(folder_path) )       
    
    ! Initialize compiler internal clock 
    call system_clock(count_rate=c_rate) ; rate = real(c_rate)
     
    open(1221, file=trim(file_path) // ".plt" )
    write(1221,'(9999A25)') " nu/v ", " E_r/v ", &
                           " Gamma_11 ", " Gamma_31 ", &
                           " Gamma_13 ", " Gamma_33 ", &
                           " Gamma_33_Spitzer ", " t_clock [s] "
    do j = 1, 5, 4!6, 5!10!, 5 !N_E_rho
       if( j < 10 ) write(E_rho_char,'(I1)') j 
       if( 10 <= j .and. j < 100 ) write(E_rho_char,'(I2)') j 
       open(2221, file=trim(file_path)// "_" // trim(E_rho_char) // ".plt" )
       
       write(2221,'(9999A25)') " nu/v ", " E_r/v ", &
                              " Gamma_11 ", " Gamma_31 ", &
                              " Gamma_13 ", " Gamma_33 ", &
                              " Gamma_33_Spitzer ", " t_clock [s] "
       do i = 1, N_nu
          
          write(*,*) " nu, E_rho = ", nu(i), E_rho(j)
          call system_clock(c0)
          call Solve_BTD_DKE_Legendre_Orthonormal( N_theta, N_zeta, N_xi,          &
                                       nu(i), E_rho(j), Gamma(:,:,i,j), &
                                       Gamma_33_Spitzer(i,j) )
          call system_clock(c1)
          t0 = c0 / rate ; t1 = c1 / rate ! Instants of time in seconds
          write(*,'(9999e)') nu(i), E_rho(j), Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0
          write(1221,'(9999e)') nu(i), E_rho(j), Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0  
          write(2221,'(9999e)') nu(i), E_rho(j), Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0
          
       end do
       close(2221)
    end do
    close(1221)
     

  end subroutine
   
   
  
  
  ! *** Computes a monoenergetic coefficients database reading the input file
  ! "monkes_input.parameters"
  subroutine Monoenergetic_Database_Input_Orthonormal
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
     
     call Monoenergetic_Database_Scan( nu(1:N_nu), E_r(1:N_E_r),       &
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
  subroutine Monoenergetic_Database_Scan( nu, E_r, N_theta, N_zeta, N_xi, DD, DD_33_Sp  )
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
     file_path = "monkes_Monoenergetic_Database_Orthonormal.dat"     
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
                    call Solve_BTD_DKE_Legendre_Orthonormal( N_theta(ii),          &
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
  
   
  subroutine DKE_zeta_Convergence_Example_1_Orthonormal
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 28, M_theta = 9, M_xi = 1
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma(3,3,N_nu,N_E_rho, M_zeta), Gamma_33_Spitzer(N_nu,N_E_rho, M_zeta) 
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char, N_xi_char 
    character(len=2000) :: folder_path, file_path
    
    integer :: c0, c1, c_rate
    
    ! Initialize compiler internal clock 
    call system_clock(count_rate=c_rate)  
    
    ! *** Vector of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [1d-6, 3d-6, 1d-5, 3d-5, 1d-4]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! *** Vectors of theta, zeta and xi sizes
    N_theta = [(  11 + 4*i, i=1,M_theta )] 
    N_zeta  = [( 11 + 4*j, j=1,M_zeta  )] 
    N_xi    = 200!+[20, 40, 60, 80, 100] 
    
    ! *** Create folder for output
    N_xi_char = int_2_str(N_xi(1))  
    folder_path = "results/DKE_zeta_Convergence_Example_Orthonormal_Nxi_" // trim(N_xi_char) 
    call execute_command_line( "mkdir " // trim(folder_path) )    
        
    do i_E_rho = 5, 5! 
       do i_nu = 3, 3!N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Monoenergetic_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,'(9999A25)') " nu/v ", " E_r/v ", &
                                 " N_theta ", " N_zeta ", " N_xi ", &
                                 " Gamma_11 ", " Gamma_31 ", &
                                 " Gamma_13 ", " Gamma_33 ", &
                                 " Gamma_33_Spitzer ", " t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call system_clock(c0)
                   call Solve_BTD_DKE_Legendre_Orthonormal( N_theta(i), N_zeta(j), N_xi(k), &
                                                nu(i_nu), E_rho(i_E_rho), &
                                                Gamma(:,:,i_nu,i_E_rho,j), &
                                                Gamma_33_Spitzer(i_nu,i_E_rho,j) )
                   call system_clock(c1)
                   
                   t_clock(i_nu,i_E_rho,j) = (c1-c0)/real(c_rate)
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma(1,1,i_nu,i_E_rho,j), Gamma(3,1,i_nu,i_E_rho,j), &
                                       Gamma(1,3,i_nu,i_E_rho,j), Gamma(3,3,i_nu,i_E_rho,j), &
                                       Gamma_33_Spitzer(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma(1,1,i_nu,i_E_rho,j), Gamma(3,1,i_nu,i_E_rho,j), &
                                       Gamma(1,3,i_nu,i_E_rho,j), Gamma(3,3,i_nu,i_E_rho,j), &
                                       Gamma_33_Spitzer(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(221)
          
       end do
    end do
    
  end subroutine
  
 
  function int_2_str(int) result(str)
    integer, intent(in) :: int
    character(len=4) :: str
  
    if( int < 10 ) write(str, '(I1)') int
    if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
    if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
    if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int
  
  end function 




end module