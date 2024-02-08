module API_Example_DKE_BTD_Solution_Legendre

  use DKE_BTD_Solution_Legendre
  use Finite_differences
  implicit none
  
  private
  
  public :: Monoenergetic_Database_Input
  public :: Monoenergetic_Database_Input_NEW
  public :: Monoenergetic_Database_Example_1
  public :: DKE_zeta_Convergence_Example_1
  public :: Monoenergetic_convergence_xi
  
  
  
  
  public :: Test_Monoenergetic_lambda_function
  
  
  contains  
  
  
  subroutine Test_Monoenergetic_lambda_function
     integer, parameter :: N_theta=39, N_zeta= 129, N_xi = 180, M_xi = 60
     integer, parameter :: N_lambda = 300
     real, parameter :: nu = 1e-5, E_r=0*1e-3
     
     real :: lambda(0:N_lambda), lambda_c 
     real :: Gamma(3,3,0:N_lambda), Gamma_33_Spitzer 
     
     integer :: k
     
     call Monoenergetic_lambda_function( N_lambda, N_theta, N_zeta, &
                                         N_xi, M_xi, nu, E_r,       &
                                         lambda, lambda_c, &
                                         Gamma, Gamma_33_Spitzer )     
     
     
     open(21, file="Test_Monoenergetic_lambda_function.dat")
     write(21, '(9999A25)') " lambda ", "lambda_c", " Gamma(1,1,k) ", " Gamma(3,1,k) ", " Gamma(1,3,k) ", " Gamma(3,3,k) "
     do k = 0, N_lambda
        write(12345, *) lambda(k), lambda_c, Gamma(1,1,k), Gamma(3,1,k), Gamma(1,3,k), Gamma(3,3,k)    
     end do        
     close(21)
     
  
  end subroutine 
  
  
  ! *** Computes a monoenergetic coefficients database reading the input file
  ! "monkes_input.parameters"
  subroutine Monoenergetic_Database_Input
     integer, parameter :: N_nu_max = 500, N_E_r_max = 500
     integer, parameter :: M_theta_max = 500, M_zeta_max = 500, M_xi_max = 500
     real :: nu(N_nu_max), E_r(N_E_r_max)
     integer :: N_theta(M_theta_max), N_zeta(M_zeta_max) , N_xi(M_xi_max) 
     integer :: N_nu, N_E_r, M_theta, M_zeta, M_xi, ierr 
     namelist /parameters/N_theta, N_zeta, N_xi, nu, E_r  
     
     ! *** Read input parameters for monoenergetic database from
     ! "monkes_input.parameters" file.
     N_theta = -1 ; N_zeta = -1 ; N_xi = -1 ; nu = -1d14 ; E_r = -1d14 
     open(1, file= "monkes_input.parameters", status="old") 
     read(1, nml=parameters, iostat=ierr)     
     close(1)     
     
     ! Count number of collisionalities and radial electric field to be
     ! included in the database
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
  
  
  ! *** Computes a monoenergetic coefficients database reading the input file
  ! "monkes_input.parameters"
  subroutine Monoenergetic_Database_Input_NEW
     integer, parameter :: N_nu_max = 500, N_E_r_max = 500
     integer, parameter :: M_theta_max = 500, M_zeta_max = 500, M_xi_max = 500
     real :: nu(N_nu_max), E_r(N_E_r_max)
     integer :: N_theta(M_theta_max), N_zeta(M_zeta_max) , N_xi(M_xi_max) 
     integer :: N_nu, N_E_r, M_theta, M_zeta, M_xi, ierr 
     namelist /parameters/N_theta, N_zeta, N_xi, nu, E_r       
     namelist /parameters_DF/N_xi_DF, N_lambda       
     integer, parameter :: M_xi_DF_max = 500, M_lambda_max = 500 
     integer :: N_xi_DF(M_xi_DF_max), N_lambda(M_lambda_max)  
     integer :: M_xi_DF, M_lambda   
     logical :: Monoenergetic_lambda
     integer :: i, j, k, ii, jj, kk, iii, kkk, c0, c1, rate
     real :: t_clock, t_cpu, t0, t1, lambda_c
     real, allocatable :: F1(:,:,:), F3(:,:,:), Gamma_ij(:,:,:), lambda(:)
     real, allocatable :: dGamma_ij_dlambda(:,:,:) 
     real, allocatable :: D(:,:,:,:,:,:,:), D_33_Sp(:,:,:,:,:)  
     character(len=500) :: file_path
     
     call system_clock(count_rate=rate) ! Get click rate
     
     ! *** Read input parameters for monoenergetic database from
     ! "monkes_input.parameters" file.
     N_theta = -1 ; N_zeta = -1 ; N_xi = -1 ; nu = -1d14 ; E_r = -1d14 
     open(1, file= "monkes_input.parameters", status="old") 
     read(1, nml=parameters, iostat=ierr)     
     close(1)     
     
     ! Count number of collisionalities and radial electric field to be
     ! included in the database
     M_theta = count(N_theta > 0)  
     M_zeta = count(N_zeta > 0)  
     M_xi = count(N_xi > 1)  
     N_nu = count(nu > 0) ; N_E_r = count(E_r /= -1d14)
     
     allocate( D( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
     allocate( D_33_Sp(  N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
          
     ! *** Read input parameters for distribution function and lambda from
     ! dependence of the monoenergetic coefficients from "monkes_input.parameters" file.
     N_xi_DF = -1 ; N_lambda = -1
     open(1, file= "monkes_input.parameters", status="old") 
     read(1, nml=parameters_DF, iostat=ierr)     
     close(1)
     Monoenergetic_lambda = (ierr == 0)
     
     if ( Monoenergetic_lambda ) then
       M_xi_DF  = count(N_xi_DF > 0)  
       M_lambda = count(N_lambda > 0)
     else
       M_xi_DF = 1 ; M_lambda = 1 
       N_xi_DF = 2 
     end if
     
     write(*,*) ierr, Monoenergetic_lambda
     write(*,*) N_theta(1:M_theta), N_zeta(1:M_zeta), N_xi(1:M_xi), nu(1:N_nu), E_r(1:N_E_r) 
     write(*,*) N_xi_DF(1:M_xi_DF), N_lambda(1:M_lambda)
     
     
     where( mod(N_theta(1:M_theta),2) == 0 )  N_theta(1:M_theta) = N_theta(1:M_theta) + 1
     where( mod(N_zeta(1:M_zeta),2) == 0 )    N_zeta(1:M_zeta) = N_zeta(1:M_zeta) + 1
 
     ! Location for output
     file_path = "monkes_Monoenergetic_Database.dat"     
     ! Open output file and write header
     open(21, file=trim(file_path))
     write(21,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_33_Spitzer ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     ! OPEN (if necessary) monkes_Monoenergetic_lambda.dat
     if( Monoenergetic_lambda ) then     
     open(31, file=trim("monkes_Monoenergetic_lambda.dat"))
     write(31,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_33_Spitzer ",   &
                             " M_xi ", " lambda ", " lambda_c ", & 
                             " d_11 ", " d_31 ", &
                             " d_13 ", " d_33 ", &
                             " d d_11 / dlambda ", " d d_31 / dlambda ", &
                             " d d_13 / dlambda ", " d d_33 / dlambda  "
     endif     
     do j = 1, N_E_r ! Loop electric field value
        do i = 1, N_nu ! Loop collisionality value   
           do kk = 1, M_xi ! Loop number of Legendre modes
              do ii = 1, M_theta ! Loop number of theta points
                 do jj = 1, M_zeta ! Loop number of zeta points                  
                    do k = 1, M_xi_DF ! Loop on number of modes extracted of the distribution function
                       do iii = 1, M_lambda

                          call system_clock(c0) ; call cpu_time(t0)
                          allocate( F1(0:N_theta(ii)-1, 0:N_zeta(jj)-1, 0:N_xi_DF(k) ) )
                          allocate( F3(0:N_theta(ii)-1, 0:N_zeta(jj)-1, 0:N_xi_DF(k) ) )

                          ! Solve the DKE and extract N_xi_DF(k)+1 Legendre modes
                          call Solve_BTD_DKE_Legendre_DF( N_theta(ii),          &
                                                       N_zeta(jj),           &
                                                       N_xi(kk),             &
                                                       nu(i), E_r(j),        &
                                                       D(:,:,i,j,ii,jj,kk),  &
                                                       D_33_Sp(i,j,ii,jj,kk), &
                                                       N_xi_DF(k), F1, F3 )
                          call system_clock(c1) ; call cpu_time(t1)                          
                          ! Wall-clock time in seconds
                          t_clock = ( c1 - c0 ) / rate
                          
                          ! Writing results on "monkes_Monoenergetic_Database.plt"
                          write(21,*) nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_33_Sp(i,j,ii,jj,kk), &
                                        t_clock, t1-t0 


                         if( Monoenergetic_lambda ) then

                            allocate( lambda(0:N_lambda(iii)), Gamma_ij(3,3,0:N_lambda(iii)) )
                            allocate( dGamma_ij_dlambda(3,3,0:N_lambda(iii)) )
                            ! Call the routine that computs Dij as lambda functions
                            call Monoenergetic_lambda_function_NEW( N_lambda(iii), &
                                                                    N_theta(ii), N_zeta(jj), &
                                                                    N_xi(kk), N_xi_DF(k), & 
                                                                    F1, F3,            &
                                                                    lambda, lambda_c, &
                                                                    Gamma_ij )


                            call Grid_initialization( "lambda", lambda, 2 )                            
                            
                            dGamma_ij_dlambda(1,1,:) = Derivative( "lambda", Gamma_ij(1,1,:), 1 )
                            dGamma_ij_dlambda(3,1,:) = Derivative( "lambda", Gamma_ij(3,1,:), 1 )
                            dGamma_ij_dlambda(1,3,:) = Derivative( "lambda", Gamma_ij(1,3,:), 1 )
                            dGamma_ij_dlambda(3,3,:) = Derivative( "lambda", Gamma_ij(3,3,:), 1 )
                            
                            ! Write Dij(lambda) in monkes_Monoenergetic_lambda.dat
                            do kkk = 0, N_lambda(iii)
                             
                               write(31,*) nu(i), E_r(j), &   
                                                   real(N_theta(ii)), &
                                                   real(N_zeta(jj)), &
                                                   real(N_xi(kk)), &
                                                   D(1,1,i,j,ii,jj,kk), &
                                                   D(3,1,i,j,ii,jj,kk), &
                                                   D(1,3,i,j,ii,jj,kk), &
                                                   D(3,3,i,j,ii,jj,kk), &
                                                   D_33_Sp(i,j,ii,jj,kk), &
                                                   real(N_xi_DF(k)),      & 
                                                   lambda(kkk), lambda_c, & 
                                                   Gamma_ij(1,1,kkk), &
                                                   Gamma_ij(3,1,kkk), &
                                                   Gamma_ij(1,3,kkk), &
                                                   Gamma_ij(3,3,kkk), & 
                                                   dGamma_ij_dlambda(1,1,kkk), &
                                                   dGamma_ij_dlambda(3,1,kkk), &
                                                   dGamma_ij_dlambda(1,3,kkk), &
                                                   dGamma_ij_dlambda(3,3,kkk)
                            end do
                            
                            deallocate( lambda, Gamma_ij, dGamma_ij_dlambda )
                          end if                       
                          deallocate( F1, F3 )
                            
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do    
     close(21) ! close monkes_Monoenergetic_database.dat
     ! CLOSE (if necessary) monkes_Monoenergetic_lambda.dat 
     if( Monoenergetic_lambda ) close(31)
     
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
     file_path = "monkes_Monoenergetic_Database.dat"     
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
                    call Solve_BTD_DKE_Legendre( N_theta(ii),          &
                                                 N_zeta(jj),           &
                                                 N_xi(kk),             &
                                                 nu(i), E_r(j),        &
                                                 D(:,:,i,j,ii,jj,kk),  &
                                                 D_33_Sp(i,j,ii,jj,kk) )                 
                    call system_clock(c1) ; call cpu_time(t1)
                    
                    ! Wall-clock time in seconds
                    t_clock = ( c1 - c0 ) / rate   
                                  
                    ! Writing results on "monkes_Monoenergetic_Database.plt"
                    write(21,*) nu(i), E_r(j), &
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
   
  ! *** Computes the monoenergetic coefficients reading the input file
  ! "input.monoenergetic"
  subroutine Monoenergetic_Database_Input_OLDDD
     integer, parameter :: N_nu_max = 500, N_E_r_max = 500
     integer :: N_theta, N_zeta, N_xi, N_nu, N_E_r, i, j, ierr 
     real :: nu(N_nu_max), E_r(N_E_r_max)
     real, allocatable :: Gamma(:,:,:,:), Gamma_33_Spitzer(:,:) 
     namelist /parameters/N_theta, N_zeta, N_xi, nu, E_r     
     character(len=500) :: folder_path, file_path
     character(len=4) E_r_char, N_theta_char, N_zeta_char, N_xi_char
     real :: t0, t1, rate 
     integer :: c0, c1, c_rate  
     
     ! *** Read input parameters from "monkes.parameters" file
     nu = -1 ; E_r = -1d14 
     open(1, file= "input_monkes.parameters", status="old") ; read(1, nml=parameters, iostat=ierr)     
     close(1)  
     
     ! Count number of collisionalities and radial electric field to be
     ! includede in the database
     N_nu = count(nu > 0) ; N_E_r = count(E_r /= -1d14)      
     allocate( Gamma(3,3,N_nu,N_E_r), Gamma_33_Spitzer(N_nu,N_E_r) )
     
     write(*,*) " N_theta ",       N_theta 
     write(*,*) " N_zeta ",        N_zeta 
     write(*,*) " N_xi ",          N_xi 
     write(*,*) " N_nu ",          N_nu 
     write(*,*) " N_E_r ",         N_E_r 
     write(*,*) " nu(1:N_nu) ",    nu(1:N_nu) 
     write(*,*) " E_r(1:N_E_r) ",  E_r(1:N_E_r) 
     
     ! Initialize compiler internal clock for computing wall-clock time
     call system_clock(count_rate=c_rate) ; rate = real(c_rate)
     
     ! Naming the output files
     N_theta_char = int_2_str(N_theta)
     N_zeta_char = int_2_str(N_zeta)
     N_xi_char = int_2_str(N_xi)
      
     folder_path = "results/Monoenergetic_Database_"&
                // trim(N_theta_char) //"_" &
                // trim(N_zeta_char) //"_"// trim(N_xi_char)
    
     file_path = trim(folder_path) // "/Monoenergetic"
     call execute_command_line( "mkdir "//trim(folder_path) )    
     
     ! Writing header on "Monoenergetic.plt"
     open(21, file=trim(file_path) // ".plt" )
     write(21,'(9999A25)') " nu/v ", " E_r/v ", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " Gamma_11 ", " Gamma_31 ", &
                             " Gamma_13 ", " Gamma_33 ", &
                             " Gamma_33_Spitzer ", " t_clock [s] "
     do j = 1, N_E_r
     
        E_r_char = int_2_str(j)
        open(31, file=trim(file_path) //"_"// trim(E_r_char)//".plt" )          
        ! Writing header on "Monoenergetic_j.plt"
        write(31,'(9999A25)') " nu/v ", " E_r/v ", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " Gamma_11 ", " Gamma_31 ", &
                             " Gamma_13 ", " Gamma_33 ", &
                             " Gamma_33_Spitzer ", " t_clock [s] "
        do i = 1, N_nu
           call system_clock(c0)
           call Solve_BTD_DKE_Legendre( N_theta, N_zeta, N_xi,          &
                                        nu(i), E_r(j), Gamma(:,:,i,j), &
                                        Gamma_33_Spitzer(i,j) )
           call system_clock(c1)
           t0 = c0 / rate ; t1 = c1 / rate ! Instants of time in seconds
                                    
           ! Writing results on "Monoenergetic.plt"
           write(21,*) nu(i), E_r(j), &
                              real(N_theta), real(N_zeta), real(N_xi), &
                              Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), &
                              Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0                              
                       
           ! Writing results on "Monoenergetic_j.plt"       
           write(31,*) nu(i), E_r(j), &
                              real(N_theta), real(N_zeta), real(N_xi), &
                              Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), &
                              Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0
           
           ! Writing on output file
           write(*,'(9999A25)') " nu/v ", " E_r/v ", &
                                " N_theta ", " N_zeta ", " N_xi ", &
                                " Gamma_11 ", " Gamma_31 ", &
                                " Gamma_13 ", " Gamma_33 ", &
                                " Gamma_33_Spitzer ", " t_clock [s] "
           write(*,*) nu(i), E_r(j), &
                              real(N_theta), real(N_zeta), real(N_xi), &
                              Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), &
                              Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0
        end do
        close(31) ! Close "Monoenergetic_j.plt"
        write(21,*) 
     end do  
     close(21)    ! Close "Monoenergetic.plt"
     
  end subroutine

  subroutine Monoenergetic_convergence_xi
    integer, parameter :: N_nu = 18, N_E_rho = 10 
    integer, parameter :: M_zeta = 1, M_theta = 1, M_xi = 20
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma(3,3,N_nu,N_E_rho, M_zeta), Gamma_33_Spitzer(N_nu,N_E_rho, M_zeta) 
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char, N_theta_char, N_zeta_char 
    character(len=2000) :: folder_path, file_path
    
    integer :: c0, c1, c_rate
    
    ! Initialize compiler internal clock 
    call system_clock(count_rate=c_rate)  
    
    ! Vectors of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! Vectors of theta, zeta and xi sizes
    N_theta = [(  29 + 4*i, i=1,M_theta )] 
    N_zeta  = 38+[( 11 + 4*j, j=1,M_zeta  )] 
    N_xi    = [( 10 + 10*j, j=1,M_xi  )] 
    
    ! Create folder for output
    N_theta_char = int_2_str(N_theta(1))
    N_zeta_char = int_2_str(N_zeta(1))  
    folder_path = "results/Monoenergetic_convergence_Legendre"
    call execute_command_line( "mkdir " // trim(folder_path) )    
        
    do i_E_rho = 1, 1! 
       do i_nu = 3, 3!N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of electric field an output file
          file_path = trim(folder_path) // "/Monoenergetic_convergence_xi_N_theta_" // &
                      trim(N_theta_char) // "_N_zeta_" // trim(N_zeta_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,'(9999A25)') " nu/v ", " E_r/v ", &
                                 " N_theta ", " N_zeta ", " N_xi ", &
                                 " Gamma_11 ", " Gamma_31 ", &
                                 " Gamma_13 ", " Gamma_33 ", &
                                 " Gamma_33_Spitzer ", " t_clock [s] "
          do j = 1, M_zeta 
             do i = 1, M_theta 
                
                ! Write xi convergence plot.     
                do k = 1, M_xi
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call system_clock(c0)
                   call Solve_BTD_DKE_Legendre( N_theta(i), N_zeta(j), N_xi(k), &
                                                nu(i_nu), E_rho(i_E_rho), &
                                                Gamma(:,:,i_nu,i_E_rho,j), &
                                                Gamma_33_Spitzer(i_nu,i_E_rho,j) )
                   call system_clock(c1)
                   
                   t_clock(i_nu,i_E_rho,j) = (c1-c0)/real(c_rate)
                   write(*,*) nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma(1,1,i_nu,i_E_rho,j), Gamma(3,1,i_nu,i_E_rho,j), &
                                       Gamma(1,3,i_nu,i_E_rho,j), Gamma(3,3,i_nu,i_E_rho,j), &
                                       Gamma_33_Spitzer(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,*) nu(i_nu), E_rho(i_E_rho), &
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
  ! *** Example on how to compute a monoenergetic database
  ! using MONKES. 
  subroutine Monoenergetic_Database_Example_1   
    integer, parameter :: N_theta = 15, N_zeta = 127, N_xi = 180
    integer, parameter :: N_nu = 36, N_E_rho = 10
    real, parameter :: nu_min = 1d-5 , nu_max = 3e2!300 
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma(3,3,N_nu,N_E_rho), Gamma_33_Spitzer(N_nu,N_E_rho) 
    !real :: nu_factor, E_rho_factor, D11_factor, D31_factor
    character(len=500) :: folder_path, file_path
    character(len=4) E_rho_char, N_theta_char, N_zeta_char, N_xi_char
    integer :: i, j
     
    real :: t0, t1, rate 
    integer :: c0, c1, c_rate   
      
    ! *** Vector of collisionality
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu  
    
    ! Vector of radial electric field E_rho / v where rho = r/a.
    E_rho = [ 0d0, 1d-5, 3d-5, 1d-4, 3d-4, 1d-3, 3d-3, 1d-2, 3d-2, 1d-1 ]

    ! Naming the output files
    N_theta_char = int_2_str(N_theta)
    N_zeta_char = int_2_str(N_zeta)
    N_xi_char = int_2_str(N_xi)
    
    folder_path = "results/Monoenergetic_Database_Example_"&
                // trim(N_theta_char) //"_" &
                // trim(N_zeta_char) //"_"// trim(N_xi_char)
    
    file_path = trim(folder_path) // "/Monoenergetic"
    call execute_command_line( "mkdir "//trim(folder_path) )       
    
    ! Initialize compiler internal clock for computing wall-clock time
    call system_clock(count_rate=c_rate) ; rate = real(c_rate)
     
    open(1221, file=trim(file_path) // ".plt" )
    write(1221,'(9999A25)') " nu/v ", " E_r/v ", &
                           " Gamma_11 ", " Gamma_31 ", &
                           " Gamma_13 ", " Gamma_33 ", &
                           " Gamma_33_Spitzer ", " t_clock [s] ", &
                           " N_theta ", " N_zeta ", " N_xi "
    do j = 1, 6, 5!10!, 5 !N_E_rho
       if( j < 10 ) write(E_rho_char,'(I1)') j 
       if( 10 <= j .and. j < 100 ) write(E_rho_char,'(I2)') j 
       open(2221, file=trim(file_path)// "_" // trim(E_rho_char) // ".plt" )
       
       write(2221,'(9999A25)') " nu/v ", " E_r/v ", &
                              " Gamma_11 ", " Gamma_31 ", &
                              " Gamma_13 ", " Gamma_33 ", &
                              " Gamma_33_Spitzer ", " t_clock [s] ", &
                              " N_theta ", " N_zeta ", " N_xi "
       do i = 1, N_nu
          
          write(*,*) " nu, E_rho = ", nu(i), E_rho(j)
          call system_clock(c0)
          call Solve_BTD_DKE_Legendre( N_theta, N_zeta, N_xi,          &
                                       nu(i), E_rho(j), Gamma(:,:,i,j), &
                                       Gamma_33_Spitzer(i,j) )
          call system_clock(c1)
          t0 = c0 / rate ; t1 = c1 / rate ! Instants of time in seconds
          write(*,*) nu(i), E_rho(j), Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0, &
                            real(N_theta) ,  real(N_zeta) ,  real(N_xi)  
          write(1221,*) nu(i), E_rho(j), Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0, &
                            real(N_theta) ,  real(N_zeta) ,  real(N_xi)  
          write(2221,*) nu(i), E_rho(j), Gamma(1,1,i,j), Gamma(3,1,i,j), Gamma(1,3,i,j), Gamma(3,3,i,j), Gamma_33_Spitzer(i,j), t1-t0, &
                            real(N_theta) ,  real(N_zeta) ,  real(N_xi)  
          
       end do
       close(2221)
    end do
    close(1221)
     
  end subroutine

  subroutine DKE_zeta_Convergence_Example_1
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 40, M_theta = 6, M_xi = 1
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
    
    ! Vectors of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [ 1d-6, 3d-6, 1d-5, 3d-5, 1d-4 ]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! Vectors of theta, zeta and xi sizes
    N_theta = 11 + [( 4*i, i=1,M_theta )] 
    N_zeta  = 11 + [( 4*j, j=1,M_zeta  )] 
    N_xi    = 200 !+ [20, 40, 60, 80, 100] 
    
    ! Create folder for output
    N_xi_char = int_2_str(N_xi(1))  
    folder_path = "results/DKE_zeta_Convergence_Example_Nxi_" // trim(N_xi_char) 
    call execute_command_line( "mkdir " // trim(folder_path) )    
        
    do i_E_rho = 5, 5!5, 4! 
       do i_nu = 3, 3!N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of electric field an output file
          file_path = trim(folder_path) // "/Monoenergetic_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,'(9999A25)') " nu/v ", " E_r/v ", &
                                 " N_theta ", " N_zeta ", " N_xi ", &
                                 " D_11 ", " D_31 ", &
                                 " D_13 ", " D_33 ", &
                                 " D_33_Spitzer ", " t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call system_clock(c0)
                   call Solve_BTD_DKE_Legendre( N_theta(i), N_zeta(j), N_xi(k), &
                                                nu(i_nu), E_rho(i_E_rho), &
                                                Gamma(:,:,i_nu,i_E_rho,j), &
                                                Gamma_33_Spitzer(i_nu,i_E_rho,j) )
                   call system_clock(c1)
                   
                   t_clock(i_nu,i_E_rho,j) = (c1-c0)/real(c_rate)
                   write(*,*) nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma(1,1,i_nu,i_E_rho,j), Gamma(3,1,i_nu,i_E_rho,j), &
                                       Gamma(1,3,i_nu,i_E_rho,j), Gamma(3,3,i_nu,i_E_rho,j), &
                                       Gamma_33_Spitzer(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,*) nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma(1,1,i_nu,i_E_rho,j), Gamma(3,1,i_nu,i_E_rho,j), &
                                       Gamma(1,3,i_nu,i_E_rho,j), Gamma(3,3,i_nu,i_E_rho,j), &
                                       Gamma_33_Spitzer(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                write(221,*)
             end do
             write(221,*)
          end do           
          close(221)
          
       end do
    end do
    
  end subroutine
  
  ! Takes a positive integer "int" between 0 and 9999 and stores
  ! its value in a string "str". 
  function int_2_str(int) result(str)
    integer, intent(in) :: int
    character(len=4) :: str
  
    if( int < 10 ) write(str, '(I1)') int
    if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
    if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
    if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int
  
  end function 




end module
