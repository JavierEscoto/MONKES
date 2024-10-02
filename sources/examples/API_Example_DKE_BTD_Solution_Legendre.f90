module API_Example_DKE_BTD_Solution_Legendre

  use DKE_BTD_Solution_Legendre
  use Finite_differences
  use Barycentric_grids
  
  implicit none
  
  private
  
  public :: Monoenergetic_Database_Input
  public :: Monoenergetic_Database_Maxwell_points
  
  
  ! TO BE ERASED:
  
  public :: Test_Maxwell_points
  
  contains  
  
  ! *** Computes a monoenergetic coefficients database reading the input file
  ! "monkes_input.parameters"
  subroutine Monoenergetic_Database_Input
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
     logical :: Monoenergetic_lambda, Monoenergetic_theta_zeta = .false.
     integer :: i, j, k, ii, jj, kk, iii, kkk, c0, c1, rate, i_theta, i_zeta
     real :: t_clock, t_cpu, t0, t1, lambda_c
     real, allocatable :: F1(:,:,:), F3(:,:,:), Gamma_ij(:,:,:), lambda(:)
     real, allocatable :: dGamma_ij_dlambda(:,:,:) 
     real, allocatable :: D(:,:,:,:,:,:,:), D_33_Sp(:,:,:,:,:)  
     real, allocatable :: Dij_FS(:,:,:,:)
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
      
     
     write(*,*) " *** Performing scan in collisionality and radial electric field "
     write(*,*) " nu/v [m^-1] "
     write(*,*) nu(1:N_nu)
     write(*,*)
     write(*,*) " E_r/v [kV s /m^2] "
     write(*,*) E_r(1:N_E_r) 
     write(*,*)
     write(*,*) " *** Scan done using the resolutions "
     write(*,*) " N_theta (# points for poloidal angle) " 
     write(*,*)   N_theta(1:M_theta)   
     write(*,*) " N_zeta (# points for toroidal angle) " 
     write(*,*)   N_zeta(1:M_zeta)   
     write(*,*) " N_xi (# Legendre modes in pitch-angle cosine) " 
     write(*,*)   N_xi(1:M_xi)  
     write(*,*)
     if( Monoenergetic_lambda ) then
       write(*,*) " *** Extracting dependence of the monoenergetic coefficients "
       write(*,*) " on lambda for each case of the scan"     
       write(*,*) " N_xi_DF " 
       write(*,*) N_xi_DF(1:M_xi_DF) 
       write(*,*) " N_lambda " 
       write(*,*) N_lambda(1:M_lambda)
       write(*,*)
     end if   
     
     where( mod(N_theta(1:M_theta),2) == 0 )  N_theta(1:M_theta) = N_theta(1:M_theta) + 1
     where( mod(N_zeta(1:M_zeta),2) == 0 )    N_zeta(1:M_zeta) = N_zeta(1:M_zeta) + 1     
     
     write(*,*) " *** Monoenergetic Database " 
     write(*,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_33_Spitzer ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     
     
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
     if( Monoenergetic_theta_zeta ) then     
       open(41, file=trim("monkes_Monoenergetic_theta_zeta.dat"))
       write(41,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                               " N_theta ", " N_zeta ", " N_xi ", &
                               " D_11 ", " D_31 ", &
                               " D_13 ", " D_33 ", &
                               " D_33_Spitzer ",   &
                               " theta ", " zeta ", & 
                               " d_11 ", " d_31 ", &
                               " d_13 ", " d_33 " 
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
                          t_clock = ( c1 - c0 )*1d0 / rate
                          
                          ! Writing results on terminal
                          write(*,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_33_Sp(i,j,ii,jj,kk), &
                                        t_clock, t1-t0 
                          
                          ! Writing results on "monkes_Monoenergetic_Database.dat"
                          write(21,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_33_Sp(i,j,ii,jj,kk), &
                                        t_clock, t1-t0 
                          flush(21)

                          ! Call the routine that computes Dij as (theta,zeta) functions                             
                          if( Monoenergetic_theta_zeta ) then     
                            
                            if( allocated(Dij_FS) ) deallocate(Dij_FS) 
                            allocate( Dij_FS(0:N_theta(ii),0:N_zeta(jj),3,3) )                        
                            call Monenergetic_theta_zeta_function( N_theta(ii), N_zeta(jj), F1(:,:,0:2), F3(:,:,0:2), Dij_FS )
                            
                            do i_zeta = 0, N_zeta(jj)
                               do i_theta = 0, N_theta(ii)
                                  write(41,'(9999e25.16)') nu(i), E_r(j),         &   
                                                    real(N_theta(ii)),            &
                                                    real(N_zeta(jj)),             &
                                                    real(N_xi(kk)),               &
                                                    D(1,1,i,j,ii,jj,kk),          &
                                                    D(3,1,i,j,ii,jj,kk),          &
                                                    D(1,3,i,j,ii,jj,kk),          &
                                                    D(3,3,i,j,ii,jj,kk),          &
                                                    D_33_Sp(i,j,ii,jj,kk),        &
                                                    theta(i_theta),               & 
                                                    zeta(i_zeta),                 & 
                                                    Dij_FS(i_theta, i_zeta, 1,1), &
                                                    Dij_FS(i_theta, i_zeta, 3,1), &
                                                    Dij_FS(i_theta, i_zeta, 1,3), &
                                                    Dij_FS(i_theta, i_zeta, 3,3) 
                                  flush(41)
                               end do
                            end do 
                          endif  
                          
                               
                          if( Monoenergetic_lambda ) then
					   
                             allocate( lambda(0:N_lambda(iii)), Gamma_ij(3,3,0:N_lambda(iii)) )
                             allocate( dGamma_ij_dlambda(3,3,0:N_lambda(iii)) )
                             ! Call the routine that computes Dij as lambda functions
                             call Monoenergetic_lambda_function( N_lambda(iii), &
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
                              
                                write(31,'(9999e25.16)') nu(i), E_r(j), &   
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
                                flush(31)
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
     ! CLOSE (if necessary) monkes_Monoenergetic_theta_zeta.dat 
     if( Monoenergetic_theta_zeta )  close(41)
     
  end subroutine
  
  
  ! *** Compute a mononergetic database for nu(v)/v correspondent to
  ! Maxwell points.
  subroutine Monoenergetic_Database_Maxwell_points
     integer, parameter :: Nx = 100 
     real, parameter :: pi = acos(-1d0) 
     real :: x(0:Nx), w(0:Nx), l(0:Nx+1), d(0:Nx) 
     integer :: k 
     
     
     call Maxwell_polynomials( Nx, l, d, x, w )
     
     write(*,*) " l(k), d(k) "
     do k = 0, Nx
        write(*,*) l(k), d(k) 
     end do
     write(*,*) " x(k), w(k) "
     do k = 0, Nx
        write(*,*) x(k), w(k)
     end do
     write(*,*)
     
     do k = 0, Nx
         
        call Maxwell_polynomials( k, l(0:k), d(0:k), x(0:k), w(0:k) )  
        
        write(*,*) k , dot_product( sin(x(0:k)**2), w(0:k) ), 0.5* sqrt( pi /sqrt(2d0) ) * sin(pi/8)
        write(*,*) k , dot_product( sin(x(0:k))+2*x(0:k)*cos(x(0:k)), w(0:k) ), 1d0
        write(*,*) k , dot_product( 1d0/(1+x(0:k))-2*x(0:k)*log(1+x(0:k)), w(0:k) ), 0d0
        write(*,*) k , dot_product( 0.5/sqrt( 1d-2 + x(0:k) )-2*x(0:k)*sqrt(1d-2+ x(0:k) ), w(0:k) ), -sqrt(1d-2)
         
     end do
  
  end subroutine 


  
  ! *** Test to check up until which Nx, using Maxwell polynomials works. 
  ! For smooth functions, convergence is excellent.
  subroutine Test_Maxwell_points
     integer, parameter :: Nx = 200 
     real, parameter :: pi = acos(-1d0) 
     real :: x(0:Nx), w(0:Nx), l(0:Nx+1), d(0:Nx) 
     integer :: k 
     
     
     !call Maxwell_polynomials( Nx, l, d, x, w )
     !
     !write(*,*) " l(k), d(k) "
     !do k = 0, Nx
     !   write(*,*) l(k), d(k) 
     !end do
     !write(*,*) " x(k), w(k) "
     !do k = 0, Nx
     !   write(*,*) x(k), w(k)
     !end do
     !write(*,*)
     
     do k = 0, Nx
         
        call Maxwell_polynomials( k, l(0:k), d(0:k), x(0:k), w(0:k) )  
        
        write(*,*) k , dot_product( sin(x(0:k)**2), w(0:k) ), 0.5* sqrt( pi /sqrt(2d0) ) * sin(pi/8)
        write(*,*) k , dot_product( sin(x(0:k))+2*x(0:k)*cos(x(0:k)), w(0:k) ), 1d0
        write(*,*) k , dot_product( 1d0/(1+x(0:k))-2*x(0:k)*log(1+x(0:k)), w(0:k) ), 0d0
        write(*,*) k , dot_product( 0.5/sqrt( 1d-2 + x(0:k) )-2*x(0:k)*sqrt(1d-2+ x(0:k) ), w(0:k) ), -sqrt(1d-2)
         
     end do
  
  end subroutine 


end module
