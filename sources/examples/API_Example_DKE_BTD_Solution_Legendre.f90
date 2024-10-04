module API_Example_DKE_BTD_Solution_Legendre

  use DKE_BTD_Solution_Legendre
  use Finite_differences
  use Barycentric_grids
  use Magnetic_configuration
  
  implicit none
  
  private
  
  public :: Monoenergetic_Database_Input 
  public :: Neoclassical_fluxes_Lorentz

  
  
  
  ! TO BE ERASED:  
  public :: Test_Maxwell_points
  public :: Test_Adjoint_derivatives
  
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
           do kk = 1, M_xi ! Loop number of Legendre modes used for solving DKE
              do ii = 1, M_theta ! Loop number of theta points used for solving DKE
                 do jj = 1, M_zeta ! Loop number of zeta points used for solving DKE                  
                    do k = 1, M_xi_DF ! Loop on number of modes extracted of the distribution function
                       do iii = 1, M_lambda

                          call system_clock(c0) ; call cpu_time(t0)
                          allocate( F1(0:N_theta(ii)-1, 0:N_zeta(jj)-1, 0:N_xi_DF(k) ) )
                          allocate( F3(0:N_theta(ii)-1, 0:N_zeta(jj)-1, 0:N_xi_DF(k) ) )

                          ! Solve the DKE and extract N_xi_DF(k)+1 Legendre modes
                          ! TO DO? call Solve_BTD_DKE_Legendre_DF_Matrix_free( N_theta(ii),          &
                          call Solve_BTD_DKE_Legendre_DF( N_theta(ii),           & 
                                                          N_zeta(jj),            &
                                                          N_xi(kk),              &
                                                          nu(i), E_r(j),         &
                                                          D(:,:,i,j,ii,jj,kk),   &
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
  
   
  
  ! *** Test to check up until which Nx, using Maxwell polynomials works. 
  ! For smooth functions, convergence is excellent.
  subroutine Test_Maxwell_points
     integer, parameter :: Nx = 200 
     real, parameter :: pi = acos(-1d0), a = 1
     real :: x(0:Nx), w(0:Nx), l(0:Nx+1), d(0:Nx)
     integer :: k, kk 
     
     
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
     open(1211, file="Test_Maxwell_points.dat")
     do k = 0, Nx
         
        call Maxwell_polynomials( k, l(0:k), d(0:k), x(0:k), w(0:k) )  
        
        write(1211,*) " Nx ", k
        write(1211,'(9999A25)') " k ", "l", "d", "x", "w"
        do kk = 0, k
           write(1211,'(9999e25.16)') 1d0*kk , l(kk), d(kk), x(kk), w(kk)       
        end do  
        
        
        write(1211,'(9999A25)') " k ", "sin(x**2)exp(-x**2)", " Exact solution "
        write(1211,'(9999e25.16)') 1d0*k , dot_product( sin(x(0:k)**2), w(0:k) ), 0.5* sqrt( pi /sqrt(2d0) ) * sin(pi/8)
        write(1211,'(9999A25)') " k ", "sin(x)+2xcos(x)", " Exact solution "
        write(1211,'(9999e25.16)') 1d0*k , dot_product( sin(x(0:k))+2*x(0:k)*cos(x(0:k)), w(0:k) ), 1d0
        write(1211,'(9999A25)') " k ", "1/(1+x)-2xlog(1+x)", " Exact solution "
        write(1211,'(9999e25.16)') 1d0*k , dot_product( 1d0/(1+x(0:k))-2*x(0:k)*log(1+x(0:k)), w(0:k) ), 0d0
        write(1211,'(9999A25)') " k ", "0.5/sqrt(0.01+x)-2xsqrt(0.01+x)", " Exact solution "
        write(1211,'(9999e25.16)') 1d0*k , dot_product( 0.5/sqrt( 1d-2 + x(0:k) )-2*x(0:k)*sqrt(1d-2+ x(0:k) ), w(0:k) ), -sqrt(1d-2)
        write(1211,'(9999A25)') " k ", "2axcos(ax^2)-2x(2+sin(ax^2))", " Exact solution "
        write(1211,'(9999e25.16)') 1d0*k , dot_product( 2*a*x(0:k)*cos(a*x(0:k)**2) - 2*x(0:k)*(2+sin(a*x(0:k)**2)), w(0:k) ), -2d0
        write(1211,*)
     end do
     close(1211)
  end subroutine 
  
  ! Reads the input monkes_flux.input to compute the neoclassical fluxes
  ! in the Lorentz limit
  subroutine Neoclassical_fluxes_Lorentz !WORK IN PROGRESS
     real, parameter :: pi = acos(-1d0), e = 1.6022d-19, c = 2.9979246d8
     real, parameter :: eps_0 = 1d7/(4*pi*c**2), c_units = 1d19*pi*e**2/(4*pi*eps_0)**2  
     integer, parameter :: N_Er_max = 500, N_species_max = 40
     namelist /parameters/N_theta, N_zeta, N_xi, N_x, E_r, Z, m, n, T, n_prime, T_prime   
     
     ! Density n [10^19 m^3], Atomic number Z [AU], Temperature T [eV], Mass m [u]
     real :: Z(N_species_max), m(N_species_max), n(N_species_max), T(N_species_max) 
     real :: n_prime(N_species_max), T_prime(N_species_max)     ! dn/dr and dT/dr (r is in metres)                       
     real :: E_r(N_Er_max)                                      ! Electric field in kV/m
     
     real, allocatable :: v_th(:)        ! Thermal velocity       [m/s]:  dimension(N_species)
     real, allocatable :: nu(:,:)        ! Collisionality         [m^-1]: dimension(N_species,0:N_x)
     real, allocatable :: A(:,:,:)       ! Thermodynamical forces [m^-1]: dimension(3,N_Er,N_species)     
     real, allocatable :: L_ij(:,:,:,:)  ! Onsager matrix:                dimension(3,3,N_Er,N_species)         
     real, allocatable :: D_ij(:,:,:,:,:)! Monoenergetic Onsager matrix:  dimension(3,3,N_Er,N_species,0:N_x)       
     real, allocatable :: D_33_Sp(:,:,:) ! Spitzer D33:  dimension(N_Er,N_species,0:N_x)       
     real, allocatable :: PartFlux(:,:),  &    ! Radial flux of particles Gamma_a^r in 1/(m^2 s)
                          HeatFlux(:,:),  &    ! Radial heat flux  Q_a^r/T in 1/(m^2 s)
                          ParFlow(:,:),   &    ! Parallel flow <n_a V_|| B>/B_00 in 1e19/(m^2 s)
                          ParCurrent(:,:)      ! Parallel current Z_a<n_a e V_|| B>/B_00 in A/m^2 
     
     real, allocatable :: x(:), l(:), d(:), w(:), F1(:,:,:), F3(:,:,:), t_clock(:), t_cpu(:) 
     real, allocatable :: C_ij(:,:,:,:) ! Factor for integrating monoenergetic Onsager matrix:     dimension(3,3,0:N_x,N_species)    
     real :: hat_E_r, t0, t1 
     integer :: N_theta, N_zeta, N_xi, N_x, N_Er, N_species, ierr
     integer :: i, j, k, ii, jj, c0, c1, rate
     
     ! *** Read parameters 
     E_r = -1d14 ; n = -1d14 ; Z = -1d14 ; T = -1d14 ; m = -1d14 
     open(1, file= "monkes_input.fluxes", status="old") 
     read(1, nml=parameters, iostat=ierr) ; close(1)
      
     if( mod(N_theta,2) == 0 ) N_theta = N_theta+1
     if( mod(N_zeta,2) == 0  ) N_zeta  = N_zeta+1
     N_Er = count(E_r /= -1d14)
     N_species = count(n /= -1d14)
     
     write(*,*) " Plasma species. Z is atomic number. m is mass in atomic mass units. "
     write(*,*) " n is density in 10^19/m^3. T is temperature in eV. "
     write(*,*) " n_prime is density gradient in 1/m.  "
     write(*,*) " T_prime is temperature gradient in 1/m.  "
     write(*,'(9999A25)') "Z", "m [u]", "n [10^19 m^-3]", "T [eV]", "n_prime [m^-1]", "T_prime [m^-1]" 
     do j = 1, N_species
        write(*,'(9999e25.16)') Z(j), m(j), n(j), T(j), n_prime(j), T_prime(j)     
     end do 
     
     ! Allocate arrays for quantities
     allocate( v_th(N_species), A(3,2*N_Er,N_species) )    
     ! Thermal velocity in m/s
     v_th = 1.389165e4* sqrt(T/m)    
     
     ! Thermodynamical forces in m^-1
     do j = 1, N_species     
        ! Cases with -Er
        do i = 1, N_Er
           A(1,i,j) = n_prime(j)/n(j) - 3*(T_prime(j)/T(j))/2 + Z(j)*E_r(i)*1d3/T(j)
           A(2,i,j) = T_prime(j)/T(j) 
           A(3,i,j) = 0            ! TO DO: read also <E_||B>   
        end do
        
        ! Cases with +Er
        do i = 1, N_Er
           A(1,i+N_Er,j) = n_prime(j)/n(j) - 3*(T_prime(j)/T(j))/2 - Z(j)*E_r(i)*1d3/T(j)
           A(2,i+N_Er,j) = T_prime(j)/T(j) 
           A(3,i+N_Er,j) = 0            ! TO DO: read also <E_||B>   
        end do
     end do     
     
     ! *** Grid in normalized speed x 
     allocate( x(0:N_x), l(0:N_x), d(0:N_x), w(0:N_x) )
     call Maxwell_polynomials( N_x, l, d, x, w ) 
     
     ! *** Compute collisionality for each species at each point of the
     ! grid in x = v/v_th
     allocate( nu(N_species, 0:N_x)  )
     do k = 0, N_x
        nu(:,k) = Collisionality_a( Z(1:N_species), n(1:N_species), &
                                    T(1:N_species), m(1:N_species), x(k) ) 
     end do 
     call Write_monkes_speed_grid 
     
     ! *** Factor for integrating the monoenergetic coefficients
     allocate( C_ij(3,3,0:N_x,N_species) )
     do j = 1, N_species 
        C_ij(1,1,:,j) = ( 2*n(j)/sqrt(pi) ) * 2* 4.552878* sqrt( ( T(j)*1d-3 )**3 * m(j)  ) / Z(j)**2  * x**2 * (-1)* x**3
        C_ij(1,2,:,j) = x**2 * C_ij(1,1,:,j) ; C_ij(2,1,:,j) = C_ij(1,2,:,j) 
        C_ij(2,2,:,j) = x**2 * C_ij(1,2,:,j)  
 
        C_ij(1,3,:,j) = -( 2*n(j)/sqrt(pi) ) *( 2*T(j)/Z(j) ) * x**2 * x**2 
        C_ij(2,3,:,j) = x**2 * C_ij(1,3,:,j)

        C_ij(3,1,:,j) = -C_ij(1,3,:,j)
        C_ij(3,2,:,j) = x**2 * C_ij(3,1,:,j)

        C_ij(3,3,:,j) = ( 2d0*n(j)/sqrt(pi) )* v_th(j) * x**2 * x          
     end do
     
     ! *** Computation of Onsager matrix and fluxes       
     call system_clock(count_rate=rate) ! Get click rate for timing
     allocate( t_clock(0:N_x), t_cpu(0:N_x)  ) 
     allocate( L_ij(3,3,N_Er,N_species), D_ij(3,3,N_Er,N_species,0:N_x), D_33_Sp(N_Er,N_species,0:N_x)  ) 
     allocate( PartFlux(2*N_Er,N_species), HeatFlux(2*N_Er,N_species), ParFlow(2*N_Er,N_species), ParCurrent(2*N_Er,N_species) ) 
     allocate( F1(0:N_theta-1, 0:N_zeta-1,0:2), F3(0:N_theta-1, 0:N_zeta-1,0:2) )
     
     open(111,file="monkes_Monoenergetic_Database.dat")
     write(111,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]",        &
                             " N_theta ", " N_zeta ", " N_xi ",          &
                             " D_11 ", " D_31 ",                         &
                             " D_13 ", " D_33 ",                         &
                             " D_33_Spitzer ",                           &
                             " E_r [kV /m] ",                            &
                             " Z ", " m ", " n ", " T ",                 &
                             " Wall-clock time [s] ",                    &
                             " CPU time [s] "
     open(112, file="monkes_Onsager_Matrix.dat")
     write(112,'(9999A25)')  "E_r [kV /m]",                             &
                             "N_theta ", " N_zeta ", " N_xi ", " N_x ", &
                             "L_11 [T^2 m/s]", "L_21 [T^2 m/s]", "L_31 [1/(m s)]",  &
                             "L_12 [T^2 m/s]", "L_22 [T^2 m/s]", "L_32 [1/(m s)]",  &
                             "L_13 [1/(m s)]", "L_23 [1/(m s)]", "L_33 [1/(m s)]",  &
                             "Z", "m", "n", "T",                                    &
                             "Wall-clock time [s]",                                 &
                             "CPU time [s]"
     open(113, file="monkes_fluxes.dat")
     write(113,'(9999A25)') "s", "E_r [kV/m]", "Z", "m [u]", "n [10^19 1/m^3]", &
                            "T [eV]", "dn/dr [1/m]", "dT/dr [1/m]",             &
                            "PartFlux [1/(m^2 s)]" , "HeatFlux [1/(m^2 s)]",    &
                            "ParFlow [1/(m^2 s)]"  , "ParCurrent [A/m^2]",      & 
                            "Wall-clock time [s]"  , "CPU time [s]"
     
     open(114, file="monkes_ambipolarity.dat")
     write(114,'(9999A25)') "s", "E_r [kV/m]", "RadCurrent [1/(m^2 s)]" 
     
     do i = 1, N_Er
        do j = 1, N_species
           
           ! Compute the monoenergetic database for each radial electric field, species and collisionality
           do k = N_x, 0, -1           
           
              ! Solve the DKE and extract 3 Legendre modes (N_xi_DF=2) of the solution  
              hat_E_r = ( E_r(i)*1d+3 )/( x(k)*v_th(j) ) !
              
              call system_clock(c0) ; call cpu_time(t0)
              call Solve_BTD_DKE_Legendre_DF( N_theta, N_zeta, N_xi,    &
                                              nu(j,k), hat_E_r,           &
                                              D_ij(:,:,i,j,k),          &
                                              D_33_Sp(i,j,k),           &
                                              2, F1, F3 )
              call system_clock(c1) ; call cpu_time(t1)
              t_clock(k) = ( c1 - c0 )*1d0 / rate
              t_cpu(k)  = t1-t0

              write(111,'(9999e25.16)') nu(j,k), hat_E_r, &
                   real(N_theta), real(N_zeta), real(N_xi), &
                   D_ij(1,1,i,j,k), D_ij(3,1,i,j,k), &
                   D_ij(1,3,i,j,k), D_ij(3,3,i,j,k),  D_33_Sp(i,j,k), &
                   E_r(i), Z(j), m(j), n(j), T(j), &
                   t_clock(k), t_cpu(k)                   
           end do
           
           ! Compute thermal coefficients for each species and radial electric field            
           do jj = 1, 3 ; do ii = 1, 3
              L_ij(ii,jj,i,j) = dot_product( D_ij(ii,jj,i,j,:), C_ij(ii,jj,:,j)* w )          
           end do ; end do 
            
           write(112,'(9999e25.16)') E_r(i), real(N_theta), real(N_zeta), &
                                     real(N_xi), real(N_x), &
                                     L_ij(:,1,i,j), L_ij(:,2,i,j), L_ij(:,3,i,j), &
                                     Z(j), m(j), n(j), T(j), &
                                     sum(t_clock), sum(t_cpu)
                                     
           ! Compute neoclassical fluxes using Onsager matrix for -Er           
           PartFlux(i,j)   = dot_product( L_ij(1,:,i,j), A(:,i,j) )
           HeatFlux(i,j)   = dot_product( L_ij(2,:,i,j), A(:,i,j) )
           ParFlow(i,j)    = dot_product( L_ij(3,:,i,j), A(:,i,j) )           
           ParCurrent(i,j) = (e*1d19) * ParFlow(i,j)
                                     
           write(113,'(9999e25.16)') s, -E_r(i), Z(j), m(j), n(j), T(j), & 
                                     n_prime(j), T_prime(j),  &
                                     PartFlux(i,j), HeatFlux(i,j), &
                                     ParFlow(i,j) , ParCurrent(i,j), & 
                                     sum(t_clock), sum(t_cpu)
                                     
           ! Compute neoclassical fluxes using Onsager matrix for +Er           
           PartFlux(i+N_Er,j)   = dot_product( L_ij(1,:,i,j), A(:,i+N_Er,j) )
           HeatFlux(i+N_Er,j)   = dot_product( L_ij(2,:,i,j), A(:,i+N_Er,j) )
           ParFlow(i+N_Er,j)    = dot_product( L_ij(3,:,i,j), A(:,i+N_Er,j) )            
           ParCurrent(i+N_Er,j) = (e*1d19) * ParFlow(i+N_Er,j)
           
           
           write(113,'(9999e25.16)') s, E_r(i), Z(j), m(j), n(j), T(j), & 
                                     n_prime(j), T_prime(j),  &
                                     PartFlux(i+N_Er,j), HeatFlux(i+N_Er,j), &
                                     ParFlow(i+N_Er,j) , ParCurrent(i+N_Er,j), & 
                                     sum(t_clock), sum(t_cpu)
           
        end do 
        
        write(114,'(9999e25.16)') s, -E_r(i), sum( Z(1:N_species)*PartFlux(i,:) )
        write(114,'(9999e25.16)') s,  E_r(i), sum( Z(1:N_species)*PartFlux(i+N_Er,:) )
        
     end do 
     close(111) ; close(112) ; close(113); close(114)
     
     
     contains
     
     subroutine Write_monkes_speed_grid     
        integer :: j, k 
        
        open(121,file="monkes_speed_grid.dat")  
        do j = 1, N_species 
           write(121,*) " Species ", j
           write(121,'(9999A25)')     "n(j)", "Z(j)", "T(j)", "m(j)", "v_th(j)"
           write(121,'(9999e25.16)')   n(j) ,  Z(j) ,  T(j) ,  m(j) ,  v_th(j) 
           write(121,'(9999A25)') " x ", "nu [m^-1]", "weight"
           do k = 0, N_x                
              write(121,'(9999e25.16)') x(k), nu(j,k), w(k)
           end do 
        end do 
        close(121)
     
     end subroutine 
     
     ! Gives the x-dependent factor in the collision frequency nu_ab(v)/v,
     ! evaluated at x = v/v_ta. Here, k_ab = v_ta / v_tb  is the ratio between
     ! the thermal velocities of species "a" and "b". 
     ! For self-collisions k_ab=1
     real function nu_function(x, k_ab) result(F)
         real, intent(in) :: x, k_ab 
         
         real ::  y, G
         
         y = k_ab * x ! x_b
         G = ( erf(y)/2 - y*exp(-y**2)/sqrt(pi) )/ y**2 ! Chandrasekhar function       
         F = ( erf(y) - G )/ x**4
         
     end function 
     
     ! Given the atomic numbers ( Z(1), Z(2) ), densities ( n(1), n(2) ),
     ! temperatures ( T(1), T(2) ) and masses ( m(1), m(2) ), computes the Coulomb logarithm
     ! associated to collisions of species "1" with species "2" (i.e. for nu_12=nu_ab).
     ! Densities should be in 10^{-13} particles/m^3
     ! Temperatures should be in eV
     ! Masses can be in whatever units, the Coulomb logarithm only depends on
     ! the ratio m_a/m_b. We can also use the mass number.
     ! Formulas extracted from NRL plasma formulary (JD Huba).
     real function Coulomb_logarithm(Z, n, T, m) result(F) 
         real, intent(in) :: Z(2), n(2), T(2), m(2)
         
         real :: n_a, n_b, Z_a, Z_b, T_a, T_b, m_a, m_b, C
         
         n_a = n(1)*1e13 ; T_a = T(1) ; Z_a = Z(1) ; m_a = m(1) 
         n_b = n(2)*1e13 ; T_b = T(2) ; Z_b = Z(2) ; m_b = m(2) 
         
         if( Z_a < 0 ) then
           if( Z_b < 0 ) then    !ee
              F = 23.5 - log( sqrt(n_a)*( T_a**(-5.0/4) ) )    &
                       - sqrt( 1d-5 + ( (log(T_a) - 2)**2 )/16 )
                        
           else                    !ei and iz
              F = 24.0 - log( sqrt(n_a)/ T_a )
           end if
         else
           if( Z_b < 0 ) then    !ie and ze
              F = 24 - log( sqrt(n_b)/T_b )
           else                    !ii, iz, zi, and zz
              C = Z_a * Z_b *(m_a + m_b) / (m_a*T_b + m_b*T_a) 
              F = 23 &
                - log( C * sqrt( n_a*Z_a**2 / T_a + n_b*Z_b**2 / T_b ) )
           end if 
         end if
     
     end function 
     
     ! *** Gives the collisionality 
     ! nu^{ab}/v(x_a) =(1d-7*c**2)**2 * pi * n_b * log Lambda_ab * (2 *e_a * e_b / m_a v_tha^2 )^2  
     !                  * (erf(x_b) - G(x_b))/x_a^4
     ! in [m^-1]
     ! Temperatures are in [eV]. Densities are in [10^19 m^3]. 
     ! In these units v_th[m/s] =sqrt(2k[J/K] T[K]/m[kg])=~ 1.389165e4* sqrt(T [eV]/m [amu]).
     ! The collisionality is given by nu = nu_0 * nu_function(x, k_ab)
     ! where nu_0 = c_units * Z_a^2 Z_b^2 * n_b[10^19 m^-3] / T_a[eV]**2. Here,
     ! c_units = 1d19*pi*e**2/(4*pi*eps_0)**2 
     real function Collisionality_ab( Z, n, T, m, x ) result(nu) 
         real, intent(in) :: Z(2), n(2), T(2), m(2), x
          
         real :: log_lambda, k_ab, nu_0     
         
         ! Coulomb logarithm, nu_0 and k_ab=v_ta / v_tb
         log_lambda = Coulomb_logarithm(Z, n, T, m)
         nu_0 = c_units * n(2) * Z(1)**2 * Z(2)**2 * log_lambda / T(1)**2 
         
         ! Ratio between thermal velocities v_th(1)/v_th(2)
         k_ab = sqrt( T(1)/T(2) ) * sqrt( m(2) / m(1) )
          
         ! Collision frequency
         nu = nu_0 * nu_function(x, k_ab)       
         
     end function 
     
     ! *** Given the set of N_species characterized by their atomic number Z(N_species), 
     ! density n(N_species), temperature T(N_species) and mass m(N_species), 
     ! computes the collisionality nu^a/v(x) correspondent to a value
     ! of normalized speed x=v/v_ta.
     function Collisionality_a( Z, n, T, m, x ) result(nu) 
         real, intent(in) :: Z(:), n(size(Z)), T(size(Z)), m(size(Z)), x
         real :: nu( size(Z) )
         
         integer :: i, j, N_species
         
         N_species = size(Z)
         nu = 0 ; 
         do i = 1, N_species
            do j = 1, N_species
               nu(i) = nu(i) + Collisionality_ab( [Z(i),Z(j)], [n(i),n(j)], &
                                                  [T(i),T(j)], [m(i),m(j)], x ) 
            end do 
         end do 
         
     end function
              
   end subroutine


  

   ! *** Computes a monoenergetic coefficients database reading the input file
   ! "monkes_input.parameters" and makes a test derivative of the monoenergetic
   ! coefficients along nu/v. This derivative is particularly simple to compute as
   ! < g_i, Lorentz( f_j) > = -sum_{k>0} (k*(k+1)/(2k+1)) < g_i^(k) * f_j^(k) > where
   ! <f,g> = < int fg dxi  >is the inner product and < f > the FSA. 
   subroutine Adjoint_derivatives_Input
     real, parameter :: pi = acos(-1d0) 
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
     real, allocatable :: G1(:,:,:), G3(:,:,:)
     real, allocatable :: dGamma_ij_dlambda(:,:,:) 
     real, allocatable :: D(:,:,:,:,:,:,:), D_33_Sp(:,:,:,:,:)  
     real, allocatable :: D_nu(:,:,:,:,:,:,:), D_Er(:,:,:,:,:,:,:)
     real, allocatable :: D_adj(:,:,:,:,:,:,:), D_33_Sp_adj(:,:,:,:,:)  
     real, allocatable :: Dij_FS(:,:,:,:)
     real, allocatable :: D_mn(:,:,:,:,:,:,:,:)
     integer, parameter :: N_mn = 1 ! Number of derivatives w.r.t B_mn
     integer :: mn_modes(1,2)       ! Set of modes for derivatives
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
     allocate( D_nu( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
     allocate( D_Er( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
     allocate( D_mn( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi, N_mn ) )

     
     allocate( D_adj( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
     allocate( D_33_Sp_adj(  N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
          
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
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database_Adjoint.dat"     
     ! Open output file and write header
     open(121, file=trim(file_path))
     write(121,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_33_Spitzer ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database_nu_derivatives.dat"     
     ! Open output file and write header
     open(221, file=trim(file_path))
     write(221,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_11_nu ", " D_31_nu ", &
                             " D_13_nu ", " D_33_nu ", &
                             " D_33_Spitzer_nu ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database_Er_derivatives.dat"     
     ! Open output file and write header
     open(222, file=trim(file_path))
     write(222,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_11_Er ", " D_31_Er ", &
                             " D_13_Er ", " D_33_Er ", &
                             " D_33_Spitzer_nu ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database_Bmn_derivatives.dat"     
     ! Open output file and write header
     open(223, file=trim(file_path))
     write(223,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_11_mn ", " D_31_mn ", &
                             " D_13_mn ", " D_33_mn ", & 
                             " m ", " n ", " B_mn ",   & 
                             " Wall-clock time [s] ",  &
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
                          allocate( G1(0:N_theta(ii)-1, 0:N_zeta(jj)-1, 0:N_xi_DF(k) ) )
                          allocate( G3(0:N_theta(ii)-1, 0:N_zeta(jj)-1, 0:N_xi_DF(k) ) )

                          ! Solve the DKE and extract N_xi_DF(k)+1 Legendre modes
                          call Compute_Monoenergetic_derivatives_Adjoint( N_theta(ii),          &
                                                                          N_zeta(jj),           &
                                                                          N_xi(kk),             &
                                                                          N_xi_DF(k),           & 
                                                                          nu(i), E_r(j),        &        
                                                                          D(:,:,i,j,ii,jj,kk),  &
                                                                          D_33_Sp(i,j,ii,jj,kk),&
                                                                          D_adj(:,:,i,j,ii,jj,kk),  &
                                                                          D_33_Sp_adj(i,j,ii,jj,kk),&
                                                                          D_nu(:,:,i,j,ii,jj,kk),   & ! nu derivatives
                                                                          D_Er(:,:,i,j,ii,jj,kk),   & ! E_r derivatives
                                                                          mn_modes,                 & ! List of modes for B_mn derivatives
                                                                          D_mn(:,:,i,j,ii,jj,kk,:), & ! B_mn derivatives
                                                                          F1, F3, G1, G3 )
 
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
                                        t_clock/2, (t1-t0)/2
                          
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
                                        t_clock/2, (t1-t0)/2
                          flush(21)

                                                    
                          ! Writing results on "monkes_Monoenergetic_Database_Adjoint.dat"
                          write(121,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D_adj(1,1,i,j,ii,jj,kk), &
                                        D_adj(3,1,i,j,ii,jj,kk), &
                                        D_adj(1,3,i,j,ii,jj,kk), &
                                        D_adj(3,3,i,j,ii,jj,kk), &
                                        D_33_Sp_adj(i,j,ii,jj,kk), &
                                        t_clock/2, (t1-t0)/2
                          flush(121)

                          
                          
                          ! Writing results on "monkes_Monoenergetic_Database_nu_Derivatives.dat"
                          write(221,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_nu(1,1,i,j,ii,jj,kk), &
                                        D_nu(3,1,i,j,ii,jj,kk), &
                                        D_nu(1,3,i,j,ii,jj,kk), &
                                        D_nu(3,3,i,j,ii,jj,kk), &
                                        -D_33_Sp(i,j,ii,jj,kk)/nu(i), & 
                                        t_clock, (t1-t0)
                          flush(221)
                          
                          ! Writing results on "monkes_Monoenergetic_Database_Er_Derivatives.dat"
                          write(222,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_Er(1,1,i,j,ii,jj,kk), &
                                        D_Er(3,1,i,j,ii,jj,kk), &
                                        D_Er(1,3,i,j,ii,jj,kk), &
                                        D_Er(3,3,i,j,ii,jj,kk), & 
                                        0d0, &
                                        t_clock, (t1-t0)
                          flush(222)
                          
                          ! Writing results on "monkes_Monoenergetic_Database_Er_Derivatives.dat"
                          write(223,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_mn(1,1,i,j,ii,jj,kk,1), &
                                        D_mn(3,1,i,j,ii,jj,kk,1), &
                                        D_mn(1,3,i,j,ii,jj,kk,1), &
                                        D_mn(3,3,i,j,ii,jj,kk,1), & 
                                        0d0, &
                                        t_clock, (t1-t0)
                          flush(223)

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
                          deallocate( G1, G3 )
                            
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do    
     close(21) ! close monkes_Monoenergetic_database.dat
     close(121) ! close monkes_Monoenergetic_database_Adjoint.dat
     close(221) ! close monkes_Monoenergetic_database_nu_derivative.dat
     close(222) ! close monkes_Monoenergetic_database_Er_derivative.dat
     close(223) ! close monkes_Monoenergetic_database_Bmn_derivative.dat
     ! CLOSE (if necessary) monkes_Monoenergetic_lambda.dat 
     if( Monoenergetic_lambda ) close(31)
     ! CLOSE (if necessary) monkes_Monoenergetic_theta_zeta.dat 
     if( Monoenergetic_theta_zeta )  close(41)
     
  end subroutine






   ! *** Computes a monoenergetic coefficients database reading the input file
   ! "monkes_input.parameters" and makes a test derivative of the monoenergetic
   ! coefficients along nu/v, Er and some B_mn
   subroutine Test_Adjoint_derivatives
     real, parameter :: pi = acos(-1d0) 
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
     real, allocatable :: G1(:,:,:), G3(:,:,:)
     real, allocatable :: dGamma_ij_dlambda(:,:,:) 
     real, allocatable :: D(:,:,:,:,:,:,:), D_33_Sp(:,:,:,:,:)  
     real, allocatable :: D_nu(:,:,:,:,:,:,:), D_Er(:,:,:,:,:,:,:)
     real, allocatable :: D_adj(:,:,:,:,:,:,:), D_33_Sp_adj(:,:,:,:,:)  
     real, allocatable :: Dij_FS(:,:,:,:)
     real, allocatable :: D_mn(:,:,:,:,:,:,:,:)
     integer, parameter :: N_mn = 1 ! Number of derivatives w.r.t B_mn
     integer :: mn_modes(1,2)       ! Set of modes for derivatives
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
     allocate( D_nu( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
     allocate( D_Er( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
     allocate( D_mn( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi, N_mn ) )

     
     allocate( D_adj( 3, 3, N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
     allocate( D_33_Sp_adj(  N_nu, N_E_r, M_theta, M_zeta, M_xi ) )
          
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
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database_Adjoint.dat"     
     ! Open output file and write header
     open(121, file=trim(file_path))
     write(121,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_33_Spitzer ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database_nu_derivatives.dat"     
     ! Open output file and write header
     open(221, file=trim(file_path))
     write(221,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_11_nu ", " D_31_nu ", &
                             " D_13_nu ", " D_33_nu ", &
                             " D_33_Spitzer_nu ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database_Er_derivatives.dat"     
     ! Open output file and write header
     open(222, file=trim(file_path))
     write(222,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_11_Er ", " D_31_Er ", &
                             " D_13_Er ", " D_33_Er ", &
                             " D_33_Spitzer_nu ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database_Bmn_derivatives.dat"     
     ! Open output file and write header
     open(223, file=trim(file_path))
     write(223,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_11_mn ", " D_31_mn ", &
                             " D_13_mn ", " D_33_mn ", & 
                             " m ", " n ", " B_mn ",   & 
                             " Wall-clock time [s] ",  &
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
                          allocate( G1(0:N_theta(ii)-1, 0:N_zeta(jj)-1, 0:N_xi_DF(k) ) )
                          allocate( G3(0:N_theta(ii)-1, 0:N_zeta(jj)-1, 0:N_xi_DF(k) ) )

                          ! Solve the DKE and extract N_xi_DF(k)+1 Legendre modes
                          call Compute_Monoenergetic_derivatives_Adjoint( N_theta(ii),          &
                                                                          N_zeta(jj),           &
                                                                          N_xi(kk),             &
                                                                          N_xi_DF(k),           & 
                                                                          nu(i), E_r(j),        &        
                                                                          D(:,:,i,j,ii,jj,kk),  &
                                                                          D_33_Sp(i,j,ii,jj,kk),&
                                                                          D_adj(:,:,i,j,ii,jj,kk),  &
                                                                          D_33_Sp_adj(i,j,ii,jj,kk),&
                                                                          D_nu(:,:,i,j,ii,jj,kk),   & ! nu derivatives
                                                                          D_Er(:,:,i,j,ii,jj,kk),   & ! E_r derivatives
                                                                          mn_modes,                 & ! List of modes for B_mn derivatives
                                                                          D_mn(:,:,i,j,ii,jj,kk,:), & ! B_mn derivatives
                                                                          F1, F3, G1, G3 )
 
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
                                        t_clock/2, (t1-t0)/2
                          
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
                                        t_clock/2, (t1-t0)/2
                          flush(21)

                                                    
                          ! Writing results on "monkes_Monoenergetic_Database_Adjoint.dat"
                          write(121,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D_adj(1,1,i,j,ii,jj,kk), &
                                        D_adj(3,1,i,j,ii,jj,kk), &
                                        D_adj(1,3,i,j,ii,jj,kk), &
                                        D_adj(3,3,i,j,ii,jj,kk), &
                                        D_33_Sp_adj(i,j,ii,jj,kk), &
                                        t_clock/2, (t1-t0)/2
                          flush(121)

                          
                          
                          ! Writing results on "monkes_Monoenergetic_Database_nu_Derivatives.dat"
                          write(221,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_nu(1,1,i,j,ii,jj,kk), &
                                        D_nu(3,1,i,j,ii,jj,kk), &
                                        D_nu(1,3,i,j,ii,jj,kk), &
                                        D_nu(3,3,i,j,ii,jj,kk), &
                                        -D_33_Sp(i,j,ii,jj,kk)/nu(i), & 
                                        t_clock, (t1-t0)
                          flush(221)
                          
                          ! Writing results on "monkes_Monoenergetic_Database_Er_Derivatives.dat"
                          write(222,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_Er(1,1,i,j,ii,jj,kk), &
                                        D_Er(3,1,i,j,ii,jj,kk), &
                                        D_Er(1,3,i,j,ii,jj,kk), &
                                        D_Er(3,3,i,j,ii,jj,kk), & 
                                        0d0, &
                                        t_clock, (t1-t0)
                          flush(222)
                          
                          ! Writing results on "monkes_Monoenergetic_Database_Er_Derivatives.dat"
                          write(223,'(9999e25.16)') nu(i), E_r(j), &                               
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_mn(1,1,i,j,ii,jj,kk,1), &
                                        D_mn(3,1,i,j,ii,jj,kk,1), &
                                        D_mn(1,3,i,j,ii,jj,kk,1), &
                                        D_mn(3,3,i,j,ii,jj,kk,1), & 
                                        0d0, &
                                        t_clock, (t1-t0)
                          flush(223)

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
                          deallocate( G1, G3 )
                            
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do    
     close(21) ! close monkes_Monoenergetic_database.dat
     close(121) ! close monkes_Monoenergetic_database_Adjoint.dat
     close(221) ! close monkes_Monoenergetic_database_nu_derivative.dat
     close(222) ! close monkes_Monoenergetic_database_Er_derivative.dat
     close(223) ! close monkes_Monoenergetic_database_Bmn_derivative.dat
     ! CLOSE (if necessary) monkes_Monoenergetic_lambda.dat 
     if( Monoenergetic_lambda ) close(31)
     ! CLOSE (if necessary) monkes_Monoenergetic_theta_zeta.dat 
     if( Monoenergetic_theta_zeta )  close(41)
     
  end subroutine






   
end module
