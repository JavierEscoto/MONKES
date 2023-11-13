module API_Example_DKE_BTD_Solution_Legendre

  use DKE_BTD_Solution_Legendre
  
  implicit none
  
  private
  
  public :: Monoenergetic_Database_Example_1
  public :: Monoenergetic_Database_Example_Reusing_1
  public :: DKE_Spatial_Convergence_Example_1
  public :: DKE_zeta_Convergence_Example_1, DKE_zeta_Convergence_Example_2
  public :: DKE_zeta_Convergence_Example_Reusing_1
  
  public :: DKE_zeta_Convergence_Example_Low_resolution
  public :: DKE_zeta_Convergence_Example_Medium_resolution
  public :: DKE_zeta_Convergence_Example_Big_Resolution
  public :: DKE_zeta_Convergence_Example_Big_Resolution_2
  public :: DKE_zeta_Convergence_Scan_Example
  
  
  contains
  
   
  
  subroutine Monoenergetic_Database_Example_1 
    integer, parameter :: N_theta = 27, N_zeta = 55, N_xi = 140
    integer, parameter :: N_nu = 36, N_E_rho = 10
    real, parameter :: nu_min = 1d-5 , nu_max = 3e2!300 
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho), Gamma_31(N_nu, N_E_rho)
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
    
    if( N_theta < 10 )  write(N_theta_char,'(I1)') N_theta 
    if( 10 <= N_theta .and. N_theta < 100 )  write(N_theta_char,'(I2)') N_theta
    if( 100 <= N_theta .and. N_theta < 1000 )  write(N_theta_char,'(I3)') N_theta
    
    if( N_zeta < 10 )  write(N_zeta_char,'(I1)') N_zeta 
    if( 10 <= N_zeta .and. N_zeta < 100 )  write(N_zeta_char,'(I2)') N_zeta
    if( 100 <= N_zeta .and. N_zeta < 1000 )  write(N_zeta_char,'(I3)') N_zeta
    
    if( N_xi < 10 )  write(N_xi_char,'(I1)') N_xi 
    if( 10 <= N_xi .and. N_xi < 100 )  write(N_xi_char,'(I2)') N_xi
    if( 100 <= N_xi .and. N_xi < 1000 )  write(N_xi_char,'(I3)') N_xi
    
    folder_path = "results/Monoenergetic_Database_Example_"&
             // trim(N_theta_char) //"_"// trim(N_zeta_char) //"_"// trim(N_xi_char)
    
    file_path = trim(folder_path) // "/Gamma_11_Gamma_31"
    call execute_command_line( "mkdir "//trim(folder_path) )
    
    
    
    ! Initialize compiler internal clock 
    call system_clock(count_rate=c_rate) ; rate = real(c_rate)
    
    open(1221, file=trim(file_path) // ".plt" )
    write(1221,*) " N_nu, N_E_rho  ", N_nu, N_E_rho
    write(1221,*) " nu, E_rho, Gamma_11, Gamma_31 "
    do j =1, 5, 4
       if( j < 10 ) write(E_rho_char,'(I1)') j 
       if( 10 <= j .and. j < 100 ) write(E_rho_char,'(I2)') j 
       open(2221, file=trim(file_path)// "_" // trim(E_rho_char) // ".plt" )
       write(2221,*) " nu ", " E_rho ", " Gamma_11 " , " Gamma_31 ", " t_clock "
       do i = 1, N_nu
           
          call system_clock(c0)
          call Compute_Monoenergetic_Database( N_theta, N_zeta, N_xi, &
                  nu(i:i), E_rho(j:j), Gamma_11(i:i,j:j), Gamma_31(i:i,j:j) )
          call system_clock(c1)
          t0 = c0 / rate ; t1 = c1 / rate ! Instants of time in seconds
          write(1221,'(9999e)') nu(i), E_rho(j), Gamma_11(i,j), Gamma_31(i,j), t1-t0  
          write(2221,'(9999e)') nu(i), E_rho(j), Gamma_11(i,j), Gamma_31(i,j), t1-t0  
       end do
       close(2221)
    end do
    close(1221)
     

  end subroutine

  subroutine Monoenergetic_Database_Example_Reusing_1
    integer, parameter :: N_theta = 15, N_zeta = 125, N_xi = 120
    integer, parameter :: N_nu = 36, N_E_rho = 10
    real, parameter :: nu_min = 1d-5 , nu_max = 3e2!300 
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho), Gamma_31(N_nu, N_E_rho)
    real :: Gamma_13(N_nu,N_E_rho), Gamma_33(N_nu, N_E_rho)
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
    
    if( N_theta < 10 )  write(N_theta_char,'(I1)') N_theta 
    if( 10 <= N_theta .and. N_theta < 100 )  write(N_theta_char,'(I2)') N_theta
    if( 100 <= N_theta .and. N_theta < 1000 )  write(N_theta_char,'(I3)') N_theta
    
    if( N_zeta < 10 )  write(N_zeta_char,'(I1)') N_zeta 
    if( 10 <= N_zeta .and. N_zeta < 100 )  write(N_zeta_char,'(I2)') N_zeta
    if( 100 <= N_zeta .and. N_zeta < 1000 )  write(N_zeta_char,'(I3)') N_zeta
    
    if( N_xi < 10 )  write(N_xi_char,'(I1)') N_xi 
    if( 10 <= N_xi .and. N_xi < 100 )  write(N_xi_char,'(I2)') N_xi
    if( 100 <= N_xi .and. N_xi < 1000 )  write(N_xi_char,'(I3)') N_xi
    
    folder_path = "results/Monoenergetic_Database_Example_Reusing_NEW_"&
    !folder_path = "results/Monoenergetic_Database_Example_Reusing_"&
             // trim(N_theta_char) //"_"// trim(N_zeta_char) //"_"// trim(N_xi_char)
    
    file_path = trim(folder_path) // "/Gamma_11_Gamma_31"
    call execute_command_line( "mkdir "//trim(folder_path) )   
    
    
    ! Initialize compiler internal clock 
    call system_clock(count_rate=c_rate) ; rate = real(c_rate)
    
    open(1221, file=trim(file_path) // ".plt" )
    write(1221,*) " N_nu, N_E_rho  ", N_nu, N_E_rho
    write(1221,*) " nu, E_rho, Gamma_11, Gamma_31, Gamma_13, Gamma_33, t_clock "
    do j = 1, 5, 4!, !5!N_E_rho
       if( j < 10 ) write(E_rho_char,'(I1)') j 
       if( 10 <= j .and. j < 100 ) write(E_rho_char,'(I2)') j 
       open(2221, file=trim(file_path)// "_" // trim(E_rho_char) // ".plt" )
       write(2221,*) " nu, E_rho, Gamma_11, Gamma_31, Gamma_13, Gamma_33, t_clock "
       do i = 1, N_nu
          
          write(*,*) " nu, E_rho = ", nu(i), E_rho(j)
          call system_clock(c0)
          call Solve_BTD_DKE_Legendre_Reusing( N_theta, N_zeta, N_xi, &
                  nu(i), E_rho(j), Gamma_11(i,j), Gamma_31(i,j), Gamma_13(i,j), Gamma_33(i,j) )
          call system_clock(c1)
          t0 = c0 / rate ; t1 = c1 / rate ! Instants of time in seconds
          write(1221,'(9999e)') nu(i), E_rho(j), Gamma_11(i,j), Gamma_31(i,j), Gamma_13(i,j), Gamma_33(i,j), t1-t0  
          write(2221,'(9999e)') nu(i), E_rho(j), Gamma_11(i,j), Gamma_31(i,j), Gamma_13(i,j), Gamma_33(i,j), t1-t0  
          
       end do
       close(2221)
    end do
    close(1221)
     

  end subroutine

  subroutine DKE_Spatial_Convergence_Example_1
    integer, parameter :: M_theta = 1, M_zeta = 4
    real, parameter ::    nu     =  1d-6,  E_rho      =  1d-1
    character(len=20) :: nu_char = "1e-6", E_rho_char = "1d-3"
    integer, parameter :: N_xi = 140 
    integer :: N_theta(M_theta), N_zeta(M_zeta)
    real :: Gamma_11(M_theta, M_zeta), Gamma_31(M_theta, M_zeta)
    real :: t_clock(M_theta, M_zeta), residual(M_theta, M_zeta)
    real :: t0, t1, rate
    integer :: c0, c1, c_rate
    character(len=500) :: folder_path, file_path
    character(len=20) N_theta_char, N_zeta_char, N_xi_char
    integer :: i, j

    N_theta = [ 31 ] ; N_zeta = [ 31, 41, 51, 61 ]
    write(   nu_char,'(e9.3)') nu
    write(E_rho_char,'(e9.3)') E_rho
    N_xi_char = int_2_str(N_xi)
    
    folder_path = "results/DKE_Spatial_Convergence_Example_nu_"&
             // trim(nu_char) //"_Erho_"// trim(E_rho_char) //"_"// trim(N_xi_char)
    call execute_command_line( "mkdir "//trim(folder_path) )
    
    file_path = trim(folder_path) // "/Gamma_11_Gamma_31"

    call system_clock(count_rate=c_rate) ; rate = c_rate
    open(1221, file=trim(file_path) // ".plt" )
    write(1221,*) " nu, E_rho, N_theta, N_zeta, Gamma_11, Gamma_31 "
    do i = 1, M_theta
       N_theta_char = int_2_str(N_theta(j))
       do j = 1, M_zeta
          open(2221, file=trim(file_path)// "_Convergence_N_theta_" // trim(N_theta_char) // ".plt" )
          write(2221,*) " nu, E_rho, N_theta, N_zeta, Gamma_11, Gamma_31 "
          
          call system_clock(c0)
          call Compute_Monoenergetic_Database( N_theta(i), N_zeta(j), N_xi, &
                          [nu], [E_rho], Gamma_11(i:i,j:j), Gamma_31(i:i,j:j) )
          call system_clock(c1)
           
          t0 = c0 / rate ; t1 = c1 / rate
          t_clock(i, j) = t1 - t0
          write(*,*) "elapsed time = ", t_clock(i, j)  
          write(1221,'(9999e)')  nu, E_rho, real(N_theta(i)), real(N_zeta(j)), Gamma_11(i,j), Gamma_31(i,j), t_clock(i, j)  
          write(2221,'(9999e)')  nu, E_rho, real(N_theta(i)), real(N_zeta(j)), Gamma_11(i,j), Gamma_31(i,j), t_clock(i, j)
       end do
       close(2221)
    end do
    close(1221)
    

  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function

  end subroutine

  subroutine DKE_zeta_Convergence_Example_1
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 30, M_theta = 10, M_xi = 1
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho, M_zeta), Gamma_31(N_nu,N_E_rho, M_zeta)
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char, N_xi_char 
    character(len=2000) :: folder_path, file_path
    
    ! *** Vector of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [1d-6, 3d-6, 1d-5, 3d-5, 1d-4]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! *** Vectors of theta, zeta and xi sizes
    N_theta = [(  3 +  2*i, i=1,M_theta )]
    N_zeta  = [(  3 +  2*j, j=1,M_zeta  )]
    N_xi    = [100] 
    
    ! *** Create folder for output
    N_xi_char = int_2_str(N_xi(1))  
    folder_path = "results/DKE_zeta_Convergence_Example_Nxi_" // trim(N_xi_char) 
    call execute_command_line( "mkdir " // trim(folder_path) )    
    
    !N_theta_char = int_2_str(N_theta)
    
    do i_E_rho = 1, 1!N_E_rho
       do i_nu = 3, 3!N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Gamma_11_Gamma_31_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,*) " nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call DKE_zeta_Convergence( nu(i_nu), E_rho(i_E_rho), N_theta(i), N_zeta(j:j), N_xi(k), &
                        Gamma_11(i_nu,i_E_rho,j:j), Gamma_31(i_nu,i_E_rho,j:j), t_clock(i_nu,i_E_rho,j:j) )
                     
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(221)
          
       end do
    end do
    
  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function
    
  end subroutine


  
  subroutine DKE_zeta_Convergence_Example_Reusing_1
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 13, M_theta = 5, M_xi = 1
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho, M_zeta), Gamma_31(N_nu,N_E_rho, M_zeta)
    real :: Gamma_13(N_nu,N_E_rho, M_zeta), Gamma_33(N_nu,N_E_rho, M_zeta)
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
    N_zeta  = [( 107 + 4*j, j=1,M_zeta  )]
    N_xi    = [160] 
    
    ! *** Create folder for output
    N_xi_char = int_2_str(N_xi(1))  
    folder_path = "results/DKE_zeta_Convergence_Example_Reusing_Nxi_" // trim(N_xi_char) 
    call execute_command_line( "mkdir " // trim(folder_path) )    
    
    !N_theta_char = int_2_str(N_theta)
    
    do i_E_rho = 6, 6!N_E_rho
       do i_nu = 3, 3!N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Gamma_11_Gamma_31_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,*) " nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, Gamma_13, Gamma_33, t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call system_clock(c0)
                   call Solve_BTD_DKE_Legendre_Reusing( N_theta(i), N_zeta(j), N_xi(k), &
                        nu(i_nu), E_rho(i_E_rho), &
                        Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                        Gamma_13(i_nu,i_E_rho,j), Gamma_33(i_nu,i_E_rho,j) )
                   call system_clock(c1)
                   
                   t_clock(i_nu,i_E_rho,j) = (c1-c0)/real(c_rate)
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       Gamma_13(i_nu,i_E_rho,j), Gamma_33(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       Gamma_13(i_nu,i_E_rho,j), Gamma_33(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(221)
          
       end do
    end do
    
  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function
    
  end subroutine
  
  


  subroutine DKE_zeta_Convergence_Example_2
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 17, M_theta = 10, M_xi = 3
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho, M_zeta), Gamma_31(N_nu,N_E_rho, M_zeta)
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char !,  N_theta_char, N_xi_char 
    character(len=2000) :: folder_path, file_path
    
    ! *** Vector of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [1d-6, 3d-6, 1d-5, 3d-5, 1d-4]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! *** Vectors of theta, zeta and xi sizes
    N_theta = [(  17 +  2*i, i=1,M_theta)]
    N_zeta  = [(  61 +  2*j, j=1,M_zeta)]
    N_xi    = [ 150, 160, 170 ]
    
    ! *** Create folder for output
    folder_path = "results/DKE_zeta_Convergence_Example_2"
    call execute_command_line( "mkdir " // trim(folder_path) )    
    
    !N_theta_char = int_2_str(N_theta)
    !N_xi_char = int_2_str(N_xi)  
    
    do i_E_rho = 5, N_E_rho
       do i_nu = 3, N_nu
       
          write(   nu_char,'(e9.3)')       nu(i_nu) ; write(*,*)    nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Gamma_11_Gamma_31_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,*) " nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call DKE_zeta_Convergence( nu(i_nu), E_rho(i_E_rho), N_theta(i), N_zeta(j:j), N_xi(k), &
                        Gamma_11(i_nu,i_E_rho,j:j), Gamma_31(i_nu,i_E_rho,j:j), t_clock(i_nu,i_E_rho,j:j) )
                     
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(221)
          
       end do
    end do
    
  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function
    
  end subroutine
  


  subroutine DKE_zeta_Convergence_Scan_Example
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 30, M_theta = 30, M_xi = 17
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho, M_zeta), Gamma_31(N_nu,N_E_rho, M_zeta)
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char !,  N_theta_char, N_xi_char 
    character(len=2000) :: folder_path, file_path
    
    ! *** Vector of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [1d-6, 3d-6, 1d-5, 3d-5, 1d-4]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! *** Vectors of theta, zeta and xi sizes
    N_theta = [(  3 +  2*i, i=1,M_theta)]
    N_zeta  = [(  3 +  2*j, j=1,M_zeta)]
    N_xi    = [( 20*k, k=1,M_xi)]
    
    ! *** Create folder for output
    folder_path = "results/DKE_zeta_Convergence_Scan"
    call execute_command_line( "mkdir " // trim(folder_path) )    
    
    !N_theta_char = int_2_str(N_theta)
    !N_xi_char = int_2_str(N_xi)  
     
    do i_E_rho = 1, 1!N_E_rho
       do i_nu = 3, N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Gamma_11_Gamma_31_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(21, file=trim(file_path) ) 
          write(21,*) " nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call DKE_zeta_Convergence( nu(i_nu), E_rho(i_E_rho), N_theta(i), N_zeta(j:j), N_xi(k), &
                        Gamma_11(i_nu,i_E_rho,j:j), Gamma_31(i_nu,i_E_rho,j:j), t_clock(i_nu,i_E_rho,j:j) )
                     
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(21,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(21)
          
       end do
    end do
    
  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function
    
  end subroutine
  
  
  
  subroutine DKE_zeta_Convergence_Example_Low_Resolution
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 25, M_theta = 25, M_xi = 8
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho, M_zeta), Gamma_31(N_nu,N_E_rho, M_zeta)
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char !,  N_theta_char, N_xi_char 
    character(len=2000) :: folder_path, file_path
    
    ! *** Vector of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [1d-6, 3d-6, 1d-5, 3d-5, 1d-4]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! *** Vectors of theta, zeta and xi sizes
    N_theta = [(  3 +  2*i, i=1,M_theta )]
    N_zeta  = [(  3 +  2*j, j=1,M_zeta  )]
    N_xi    = [(       20*k, k=1,M_xi )]
    
    ! *** Create folder for output
    folder_path = "results/DKE_zeta_Convergence_Low_Resolution"
    call execute_command_line( "mkdir " // trim(folder_path) )    
    
    !N_theta_char = int_2_str(N_theta)
    !N_xi_char = int_2_str(N_xi)  
    
    do i_E_rho = 1, 1!N_E_rho
       do i_nu = 3, 3!N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Gamma_11_Gamma_31_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,*) " nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call DKE_zeta_Convergence( nu(i_nu), E_rho(i_E_rho), N_theta(i), N_zeta(j:j), N_xi(k), &
                        Gamma_11(i_nu,i_E_rho,j:j), Gamma_31(i_nu,i_E_rho,j:j), t_clock(i_nu,i_E_rho,j:j) )
                     
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(221)
          
       end do
    end do
    
  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function
    
  end subroutine
  
  

  subroutine DKE_zeta_Convergence_Example_Medium_Resolution
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 25, M_theta = 25, M_xi = 5
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho, M_zeta), Gamma_31(N_nu,N_E_rho, M_zeta)
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char !,  N_theta_char, N_xi_char 
    character(len=2000) :: folder_path, file_path
    
    ! *** Vector of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [1d-6, 3d-6, 1d-5, 3d-5, 1d-4]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! *** Vectors of theta, zeta and xi sizes
    N_theta = [(  3 +  2*i, i=1,M_theta )]
    N_zeta  = [(  3 +  2*j, j=1,M_zeta  )]
    N_xi    = [( 100 + 20*k, k=1,M_xi)]
    
    ! *** Create folder for output
    folder_path = "results/DKE_zeta_Convergence_Medium_Resolution"
    call execute_command_line( "mkdir " // trim(folder_path) )    
    
    !N_theta_char = int_2_str(N_theta)
    !N_xi_char = int_2_str(N_xi)  
    
    do i_E_rho = 1, 1!N_E_rho
       do i_nu = 3, 3!N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Gamma_11_Gamma_31_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,*) " nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call DKE_zeta_Convergence( nu(i_nu), E_rho(i_E_rho), N_theta(i), N_zeta(j:j), N_xi(k), &
                        Gamma_11(i_nu,i_E_rho,j:j), Gamma_31(i_nu,i_E_rho,j:j), t_clock(i_nu,i_E_rho,j:j) )
                     
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(221)
          
       end do
    end do
    
  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function
    
  end subroutine
  
  

  subroutine DKE_zeta_Convergence_Example_Big_Resolution
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 15, M_theta = 15, M_xi = 6
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho, M_zeta), Gamma_31(N_nu,N_E_rho, M_zeta)
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char !,  N_theta_char, N_xi_char 
    character(len=2000) :: folder_path, file_path
    
    ! *** Vector of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [1d-6, 3d-6, 1d-5, 3d-5, 1d-4]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! *** Vectors of theta, zeta and xi sizes
    N_theta = [(  17 +  2*i, i=1,M_theta)]
    N_zeta  = [(  61 +  2*j, j=1,M_zeta)]
    N_xi    = [( 80 + 20*k, k=1,M_xi)]
    
    ! *** Create folder for output
    folder_path = "results/DKE_zeta_Convergence_Big_Resolution"
    call execute_command_line( "mkdir " // trim(folder_path) )    
    
    !N_theta_char = int_2_str(N_theta)
    !N_xi_char = int_2_str(N_xi)  
    
    do i_E_rho = 5, N_E_rho
       do i_nu = 3, N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Gamma_11_Gamma_31_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,*) " nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call DKE_zeta_Convergence( nu(i_nu), E_rho(i_E_rho), N_theta(i), N_zeta(j:j), N_xi(k), &
                        Gamma_11(i_nu,i_E_rho,j:j), Gamma_31(i_nu,i_E_rho,j:j), t_clock(i_nu,i_E_rho,j:j) )
                     
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(221)
          
       end do
    end do
    
  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function
    
  end subroutine
  


  subroutine DKE_zeta_Convergence_Example_Big_Resolution_2
    integer, parameter :: N_nu = 5, N_E_rho = 10 
    integer, parameter :: M_zeta = 15, M_theta = 15, M_xi = 6
    real, parameter :: nu_min = 1d-6, nu_max = 1d-4
    integer :: N_theta(M_theta), N_zeta(M_zeta), N_xi(M_xi)
    real :: nu(N_nu), E_rho(N_E_rho)
    real :: Gamma_11(N_nu,N_E_rho, M_zeta), Gamma_31(N_nu,N_E_rho, M_zeta)
    real :: t_clock(N_nu,N_E_rho, M_zeta) 
    
    integer :: i, j, k, i_nu, i_E_rho

    character(len=20) :: nu_char, E_rho_char !,  N_theta_char, N_xi_char 
    character(len=2000) :: folder_path, file_path
    
    ! *** Vector of collisionality and radial electric field
    nu = [( log10(nu_min) + (i-1)*(log10(nu_max)-log10(nu_min))/(N_nu-1),i=1,N_nu)]      
    nu = 10**nu
    nu = [1d-6, 3d-6, 1d-5, 3d-5, 1d-4]
    E_rho = [ 0d0, 0.1e-4, 0.3e-4, 0.1e-3, 0.3e-3, 0.1e-2, 0.3e-2, 0.1e-1, 0.3e-1, 0.1 ]

    ! *** Vectors of theta, zeta and xi sizes
    N_theta = [(  17 +  2*i, i=1,M_theta)]
    N_zeta  = [(  61 +  2*j, j=1,M_zeta)]
    N_xi    = [( 120 + 20*k, k=1,M_xi)]
    
    ! *** Create folder for output
    folder_path = "results/DKE_zeta_Convergence_Big_Resolution_2"
    call execute_command_line( "mkdir " // trim(folder_path) )    
    
    !N_theta_char = int_2_str(N_theta)
    !N_xi_char = int_2_str(N_xi)  
    
    do i_E_rho = 5, N_E_rho
       do i_nu = 3, N_nu
       
          write(   nu_char,'(e9.3)')    nu(i_nu) ; write(*,*) nu_char
          write(E_rho_char,'(e9.3)') E_rho(i_E_rho) ; write(*,*) E_rho_char
          
          ! For each value of collisionality and electric field an output file
          file_path = trim(folder_path) // "/Gamma_11_Gamma_31_nu_" // &
                      trim(nu_char) // "_E_rho_" // trim(E_rho_char)  // ".plt" 
          
          open(221, file=trim(file_path) ) 
          write(221,*) " nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, t_clock [s] "
          do k = 1, M_xi
             do i = 1, M_theta 
                
                ! Write zeta convergence plot.     
                do j = 1, M_zeta
                
                   ! Calculate zeta convergence curves for Gamma_11 and Gamma_31.
                   call DKE_zeta_Convergence( nu(i_nu), E_rho(i_E_rho), N_theta(i), N_zeta(j:j), N_xi(k), &
                        Gamma_11(i_nu,i_E_rho,j:j), Gamma_31(i_nu,i_E_rho,j:j), t_clock(i_nu,i_E_rho,j:j) )
                     
                   write(*,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                   write(221,'(9999e)') nu(i_nu), E_rho(i_E_rho), &
                                       real(N_theta(i)), real(N_zeta(j)), real(N_xi(k)), &
                                       Gamma_11(i_nu,i_E_rho,j), Gamma_31(i_nu,i_E_rho,j), &
                                       t_clock(i_nu,i_E_rho,j)               
                end do 
                
             end do
          end do           
          close(221)
          
       end do
    end do
    
  contains

    function int_2_str(int) result(str)
      integer, intent(in) :: int
      character(len=4) :: str

      if( int < 10 ) write(str, '(I1)') int
      if( 10 <= int .and. int < 100 ) write(str, '(I2)') int
      if( 100 <= int .and. int < 1000 ) write(str, '(I3)') int
      if( 1000 <= int .and. int < 10000 ) write(str, '(I4)') int

    end function
    
  end subroutine
  






end module
