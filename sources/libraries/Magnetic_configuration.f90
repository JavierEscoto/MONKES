module Magnetic_configuration
   use netcdf
   use Lagrange_interpolation
   
   implicit none
   
   private
   
   public :: read_ddkes2data, read_boozer_xform_output
   public :: read_VMEC_output, read_MONKES_input_configuration
   public :: read_boozer_txt
   public :: Select_Surface
   public :: Beidler_NF_2011_Normalization

   
   public :: Np, iota, psi_p, chi_p, B_theta, B_zeta, B00
   public :: Aspect_ratio, Major_Radius, Minor_Radius, s
   
   
   ! *** Functions for magnetic field  
   public :: Boozer_coordinates
   public :: Magnetic_Field_Boozer
   public :: Magnetic_Field_Boozer_Derivative_theta
   public :: Magnetic_Field_Boozer_Derivative_zeta
   
   ! *** Functions for the contravariant and covariant components of the magnetic field 
   public :: Magnetic_Field
   public :: Magnetic_Field_Dv_theta, Magnetic_Field_Dv_zeta
   public :: Magnetic_Field_Pol_Contr, Magnetic_Field_Tor_Contr
   public :: Magnetic_Field_Pol_Cov, Magnetic_Field_Tor_Cov
   public :: Jacobian
   
   logical, save :: Stellarator_symmetry = .true. 
   logical, save :: Boozer_coordinates = .false. 
   integer, save :: Np, N_modes, N_modes_nyq    
   real, save :: iota, psi_p, chi_p, B_theta, B_zeta, B00
   integer, save, allocatable :: mn(:,:), mn_nyq(:,:)  
   real, save, allocatable :: B_modes(:)   
   real, save, allocatable :: B_mnc(:), B_mns(:), g_mnc(:), g_mns(:)   
   real, save, allocatable :: Bu_mnc(:), Bu_mns(:), Bv_mnc(:), Bv_mns(:)    ! Contravariant
   real, save, allocatable :: B_u_mnc(:), B_u_mns(:), B_v_mnc(:), B_v_mns(:)! Covariant  
   real, save, allocatable :: lambda_mns(:)  ! Sine Fourier modes of the lambda function from VMEC
   real, save, allocatable :: r_mnc(:), z_mns(:)  ! Cosine Fourier modes of the flux surface representation in cylindrical coordinates 
   
   real, save :: Aspect_ratio, Major_Radius, Minor_Radius, s
   namelist /surface/s
   
   contains
   
     ! *** Selects the value of s=(r/a)**2=(psi(s)/psi(s=1)) from the input
     ! file "input.surface". By default it selects the middle of the plasma.
     ! Here, psi is the toroidal flux (divided by 2*pi) of the magnetic field 
     ! and s=1 corresponds to the last closed flux-surface. 
   subroutine Select_Surface     
     integer :: ierr 
     
     open(21, file= "monkes_input.surface", status="old", iostat=ierr) 
     if( ierr == 0 ) read(21, nml=surface)     
     close(21) 
     
     if( ierr /= 0 ) s = 0.25    
     
   end subroutine
   
   subroutine read_MONKES_input_configuration(file_path) 
     character(len=*), intent(in) :: file_path
     
     integer, parameter :: M_max = 500, N_max = 500 
     integer :: ierr, i, j, k 
     real :: Bmnc(-N_max:N_max,-M_max:M_max), Bmns(-N_max:N_max,-M_max:M_max)
     
     namelist /magnetic_configuration/Np, iota, B_theta, B_zeta, Bmnc!, Bmns
     
     Stellarator_symmetry = .true. 
     ! *** Read ddkes2.data parameters and magnetic field modes in DKES
     ! format.
     Bmnc = 0    
     open(1, file= file_path, status="old") ; read(1, nml=magnetic_configuration, iostat=ierr)     
     close(1)      
     write(*,*) iota, Np, B_theta, B_zeta 
     ! *** Storing flux-surface parameters  
     B00 = Bmnc(0,0)   ; psi_p = 1 ; chi_p = 1
     
     ! *** Store the magnetic field strength modes in global variables
     N_modes = count( Bmnc /= 0 )   
     if( allocated(B_mnc) ) deallocate(B_mnc, mn)   
     allocate( B_mnc(N_modes), mn(N_modes,2) ) 
     
     k = 1
     do i = -M_max, M_max
        do j = -N_max, N_max
        
           ! Store the non zero modes
           if( Bmnc(i,j) /= 0 ) then 
             mn(k,:) = [ i, j ] 
             B_mnc(k) = Bmnc(i,j)    
                   
             k = k + 1
           end if
        
        end do
     end do      
     
     Aspect_ratio = 1 ; Minor_Radius= 0 ; Major_Radius = 0
     call Write_Magnetic_Configuration
     Boozer_coordinates = .true. 
   end subroutine
      
      
   ! *** Read the DKES input file "ddkes2.data" stored in "file_path"
   ! it requires that the input file gives a namelist called "datain"
   ! It does not have outputs but initializes
   !  
   !  - Np :   number of toroidal periods of the magnetic field
   !  - iota : rotational transform of the magnetic field lines
   !  - psi_p : radial derivative of toroidal flux of the magnetic field
   !  - chi_p : radial derivative of poloidal flux of the magnetic field
   !  - B_theta : covariant poloidal component of the magnetic field (in Boozer angles)
   !  - B_zeta : covariant toroidal component of the magnetic field (in Boozer angles)
   !  - N_modes : total number of Fourier modes of the magnetic field strength
   !              in Boozer angles. 
   !  - B_modes(1:N_modes) : Fourier modes of the magnetic field strength
   !              in Boozer angles. 
   !  - mn(1:N_modes,1) : poloidal mode number "m" of the magnetic field strength
   !              in Boozer angles. 
   !  - mn(1:N_modes,2) : toroidal mode number "n" of the magnetic field strength
   !              in Boozer angles. 
   !
   ! The magnetic field strength is represented as
   !   B(theta,zeta) = \sum_{m,n} B_mn cos( m * theta + n * Np * zeta )
   ! 
   subroutine read_ddkes2data(file_path) 
     character(len=*), intent(in) :: file_path
     
     integer, parameter :: M_max = 500, N_max = 500 
     integer :: ierr, i, j, k
     integer :: nzperiod, lalpha, lfout, nrun, mpolb, ntorb, ibbi     
     real :: cmul, efield, psip, chip, btheta, bzeta, borbi(-N_max:N_max,-M_max:M_max)
     
     namelist /datain/nzperiod, lalpha, lfout, nrun, mpolb, ntorb, ibbi, nrun, &  
              cmul, efield, psip, chip, btheta, bzeta, borbi
     
     ! *** Read ddkes2.data parameters and magnetic field modes in DKES
     ! format.
     borbi = 0    
     open(1, file= file_path, status="old") ; read(1, nml=datain, iostat=ierr)     
     close(1)      
     
     ! *** Storing flux-surface parameters 
     Np = nzperiod ; psi_p = psip ; chi_p = -chip ; iota = chi_p / psi_p
     B_theta = btheta ; B_zeta = bzeta ; B00 = borbi(0,0)
     
     ! *** Store the magnetic field strength modes in global variables
     N_modes = count(borbi /= 0)   
     if( allocated(B_mnc) ) deallocate(B_mnc, mn)   
     allocate( B_mnc(N_modes), mn(N_modes,2) ) 
     
     k = 1
     do i = -M_max, M_max
        do j = -N_max, N_max
        
           ! Store the non zero modes
           if( borbi(j,i) /= 0 ) then 
             mn(k,:) = [ i, j ] 
             B_mnc(k) = borbi(j,i)           
             k = k + 1
           end if
        
        end do
     end do     
     
     ! *** Estimate s, major and minor radius from large-aspect-ratio limit
     Major_Radius = abs(B_zeta/ B00) ; Minor_Radius = abs(psi_p / B00)
     Aspect_ratio = Major_Radius / Minor_Radius 

     write(*,*) " ! *** Approximating Minor_Radius = abs(psi_p / B00)  " 
     write(*,*) " ! *** Approximating Major_Radius = abs(B_zeta/ B00)  " 
       
     
     ! *** Write magnetic configuration in output file
     call Write_Magnetic_Configuration
     Boozer_coordinates = .true.
   end subroutine 
   
   ! *** Writes on the output file the magnetic configuration parameters
   ! that are involved in the monoenergetic calculation or normalization
   ! of quantities. 
   subroutine Write_Magnetic_Configuration
     integer :: k 
     
     write(*,*) " *** Flux-surface parameters "
     write(*,*) "     s = ", s
     write(*,*) "     Number of periods = ", Np
     write(*,*) "     psi_p = ", psi_p
     write(*,*) "     chi_p = ", chi_p
     write(*,*) "     iota = ", iota
     write(*,*) "     B00 = ", B00
     write(*,*) "     B_theta = ", B_theta
     write(*,*) "     B_zeta = ", B_zeta   
     write(*,*) "     B00/(B_zeta+iota*B_theta) = ", B00/(B_zeta+iota*B_theta)  
     write(*,*)
     
     write(*,*) " *** Stellarator parameters "
     write(*,*) "     Aspect_ratio = ", Aspect_ratio
     write(*,*) "     Minor_Radius = ", Minor_Radius
     write(*,*) "     Major_Radius = ", Major_Radius ; write(*,*)

     call Beidler_NF_2011_Normalization
     
     write(*,*) " *** Magnetic field strength Fourier modes "
     write(*,'(9999A25)') "m", "n", "B_mn", "B_mn/B00"
     do k = 1, N_modes 
        write(*,'(I25,I25,e25.16,e25.16)') mn(k,1), mn(k,2), B_mnc(k), B_mnc(k)/B00           
     end do; write(*,*)

     open(21,file="monkes_Magnetic_configuration.dat")
     write(21,*) " *** Flux-surface parameters "
     write(21,*) "     s = ", s
     write(21,*) "     Number of periods = ", Np
     write(21,*) "     psi_p = ", psi_p
     write(21,*) "     chi_p = ", chi_p
     write(21,*) "     iota = ", iota
     write(21,*) "     B00 = ", B00
     write(21,*) "     B_theta = ", B_theta
     write(21,*) "     B_zeta = ", B_zeta  
     write(21,*) "     B00/(B_zeta+iota*B_theta) = ", B00/(B_zeta+iota*B_theta) 
     
     write(21,*) " *** Stellarator parameters "
     write(21,*) "     Aspect_ratio = ", Aspect_ratio
     write(21,*) "     Minor_Radius = ", Minor_Radius
     write(21,*) "     Major_Radius = ", Major_Radius ; write(21,*)
     
     write(21,*) " *** Magnetic field strength Fourier modes "
     write(21,'(9999A25)') "m", "n", "B_mn", "B_mn/B00"
     do k = 1, N_modes
        write(21,'(I25,I25,e25.16,e25.16)') mn(k,:), B_mnc(k), B_mnc(k)/B00        
     end do 
     close(21)
     
   end subroutine 
   
   
   real function Jacobian( theta, zeta ) result(g)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
       
      g = 0
      do k = 1, N_modes 
         m = mn(k,1) ; n = mn(k,2)          
         g = g + g_mnc(k) * cos( m * theta + n *Np * zeta )
         if( .not. Stellarator_symmetry ) & 
         g = g + g_mns(k) * sin( m * theta + n *Np * zeta )
      end do
              
   end function
   
   real function Magnetic_Field_Pol_Contr( theta, zeta ) result(Bu)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
      real :: g 
      
      Bu = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         Bu = Bu + Bu_mnc(k) * cos( m * theta + n *Np * zeta )
         if( .not. Stellarator_symmetry ) & 
         Bu = Bu + Bu_mns(k) * sin( m * theta + n *Np * zeta )
      end do        
      
   end function
   
   real function Magnetic_Field_Tor_Contr( theta, zeta ) result(Bv)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
      real :: g 
      
      Bv = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         Bv = Bv + Bv_mnc(k) * cos( m * theta + n *Np * zeta )
         if( .not. Stellarator_symmetry ) & 
         Bv = Bv + Bv_mns(k) * sin( m * theta + n *Np * zeta )
      end do 
      
   end function
   
   real function Magnetic_Field_Pol_Cov( theta, zeta ) result(B_u)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
       
      B_u = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         B_u = B_u + B_u_mnc(k) * cos( m * theta + n *Np * zeta )
         if( .not. Stellarator_symmetry ) & 
         B_u = B_u + B_u_mns(k) * sin( m * theta + n *Np * zeta )
      end do 
      
   end function
   
   real function Magnetic_Field_Tor_Cov( theta, zeta ) result(B_v)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
       
      B_v = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         B_v = B_v + B_v_mnc(k) * cos( m * theta + n *Np * zeta )
         if( .not. Stellarator_symmetry ) & 
         B_v = B_v + B_v_mns(k) * sin( m * theta + n *Np * zeta )
      end do
         
   end function
   
   real function Magnetic_Field( theta, zeta ) result(B)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
      
      B = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         B = B + B_mnc(k) * cos( m * theta + n *Np * zeta )
         if( .not. Stellarator_symmetry ) & 
         B = B + B_mns(k) * sin( m * theta + n *Np * zeta )
      end do
      
   end function
   
   real function Magnetic_Field_Dv_theta( theta, zeta ) result(B)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
      
      B = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         B = B - m * B_mnc(k) * sin( m * theta + n *Np * zeta )
         if( .not. Stellarator_symmetry ) & 
         B = B + m * B_mns(k) * cos( m * theta + n *Np * zeta )
      end do
      
   end function
   
   real function Magnetic_Field_Dv_zeta( theta, zeta ) result(B)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
      
      B = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         B = B - n*Np * B_mnc(k) * sin( m * theta + n *Np * zeta )
         if( .not. Stellarator_symmetry ) & 
         B = B + n*Np * B_mns(k) * cos( m * theta + n *Np * zeta )
      end do
      
   end function
   
   real function Magnetic_Field_Boozer( theta, zeta ) result(B)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
      
      B = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         B = B + B_modes(k) * cos( m * theta + n *Np * zeta )
      end do
      
   end function
   
   real function Magnetic_Field_Boozer_Derivative_theta( theta, zeta ) result(B)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
      
      B = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         B = B - m * B_modes(k) * sin( m * theta + n * Np * zeta )
      end do
      
   end function
   
   real function Magnetic_Field_Boozer_Derivative_zeta( theta, zeta ) result(B)
      real, intent(in) :: theta, zeta
      
      integer :: m, n, k 
      
      B = 0
      do k = 1, N_modes
         m = mn(k,1) ; n = mn(k,2)          
         B = B - n * Np * B_modes(k) * sin( m * theta + n *Np * zeta )
      end do
      
   end function

  
   
  ! *** Computes the normalization factor to change K_{ij}*\hat{D}_ij
  ! from [Escoto, NF (2024)] to the normalized monoenergetic coefficients
  ! D_ij^* in  [Beidler, NF (2011)]. It also gives the factor to change
  ! from nu/v to adimensional collisionality nu^* and v_E^ = E_r/(B00*v). 
  subroutine Beidler_NF_2011_Normalization( nu_f, E_r_f, D11_f, D31_f )
    real, optional, intent(out) :: nu_f, E_r_f, D11_f, D31_f  
    real, parameter :: pi = acos(-1d0) 
    real:: nu_factor, E_r_factor, D11_factor, D31_factor !, intent(out) 
     
    real :: eps_t, fc 
    
    nu_factor = Major_Radius / abs(iota)
    E_r_factor = 1d0 /( B00  ) 
    D11_factor = 8 * Major_radius * B00**2 * abs(iota) / pi
    
    eps_t = Minor_Radius / Major_Radius
    fc = 1 - 1.46* sqrt(eps_t) 
    D31_factor = 1.5 * abs(iota) * eps_t * B00 / (1-fc)

    if( present(nu_f) )  nu_f = nu_factor
    if( present(E_r_f) ) E_r_f = E_r_factor
    if( present(D11_f) ) D11_f = D11_factor
    if( present(D31_f) ) D31_f = D31_factor

    
    write(*,*) " *** Monoenergetic coefficients normalization parameters "
    write(*,*) "     eps_t = Minor_Radius / Major_Radius = ",  eps_t
    write(*,*) "     fc = 1 - 1.46* sqrt(eps_t) = ",  fc   ; write(*,*)

    write(*,*) " *** Factors relating nu/v and E_r/v to "
    write(*,*) "     nu^*  = nu * Major_radius / (v*iota) = nu_factor * nu /v  and "
    write(*,*) "     v_E^* =  E_r/(B00*v) = E_r_factor * E_rho"
    write(*,*) "     nu_factor = ", nu_factor
    write(*,*) "     E_r_factor = ", E_r_factor ; write(*,*)

    
    write(*,*) " *** Factors relating K_{ij}*\hat{D}_ij to Dij^* from [Beidler, NF (2011)] "
    write(*,*) "     Dij^* = K_{ij}*\hat{D}_ij * Dij_factor "
    write(*,*) "     D11_factor = 8 * Major_radius * B00**2 * abs(iota) / pi ", &
                     D11_factor
    write(*,*) "     D31_factor = 1.5 * Minor_radius * abs(iota) * eps_t * B00 / (1-fc) ", &
                     D31_factor
    write(*,*) 
  end subroutine 
  
    
  subroutine read_boozer_txt(s0)
    real, intent(in) :: s0
    real, parameter :: pi = acos(-1d0), mu0_o_2pi = 2d-7
    
    integer, parameter :: N_file_header = 7, &
                          N_parameters_header = 1, &
                          N_fs_parameters_header = 2, &
                          N_modes_header = 1 
    integer :: i, j, k, N_modes_txt
    character(len=100) :: file_header, parameters_header, modes_header
    integer, allocatable :: mn_v(:,:)
    real, allocatable :: r_mn_v(:,:), z_mn_v(:,:), dphi_mn_v(:,:)
    real, allocatable :: B_mn_v(:,:), s_v(:), iota_v(:), B_theta_v(:), B_zeta_v(:)  
    real, allocatable :: B_mnc_txt(:), mn_txt(:,:) 
    integer :: m_max, n_max, Ns
    real :: tor_flux, a, R, avol, rvol, Vol, dp_ds, sqrt_g_00
    real :: B_mnc_min
    
    ! Read header of version, year, author...
    open(1, file="boozer.txt")
    do i = 1, N_file_header 
       read(1,'(A)') file_header 
       write(*,*) file_header
    end do
    
    ! *** Read m_max, n_max, #surfaces, #field-periods, toroidal flux,
    ! minor radius, major radius, avol, rvol, Vol
    do i = 1, N_parameters_header ! Read header 
       read(1,'(A)') parameters_header  
       write(*,*) parameters_header
    end do
    read(1,*)  m_max, n_max, Ns, Np, tor_flux, a, R, avol, rvol, Vol
    !write(*,*) m_max, n_max, Ns, Np, tor_flux, a, R, avol, rvol, Vol
    
    ! *** Read flux-surface parameters and Fourier modes for each surface
    N_modes_txt = n_max+1 + m_max*(2*n_max+1) ! Number of Fourier modes in the txt
    allocate( s_v(Ns), iota_v(Ns), B_theta_v(Ns), B_zeta_v(Ns) )    
    allocate( mn_v(N_modes_txt,2), B_mn_v(N_modes_txt, Ns) )
    allocate( r_mn_v(N_modes_txt, Ns), z_mn_v(N_modes_txt, Ns), dphi_mn_v(N_modes_txt, Ns) )
    do j = 1, Ns ! loop on flux-surfaces
    
       do i = 1, N_fs_parameters_header  ! Read header 
          read(1,'(A)') parameters_header   ! write(*,*) parameters_header
       end do
       read(1,*) s_v(j), iota_v(j), B_zeta_v(j), B_theta_v(j), dp_ds, sqrt_g_00
       !write(*,*) s_v(j), iota_v(j), B_zeta_v(j), B_theta_v(j), dp_ds, sqrt_g_00
       
       ! *** Read Fourier modes 
       do i = 1, N_modes_header  ! Read header 
          read(1,'(A)') modes_header  ! write(*,*) modes_header 
       end do
       
       do k = 1, N_modes_txt  ! Loop on Fourier modes         
          read(1,*)  mn_v(k,1), mn_v(k,2), r_mn_v(k, j), z_mn_v(k, j), dphi_mn_v(k, j), B_mn_v(k, j)
          !write(*,*) mn_v(k,1), mn_v(k,2), r_mn_v(k, j), z_mn_v(k, j), dphi_mn_v(k, j), B_mn_v(k, j)       
       end do 
       
    end do
    close(1) 

    ! Interpolate quantities to given surface s0 and adapt normalizations
    ! and definitions.
    Minor_Radius = a ; Major_Radius = R ; Aspect_ratio = R/a
    iota  = -Linear_Interpolation( s0, s_v, iota_v )  
    B_theta  = Linear_Interpolation( s0, s_v, B_theta_v )  
    B_zeta  = Linear_Interpolation( s0, s_v, B_zeta_v )  
    B_theta=-B_theta * Np* mu0_o_2pi
    B_zeta=B_zeta * Np * mu0_o_2pi
    tor_flux = tor_flux/(2*pi)
    psi_p = 2*tor_flux*sqrt(s0)/a
    chi_p = iota * psi_p
    
    allocate( B_mnc_txt(N_modes_txt), mn_txt(N_modes_txt,2) )
    
    do k = 1, N_modes_txt 
       mn_txt(k,:) = mn_v(k,:) 
       B_mnc_txt(k) = Linear_Interpolation( s0, s_v, B_mn_v(k, :) )  
       if( mn_txt(k,1) == 0 .and. mn_txt(k,2) == 0 ) B00 = B_mnc_txt(k) 
    end do  

    ! Filter non relevant Fourier modes if necessary
    if( N_modes_txt <= 150 ) B_mnc_min = 0
    if( N_modes_txt >  150 ) B_mnc_min = 1d-5 * B00
    
    N_modes = count( abs(B_mnc_txt) >= B_mnc_min ) 
    allocate( B_mnc(N_modes), mn(N_modes,2) )
    B_mnc = pack( B_mnc_txt, abs(B_mnc_txt) >= B_mnc_min ) 
    mn(:,1) = pack( mn_txt(:,1), abs(B_mnc_txt) >= B_mnc_min ) 
    mn(:,2) = pack( mn_txt(:,2), abs(B_mnc_txt) >= B_mnc_min ) 
    call Write_Magnetic_Configuration
    Boozer_coordinates = .true. ! Set logical for Boozer coordinates as true
  
  end subroutine 


  subroutine read_boozer_xform_output(s0)
    real, intent(in) :: s0
    
    real, parameter :: pi = acos(-1d0)
    integer :: Ns_b, Ns ! Ns_b: Number of surfaces with Boozer, Ns: Total number of surfaces
    integer :: ncid, ierr_netcdf ! integer for knowing if we are reading correctly the file
    integer :: mnboz_b ! Number of (Boozer) Fourier modes in the boozer xform output
    real, allocatable :: s_b(:), ss(:), ss_half(:) ! Surfaces grid
    real, allocatable :: iota_s(:), psi_s(:), psi_p_s(:), B_theta_s(:), B_zeta_s(:)   
    real, allocatable :: B_mnc_s(:,:),  B_mns_s(:,:), BB_mnc(:), BB_mns(:)  ! Cosine and Sine Fourier amplitudes
    real, allocatable :: r_mnc_s(:,:),  r_mns_s(:,:) 
    integer, allocatable :: i_b(:), mn_s(:,:)  
    real :: B_mn_min = 5d-6   
    logical, allocatable :: bigger_earth_field(:)
    real :: psi_LCFS
    integer :: sign_lh = 1!, q_Ns_b, q_Ns      
    
    ! Open BOOZER_XFORM output file "boozmn.nc"
    ierr_netcdf = nf90_open('boozmn.nc', nf90_nowrite, ncid )
    
    call read_scalar_quantities      
    call read_profiles     
    call read_magnetic_field_modes    
    
    ! *** Global module quantities for the flux-surface: rho = r/a
    !     B_theta(rho), B_zeta(rho), psi'(rho), chi'(rho) and iota(rho).
    ! If the reference system is left-handed, change the sign of 
    ! B_theta, chi_p and iota to make it right-handed (change theta by -theta).   
    
    psi_LCFS = psi_s(Ns) / (2 * pi)
    psi_p = 2 * psi_LCFS * sqrt(s0) / Minor_Radius   
    
    sign_lh = -1 ! Sign to change from left-handed to right-handed coordinates 
    iota  = sign_lh * Linear_Interpolation( s0, ss_half(2:Ns), iota_s(2:Ns) )  
    chi_p = sign_lh * iota * psi_p 
       
    B_theta = sign_lh * Linear_Interpolation(s0, ss_half(2:Ns), B_theta_s(2:Ns) )  
    B_zeta  =           Linear_Interpolation(s0, ss_half(2:Ns),  B_zeta_s(2:Ns) ) 
    s = s0	
    
    ! Close BOOZER_XFORM output file "boozmn.nc"
    ierr_netcdf = nf90_close(ncid)     
         
    call Write_Magnetic_Configuration
    Boozer_coordinates = .true. 
    
    contains 
    
    ! Reads the number of field periods, the radial grid sizes 
    subroutine read_scalar_quantities
       integer :: rhid, idimid
       integer jsize 
       
       ! Read No. Field periods and write it on global "Np"
       ierr_netcdf = nf90_inq_varid( ncid, 'nfp_b', rhid )  
       ierr_netcdf = nf90_get_var( ncid, rhid, Np )      
       write(*,*) " Number of field periods ", Np 
       
       ! Read aspect ratio and write it on global "Aspect_ratio"
       ierr_netcdf = nf90_inq_varid( ncid, 'aspect_b', rhid )   
       ierr_netcdf = nf90_get_var( ncid, rhid, Aspect_ratio ) 
            
       ! Read total number of surfaces in original VMEC file
       ierr_netcdf = nf90_inq_varid( ncid, 'ns_b', rhid )    
       ierr_netcdf = nf90_get_var( ncid, rhid, Ns )
       write(*,*) " Number of total surfaces in VMEC file ", Ns           
        
       ! Read number of surfaces of BOOZER_XFORM
       ierr_netcdf = nf90_inq_dimid( ncid, 'pack_rad', idimid) 
       ierr_netcdf = nf90_inquire_dimension( ncid, idimid, len=Ns_b )
       write(*,*) " Number of surfaces of BOOZER_XFORM ", Ns_b    
       
       ! Set orders of interpolation 
      !  q_Ns_b = minval( [ Ns_b-1,  2] )
      !  q_Ns = minval( [ Ns-2,  2 ] )
       
       ! Read maximum number of (Boozer) Fourier modes 
       ierr_netcdf = nf90_inq_varid( ncid, 'mnboz_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, mnboz_b)
       write(*,*) " Number of (Boozer) Fourier modes  of BOOZER_XFORM ", mnboz_b
                
    end subroutine 
    
    subroutine read_profiles
       integer :: rhid, idimid
       integer :: i
       
       real, allocatable :: rmnc_b(:,:), rmns_b(:,:)
       
       ! Read index vector for surface position in which BOOZER_XFORM 
       ! has done calculations
       allocate( i_b(Ns_b) )  
       ierr_netcdf = nf90_inq_varid( ncid, 'jlist', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, i_b )   
       
       ! Vector of of 2pi*psi(s) (phi_b=2pi*psi) in VMEC radial positions
       allocate( psi_s(Ns) ) ; psi_s(1) = 0        
       ierr_netcdf = nf90_inq_varid( ncid, 'phi_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, psi_s )   
       if( maxval(abs(psi_s))==minval(abs(psi_s)) ) psi_s = [((i-1)*1d0/(Ns-1), i=1, Ns)]! Done because python's booz xform sometimes does not calculate phi
       
       ! Radial grids 
       allocate( ss(Ns), ss_half(2:Ns), s_b(Ns_b) )        
       ss = psi_s/psi_s(Ns)! Vector of VMEC radial positions 
       ss_half =( ss(1:Ns-1) + ss(2:Ns))/2
       s_b = (ss(i_b) + ss(i_b-1))/2 ! Vector of surfaces in which Boozer transform has been done
              
       ! Vector of rotational transform in VMEC radial positions
       allocate( iota_s(1:Ns) ) ; iota_s(1) = 0
       ierr_netcdf = nf90_inq_varid(ncid,'iota_b',rhid) 
       ierr_netcdf = nf90_get_var( ncid, rhid, iota_s )   
                     
                     
       ! Vector of derivative along s of psi'(s) (tor. flux=psi/2pi) in VMEC radial positions
       allocate( psi_p_s(Ns) ) ; psi_p_s(1) = 0        
       ierr_netcdf = nf90_inq_varid( ncid, 'phip_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, psi_p_s )   
       if( maxval(abs(psi_p_s))== 0 ) psi_p_s = 1
       
       ! Poloidal current B_zeta 
       allocate( B_theta_s(Ns) ) ; B_theta_s(1) = 0   
       ierr_netcdf = nf90_inq_varid( ncid, 'buco_b', rhid )
       ierr_netcdf = nf90_get_var( ncid, rhid, B_theta_s )
       
       ! Toroidal current B_zeta    
       allocate( B_zeta_s(Ns) ) ; B_zeta_s(1) = 0    
       ierr_netcdf = nf90_inq_varid( ncid, 'bvco_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, B_zeta_s )
       
       open(21, file="monkes_boozerxform_profile.dat")       
       write(21,'(9999A25)') "s", "iota", "2pi*psi(s)", "psi'(s)", "B_theta_s(s)", "B_zeta_s(s)" 
       do i = 2, Ns 
          write(21,'(9999e25.16)') ss_half(i), iota_s(i), ss_half(i)*psi_s(Ns), psi_p_s(i), B_theta_s(i), B_zeta_s(i) 
       end do         
       close(21) 
       
    end subroutine 
    
    subroutine read_magnetic_field_modes
       integer :: rhid, idimid
       integer :: i, j
       
       ! Poloidal (m) and toroidal (n) mode numbers  
       allocate( mn_s(mnboz_b,2) )   
       ierr_netcdf = nf90_inq_varid( ncid, 'ixm_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, mn_s(:,1) )
					 
       ierr_netcdf = nf90_inq_varid( ncid, 'ixn_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, mn_s(:,2) )
       mn_s(:,2) = mn_s(:,2) / Np ! Projecting to one period
       
       ! B modes B = B_mnc_s * cos + B_mns_s * sin
       allocate( B_mnc_s(mnboz_b, Ns_b) )  
       ierr_netcdf = nf90_inq_varid( ncid, 'bmnc_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, B_mnc_s )   
          
       ! *** Non Stellarator-symmetric modes (if needed)       
       ! Detect stellarator symmetry 
       ierr_netcdf = nf90_inq_varid( ncid, 'bmns_b', rhid )  
       Stellarator_symmetry = ( ierr_netcdf /= 0 )
       if( .not. Stellarator_symmetry ) then
         allocate( B_mns_s(mnboz_b, Ns_b) )          
         ierr_netcdf = nf90_inq_varid( ncid, 'bmns_b', rhid ) 
         ierr_netcdf = nf90_get_var( ncid, rhid, B_mns_s )       
         
         allocate( BB_mns(mnboz_b) )   
         do i = 1, mnboz_b                             
            BB_mns(i) = Linear_Interpolation(s0, s_b(1:Ns_b), B_mns_s(i,1:Ns_b) )                      
         end do  
         
       end if 
          
       ! Surface
       allocate( r_mnc_s(mnboz_b, Ns_b), r_mns_s(mnboz_b, Ns_b) )  
       ierr_netcdf = nf90_inq_varid( ncid, 'rmnc_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, r_mnc_s )
       ierr_netcdf = nf90_inq_varid( ncid, 'rmns_b', rhid ) 
       ierr_netcdf = nf90_get_var( ncid, rhid, r_mns_s ) 
       
       ! Major and minor radius of the device      
       do i = 1, mnboz_b 
          if( mn_s(i,1) == 0 .and. mn_s(i,2) == 0 ) &
            Major_Radius = r_mnc_s(i,2)       
       end do 
       Minor_Radius = Major_Radius / Aspect_ratio     
         
       ! *** Interpolated Fourier modes of the magnetic field strength at s0 and B_00
       allocate( BB_mnc(mnboz_b), bigger_earth_field(mnboz_b) )   
       do i = 1, mnboz_b     
          BB_mnc(i) = Linear_Interpolation(s0, s_b(1:Ns_b), B_mnc_s(i,1:Ns_b) )   
                     
          if( mn_s(i,1) == 0 .and. mn_s(i,2) == 0 ) &
            B00 = BB_mnc(i)  
       end do        
        
       ! *** Truncate magnetic field modes 
       B_mn_min = 1d-5 * B00
       bigger_earth_field = abs(BB_mnc)  > B_mn_min
       N_modes = count( bigger_earth_field )   
       write(*,*) " Truncating modes of B smaller than ", B_mn_min, " Tesla " 
			   
       if( allocated(B_mnc) ) deallocate( B_mnc, mn ) 
       allocate( B_mnc(N_modes), mn(N_modes,2) )
       B_mnc = pack( BB_mnc, bigger_earth_field )
       mn(:,1) = pack( mn_s(:,1), bigger_earth_field )
       mn(:,2) = pack( mn_s(:,2), bigger_earth_field )	
       
       if( .not. Stellarator_symmetry ) then
         if( allocated(B_mns) ) deallocate( B_mns )
         B_mns = pack( BB_mns, bigger_earth_field )       
       end if 
	 
	   open(21,file="boozxform_B_modes.plt")
	   write(21,*) " Number of surfaces = ", Ns_b
	   do j = 1, Ns_b
	      write(21,*) j, " Surface s= ", s_b(j)
	      do i = 1, mnboz_b
	         if( bigger_earth_field(i) .and. .not. Stellarator_symmetry ) & 
	         write(21,'(9999e25.16)') 1d0*j, 1d0*mn_s(i,:), B_mnc_s(i,j), B_mns_s(i,j)	
	         if( bigger_earth_field(i) .and. Stellarator_symmetry ) & 
	         write(21,'(9999e25.16)') 1d0*j, 1d0*mn_s(i,:), B_mnc_s(i,j) 	      
	      end do
	   end do 
	   close(21)	
	    	  
    end subroutine 
    
  end subroutine  

  
  subroutine read_VMEC_output(s0)
    real, intent(in) :: s0
    
    real, parameter :: pi = acos(-1d0)
    integer :: Ns ! Ns: Total number of surfaces
    integer :: ncid, ierr_netcdf ! integer for knowing if we are reading correctly the file
    integer :: mn_b ! Number of (Boozer) Fourier modes
    real, allocatable :: ss(:), ss_half(:) ! Surfaces grids
    real, allocatable :: iota_s(:), psi_s(:), psi_p_s(:), B_theta_s(:), B_zeta_s(:)   
    real, allocatable :: g_mnc_s(:,:),  g_mns_s(:,:)  ! Cosine and Sine Fourier amplitudes
    real, allocatable :: B_mnc_s(:,:),  B_mns_s(:,:)  ! Cosine and Sine Fourier amplitudes
    real, allocatable :: B_u_mnc_s(:,:),  B_v_mnc_s(:,:) 
    real, allocatable :: B_u_mns_s(:,:),  B_v_mns_s(:,:)  
    real, allocatable :: Bu_mnc_s(:,:),  Bv_mnc_s(:,:) 
    real, allocatable :: Bu_mns_s(:,:),  Bv_mns_s(:,:)  
    real, allocatable :: r_mnc_s(:,:),  r_mns_s(:,:) ! Modes of r(psi,u,v)
    real, allocatable :: z_mnc_s(:,:),  z_mns_s(:,:) ! Modes of z(psi,u,v)
    real, allocatable :: lambda_mnc_s(:,:),  lambda_mns_s(:,:) ! Modes of lambda(psi,u,v)
    real, allocatable :: Tor_flux(:) ! Toroidal flux from VMEC.   
    
    integer :: mpol, ntor, rhid, idimid 
    
    ! Open VMEC output file "VMEC.nc"
    !ierr_netcdf = nf90_open('VMEC.nc', nf90_nowrite, ncid )
        
    psi_p = 1 ;  chi_p = 1 ; B_theta = 1 ; B_zeta = 1;
    Aspect_ratio = 1 ; Minor_Radius = 1 ; Major_Radius = 1 ; 
    
    call read_scalar_quantities
    call read_profiles
   !  call Write_Modes  
    call Interpolate_magnetic_configuration 
     
    
    call Write_Magnetic_Configuration
   !  call Write_quantities_on_FS
    
    contains 
    subroutine Write_Modes
      integer :: i, j, jj 


      open(1234, file="VMEC_modes.dat")
      do j = 2, Ns
         write(1234,*) " *** Surface s = ", ss(j) 
         write(1234,'(9999A25)') "m", "n", "r_mn", "z_mn"
         do i = 1, N_modes
            write(1234,'(9999e25.16)') &
            1d0*mn(i,1), 1d0*mn(i,2), r_mnc_s(i,j), z_mns_s(i,j)
         end do       
      end do       
      close(1234) 

      
      open(1234, file="VMEC_modes_nyquist.dat")
      do j = 2, Ns ; jj = j + 1         
         write(1234,*) " *** Surface s = ", ss_half(j) 
         write(1234,'(9999A25)') "m", "n", "B_mnc", "g_mnc" , "B_u_mnc", "B_v_mnc", "Bu_mnc", "Bv_mnc"   
         do i = 1, N_modes_nyq
            write(1234,'(9999e25.16)') &
            1d0*mn_nyq(i,1), 1d0*mn_nyq(i,2), B_mnc_s(i,j), g_mnc_s(i,j) , &
            B_u_mnc_s(i,j), B_v_mnc_s(i,j), Bu_mnc_s(i,j), Bv_mnc_s(i,j)  
         end do       
      end do       
      close(1234) 
    
    
    end subroutine 

    subroutine Write_quantities_on_FS
       integer, parameter :: N_theta = 40, N_zeta = 40 
       integer :: i, j
       real :: theta(0:N_theta), zeta(0:N_zeta), B(0:N_theta, 0:N_zeta), g(0:N_theta, 0:N_zeta)
       real :: lambda(0:N_theta, 0:N_zeta), lambda_theta(0:N_theta, 0:N_zeta), lambda_zeta(0:N_theta, 0:N_zeta)
       real :: r(0:N_theta, 0:N_zeta), z(0:N_theta, 0:N_zeta)
       
       theta = [(i*2*pi/N_theta, i=0, N_theta)]
       zeta  = [(j*2*pi/(N_zeta*Np), j=0, N_zeta)]

       do j = 0, N_zeta
          do i = 0, N_theta 
             B(i,j) =  Fourier_expansion_cos_sin( B_mnc, 0*B_mnc, theta(i), zeta(j) )
             g(i,j) = Fourier_expansion_cos_sin( g_mnc, 0*g_mnc, theta(i), zeta(j) )

             
             lambda(i,j)       = Fourier_expansion_cos_sin( 0*lambda_mns, lambda_mns, theta(i), zeta(j) )
             lambda_theta(i,j) = Fourier_expansion_cos_sin_DV_theta( 0*lambda_mns, lambda_mns, theta(i), zeta(j) )
             lambda_zeta(i,j)  = Fourier_expansion_cos_sin_DV_zeta( 0*lambda_mns, lambda_mns, theta(i), zeta(j) )  
             r(i,j) = Fourier_expansion_cos_sin( r_mnc, 0*r_mnc, theta(i), zeta(j) ) 
             z(i,j) = Fourier_expansion_cos_sin( 0*z_mns, z_mns, theta(i), zeta(j) ) 
          end do 
       end do
       
       open(1234, file="Quantities_on_FS.dat")
       write(1234,'(9999A25)') "theta", "zeta", "B", "Jacobian", "lambda", "lambda_theta", "lambda_zeta"  
       do j = 0, N_zeta
          do i = 0, N_theta             
             write(1234,'(9999e25.16)') theta(i), zeta(j), B(i,j), g(i,j), lambda(i,j), lambda_theta(i,j), lambda_zeta(i,j) 
          end do 
       end do
      close(1234)
      
      open(1234, file="Cross_sections.dat")
      write(1234,'(9999A25)') "r", "z", "zeta / (2pi/N_fp)"
      do j = 0, N_zeta
         write(1234,*) "zeta / (2pi/N_fp)= ", zeta(j)/(2*pi/Np)
         do i = 0, N_theta             
            write(1234,'(9999e25.16)') r(i,j), z(i,j), zeta(j)/(2*pi/Np)
         end do  
      end do
      close(1234)

      open(1234, file="Flux_surface.dat")
      write(1234,'(9999A25)') "x","y", "z"
      do j = 0, N_zeta 
         do i = 0, N_theta             
            write(1234,'(9999e25.16)') r(i,j) * cos(zeta(j)), r(i,j) * sin(zeta(j)), z(i,j)
         end do  
      end do
      close(1234)



    end subroutine 

    subroutine Interpolate_magnetic_configuration
       integer :: i
       real :: ds_dr, dpsi_ds, psi_LCFS

       ! Rotational transform and toroidal and poloidal flux derivatives at s0
       iota = Interpolated_value( s0, ss(2:Ns), iota_s(2:Ns), 2 ) 
       psi_LCFS = Tor_flux(Ns) / (2 * pi) ! dpsi/ds = psi_LCFS
       dpsi_ds = psi_LCFS
       ds_dr = 2 * sqrt(s0) / Minor_Radius 
       
       psi_p = psi_LCFS * ds_dr
       chi_p = iota * psi_p

       ! *** Stellarator symmetric modes
       allocate( B_mnc(N_modes), g_mnc(N_modes) ) 
       allocate( B_u_mnc(N_modes), B_v_mnc(N_modes) ) 
       allocate( Bu_mnc(N_modes), Bv_mnc(N_modes) )  
       do i = 1, N_modes
           
          ! Cosine Fourier modes of magnetic field strength
          B_mnc(i) = Interpolated_value( s0, ss_half(2:Ns), B_mnc_s(i,2:Ns), 2 )   
          if( mn(i,1)==0 .and. mn(i,2)==0 ) B00 = B_mnc(i)            
          
          ! Cosine Fourier modes of covariant components of the Jacobian
          ! It must be corrected to change from coordinate s=psi/psi_LCFS to psi to be consistent with the DKE formulation for Er /= 0. 
          g_mnc(i) = - Interpolated_value( s0, ss_half(2:Ns), g_mnc_s(i,2:Ns), 2 ) / dpsi_ds 
          
          ! Cosine Fourier modes of covariant components of B
          B_u_mnc(i) =  Interpolated_value( s0, ss_half(2:Ns), B_u_mnc_s(i,2:Ns), 2 )   
          B_v_mnc(i) = -Interpolated_value( s0, ss_half(2:Ns), B_v_mnc_s(i,2:Ns), 2 )  
          if( mn(i,1)==0 .and. mn(i,2)==0 ) B_theta = B_u_mnc(i)  
          if( mn(i,1)==0 .and. mn(i,2)==0 ) B_zeta = B_v_mnc(i)
          
          ! Cosine Fourier modes of contravariant components of B
          Bu_mnc(i) =  Interpolated_value( s0, ss_half(2:Ns), Bu_mnc_s(i,2:Ns), 2 ) 
          Bv_mnc(i) = - Interpolated_value( s0, ss_half(2:Ns), Bv_mnc_s(i,2:Ns), 2 )
             
       end do  
       
       ! *** Non-Stellarator symmetric modes 
       if(.not. Stellarator_symmetry) then  
         allocate( B_mns(N_modes), g_mns(N_modes) ) 
         allocate( B_u_mns(N_modes), B_v_mns(N_modes) ) 
         allocate( Bu_mns(N_modes), Bv_mns(N_modes) ) 
         do i = 1, N_modes
            ! Sine Fourier modes of magnetic field strength
            B_mns(i) = Interpolated_value( s0, ss_half(2:Ns), B_mns_s(i,2:Ns), 2 )
            
            ! Sine Fourier modes of covariant components of the Jacobian
            g_mns(i) = Interpolated_value( s0, ss_half(2:Ns), g_mns_s(i,2:Ns), 2 ) 
          
            ! Sine Fourier modes of covariant components of B
            B_u_mns(i) = Interpolated_value( s0, ss_half(2:Ns), B_u_mns_s(i,2:Ns), 2 )
            B_v_mns(i) = Interpolated_value( s0, ss_half(2:Ns), B_v_mns_s(i,2:Ns), 2 )
          
            ! Sine Fourier modes of contravariant components of B
            Bu_mns(i) = Interpolated_value( s0, ss_half(2:Ns), Bu_mns_s(i,2:Ns), 2 )
            Bv_mns(i) = Interpolated_value( s0, ss_half(2:Ns), Bv_mns_s(i,2:Ns), 2 )
            
         end do 
       end if
    
    end subroutine 
    
    ! Reads the number of field periods, the radial grid sizes 
    subroutine read_scalar_quantities 
       
       ! Read No. Field periods and write it on global "Np" 
       call read_netcdf_integer('nfp', Np)
       write(*,*) " Number of field periods ", Np 
        
       ! Read aspect ratio and write it on global "Aspect_ratio" 
       call read_netcdf_scalar('aspect', Aspect_ratio)    
       write(*,*) " Aspect ratio ", Aspect_ratio 
       call read_netcdf_scalar('Aminor_p', Minor_Radius)
       write(*,*) " Minor radius ", Minor_Radius 
       call read_netcdf_scalar('Rmajor_p', Major_Radius)
       write(*,*) " Major radius ", Major_Radius
               
       ! Read number of surfaces and write it on "Ns" 
       call read_netcdf_integer('ns', Ns)
       write(*,*) " Number of total surfaces in VMEC file ", Ns                 
        

       call read_netcdf_integer('mpol', mpol)
       write(*,*) " Number of poloidal Fourier modes  of VMEC ", mpol
       
       ! Read maximum number of poloidal Fourier modes  
       call read_netcdf_integer('ntor', ntor)
       write(*,*) " Number of poloidal Fourier modes  of VMEC ", ntor
       
      ! Read number of total Fourier modes   
       call read_netcdf_integer('mnmax_nyq', N_modes)
       write(*,*) " Number of total Fourier modes  of VMEC ", N_modes
      
       call read_netcdf_scalar('Aminor_p', Minor_Radius)
       call read_netcdf_scalar('Rmajor_p', Major_Radius)

    end subroutine 
    
    subroutine read_profiles   
       integer :: i, j 

       real, allocatable :: mn_temp(:,:), mn_temp_nyq(:,:)
           
       ! *** Flux-surface label as normalized Toroidal flux s = psi/psi_LCFS 
       allocate( ss(Ns), ss_half(Ns) )        
       call read_netcdf_vector('phi', ss) ! In vmec output psi is denoted phi
       ss = ss/ss(Ns)  
       ss_half(2:Ns) = ( ss(2:Ns) + ss(1:Ns-1) ) /2 ; ss_half(1) = 0! VMEC's half mesh
       
       allocate( Tor_flux(Ns) )
       call read_netcdf_vector('phi', Tor_flux)     
 
           
       ! *** Rotational transform in full mesh
       allocate( iota_s(Ns) )        
       call read_netcdf_vector('iotaf', iota_s)               
       
       ! *** Integer map of the Fourier modes (m,n) = ( xm(i), xn(i) ) for
       ! 1 <= i <= N_modes. It is the same for all surfaces.            
       allocate( mn(N_modes,2), mn_temp(N_modes,2) )     
       call read_netcdf_vector('xm_nyq', mn_temp(:,1))
       call read_netcdf_vector('xn_nyq', mn_temp(:,2))   
       
       mn = nint( mn_temp ) ! The modes are read as real and we need to convert them to integer
       mn(:,2) = mn(:,2) / Np          
        
        
       ! *** Stellarator symmetric modes          
       ! Cosine Fourier modes of the magnetic field strength on half mesh.
       ! The first point does not count
       allocate( B_mnc_s(N_modes,Ns) )             
       call read_netcdf_array('bmnc', B_mnc_s)

       ! Cosine Fourier modes of the Jacobian on half mesh.
       ! The first point does not count
       allocate( g_mnc_s(N_modes,Ns) )        
       call read_netcdf_array('gmnc', g_mnc_s)
       
       ! Cosine Fourier modes of the poloidal B_u and toroidal B_v covariant
       ! components of the magnetic fieldon half mesh.
       ! The first point does not count
       allocate( B_u_mnc_s(N_modes,Ns), B_v_mnc_s(N_modes,Ns) )   
       call read_netcdf_array('bsubumnc', B_u_mnc_s)
       call read_netcdf_array('bsubvmnc', B_v_mnc_s)
               
       
       
       ! Cosine Fourier modes of the poloidal B_u and toroidal B_v contravariant
       ! components of the magnetic field
       allocate( Bu_mnc_s(N_modes,Ns), Bv_mnc_s(N_modes,Ns) )   
       call read_netcdf_array('bsupumnc', Bu_mnc_s)
       call read_netcdf_array('bsupvmnc', Bv_mnc_s)
 
       
       
       ! *** Non Stellarator-symmetric modes (if needed)       
       ! Detect stellarator symmetry 
        
       ierr_netcdf = nf90_inq_varid( ncid, 'bmns', rhid )  
       Stellarator_symmetry = ( ierr_netcdf /= 0 )  
       write(*,*) " Stellarator_symmetry ", Stellarator_symmetry
       
       if( .not. Stellarator_symmetry ) then 
       
         ! Sine Fourier modes of the magnetic field strength 
         allocate( B_mns_s(N_modes,Ns) )        
         ierr_netcdf = nf90_inq_varid( ncid, 'bmns', rhid )  
         ierr_netcdf = nf90_get_var( ncid, rhid, B_mns_s) 
         
         ! Sine Fourier modes of the Jacobian
         allocate( g_mns_s(N_modes,Ns) )        
         ierr_netcdf = nf90_inq_varid( ncid, 'gmns', rhid )  
         ierr_netcdf = nf90_get_var( ncid, rhid, g_mns_s)
         
         ! Sine Fourier modes of the poloidal B_u and toroidal B_v covariant
         ! components of the magnetic field
         allocate( B_u_mns_s(N_modes,Ns), B_v_mns_s(N_modes,Ns) )        
         ierr_netcdf = nf90_inq_varid( ncid, 'bsubumns', rhid )  
         ierr_netcdf = nf90_get_var( ncid, rhid, B_u_mns_s) 
            
         ierr_netcdf = nf90_inq_varid( ncid, 'bsubvmns', rhid ) 
         ierr_netcdf = nf90_get_var( ncid, rhid, B_v_mns_s) 
         
         ! *** Sine Fourier modes of the poloidal B_u and toroidal B_v contravariant
         ! components of the magnetic field
         allocate( Bu_mns_s(N_modes,Ns), Bv_mns_s(N_modes,Ns) )        
         ierr_netcdf = nf90_inq_varid( ncid, 'bsupumns', rhid )  
         ierr_netcdf = nf90_get_var( ncid, rhid, Bu_mns_s) 
            
         ierr_netcdf = nf90_inq_varid( ncid, 'bsupvmns', rhid )  
         ierr_netcdf = nf90_get_var( ncid, rhid, Bv_mns_s)   
         
       end if  
        
!~        do j = 1, Ns 
!~           write(*,*) " Magnetic field modes of Surface s = ", ss(j) 
!~           write(*,*) "m", "n", "B_mnc", "B_mns"
!~           do i = 1, N_modes           
!~              if( .not. Stellarator_symmetry ) write(*,*) mn(i,:), B_mnc_s(i,j), B_mns_s(i,j) 
!~              if( Stellarator_symmetry ) write(*,*) mn(i,:), B_mnc_s(i,j) 
!~           end do 
!~        end do                 
       
!~        do j = 1, Ns 
!~           write(*,*) " Jacobian modes of Surface s = ", ss(j) 
!~           write(*,*) "m", "n", "g_mnc", "g_mns"
!~           do i = 1, N_modes           
!~              if( .not. Stellarator_symmetry ) write(*,*) mn(i,:), g_mnc_s(i,j), g_mns_s(i,j) 
!~              if( Stellarator_symmetry ) write(*,*) mn(i,:), g_mnc_s(i,j)
!~           end do 
!~        end do        
       
!~        do j = 1, Ns 
!~           write(*,*) " Modes of the (poloidal) covariant component at Surface s = ", ss(j) 
!~           write(*,*) "m", "n", "B_u_mnc", "B_u_mns"
!~           do i = 1, N_modes           
!~              if( .not. Stellarator_symmetry ) write(*,*) mn(i,:), B_u_mnc_s(i,j), B_u_mns_s(i,j) 
!~              if( Stellarator_symmetry ) write(*,*) mn(i,:), B_u_mnc_s(i,j)
!~           end do 
!~        end do        
       
!~        do j = 1, Ns 
!~           write(*,*) " Modes of the (toroidal) covariant component at Surface s = ", ss(j) 
!~           write(*,*) "m", "n", "B_v_mnc", "B_v_mns"
!~           do i = 1, N_modes           
!~              if( .not. Stellarator_symmetry ) write(*,*) mn(i,:), B_v_mnc_s(i,j), B_v_mns_s(i,j) 
!~              if( Stellarator_symmetry ) write(*,*) mn(i,:), B_v_mnc_s(i,j)  
!~           end do 
!~        end do        
       
!~        do j = 1, Ns 
!~           write(*,*) " Modes of the (poloidal) contravariant component at Surface s = ", ss(j) 
!~           write(*,*) "m", "n", "B_u_mnc", "B_u_mns"
!~           do i = 1, N_modes           
!~              if( .not. Stellarator_symmetry ) write(*,*) mn(i,:), Bu_mnc_s(i,j), Bu_mns_s(i,j) 
!~              if( Stellarator_symmetry ) write(*,*) mn(i,:), Bu_mnc_s(i,j) 
!~           end do 
!~        end do        
       
!~        do j = 1, Ns 
!~           write(*,*) " Modes of the (toroidal) contravariant component at Surface s = ", ss(j) 
!~           write(*,*) "m", "n", "B_v_mnc", "B_v_mns"
!~           do i = 1, N_modes           
!~              if( .not. Stellarator_symmetry ) write(*,*) mn(i,:), Bv_mnc_s(i,j), Bv_mns_s(i,j) 
!~              if( Stellarator_symmetry ) write(*,*) mn(i,:), Bv_mnc_s(i,j) 
!~           end do 
!~        end do        
       
    end subroutine 


  end subroutine
  
  
    real function Fourier_expansion_cos_sin( F_mn_c, F_mn_s, theta, zeta ) result(F)
       real, intent(in) :: F_mn_c(:), F_mn_s(:), theta, zeta

       integer :: i, m, n, N_mn, mn_index(size(F_mn_c),2)
       
       N_mn = size(F_mn_c) ! Choose if nyquist modes or not
       if( N_mn == N_modes )     mn_index = mn
       if( N_mn == N_modes_nyq ) mn_index = mn_nyq

       F = 0 
       do i = 1, N_mn ; 
         
         m = mn_index(i,1) ; n = mn_index(i,2)
          F = F &
            + F_mn_c(i) * cos( m*theta - Np*n*zeta ) &
            + F_mn_s(i) * sin( m*theta - Np*n*zeta ) 
       end do 
    end function
   
    real function Fourier_expansion_cos_sin_DV_theta( F_mn_c, F_mn_s, theta, zeta ) result(F)
       real, intent(in) :: F_mn_c(:), F_mn_s(:), theta, zeta
       integer :: i, m, n, N_mn, mn_index(size(F_mn_c),2)
       
       N_mn = size(F_mn_c) ! Choose if nyquist modes or not
       if( N_mn == N_modes )     mn_index = mn
       if( N_mn == N_modes_nyq ) mn_index = mn_nyq

       F = 0 
       do i = 1, N_mn ;  m = mn_index(i,1) ; n = mn_index(i,2)

          F = F &
            - m* F_mn_c(i) * sin( m*theta - Np*n*zeta ) &
            + m* F_mn_s(i) * cos( m*theta - Np*n*zeta ) 
       end do 
    end function

    
    real function Fourier_expansion_cos_sin_DV_zeta( F_mn_c, F_mn_s, theta, zeta ) result(F)
       real, intent(in) :: F_mn_c(:), F_mn_s(:), theta, zeta
       integer :: i, m, n, N_mn, mn_index(size(F_mn_c),2)
       
       N_mn = size(F_mn_c) ! Choose if nyquist modes or not
       if( N_mn == N_modes )     mn_index = mn
       if( N_mn == N_modes_nyq ) mn_index = mn_nyq

       F = 0 
       do i = 1, N_mn ;  m = mn_index(i,1) ; n = mn_index(i,2)
          F = F &
            + Np*n* F_mn_c(i) * sin( m*theta - Np*n*zeta ) &
            - Np*n* F_mn_s(i) * cos( m*theta - Np*n*zeta ) 
       end do 
    end function

    
    real function Fourier_expansion_cos_sin_DV2_theta( F_mn_c, F_mn_s, theta, zeta ) result(F)
       real, intent(in) :: F_mn_c(:), F_mn_s(:), theta, zeta
       integer :: i, m, n, N_mn, mn_index(size(F_mn_c),2)
       
       N_mn = size(F_mn_c) ! Choose if nyquist modes or not
       if( N_mn == N_modes )     mn_index = mn
       if( N_mn == N_modes_nyq ) mn_index = mn_nyq

       F = 0 
       do i = 1, N_mn ;  m = mn_index(i,1) ; n = mn_index(i,2)
          F = F &
            - m*m* F_mn_c(i) * cos( m*theta - Np*n*zeta ) &
            - m*m* F_mn_s(i) * sin( m*theta - Np*n*zeta ) 
       end do 
    end function

    
    
    real function Fourier_expansion_cos_sin_DV2_theta_zeta( F_mn_c, F_mn_s, theta, zeta ) result(F)
       real, intent(in) :: F_mn_c(:), F_mn_s(:), theta, zeta
       integer :: i, m, n, N_mn, mn_index(size(F_mn_c),2)
       
       N_mn = size(F_mn_c) ! Choose if nyquist modes or not
       if( N_mn == N_modes )     mn_index = mn
       if( N_mn == N_modes_nyq ) mn_index = mn_nyq

       F = 0 
       do i = 1, N_mn ;  m = mn_index(i,1) ; n = mn_index(i,2)
          F = F &
            + m*Np*n* F_mn_c(i) * cos( m*theta - Np*n*zeta ) &
            + m*Np*n* F_mn_s(i) * sin( m*theta - Np*n*zeta ) 
       end do 
    end function



    
    real function Fourier_expansion_cos_sin_DV2_zeta( F_mn_c, F_mn_s, theta, zeta ) result(F)
       real, intent(in) :: F_mn_c(:), F_mn_s(:), theta, zeta
       integer :: i, m, n, N_mn, mn_index(size(F_mn_c),2)
       
       N_mn = size(F_mn_c) ! Choose if nyquist modes or not
       if( N_mn == N_modes )     mn_index = mn
       if( N_mn == N_modes_nyq ) mn_index = mn_nyq

       F = 0 
       do i = 1, N_mn ;  m = mn_index(i,1) ; n = mn_index(i,2)
          F = F &
            - Np*n*Np*n* F_mn_c(i) * cos( m*theta - Np*n*zeta ) &
            - Np*n*Np*n* F_mn_s(i) * sin( m*theta - Np*n*zeta ) 
       end do 
    end function


  subroutine read_netcdf_array(tag, variable)
    character(len=*), intent(in) :: tag
    real, intent(out) :: variable(:, :)
    
    integer :: rhid, ierr_netcdf, ncid

    
    ierr_netcdf = nf90_open('VMEC.nc', nf90_nowrite, ncid )
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_inq_varid( ncid, trim(tag), rhid )  
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_get_var( ncid, rhid, variable )      
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_close( ncid )
 

  end subroutine

  
  subroutine read_netcdf_vector(tag, variable)
    character(len=*), intent(in) :: tag
    real, intent(out) :: variable(:)
    
    integer :: rhid, ierr_netcdf, ncid

    ierr_netcdf = nf90_open('VMEC.nc', nf90_nowrite, ncid )
    if(ierr_netcdf/=0) write(*,*) " Error opening ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_inq_varid( ncid, tag, rhid )  
    if(ierr_netcdf/=0) write(*,*) " Error inquiring ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_get_var( ncid, rhid, variable )      
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_close( ncid )

  end subroutine

  
  subroutine read_netcdf_integer_vector(tag, variable)
    character(len=*), intent(in) :: tag
    integer, intent(out) :: variable(:)
    
    integer :: rhid, ierr_netcdf, ncid

    ierr_netcdf = nf90_open('VMEC.nc', nf90_nowrite, ncid )
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_inq_varid( ncid, tag, rhid )  
    if(ierr_netcdf/=0) write(*,*) " Error inquiring ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_get_var( ncid, rhid, variable )      
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_close( ncid )

  end subroutine

  
  subroutine read_netcdf_integer(tag, variable)
    character(len=*), intent(in) :: tag
    integer, intent(out) :: variable
    
    integer :: rhid, ierr_netcdf, ncid 

    ierr_netcdf = nf90_open('VMEC.nc', nf90_nowrite, ncid )
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_inq_varid( ncid, tag, rhid )  
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_get_var( ncid, rhid, variable )      
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_close( ncid )

  end subroutine
  
  subroutine read_netcdf_scalar(tag, variable)
    character(len=*), intent(in) :: tag
    real, intent(out) :: variable 
    
    integer :: rhid, ierr_netcdf, ncid 

    ierr_netcdf = nf90_open('VMEC.nc', nf90_nowrite, ncid )
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_inq_varid( ncid, tag, rhid )  
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_get_var( ncid, rhid, variable )      
    if(ierr_netcdf/=0) write(*,*) " Error reading ", tag, " ierr_netcdf ", ierr_netcdf
    ierr_netcdf = nf90_close( ncid )

  end subroutine
  
   
  
  
end module
