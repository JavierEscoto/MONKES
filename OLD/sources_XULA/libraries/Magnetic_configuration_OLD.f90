module Magnetic_configuration
   
   use Global, only: iout!, psip, chip, Bzeta, Btheta, nzperiod , &
                         !  Nnm, mp, np, bnmc0, bnms0 
   
   use Utilities
   use Finite_differences
   
   implicit none
   
   private
   
   public :: Magnetic_Field, Magnetic_Field_Clebsch
   public :: Magnetic_Field_Derivative_alpha, Magnetic_Field_Derivative_theta
   public :: Magnetic_Field_Derivative_zeta, Magnetic_Field_Derivative_zeta_l

   public :: Magnetic_Field_Boozer_Scaled
   public :: Magnetic_Field_Derivative_theta_Boozer_Scaled
   public :: Magnetic_Field_Derivative_zeta_Boozer_Scaled
   
   public :: Set_Bounce_points, Set_Phase_space
   public :: Initialize_Magnetic_field, Initialize_Magnetic_parameters
   public :: Set_lambda_boundaries
   
   public :: chi_p, psi_p, B_theta, B_zeta, N_p, iota
   public :: B_mn, N_modes, mn
   
   public :: Initialize_Magnetic_field_OLD
   interface Magnetic_Field
      module procedure Magnetic_Field_Boozer_mn, Magnetic_Field_Clebsch_mn
      module procedure Magnetic_Field_Boozer_Scaled_mn
   end interface

   interface Magnetic_Field_Boozer_Scaled
      module procedure Magnetic_Field_Boozer_Scaled_point
   end interface

   interface Magnetic_Field_Derivative_theta_Boozer_Scaled
      module procedure Magnetic_Field_Derivative_theta_Boozer_Scaled_point
   end interface

   interface Magnetic_Field_Derivative_zeta_Boozer_Scaled
      module procedure Magnetic_Field_Derivative_zeta_Boozer_Scaled_point
   end interface 



   interface Magnetic_Field_Clebsch
       module procedure Magnetic_Field_Clebsch_point, Magnetic_Field_Clebsch_Grid       
   end interface
   
   interface Magnetic_Field_Derivative_alpha
       module procedure Magnetic_Field_Derivative_alpha_Clebsch_mn
       module procedure Magnetic_Field_Derivative_alpha_Clebsch
       module procedure Magnetic_Field_Derivative_alpha_Clebsch_Grid
   end interface
   
   interface Magnetic_Field_Derivative_theta
       module procedure Magnetic_Field_Derivative_theta_point
       module procedure Magnetic_Field_Derivative_theta_Boozer_mn
       module procedure Magnetic_Field_Derivative_theta_Boozer_Scaled_mn
   end interface
   
   interface Magnetic_Field_Derivative_zeta
       module procedure Magnetic_Field_Derivative_zeta_Boozer_mn
       module procedure Magnetic_Field_Derivative_zeta_Boozer
       module procedure Magnetic_Field_Derivative_zeta_Boozer_Grid
       module procedure Magnetic_Field_Derivative_zeta_Clebsch
       module procedure Magnetic_Field_Derivative_zeta_Clebsch_Grid
   end interface
   
   interface Magnetic_Field_Derivative_zeta_l
       module procedure Magnetic_Field_Derivative_zeta_Clebsch_mn
       module procedure Magnetic_Field_Derivative_zeta_l_point
       module procedure Magnetic_Field_Derivative_zeta_l_Grid_1D
       module procedure Magnetic_Field_Derivative_zeta_l_Grid_2D
   end interface 
   
   
    
   integer, parameter ::  Nz_dense = 1001!120001!85001!901!35001!, Nz_l_dense = 901
   real, save, allocatable :: B_mn(:) 
   integer, save, allocatable :: mn(:,:), l_ij(:,:), ij_l(:,:)  
   
   integer, save :: N_p, N_modes
   real, save :: chi_p, psi_p, B_theta, B_zeta, iota
   
   real, parameter :: lambda_B_max = 1 !- 1d-4 ! tolerance to avoid singularity at bounce points
      
   
   
   contains
  
  
  subroutine Set_Phase_space( zeta_l, lambda, B_l, Bounce_points_l, Phase_space_l, jl_l, jl_u ) 
        real, intent(in) :: zeta_l(0:), lambda(0:), B_l(0:), Bounce_points_l(0:,0:,:)
        logical, intent(out) :: Phase_space_l(0:,0:)
        integer, intent(out) :: jl_l(0:size(zeta_l)-1,0:size(lambda)-1), &
                                jl_u(0:size(zeta_l)-1,0:size(lambda)-1)
                
        real, parameter :: k_tol = 0*0.1
        real :: delta_zeta, delta_lambda, tol
        integer :: j, l, k, Nz_l, N_lambda
        
        Nz_l = size(zeta_l)-1
        N_lambda = size(lambda)-1
        delta_zeta = ( zeta_l(Nz_l) - zeta_l(0) ) / Nz_l
        delta_lambda = ( lambda(N_lambda) - lambda(0) ) / N_lambda
        tol = max( delta_zeta*k_tol, delta_lambda*minval(B_l) )
        
        do k = 0, N_lambda
           ! *** First guess for phase_space_l
           Phase_space_l(:,k) = Bounce_points_l(:,k,1) <= zeta_l &
                         .and. zeta_l <= Bounce_points_l(:,k,2)
           do l = 0, Nz_l
              call Interval_bounds( zeta_l, Bounce_points_l(l,k,:), &
                                       jl_l(l,k), jl_u(l,k) )            
           end do
        
           do j = 1, 10   
              do l = 0, Nz_l        
                 if( Phase_space_l(l,k) ) then
                   if ( l==jl_l(l,k) .and.  1 - lambda(k) * B_l(l) < tol ) then
                       Phase_space_l(l,k) = 1 - lambda(k) * B_l(l) >= tol
                   end if
                   
                   if ( l==jl_u(l,k) .and.  1 - lambda(k) * B_l(l) < tol ) then
                       Phase_space_l(l,k) = 1 - lambda(k) * B_l(l) >= tol  
                   end if
                 
                 end if
              
              end do
              
              ! *** Third guess for phase_space_l
              do l = 0, Nz_l
                 call Interval_bounds( zeta_l, Bounce_points_l(l,k,:), Phase_space_l(:,k), &
                                          jl_l(l,k), jl_u(l,k) )            
              end do     
              
           end do
        end do
  
  end subroutine
  
  subroutine Set_lambda_boundaries( alpha, zeta_l, lambda_c, lambda_max, B_l ) 
     real, intent(in) :: alpha(0:), zeta_l(0:)
     real, intent(out) :: lambda_c, lambda_max
     real, optional, intent(in) :: B_l(0:) ! Magnetic field in dense grid
     
     integer :: l, Na, Nz_l, Nz_l_dense 
     real :: delta_zeta
     real, allocatable :: zeta_l_dense(:), B_dense(:)     
     
     Na = size(alpha) - 1
     Nz_l = size(zeta_l) - 1
     Nz_l_dense = Na*Nz_dense - 1
     
     ! *** Dense grid to determine bounce points
     allocate( zeta_l_dense(0:Nz_l_dense), B_dense(0:Nz_l_dense) ) 
     delta_zeta = ( zeta_l(Nz_l) - zeta_l(0) ) / Nz_l_dense
     zeta_l_dense = [( zeta_l(0) + l*delta_zeta, l=0, Nz_l_dense )]
     
     ! *** Magnetic field in dense grid
     if( present(B_l) ) then
       B_dense = B_l
     else
       do l = 0, Nz_l_dense
          B_dense(l) = Magnetic_Field_Clebsch( alpha(0), zeta_l_dense(l) ) 
       end do
     end if
     
     ! *** 
     lambda_c = lambda_B_max / maxval(B_dense)
     lambda_max = lambda_B_max / minval(B_dense)
     
     deallocate( zeta_l_dense, B_dense )
  end subroutine
  
  
  
  
  subroutine Set_Bounce_points( alpha, zeta_l, lambda, Bounce_points_l, B_l )
      real, intent(in) :: alpha(0:), zeta_l(0:), lambda(0:)
      real, intent(out) :: Bounce_points_l(0:size(zeta_l)-1,0:size(lambda)-1,2) 
      real, optional, intent(in) :: B_l(0:) ! Magnetic field in dense grid
      
      integer :: j, l, l_max, k, Na, Nz_l, N_lambda, Nz_l_dense
      integer, allocatable :: jl_l(:), jl_u(:) 
      logical, allocatable :: Phase_space_l(:)
      logical :: no_bounce_points
      real, allocatable :: zeta_l_dense(:), B_dense(:)
      real :: delta_zeta
      
      Na = size(alpha) - 1 
      Nz_l = size(zeta_l) - 1 
      N_lambda = size(lambda) - 1 
      Nz_l_dense = Na*Nz_dense - 1
      
      ! *** Dense grid to determine bounce points   
      allocate( zeta_l_dense(0:Nz_l_dense) )
      delta_zeta = (zeta_l(Nz_l)-zeta_l(0))/Nz_l_dense
      zeta_l_dense = [( zeta_l(0) + l*delta_zeta, l=0,Nz_l_dense )]
            
      ! *** Magnetic field in dense grid      
      allocate( B_dense(0:Nz_l_dense) )
      if( present(B_l) ) then
        B_dense = B_l
      else
        
        do l = 0, Nz_l_dense
           B_dense(l) = Magnetic_Field_Clebsch( alpha(0), zeta_l_dense(l) ) 
        end do
      end if
      
      ! *** Position of maximum magnetic field along zeta_l_dense
      l_max = maxloc( B_dense, dim=1 ) - 1
      allocate( jl_l(0:Nz_l_dense), jl_u(0:Nz_l_dense), Phase_space_l(0:Nz_l_dense) )
      do k = 0, N_lambda
             
        ! *** Phase space
        Phase_space_l = 1 - lambda(k)  * B_dense >= ( 1 - lambda_B_max ) !0  
        call Grid_bounds_1D( zeta_l_dense, Phase_space_l, jl_l, jl_u )
        jl_l = jl_l - 1 
        jl_u = jl_u - 1      
        
        ! *** Bounce points determination
        Bounce_points_l(:,k,2) = 0 ! Start with absurd values for 
        Bounce_points_l(:,k,1) = 1 ! bounce points
        
          do l = 0, Nz_l_dense
             if ( Phase_space_l(l) .and. l == jl_l(l) ) then   !jl_l(l) >= 0 ) then
               do j = 0, Nz_l
               
                  if ( zeta_l(j) >= zeta_l_dense( jl_l(l) )      &
                       .and. zeta_l(j) <= zeta_l_dense( jl_u(l)) ) then 
                       
                    ! *** Bounce points using Bisection iteration 
                    Bounce_points_l(j,k,1) = zeta_l_dense( jl_l(l) ) 
                    Bounce_points_l(j,k,2) = zeta_l_dense( jl_u(l) )
                    
                    if( jl_l(l) > 0 ) & 
                    Bounce_points_l(j,k,1) =&
                    Bisection_Bounce_points( zeta_l_dense(jl_l(l)-1), Bounce_points_l(j,k,1), lambda(k) )
				                          
                    if( jl_u(l) < Nz_l_dense ) &              
                    Bounce_points_l(j,k,2) = &
                    Bisection_right_Bounce_points( zeta_l_dense(jl_u(l)+1), Bounce_points_l(j,k,2), lambda(k) )
                        
                    ! *** For orbits in which we have no bounce points
                    ! take the position of maximum l as left or right
                    ! bounce point and consider it as the boundary.
                    ! no_bounce_points =  &
                    !     ( jl_l(l) == 0 ) .and. ( jl_u(l) == Nz_l_dense ) 
                    ! if( no_bounce_points ) then
                    !    if ( zeta_l(j) < zeta_l_dense(l_max) ) then 
                    !    !Bounce_points_l(j,k,1) = zeta_l_dense(0)
                    !     Bounce_points_l(j,k,2) = zeta_l_dense(l_max)
                    !  
                    !    
                    !    
                    !    elseif ( zeta_l(j) > zeta_l_dense(l_max) ) then
                    !      Bounce_points_l(j,k,1) = zeta_l_dense(l_max)
                    !      !Bounce_points_l(j,k,2) = zeta_l_dense(Nz_l_dense)
                    !    else
                    !      Bounce_points_l(j,k,1) = 1
                    !      Bounce_points_l(j,k,2) = 0
                    !    end if
                    ! 
                    ! end if
                    
                  end if
               
               end do
             end if
          end do
          
          ! *** Second iteration to avoid casuistics in which zeta_l is
          ! between bounce points but not between zeta_l_dense(jl_l) and
          ! zeta_l_dense(jl_u) 
          
          do l = 0, Nz_l
             do j = 0, Nz_l
                if( Bounce_points_l(l,k,1) <= zeta_l(j) .and. & 
                    zeta_l(j) <= Bounce_points_l(l,k,2) ) &
                    Bounce_points_l(j,k,:) = Bounce_points_l(l,k,:)
             end do
          end do
        
     end do
     
     
      
     deallocate( zeta_l_dense, B_dense, Phase_space_l, jl_l, jl_u )
     
     
     
     contains
     
    real function Newton_Bounce_points( zeta_0, lambda ) result( zeta_b )
       real, intent(in) :: zeta_0, lambda
       
       integer, parameter :: it_max = 35000
       integer :: i
       real, parameter :: tol = 1d-15
       real :: B, dBdz, eps, zeta_1, F, zeta_b_old
       
       zeta_b = zeta_0
       eps = 1 ; i = 0
       do while ( eps > tol .and. i <= it_max ) 
                 
          B = Magnetic_Field_Clebsch( alpha(0), zeta_b )
          dBdz = Magnetic_Field_Derivative_zeta_l( alpha(0), zeta_b )
          F = 1 - lambda * B
          zeta_b = zeta_b + F/( lambda*dBdz )
          eps = abs(F)
          i = i +1 
            
       end do
    end function
      
    real function Bisection_Bounce_points( zeta_0, zeta_1, lambda ) result( zeta_b )
       real, intent(in) :: zeta_0, zeta_1, lambda
       
       integer, parameter :: it_max = 35000
       integer :: i
       real, parameter :: tol = 1d-12!1d-6 !1d-15
       real :: eps, F, zeta_b_minus, zeta_b_plus, B
       
       zeta_b_minus = zeta_0 ; zeta_b_plus = zeta_1
       eps = 1 ; i = 0
       do while ( eps > tol .and. i <= it_max ) 
       
          zeta_b = (zeta_b_plus + zeta_b_minus)/2
          B = Magnetic_Field_Clebsch( alpha(0), zeta_b )
          F  = 1 - lambda * B - ( 1 - lambda_B_max )  
          
          if( F >= 0 ) zeta_b_plus = zeta_b
          if( F < 0 ) zeta_b_minus = zeta_b
          eps = abs(F)
          i = i + 1 
          
       end do
       
       zeta_b = zeta_b_plus
       
     end function 
      
    real function Bisection_right_Bounce_points( zeta_0, zeta_1, lambda ) result( zeta_b )
       real, intent(in) :: zeta_0, zeta_1, lambda
       
       integer, parameter :: it_max = 35000
       integer :: i
       real, parameter :: tol = 1d-12! 1d-8!1d-6!!1d-15
       real :: eps, F, zeta_b_minus, zeta_b_plus, B
       
       zeta_b_minus = zeta_0 ; zeta_b_plus = zeta_1
       eps = 1 ; i = 0
       do while ( eps > tol .and. i <= it_max ) 
       
          zeta_b = (zeta_b_plus + zeta_b_minus)/2
          B = Magnetic_Field_Clebsch( alpha(0), zeta_b )
          F  = 1 - lambda * B - ( 1 - lambda_B_max )
          
          if( F >= 0 ) zeta_b_plus = zeta_b
          if( F < 0 ) zeta_b_minus = zeta_b
          eps = abs(F)
          i = i + 1 
          
       end do
       
       zeta_b = zeta_b_plus
       
     end function 
     
     function Function_Bounce_Points(zeta) result(F) 
      real, intent(in) :: zeta(:) 
      real :: F( size(zeta) ) 
      
      real :: B
      
      B = Magnetic_Field_Clebsch( alpha(0), zeta(1) ) 
      
      F(1) = 1 - lambda(k) * B 
      
     end function
  end subroutine
  
  
  subroutine Set_Bounce_points_OLD( alpha, zeta_l, lambda, Bounce_points_l, B_l )
      real, intent(in) :: alpha(0:), zeta_l(0:), lambda(0:)
      real, intent(out) :: Bounce_points_l(0:size(zeta_l)-1,0:size(lambda)-1,2) 
      real, optional, intent(in) :: B_l(0:) ! Magnetic field in dense grid
      
      integer :: j, l, l_max, k, Na, Nz_l, N_lambda, Nz_l_dense
      integer, allocatable :: jl_l(:), jl_u(:) 
      logical, allocatable :: Phase_space_l(:)
      logical :: no_bounce_points
      real, allocatable :: zeta_l_dense(:), B_dense(:)
      real :: delta_zeta
      
      Na = size(alpha) - 1 
      Nz_l = size(zeta_l) - 1 
      N_lambda = size(lambda) - 1 
      Nz_l_dense = Na*Nz_dense - 1
      
      ! *** Dense grid to determine bounce points   
      allocate( zeta_l_dense(0:Nz_l_dense) )
      delta_zeta = (zeta_l(Nz_l)-zeta_l(0))/Nz_l_dense
      zeta_l_dense = [( zeta_l(0) + l*delta_zeta, l=0,Nz_l_dense )]
      
      
      ! *** Magnetic field in dense grid      
      allocate( B_dense(0:Nz_l_dense) )
      if( present(B_l) ) then
        B_dense = B_l
      else
        
        do l = 0, Nz_l_dense
           B_dense(l) = Magnetic_Field_Clebsch( alpha(0), zeta_l_dense(l) ) 
        end do
      end if
      
      ! *** Position of maximum magnetic field along zeta_l_dense
      l_max = maxloc( B_dense, dim=1 ) - 1
      allocate( jl_l(0:Nz_l_dense), jl_u(0:Nz_l_dense), Phase_space_l(0:Nz_l_dense) )
      do k = 0, N_lambda
             
        ! *** Phase space
        Phase_space_l = 1 - lambda(k)  * B_dense >= ( 1 - lambda_B_max ) !0  
        call Grid_bounds_1D( zeta_l_dense, Phase_space_l, jl_l, jl_u )
        jl_l = jl_l - 1 
        jl_u = jl_u - 1      
        
        ! *** Bounce points determination
        Bounce_points_l(:,k,2) = 0 ! Start with absurd values for 
        Bounce_points_l(:,k,1) = 1 ! bounce points
        
          do l = 0, Nz_l_dense
             if ( Phase_space_l(l) .and. l == jl_l(l) ) then   !jl_l(l) >= 0 ) then
               do j = 0, Nz_l
               
                  if ( zeta_l(j) >= zeta_l_dense( jl_l(l) )      &
                       .and. zeta_l(j) <= zeta_l_dense( jl_u(l)) ) then 
                       
                    ! *** Bounce points using Bisection iteration 
                    Bounce_points_l(j,k,1) = zeta_l_dense( jl_l(l) ) 
                    Bounce_points_l(j,k,2) = zeta_l_dense( jl_u(l) )
                    
                    if( jl_l(l) > 0 ) & 
                    Bounce_points_l(j,k,1) =&
                    Bisection_Bounce_points( zeta_l_dense(jl_l(l)-1), Bounce_points_l(j,k,1), lambda(k) )
				                          
                    if( jl_u(l) < Nz_l_dense ) &              
                    Bounce_points_l(j,k,2) = &
                    Bisection_right_Bounce_points( zeta_l_dense(jl_u(l)+1), Bounce_points_l(j,k,2), lambda(k) )
                        
                    ! *** For orbits in which we have no bounce points
                    ! take the position of maximum l as left or right
                    ! bounce point and consider it as the boundary.
                    ! no_bounce_points =  &
                    !     ( jl_l(l) == 0 ) .and. ( jl_u(l) == Nz_l_dense ) 
                    ! if( no_bounce_points ) then
                    !    if ( zeta_l(j) < zeta_l_dense(l_max) ) then 
                    !    !Bounce_points_l(j,k,1) = zeta_l_dense(0)
                    !     Bounce_points_l(j,k,2) = zeta_l_dense(l_max)
                    !  
                    !    
                    !    
                    !    elseif ( zeta_l(j) > zeta_l_dense(l_max) ) then
                    !      Bounce_points_l(j,k,1) = zeta_l_dense(l_max)
                    !      !Bounce_points_l(j,k,2) = zeta_l_dense(Nz_l_dense)
                    !    else
                    !      Bounce_points_l(j,k,1) = 1
                    !      Bounce_points_l(j,k,2) = 0
                    !    end if
                    ! 
                    ! end if
                    
                  end if
               
               end do
             end if
          end do
          
          ! *** Second iteration to avoid casuistics in which zeta_l is
          ! between bounce points but not between zeta_l_dense(jl_l) and
          ! zeta_l_dense(jl_u) 
          
          do l = 0, Nz_l
             do j = 0, Nz_l
                if( Bounce_points_l(l,k,1) <= zeta_l(j) .and. & 
                    zeta_l(j) <= Bounce_points_l(l,k,2) ) &
                    Bounce_points_l(j,k,:) = Bounce_points_l(l,k,:)
             end do
          end do
        
     end do
     
     
      
     deallocate( zeta_l_dense, B_dense, Phase_space_l, jl_l, jl_u )
     
     
     
     contains
     
    real function Newton_Bounce_points( zeta_0, lambda ) result( zeta_b )
       real, intent(in) :: zeta_0, lambda
       
       integer, parameter :: it_max = 35000
       integer :: i
       real, parameter :: tol = 1d-15
       real :: B, dBdz, eps, zeta_1, F, zeta_b_old
       
       zeta_b = zeta_0
       eps = 1 ; i = 0
       do while ( eps > tol .and. i <= it_max ) 
                 
          B = Magnetic_Field_Clebsch( alpha(0), zeta_b )
          dBdz = Magnetic_Field_Derivative_zeta_l( alpha(0), zeta_b )
          F = 1 - lambda * B
          zeta_b = zeta_b + F/( lambda*dBdz )
          eps = abs(F)
          i = i +1 
            
       end do
    end function
      
    real function Bisection_Bounce_points( zeta_0, zeta_1, lambda ) result( zeta_b )
       real, intent(in) :: zeta_0, zeta_1, lambda
       
       integer, parameter :: it_max = 35000
       integer :: i
       real, parameter :: tol = 1d-12!1d-6 !1d-15
       real :: eps, F, zeta_b_minus, zeta_b_plus, B
       
       zeta_b_minus = zeta_0 ; zeta_b_plus = zeta_1
       eps = 1 ; i = 0
       do while ( eps > tol .and. i <= it_max ) 
       
          zeta_b = (zeta_b_plus + zeta_b_minus)/2
          B = Magnetic_Field_Clebsch( alpha(0), zeta_b )
          F  = 1 - lambda * B - ( 1 - lambda_B_max )  
          
          if( F >= 0 ) zeta_b_plus = zeta_b
          if( F < 0 ) zeta_b_minus = zeta_b
          eps = abs(F)
          i = i + 1 
          
       end do
       
       zeta_b = zeta_b_plus
       
     end function 
      
    real function Bisection_right_Bounce_points( zeta_0, zeta_1, lambda ) result( zeta_b )
       real, intent(in) :: zeta_0, zeta_1, lambda
       
       integer, parameter :: it_max = 35000
       integer :: i
       real, parameter :: tol = 1d-12! 1d-8!1d-6!!1d-15
       real :: eps, F, zeta_b_minus, zeta_b_plus, B
       
       zeta_b_minus = zeta_0 ; zeta_b_plus = zeta_1
       eps = 1 ; i = 0
       do while ( eps > tol .and. i <= it_max ) 
       
          zeta_b = (zeta_b_plus + zeta_b_minus)/2
          B = Magnetic_Field_Clebsch( alpha(0), zeta_b )
          F  = 1 - lambda * B - ( 1 - lambda_B_max )
          
          if( F >= 0 ) zeta_b_plus = zeta_b
          if( F < 0 ) zeta_b_minus = zeta_b
          eps = abs(F)
          i = i + 1 
          
       end do
       
       zeta_b = zeta_b_plus
       
     end function 
     
     function Function_Bounce_Points(zeta) result(F) 
      real, intent(in) :: zeta(:) 
      real :: F( size(zeta) ) 
      
      real :: B
      
      B = Magnetic_Field_Clebsch( alpha(0), zeta(1) ) 
      
      F(1) = 1 - lambda(k) * B 
      
     end function
  end subroutine
  
  subroutine Set_Bounce_points_OLDEST( alpha, zeta_l, lambda, Bounce_points_l, B_l )
      real, intent(in) :: alpha(0:), zeta_l(0:), lambda(0:)
      real, intent(out) :: Bounce_points_l(0:size(zeta_l)-1,0:size(lambda)-1,2) 
      real, optional, intent(in) :: B_l(0:) ! Magnetic field in dense grid
      
      integer :: j, l, k, Na, Nz_l, N_lambda, Nz_l_dense
      integer, allocatable :: jl_l(:), jl_u(:) 
      logical, allocatable :: Phase_space_l(:)
      real, allocatable :: zeta_l_dense(:), B_dense(:)
      real :: delta_zeta
      
      Na = size(alpha) - 1 
      Nz_l = size(zeta_l) - 1 
      N_lambda = size(lambda) - 1 
      Nz_l_dense = Na*Nz_dense - 1
      
      ! *** Dense grid to determine bounce points   
      allocate( zeta_l_dense(0:Nz_l_dense) )
      delta_zeta = (zeta_l(Nz_l)-zeta_l(0))/Nz_l_dense
      zeta_l_dense = [( zeta_l(0) + l*delta_zeta, l=0,Nz_l_dense )]
      
      
      ! *** Magnetic field in dense grid      
      allocate( B_dense(0:Nz_l_dense) )
      if( present(B_l) ) then
        B_dense = B_l
      else
        
        do l = 0, Nz_l_dense
           B_dense(l) = Magnetic_Field_Clebsch( alpha(0), zeta_l_dense(l) ) 
        end do
      end if
      
      allocate( jl_l(0:Nz_l_dense), jl_u(0:Nz_l_dense), Phase_space_l(0:Nz_l_dense) )
      do k = 0, N_lambda
             
        ! *** Phase space
        Phase_space_l = 1 - lambda(k)  * B_dense >= 0  
        call Grid_bounds_1D( zeta_l_dense, Phase_space_l, jl_l, jl_u )
        jl_l = jl_l - 1 
        jl_u = jl_u - 1      
        
        ! *** First guess for phase_space_l
        Phase_space_l = 1 - lambda(k)  * B_dense >= 0
        ! *** Bounce points determination
        Bounce_points_l(:,k,2) = 0 ! Start with absurd values for 
        Bounce_points_l(:,k,1) = 1 ! bounce points
        
        do l = 0, Nz_l_dense
           if ( Phase_space_l(l) .and. jl_l(l) >= 0 ) then
             do j = 0, Nz_l
             
                if ( zeta_l(j) >= zeta_l_dense( jl_l(l) )      &
                     .and. zeta_l(j) <= zeta_l_dense( jl_u(l)) ) then 
                    Bounce_points_l(j,k,1) = zeta_l_dense( jl_l(l) )
                    Bounce_points_l(j,k,2) = zeta_l_dense( jl_u(l) )
                    
                end if
             
             end do
           end if
        end do
        
     end do
      
     deallocate( zeta_l_dense, B_dense, Phase_space_l, jl_l, jl_u )
     
     
     
     contains
     
     
     function Function_Bounce_Points(zeta) result(F) 
      real, intent(in) :: zeta(:) 
      real :: F( size(zeta) ) 
      
      real :: B
      
      B = Magnetic_Field_Clebsch( alpha(0), zeta(1) ) 
      F(1) = 1 - lambda(k) * B 
      
     end function
  end subroutine
  
  subroutine Set_Bounce_points_NEW( alpha, zeta_l, lambda, Bounce_points_l, B_l )
      real, intent(in) :: alpha(0:), zeta_l(0:), lambda(0:)
      real, intent(out) :: Bounce_points_l(0:size(zeta_l)-1,0:size(lambda)-1,2) 
      real, optional, intent(in) :: B_l(0:) ! Magnetic field in dense grid
      
      integer :: j, l, k, Na, Nz_l, N_lambda, Nz_l_dense, l0 , l1
      integer, allocatable :: jl_l(:), jl_u(:) 
      logical, allocatable :: Phase_space_l(:)
      real, allocatable :: zeta_l_dense(:), B_dense(:)
      real :: delta_zeta, F_minus, F_plus
      
      Na = size(alpha) - 1 
      Nz_l = size(zeta_l) - 1 
      N_lambda = size(lambda) - 1 
      Nz_l_dense = Na*Nz_dense - 1
      
      ! *** Dense grid to determine bounce points   
      allocate( zeta_l_dense(0:Nz_l_dense) )
      delta_zeta = (zeta_l(Nz_l)-zeta_l(0))/Nz_l_dense
      zeta_l_dense = [( zeta_l(0) + l*delta_zeta, l=0,Nz_l_dense )]
      
      
      ! *** Magnetic field in dense grid      
      allocate( B_dense(0:Nz_l_dense) )
      if( present(B_l) ) then
        B_dense = B_l
      else
        
        do l = 0, Nz_l_dense
           B_dense(l) = Magnetic_Field_Clebsch( alpha(0), zeta_l_dense(l) ) 
        end do
      end if
      
      allocate( jl_l(0:Nz_l_dense), jl_u(0:Nz_l_dense), Phase_space_l(0:Nz_l_dense) )
      do k = 0, N_lambda
             
        ! *** Phase space
        !Phase_space_l = 1 - lambda(k)  * B_dense >= 0  
        
        ! *** Bounce points determination
        Bounce_points_l(:,k,2) = 0 ! Start with absurd values for 
        Bounce_points_l(:,k,1) = 1 ! bounce points
        
        do l = 0, Nz_l_dense
           if ( Phase_space_l(l) .and. jl_l(l) >= 0 ) then
             do j = 0, Nz_l
             
                if ( zeta_l(j) >= zeta_l_dense( jl_l(l) )      &
                     .and. zeta_l(j) <= zeta_l_dense( jl_u(l)) ) then 
                    Bounce_points_l(j,k,1) = zeta_l_dense( jl_l(l) )
                    Bounce_points_l(j,k,2) = zeta_l_dense( jl_u(l) )
                    
                end if
             
             end do
           end if
        end do
        
     end do
      
     deallocate( zeta_l_dense, B_dense, Phase_space_l, jl_l, jl_u )
  end subroutine
  
   
  subroutine Initialize_Magnetic_parameters(psip, chip, Bzeta, Btheta, nzperiod)
      real, intent(in) :: psip, chip, Bzeta, Btheta              
      integer, intent(in) :: nzperiod             
            
      chi_p = chip
      psi_p = psip
      B_theta = Btheta
      B_zeta = Bzeta
      N_p = nzperiod
      iota = chi_p / psi_p
      
      write(iout,*) " chi_p, psi_p, B_theta, B_zeta, N_p , iota "
      write(iout,*)   chi_p, psi_p, B_theta, B_zeta, N_p , iota  
  end subroutine
  
  subroutine Initialize_Magnetic_field( Nnm, mp, np, bnmc0, bnms0 )
      integer, intent(in) :: Nnm
      real, intent(in) :: mp(Nnm), np(Nnm), bnmc0(Nnm), bnms0(Nnm)
      integer :: k
      
     
     if ( allocated( B_mn ) ) deallocate( B_mn, mn )  
     N_modes = Nnm
     allocate( B_mn(N_modes), mn(N_modes,2) )
     
     ! *** Translate modes to nearest integer
     mn(:,1) = nint(mp) ; mn(:,2) = nint(np) 
     B_mn = bnmc0
     
     
  end subroutine
      
  subroutine Initialize_Magnetic_field_OLD(path)    
     character(len=*), intent(in) :: path
     integer :: k
     
     if ( allocated( B_mn ) ) deallocate( B_mn, mn )     
     call Read_magnetic_field( path, mn, B_mn )
     N_modes = size(B_mn) 
     
     call Read_magnetic_parameters( path, N_p, chi_p, psi_p, iota, B_theta, B_zeta )
      
      
     open(21, file="results/Magnetic_field_modes_FJEL.plt")
     do k = 1, N_modes
        write(21,*) mn(k,:), B_mn(k) 
     end do 
     close(21)
     
     !call Initialize_Magnetic_field_NEW
     
     !stop
  end subroutine 
   
  subroutine Read_magnetic_field( path, mn, B_mn ) 
     character(len=*), intent(in) :: path
     integer, allocatable, intent(out) :: mn(:,:) 
     real, allocatable, intent(out) :: B_mn(:)
     
     integer, parameter :: N_header = 1 
     integer :: k, N_lines 
     
     ! *** Number of lines in the file with magnetic field modes
     N_lines = Count_lines(trim(path)//"Magnetic_field_modes.dat")
     
     allocate( mn(N_lines-N_header,2), B_mn(N_lines-N_header) )
     open(21,file=trim(path)//"Magnetic_field_modes.dat")
     
     ! *** Read header
     do k = 1, N_header
        read(21,*)
     end do
     ! *** Read modes
     do k = 1, N_lines - N_header
        read(21,*) mn(k,2), mn(k,1), B_mn(k)       
     end do
     close(21)
     
  end subroutine
   
  subroutine Read_magnetic_parameters( path, N_p, chi_p, psi_p, iota, B_theta, B_zeta ) 
     character(len=*), intent(in) :: path
     integer, intent(out) :: N_p 
     real, intent(out) :: chi_p, psi_p, iota, B_theta, B_zeta
     
     integer, parameter :: N_header = 1 
     integer :: k, N_lines 
     
     ! *** Number of lines in the file with magnetic field modes
     N_lines = 2
     open(21,file=trim(path)//"Magnetic_parameters.dat")
     
     ! *** Read header
     do k = 1, N_header
        read(21,*)
     end do
     
     ! *** Read parameters
     do k = 1, N_lines - N_header
        read(21,*) N_p, chi_p, psi_p, B_theta, B_zeta      
     end do
     close(21)
     
     ! **** Rotational transform
     iota = -chi_p / psi_p
     
  end subroutine


  real function Magnetic_Field_Clebsch_mn( alpha, zeta, iota, Np, mn, B_mn ) result(B)
     real, intent(in) :: alpha, zeta, iota
     integer, intent(in) :: Np, mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          + B_mn(k) * cos( m * ( alpha + iota*zeta ) + n *Np * zeta )
     end do
     
  end function

  real function Magnetic_Field_Derivative_alpha_Clebsch_mn( alpha, zeta, iota, Np, mn, B_mn ) &
                                           result(dB_alpha)
     real, intent(in) :: alpha, zeta, iota
     integer, intent(in) :: Np, mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     dB_alpha = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        dB_alpha = dB_alpha &
          - m * B_mn(k) * sin( m * ( alpha + iota*zeta ) + n *Np * zeta )
     end do
     
  end function
  
  real function Magnetic_Field_Derivative_alpha_Clebsch( alpha, zeta ) &
                                           result(dB_alpha)
     real, intent(in) :: alpha, zeta
     
     integer :: m, n, k 
     
     dB_alpha = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        dB_alpha = dB_alpha &
          - m * B_mn(k) * sin( m * ( alpha + iota*zeta ) + n *N_p * zeta )
     end do
     
  end function
  
  
  function Magnetic_Field_Derivative_alpha_Clebsch_Grid( alpha, zeta ) &
                                           result(dB_alpha)
     real, intent(in) :: alpha(0:), zeta(0:) 
     real :: dB_alpha(0:size(alpha)-1,0:size(zeta)-1)
     
     integer :: i, j, Na, Nz 
     
     Na = size(alpha)-1 ; Nz = size(zeta)-1
     do j = 0, Nz
        do i = 0, Na
           dB_alpha(i,j) = Magnetic_Field_Derivative_alpha_Clebsch( alpha(i), zeta(j) )
        end do
     end do
     
  end function
  
  ! *** Derivative of magnetic field along toroidal angle following a field
  ! line, i.e. keeping alpha constant.
  real function Magnetic_Field_Derivative_zeta_Clebsch_mn( alpha, zeta, iota, Np, mn, B_mn ) &
                        result(dB_zeta)
     real, intent(in) :: alpha, zeta, iota
     integer, intent(in) :: Np, mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     dB_zeta = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        dB_zeta = dB_zeta &
          - ( n * Np + iota * m ) &
          * B_mn(k) * sin( m * ( alpha + iota*zeta ) + n *Np * zeta )
     end do
     
  end function
  
  ! *** Derivative of magnetic field along toroidal angle following a field
  ! line, i.e. keeping alpha constant.
  real function Magnetic_Field_Derivative_zeta_l_point( alpha, zeta ) &
                        result(dB_zeta)
     real, intent(in) :: alpha, zeta
     
     integer :: m, n, k  
     
     dB_zeta = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        dB_zeta = dB_zeta &
          - ( n * N_p + iota * m ) &
          * B_mn(k) * sin( m * ( alpha + iota*zeta ) + n *N_p * zeta )
     end do
     
  end function
  
  ! *** Derivative of magnetic field along toroidal angle following a field
  ! line, i.e. keeping alpha constant. It gives back a vector.
  function Magnetic_Field_Derivative_zeta_l_Grid_1D( alpha_0, zeta_l ) &
                        result(dB_zeta_l) 
     real, intent(in) :: alpha_0, zeta_l(0:) 
     real :: dB_zeta_l(0:size(zeta_l)-1)
     
     integer :: l, Nz_l 
     
     Nz_l = size(zeta_l)-1
     do l = 0, Nz_l
        dB_zeta_l(l) = &
          Magnetic_Field_Derivative_zeta_l_point( alpha_0, zeta_l(l) )
     end do
     
  end function
  
  ! *** Derivative of magnetic field along toroidal angle following a field
  ! line, i.e. keeping alpha constant. It gives back a 2D array.
  function Magnetic_Field_Derivative_zeta_l_Grid_2D( alpha, zeta ) &
                        result(dB_zeta_l) 
     real, intent(in) :: alpha(0:), zeta(0:) 
     real :: dB_zeta_l(0:size(alpha)-1,0:size(zeta)-1)
     
     integer :: i, j, Na, Nz 
     
     Na = size(alpha)-1 ; Nz = size(zeta)-1
     do j = 0, Nz
        do i = 0, Na
           dB_zeta_l(i,j) = &
           Magnetic_Field_Derivative_zeta_l( alpha(i), zeta(j) )
        end do
     end do
     
  end function

  real function Magnetic_Field_Boozer_mn( theta, zeta, Np, mn, B_mn ) result(B)
     real, intent(in) :: theta, zeta
     integer, intent(in) :: Np, mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          + B_mn(k) * cos( m * theta + n *Np * zeta )
     end do
     
   end function 


  real function Magnetic_Field_Boozer_Scaled_mn( theta, zeta, mn, B_mn ) result(B)
     real, intent(in) :: theta, zeta
     integer, intent(in) :: mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          + B_mn(k) * cos( m * theta + n * zeta )
     end do
     
  end function

  real function Magnetic_Field_Derivative_theta_Boozer_Scaled_mn( theta, zeta, mn, B_mn ) result(B)
     real, intent(in) :: theta, zeta
     integer, intent(in) :: mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          - m * B_mn(k) * sin( m * theta + n * zeta )
     end do
     
  end function 

  real function Magnetic_Field_Derivative_zeta_Boozer_Scaled_mn( theta, zeta, mn, B_mn ) result(B)
     real, intent(in) :: theta, zeta
     integer, intent(in) :: mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          - n * B_mn(k) * sin( m * theta + n * zeta )
     end do
     
  end function
  real function Magnetic_Field_Derivative_theta_Boozer_mn( theta, zeta, Np, mn, B_mn ) result(dB_theta)
     real, intent(in) :: theta, zeta
     integer, intent(in) :: Np, mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     dB_theta = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        dB_theta = dB_theta &
          - m * B_mn(k) * sin( m * theta + n *Np * zeta )
     end do
     
  end function


  real function Magnetic_Field_Derivative_zeta_Boozer_mn( theta, zeta, Np, mn, B_mn ) result(dB_zeta)
     real, intent(in) :: theta, zeta
     integer, intent(in) :: Np, mn(:,:) 
     real, intent(in) :: B_mn(:)
     
     integer :: m, n, k, N_modes 
  
     N_modes = size(B_mn) 
     
     dB_zeta = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        dB_zeta = dB_zeta &
          - n * Np * B_mn(k) * sin( m * theta + n *Np * zeta )
     end do
     
  end function



  real function Magnetic_Field_Clebsch_point( alpha, zeta ) result(B)
     real, intent(in) :: alpha, zeta 
     
     integer :: m, n, k 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          + B_mn(k) * cos( m * ( alpha + iota*zeta ) + n *N_p * zeta )
     end do
     
  end function


  real function Magnetic_Field_Boozer_point( theta, zeta ) result(B)
     real, intent(in) :: theta, zeta 
     
     integer :: m, n, k 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          + B_mn(k) * cos( m * theta + n *N_p * zeta )
     end do
     
  end function


  real function Magnetic_Field_Boozer_Scaled_point( theta, zeta ) result(B)
     real, intent(in) :: theta, zeta 
     
     integer :: m, n, k 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          + B_mn(k) * cos( m * theta + n * zeta )
     end do
     
  end function



  real function Magnetic_Field_Derivative_theta_Boozer_Scaled_point( theta, zeta ) result(B)
     real, intent(in) :: theta, zeta 
     
     integer :: m, n, k 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          - m * B_mn(k) * sin( m * theta + n * zeta )
     end do
     
  end function



  real function Magnetic_Field_Derivative_zeta_Boozer_Scaled_point( theta, zeta ) result(B)
     real, intent(in) :: theta, zeta 
     
     integer :: m, n, k 
     
     B = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        B = B &
          - n * B_mn(k) * sin( m * theta + n * zeta )
     end do
     
  end function
  
  function Magnetic_Field_Clebsch_Grid( alpha, zeta ) result(B)
     real, intent(in) :: alpha(0:), zeta(0:) 
     real :: B(0:size(alpha)-1,0:size(zeta)-1)
     
     integer :: i, j, Na, Nz 
     
     Na = size(alpha)-1 ; Nz = size(zeta)-1
     do j = 0, Nz
        do i = 0, Na
           B(i,j) = Magnetic_Field_Clebsch( alpha(i), zeta(j) )
        end do
     end do
     
  end function

  
  


  real function Magnetic_Field_Derivative_zeta_Boozer( theta, zeta ) result(dB_zeta)
     real, intent(in) :: theta, zeta
     
     integer :: m, n, k
       
     dB_zeta = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        dB_zeta = dB_zeta &
          - n * N_p * B_mn(k) * sin( m * theta + n *N_p * zeta )
     end do
     
  end function


  function Magnetic_Field_Derivative_zeta_Boozer_Grid( theta, zeta ) result(dB_zeta) 
     real, intent(in) :: theta(0:), zeta(0:) 
     real :: dB_zeta(0:size(theta)-1,0:size(zeta)-1)
     
     integer :: i, j, Nt, Nz 
     
     Nt = size(theta)-1 ; Nz = size(zeta)-1
     do j = 0, Nz
        do i = 0, Nt
           dB_zeta(i,j) = &
           Magnetic_Field_Derivative_zeta_Boozer( theta(i), zeta(j) )
        end do
     end do
     
  end function

  function Magnetic_Field_Derivative_zeta_Clebsch_Grid( alpha, zeta, iota ) result(dB_zeta) 
     real, intent(in) :: alpha(0:), zeta(0:), iota
     real :: dB_zeta(0:size(alpha)-1,0:size(zeta)-1)
     
     real :: theta
     integer :: i, j, Na, Nz 
     
     Na = size(alpha)-1 ; Nz = size(zeta)-1
     do j = 0, Nz
        do i = 0, Na
           theta = alpha(i) + iota * zeta(j) 
           dB_zeta(i,j) = &
           Magnetic_Field_Derivative_zeta_Boozer( theta, zeta(j) )
        end do
     end do
     
  end function

  real function Magnetic_Field_Derivative_zeta_Clebsch( alpha, zeta, iota ) result(dB_zeta) 
     real, intent(in) :: alpha, zeta, iota
     
     real :: theta
     
     theta = alpha + iota * zeta 
     dB_zeta = Magnetic_Field_Derivative_zeta_Boozer( theta, zeta )
     
  end function


  real function Magnetic_Field_Derivative_theta_point( theta, zeta ) result(dB_theta)
     real, intent(in) :: theta, zeta
     
     integer :: m, n, k 
     
     dB_theta = 0
     do k = 1, N_modes
        m = mn(k,1) ; n = mn(k,2) 
        
        dB_theta = dB_theta &
          - m * B_mn(k) * sin( m * theta + n *N_p * zeta )
     end do
     
  end function




end module
