module Adjoint_DKE_BTD_Solution_Legendre 
   
   
   use Magnetic_configuration
   
   use Finite_Differences
   use Barycentric_grids   
   
   use Linear_systems
   use plots   
   
   implicit none

   private

   public :: Solve_BTD_Adjoint_DKE_Legendre 
   
   real, parameter :: pi = acos(-1d0)   
      
   ! Phase space grids
   logical, save :: Grids_and_Coefficients_DKE = .false. 
   real, save, allocatable :: theta(:), zeta(:) 
   
   ! Discretization matrices for Fourier Pseudospectral
   real, save, allocatable :: D_theta(:,:), D_zeta(:,:)   
      
   ! Source term and coefficients of the drift-kinetic equation
   ! DKE in Legendre basis: 
   !   a_theta(k-1) * u_theta(k-1) + a_zeta(k-1) * u_zeta(k-1) + a(k-1) * u(k-1)
   ! + a_theta(k) * u_theta(k) + a_zeta(k) * u_zeta(k) + a(k) * u(k)
   ! + a_theta(k+1) * u_theta(k+1) + a_zeta(k+1) * u_zeta(k+1)   + a(k+1) * u(k+1)
   ! = s(l).
   real, save, allocatable :: a_theta(:,:,:,:), a_zeta(:,:,:,:), a(:,:,:,:) 
   real, save, allocatable ::  vm(:,:,:)
   real, save, allocatable :: B(:,:), g(:,:), V_prime 
   
 contains  
  
   ! Solves the monoenergetic DKE using an algorithm for inverting
   ! block tridiagonal matrices. N_xi must be greater or equal to 2
   subroutine Solve_BTD_Adjoint_DKE_Legendre( N_theta, N_zeta, N_xi, nu, E_rho, &
                                       Gamma, Gamma_33_Spitzer )
                                       
      integer, intent(in) :: N_theta, N_zeta, N_xi
      real, intent(in) :: nu, E_rho
      real, intent(out) :: Gamma(3,3), Gamma_33_Spitzer 
      
      ! Block matrices for discretization and forward elimination
      real :: L(N_theta*N_zeta,N_theta*N_zeta,1:2)   ! Lower  
      real :: D(N_theta*N_zeta,N_theta*N_zeta,0:2)   ! Diagonal  
      real :: U(N_theta*N_zeta,N_theta*N_zeta,0:2)   ! Upper
      
      real :: Lower(N_theta*N_zeta,N_theta*N_zeta,1:2)   ! Lower for Backward substitution
      real :: Delta(N_theta*N_zeta,N_theta*N_zeta,0:2)   ! Diagonal 
      
      ! Source term for backward substitution
      real, target :: Sigma_1(0:N_theta-1,0:N_zeta-1, 0:2) ! Source for D11 and D31
      real, pointer :: pSigma_1(:,:) 
      real, target :: Sigma_3(0:N_theta-1,0:N_zeta-1, 0:2) ! Source for D13 and D33
      real, pointer :: pSigma_3(:,:) 

      ! Distribution function 
      real, target :: F1(0:N_theta-1,0:N_zeta-1, 0:2), F3(0:N_theta-1,0:N_zeta-1, 0:2) 
      real, pointer :: pF1(:,:),  pF3(:,:) 

      ! Source term for original BTD system
      real, target :: S1(0:N_theta-1,0:N_zeta-1, 0:N_xi), S3(0:N_theta-1,0:N_zeta-1, 0:N_xi) 
      real, pointer :: pS1(:,:), pS3(:,:)
      
      integer :: N_fs, N  
      
      N_fs = N_theta * N_zeta
      N = N_fs * (N_xi+1) 

      ! *** Fixing coefficients and quantities required for the calculation
      call Set_Adjoint_DKE_Grids_and_coefficients( nu, E_rho/psi_p, &
                                           N_theta, N_zeta, N_xi ) 

      ! *** Bounds remapping and source term
      pF1(1:N_fs,0:2) => F1     ; pF3(1:N_fs,0:2) => F3
      pS1(1:N_fs,0:N_xi) => S1  ; pS3(1:N_fs,0:N_xi) => S3
      pSigma_1(1:N_fs,0:2) => Sigma_1 ; pSigma_3(1:N_fs,0:2) => Sigma_3
      
      S1 = -vm ; S1(0,0,0) = 0 ; Sigma_1 = S1(:,:,0:2)
      S3 = 0 ; S3(:,:,1) = - B(0:N_theta-1,0:N_zeta-1) ; Sigma_3 = S3(:,:,0:2)

      ! Solve the drift-kinetic equation with the BTD algorithm
      call Forward_elimination
      call Backward_substitution

      ! Once solved compute the monoenergetic coefficients      
      call Compute_monoenergetic_coefficients( Gamma, Gamma_33_Spitzer )      
      
    contains

      ! *** Subroutine to calculate the monoenergetic coefficients once
      ! the solution W(0:N_theta-1,0:N_zeta-1,0:N_xi) is known.
      subroutine Compute_monoenergetic_coefficients( Gamma, Gamma_33_Spitzer )
        real, intent(out) :: Gamma(3,3), Gamma_33_Spitzer
             
        real :: U1(0:N_theta,0:N_zeta,0:2), U3(0:N_theta,0:N_zeta,0:2)         
        real :: Gamma_11, Gamma_31, Gamma_13, Gamma_33     
        real :: g_Spitzer(0:N_theta,0:N_zeta) ! Spitzer distribution function (only has Legendre mode 1)
        integer :: i, j 
     
        ! *** Legendre modes 1 and 2 for radial transport and bootstrap
        ! in periodic grid with repeated point.   
        do j = 0, N_zeta
           do i = 0, N_theta
              U1(i,j,:) = F1( modulo(i,N_theta), modulo(j,N_zeta), 0:2)
              U3(i,j,:) = F3( modulo(i,N_theta), modulo(j,N_zeta), 0:2)
           end do
        end do
     
        ! *** Flux surface average        
        Gamma_11 = Integral( ["theta","zeta"], g * ( -2*vm(:,:,0) * U1(:,:,0) - 2*vm(:,:,2) * U1(:,:,2)/5 )  ) / V_prime
        Gamma_31 = Integral( ["theta","zeta"], g * 2* U1(:,:,1) * B/3 ) / V_prime 
        Gamma_13 = Integral( ["theta","zeta"], g * ( 2 * vm(:,:,0) * U3(:,:,0) + 2* vm(:,:,2) * U3(:,:,2) / 5 )  ) / V_prime
        Gamma_33 = Integral( ["theta","zeta"], g * 2 * B * U3(:,:,1) /3) / V_prime  
     
        ! ** Scaling to match DKES normalization
        Gamma_11 = Gamma_11 / psi_p**2 
        Gamma_31 = Gamma_31 / (psi_p * B00)  
        Gamma_13 = Gamma_13 / (psi_p * B00)
        Gamma_33 = Gamma_33 / B00**2           
        
        ! Onsager matrix including the relations
        Gamma(1,:) = [ Gamma_11, Gamma_11, Gamma_13 ] 
        Gamma(2,:) = [ Gamma_11, Gamma_11, Gamma_13 ] 
        Gamma(3,:) = [ Gamma_31, Gamma_31, Gamma_33 ]        
        
        ! Spitzer Gamma_33 
        g_Spitzer = - B / nu
        Gamma_33_Spitzer = Integral( ["theta","zeta"], g * 2 * B * g_Spitzer /3) / V_prime 
        Gamma_33_Spitzer = Gamma_33_Spitzer / B00**2          

      end subroutine

      ! Performs the forward elimination storing only the matrices
      ! associated at the modes k=0,1 and 2 to save memory. 
      subroutine Forward_elimination
        integer :: k
         
        integer :: c0, c1, c_rate
        real :: t0, t1,  Delta_k( N_fs, N_fs ), X(N_fs, N_fs) ! Auxiliary variables for calculations 
        
        call system_clock(count_rate=c_rate) ! Get click rate
                
        L(:,:,2) = Block_L(N_xi)
        Delta_k  = Block_D(N_xi)  ! Delta_{N_\xi}       

        ! Computation of Delta_2^{-1} and L_2
        do k = N_xi-1, 2, -1
           
           write(*,*) " Legendre Mode k = ", k                        
           call system_clock(c0)
           ! Solve Delta_{k+1} * X_{k+1} = L_{k+1} for X_{k+1}
           X = Solve_system_LU_LAPACK( Delta_k, L(:,:,2) )
           call system_clock(c1)
           write(*,*) "      Time in Solving Delta * X = L ", ( c1 - c0  ) / real(c_rate)
                              
           call system_clock(c0)
           L(:,:,2) =  Block_L(k) ! L_k for next iteration
           D(:,:,1) =  Block_D(k) ! D_k  
           U(:,:,1) =  Block_U(k) ! U_k             
           call system_clock(c1)
           write(*,*) "      Time in Defining L, D, U = ", ( c1 - c0  ) / real(c_rate)   

           ! Schur complement Delta_k for next iteration
           call system_clock(c0)
           Delta_k = D(:,:,1) - matmul_LAPACK( U(:,:,1), X )          
           call system_clock(c1)
           write(*,*) "      Time in Defining Delta = ", ( c1 - c0  ) / real(c_rate)

           ! Extract Delta_2^{-1}, L_2 for Backward substitution and U_2
           if( k == 2 ) Delta(:,:,k) = Invert_LU_LAPACK( Delta_k )
           if( k == 2 ) Lower(:,:,k) = L(:,:,2)
           if( k == 2 )     U(:,:,k) = U(:,:,1)
           
        end do

        ! Computation of the inverses of Delta_0, Delta_1
        do k = 1, 0, -1          
           write(*,*) " Legendre Mode k = ", k                        
           call system_clock(c0)
           D(:,:,k) =  Block_D(k) ! D_k 
           if( k > 0 ) Lower(:,:,k) =  Block_L(k) ! L_k
           U(:,:,k) =  Block_U(k) ! U_k           
           call system_clock(c1)
           write(*,*) "      Time in Defining L, D, U = ", ( c1 - c0  ) / real(c_rate)
           
           call system_clock(c0)
           Delta_k = D(:,:,k) &
                   - matmul_LAPACK( U(:,:,k), matmul_LAPACK( Delta(:,:,k+1), Lower(:,:,k+1) ) )          
           call system_clock(c1)
           write(*,*) "      Time in Defining Delta = ", ( c1 - c0  ) / real(c_rate)
           
           call system_clock(c0)
           Delta(:,:,k) = Invert_LU_LAPACK( Delta_k )  
           call system_clock(c1)
           write(*,*) "      Time in Inverting Delta = ", ( c1 - c0  ) / real(c_rate)

           ! Source terms for equivalent system 
           pSigma_1(:,k) = pS1(:,k) - matmul( U(:,:,k), matmul( Delta(:,:,k+1), pSigma_1(:,k+1) ) )           
           ! The source term S3 only has Legendre mode 1
           pSigma_3(:,k) = pS3(:,k) - matmul( U(:,:,k), matmul( Delta(:,:,k+1), pSigma_3(:,k+1) ) )
           
        end do
        
      end subroutine
      
      ! Solves the drift-kinetic equation once it is in lower block triangular
      ! form.
      subroutine Backward_substitution       
        integer :: k

        pF1(:,0) = matmul( Delta(:,:,0), pSigma_1(:,0) )
        pF3(:,0) = matmul( Delta(:,:,0), pSigma_3(:,0) )
        do k = 1, 2
           pF1(:,k) = matmul( Delta(:,:,k), pSigma_1(:,k)- matmul( Lower(:,:,k), pF1(:,k-1) ) )
           pF3(:,k) = matmul( Delta(:,:,k), pSigma_3(:,k)- matmul( Lower(:,:,k), pF3(:,k-1) ) ) 
        end do
        
      end subroutine

      ! *** Computes the block matrix L_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function.
      function Block_L(k) result(L)
        integer, intent(in) :: k
        real :: L(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 
        L = 0
        do j = 0, N_zeta-1   ! Loops for rows
           do i = 0, N_theta-1
              ! Row number     
              i_row = 1 + i + j * N_theta                   
              
              ! Diagonal terms 
              L(i_row,i_row) = a_theta(i,j,k,-1) * D_theta(i,i) &
                             +  a_zeta(i,j,k,-1) * D_zeta(j,j)  &
                             + a(i,j,k,-1) 
              
              ! Non-diagonal terms of derivative in theta
              do ii = 0, N_theta-1 ;  i_col = 1 + ii + j * N_theta
                 if( i_row /= i_col ) & 
                 L(i_row,i_col) = a_theta(i,j,k,-1) * D_theta(i,ii)
              end do         
              
              ! Non-diagonal terms of derivative in zeta
              do jj = 0, N_zeta-1 ;  i_col = 1 + i + jj * N_theta
                 if( i_row /= i_col ) & 
                 L(i_row,i_col) = a_zeta(i,j,k,-1) * D_zeta(j,jj)
              end do   
              
           end do
        end do
            
      end function
       
       
      ! *** Computes the block matrix D_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function.
      function Block_D(k) result(D)
        integer, intent(in) :: k
        real :: D(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 
        D = 0
        do j = 0, N_zeta-1   ! Loops for rows
           do i = 0, N_theta-1
              ! Row number     
              i_row = 1 + i + j * N_theta                   
              
              ! Diagonal terms 
              D(i_row,i_row) = a_theta(i,j,k,0) * D_theta(i,i) &
                             + a_zeta(i,j,k,0)  * D_zeta(j,j)  &
                             + a(i,j,k,0) 
              
              ! Non-diagonal terms of derivative in theta
              do ii = 0, N_theta-1 ; i_col = 1 + ii + j * N_theta
                 if( i_row /= i_col ) & 
                 D(i_row,i_col) = a_theta(i,j,k,0) * D_theta(i,ii)
              end do         
              
              ! Non-diagonal terms of derivative in zeta
              do jj = 0, N_zeta-1 ; i_col = 1 + i + jj * N_theta
                 if( i_row /= i_col ) & 
                 D(i_row,i_col) = a_zeta(i,j,k,0) * D_zeta(j,jj)
              end do   
              
           end do
        end do
        
        ! Elimination of nullspace
        if( k==0 ) D(1,1) = 1 ; if( k==0 ) D(1,2:N_fs) = 0
      end function
       
       
      ! *** Computes the block matrix U_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function.
      function Block_U(k) result(U)
        integer, intent(in) :: k
        real :: U(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 
        
        U = 0
        do j = 0, N_zeta-1   ! Loops for rows
           do i = 0, N_theta-1
              ! Row number     
              i_row = 1 + i + j * N_theta                   
              
              ! Diagonal terms 
              U(i_row,i_row) = a_theta(i,j,k,+1) * D_theta(i,i) &
                             + a_zeta(i,j,k,+1)  * D_zeta(j,j)  &
                             + a(i,j,k,+1) 
              
              ! Non-diagonal terms of derivative in theta
              do ii = 0, N_theta-1 ; i_col = 1 + ii + j * N_theta
                 if( i_row /= i_col ) & 
                 U(i_row,i_col) = a_theta(i,j,k,+1) * D_theta(i,ii)
              end do         
              
              ! Non-diagonal terms of derivative in zeta
              do jj = 0, N_zeta-1 ; i_col = 1 + i + jj * N_theta
                 if( i_row /= i_col ) & 
                 U(i_row,i_col) = a_zeta(i,j,k,+1) * D_zeta(j,jj)
              end do   
              
           end do
        end do
        
        ! Elimination of nullspace
        if( k==0 ) U(1,:) = 0
        
      end function
             
     end subroutine
 
  ! *** Initializes the grids and coefficients/magnetic configuration
  ! related quantities required for solving the DKE and posterior computation
  ! of monoenergetic coefficients. 
   subroutine Set_Adjoint_DKE_Grids_and_coefficients( nu, E_psi, &
                                         N_theta, N_zeta, N_xi ) 
      real, intent(in) :: nu, E_psi
      integer, intent(in) :: N_theta, N_zeta, N_xi
      
      logical :: Grids_DKE = .false., Coefficients_DKE = .false. 
      
      ! *** If the variables have been allocated, deallocate them.
      if( Grids_and_Coefficients_DKE ) then
       deallocate( theta, zeta ) 
       deallocate( a_theta, a_zeta, a, B, g, vm ) 
      end if        
            
      ! *** Computation of grids weights and DKE coefficients
      call Set_Phase_space_grids  
      call Set_DKE_Coefficients  
      
      ! *** Setting logical to know state of grids and coefficients
      Grids_and_Coefficients_DKE = Grids_DKE .and. Coefficients_DKE
      
      contains 
      
      subroutine Set_Phase_space_grids
         integer :: i, j, k
         
         allocate( theta(0:N_theta), zeta(0:N_zeta) )
         theta = [( i * 2d0*pi/N_theta, i = 0, N_theta)]
          zeta = [( j * 2d0*pi/(Np*N_zeta), j = 0, N_zeta)]
                  
         ! *** For integrals using trapezoidal rule
         call Grid_initialization( "theta", theta, 1 )
         call Grid_initialization( "zeta", zeta, 1 )

         ! *** Discretization matrices for Pseudospectral Fourier
         if( allocated(D_theta) ) deallocate( D_theta ) 
         if( allocated(D_zeta) )  deallocate( D_zeta ) 
         
         allocate( D_theta(0:N_theta-1,0:N_theta-1), D_zeta(0:N_zeta-1, 0:N_zeta-1) )
         do i = 0, N_theta-1
            D_theta(i,:) = Fourier_Derivative_row( theta(i), theta, 1 )
         end do 
         do j = 0, N_zeta-1
            D_zeta(j,:) = Fourier_Derivative_row( zeta(j), zeta, 1 )
         end do 
                  
         Grids_DKE = .true.
       end subroutine
       
      subroutine Set_DKE_Coefficients
        integer :: i, j, k
        real ::  B2_mean , vds_l(0:N_theta,0:N_zeta) 
        real :: dBdzeta(0:N_theta,0:N_zeta), dBdtheta(0:N_theta,0:N_zeta)
        real :: dBdz_l(0:N_theta,0:N_zeta), dldz_l(0:N_theta,0:N_zeta)                
        
        ! *** Magnetic configuration
        allocate( B(0:N_theta,0:N_zeta), g(0:N_theta,0:N_zeta) )
        do j = 0, N_zeta
           do i = 0, N_theta
              B(i,j) = Magnetic_Field_Boozer( theta(i), zeta(j) ) 
              dBdtheta(i,j) = Magnetic_Field_Boozer_Derivative_theta( theta(i), zeta(j) )
              dBdzeta(i,j) = Magnetic_Field_Boozer_Derivative_zeta( theta(i), zeta(j) )
           end do
        end do
        
        g = abs( B_zeta + iota * B_theta  ) / B**2    ! Jacobian
        dldz_l = abs( B_zeta + iota * B_theta  ) / B  ! Length differential form
        vds_l = ( B_theta * dBdzeta - B_zeta * dBdtheta ) & ! Spatial dependence of 
              / ( dldz_l * B**2 )                           ! magnetic radial drift
        dBdz_l = dBdzeta + iota * dBdtheta  ! Angular derivative of B along field lines
               
        V_prime = Integral( ["theta","zeta"], g ) 
        B2_mean = Integral( ["theta","zeta"], g * B**2 ) / V_prime 
        
        ! *** Coefficients for Adjoint DKE (V+L)f=s
        allocate( a_theta( 0:N_theta, 0:N_zeta, 0:N_xi, -1:1), &
                   a_zeta( 0:N_theta, 0:N_zeta, 0:N_xi, -1:1) )
        allocate( a( 0:N_theta, 0:N_zeta, 0:N_xi, -1:1 )  )
        do k = 0, N_xi
           a_zeta(:,:,k,-1)  = k * 1d0 / ( g*B * (2*k-1)) 
           a_zeta(:,:,k,+1)  = (k+1) * 1d0 / ( g * B * (2*k+3)) 
           a_zeta(:,:,k,0)  = E_psi * B_theta / ( g * B2_mean )
           
           a_theta(:,:,k,-1)  = k * iota / ( g*B * (2*k-1)) 
           a_theta(:,:,k,1)  = (k+1) * iota / ( g * B * (2*k+3)) 
           a_theta(:,:,k,0)  = -E_psi * B_zeta / ( g * B2_mean )            
		   
           a(:,:,k,-1) = (k*(k-1)) * ( dBdzeta + iota * dBdtheta )      &
                       / ( 2*(2*k-1) * g * B**2 ) 
           a(:,:,k,+1) = -((k+1)*(k+2)) * ( dBdzeta + iota * dBdtheta ) &
                       / ( 2*(2*k+3) * g * B**2 )   
           a(:,:,k,0)  =  (-1)*0.5*nu * k * (k+1)              
        end do            
        
        ! *** Legendre modes of the magnetic radial drift 
        !  vm(:,:,k) = vds_l * ( 2 + P_2(xi) ) / 3
        allocate( vm(0:N_theta,0:N_zeta,0:N_xi) )
        vm = 0
        vm(:,:,0) = vds_l * 2 / 3
        vm(:,:,2) = vds_l / 3
        
        Coefficients_DKE = .true.  
        
      end subroutine
      
   end subroutine 
  

end module