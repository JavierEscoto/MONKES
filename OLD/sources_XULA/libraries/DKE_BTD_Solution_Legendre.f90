module DKE_BTD_Solution_Legendre 
   
   
   use Magnetic_configuration
   
   use Finite_Differences
   use Barycentric_grids   
   
   use Linear_systems
   
   implicit none

   private

   public :: Solve_BTD_DKE_Legendre 
   public :: Solve_BTD_DKE_Legendre_DF 
   public :: Monoenergetic_lambda_function
   public :: Monoenergetic_lambda_function_NEW
   
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
   ! block tridiagonal matrices. N_xi must be greater or equal to 2.
   ! 
   ! INPUTS:
   !     N_theta: Number of points in which the poloidal angle is discretized
   !      N_zeta: Number of points in which the toroidal angle is discretized
   !        N_xi: Number of Legendre modes
   !          nu: Collision frequency in s^{-1} divided by the speed in m/s.
   !         E_r: Radial electric field in kV/m divided by the speed in m/s.
   ! 
   ! OUTPUTS:
   !       Gamma: 3x3 matrix containing the monoenergetic coefficients
   !       Gamma_33_Spitzer: Monoenergetic coefficient of Spitzer conductivity.
   ! 
   subroutine Solve_BTD_DKE_Legendre( N_theta, N_zeta, N_xi, nu, E_r, &
                                       Gamma, Gamma_33_Spitzer ) 
                                       
      integer, intent(in) :: N_theta, N_zeta, N_xi
      real, intent(in) :: nu, E_r
      real, intent(out) :: Gamma(3,3), Gamma_33_Spitzer 
      
      ! Intermediate block matrices for discretization and Forward elimination
      real :: L( N_theta*N_zeta, N_theta*N_zeta )   ! Lower  
      real :: D( N_theta*N_zeta, N_theta*N_zeta )   ! Diagonal  
      real :: U( N_theta*N_zeta, N_theta*N_zeta )   ! Upper
      
      ! Block matrices to be stored for Backward substitution
      real :: Lower( N_theta*N_zeta, N_theta*N_zeta, 1:2 )   ! Lower 
      real :: Delta( N_theta*N_zeta, N_theta*N_zeta, 0:2 )   ! Diagonal 
      real :: Upper( N_theta*N_zeta, N_theta*N_zeta, 0:2 )   ! Upper 
      
      ! Source terms for Backward substitution
      real, target :: Sigma_1(0:N_theta-1,0:N_zeta-1, 0:2) ! Source for D11 and D31
      real, target :: Sigma_3(0:N_theta-1,0:N_zeta-1, 0:2) ! Source for D13 and D33
      real, pointer :: pSigma_1(:,:), pSigma_3(:,:) 

      ! Distribution functions 
      real, target :: F1(0:N_theta-1,0:N_zeta-1, 0:2), F3(0:N_theta-1,0:N_zeta-1, 0:2) 
      real, pointer :: pF1(:,:),  pF3(:,:) 

      ! Source terms for original BTD system
      real, target :: S1(0:N_theta-1,0:N_zeta-1, 0:N_xi), S3(0:N_theta-1,0:N_zeta-1, 0:N_xi) 
      real, pointer :: pS1(:,:), pS3(:,:)
      
      integer :: N_fs, ipiv(N_theta*N_zeta,0:2)
      
      N_fs = N_theta * N_zeta  ! Number of points for discretizing the flux-surface

      ! *** Fixing coefficients and quantities required for the calculation
      call Set_DKE_Grids_and_coefficients( nu, E_r/psi_p, &
                                           N_theta, N_zeta, N_xi ) 

      ! *** Bounds remapping for distribution functions and source terms
      pF1(1:N_fs,0:2) => F1     ; pF3(1:N_fs,0:2) => F3
      pS1(1:N_fs,0:N_xi) => S1  ; pS3(1:N_fs,0:N_xi) => S3
      pSigma_1(1:N_fs,0:2) => Sigma_1 ; pSigma_3(1:N_fs,0:2) => Sigma_3
      
      S1 = -vm ; S1(0,0,0) = 0 ; Sigma_1 = S1(:,:,0:2)
      S3 = 0 ; S3(:,:,1) = B(0:N_theta-1,0:N_zeta-1) ; Sigma_3 = S3(:,:,0:2)
 
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
        Gamma_11 = Integral( ["theta","zeta "], g * ( -2*vm(:,:,0) * U1(:,:,0) - 2*vm(:,:,2) * U1(:,:,2)/5 )  ) / V_prime
        Gamma_31 = Integral( ["theta","zeta "], g * 2* U1(:,:,1) * B/3 ) / V_prime 
        Gamma_13 = Integral( ["theta","zeta "], g * ( -2 * vm(:,:,0) * U3(:,:,0) - 2* vm(:,:,2) * U3(:,:,2) / 5 )  ) / V_prime
        Gamma_33 = Integral( ["theta","zeta "], g * 2 * B * U3(:,:,1) /3) / V_prime  
     
        ! ** Scaling to match DKES normalization
        Gamma_11 = Gamma_11 / psi_p**2 
        Gamma_31 = Gamma_31 / (psi_p * B00)  
        Gamma_13 = Gamma_13 / (psi_p * B00)
        Gamma_33 = Gamma_33 / B00**2           
        
        ! Onsager matrix including the relations D_{2j} = D_{1j}, D_{i2}= D_{i1}
        Gamma(1,:) = [ Gamma_11, Gamma_11, Gamma_13 ] 
        Gamma(2,:) = [ Gamma_11, Gamma_11, Gamma_13 ] 
        Gamma(3,:) = [ Gamma_31, Gamma_31, Gamma_33 ]        
        
        ! Spitzer Gamma_33 
        g_Spitzer = B / nu
        Gamma_33_Spitzer = Integral( ["theta","zeta "], g * 2 * B * g_Spitzer /3) / V_prime 
        Gamma_33_Spitzer = Gamma_33_Spitzer / B00**2          

      end subroutine

      ! Performs the forward elimination storing only the matrices
      ! associated at the modes k=0,1 and 2 to save memory. 
      subroutine Forward_elimination
        integer :: k
         
        integer :: c0, c1, c_rate, ipiv_temp(N_fs)
        real :: t0, t1, Delta_k(N_fs, N_fs), X(N_fs, N_fs) ! Auxiliary variables for calculations 
        real :: y1(N_fs,1), y3(N_fs,1)
        
        call system_clock(count_rate=c_rate) ! Get click rate  
        ! Initialize matrices as zero
        L = 0 ; D = 0 ; U = 0 ; Delta_k = 0
        
        call Block_L(N_xi, L)        ! L_{N_\xi}  
        call Block_D(N_xi, Delta_k)  ! Delta_{N_\xi}    
        call Solve_Matrix_System_LAPACK( Delta_k, L, .false., ipiv_temp, X ) ! X_{N_\xi}     
 
        ! Computation of Delta_2^{-1} and L_2
        do k = N_xi-1, 0, -1
                                
           write(*,*) " Legendre Mode k = ", k  
           call system_clock(c0)
           call Block_L(k, L) ! L_k    
           call Block_D(k, D) ! D_k  
           call Block_U(k, U) ! U_k             
           call system_clock(c1)
           write(*,*) "      Time in Defining L, D, U = ", ( c1 - c0  ) / real(c_rate)               

           ! Schur complement Delta_k = D_k - U_k X_{k+1}
           call system_clock(c0)
           Delta_k = D       
           call dgemm_LAPACK( U, X, Delta_k, -1d0, 1d0 )      
           call system_clock(c1)
           write(*,*) "      Time in Defining Delta = ", ( c1 - c0  ) / real(c_rate)  
           
           ! Solve Delta_{k} * X_{k} = L_{k} for X_{k} for next iteration                         
           call system_clock(c0)        
           if(k>0) call Solve_Matrix_System_LAPACK( Delta_k, L, .false., ipiv_temp, X )
           call system_clock(c1)
           write(*,*) "      Time in Solving Delta * X = L ", ( c1 - c0  ) / real(c_rate) 
           
           ! Extract required matrices for sources and backward substitution
           if( k <= 2 ) then
              if(k>0) Lower(:,:,k) = L
              Delta(:,:,k) = Delta_k
              Upper(:,:,k) = U
              ipiv(:,k) = ipiv_temp      ! LU factorization information
           end if 
           
        end do
        
        ! *** Source terms for equivalent lower triangular system are given by formula: 
        !   sigma^k = s^k - Uk y^{k+1},
        !   where Delta_{k+1} y^{k+1} = sigma^{k+1}
        pSigma_1(:,2) = pS1(:,2) ! S1 has only Legendre modes 0 and 2
        pSigma_3(:,1) = pS3(:,1) ! S3 only has Legendre mode 1        
        do k = 1, 0, -1           
           ! Solve Delta_{k+1} y^{k+1} = sigma^{k+1} for y^{k+1} reusing LU factorization
           call Solve_Matrix_System_LAPACK( Delta(:,:,k+1), pSigma_1(:,k+1:k+1), .true., ipiv(:,k+1), y1 )           
           if( k == 0 ) &       
           call Solve_Matrix_System_LAPACK( Delta(:,:,k+1), pSigma_3(:,k+1:k+1), .true., ipiv(:,k+1), y3 )              
           
           ! Compute sources for Backward substitution
           pSigma_1(:,k) = pS1(:,k) - matmul( Upper(:,:,k), y1(:,1) )            
           if( k == 0 ) &  
           pSigma_3(:,k) = - matmul( Upper(:,:,k), y3(:,1) )
        end do               
        
      end subroutine
      
      ! Solves the drift-kinetic equation once it is in lower block triangular
      ! form.
      subroutine Backward_substitution      
        integer :: k
        
        real :: y1(N_fs,1), y3(N_fs,1)
        
        ! Solve Delta_0 f^{0} = sigma^{0} reusing LU factorization when possible
        call Solve_Matrix_System_LAPACK( Delta(:,:,0), pSigma_1(:,0:0), .false., ipiv(:,0), pF1(:,0:0) )
        call Solve_Matrix_System_LAPACK( Delta(:,:,0), pSigma_3(:,0:0),  .true., ipiv(:,0), pF3(:,0:0) ) 
        do k = 1, 2
           y1(:,1) = pSigma_1(:,k)- matmul( Lower(:,:,k), pF1(:,k-1) )
           y3(:,1) = pSigma_3(:,k)- matmul( Lower(:,:,k), pF3(:,k-1) )
           call Solve_Matrix_System_LAPACK( Delta(:,:,k), y1, .true., ipiv(:,k), pF1(:,k:k) )
           call Solve_Matrix_System_LAPACK( Delta(:,:,k), y3, .true., ipiv(:,k), pF3(:,k:k) ) 
        end do
        
      end subroutine
      
      ! *** Computes the block matrix L_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function. 
      ! It requires to be given a zero matrix as it fills only the required
      ! positions
      subroutine Block_L(k, L)
        integer, intent(in) :: k
        real, intent(inout) :: L(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 
        
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
            
      end subroutine
       
      ! *** Computes the block matrix D_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function. 
      ! It requires to be given a zero matrix as it fills only the required
      ! positions 
      subroutine Block_D(k, D)
        integer, intent(in) :: k 
        real, intent(inout) :: D(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 

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
      end subroutine
       
      ! *** Computes the block matrix U_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function. 
      ! It requires to be given a zero matrix as it fills only the required
      ! positions
      subroutine Block_U(k, U)
        integer, intent(in) :: k
        real, intent(inout) :: U(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 
        
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
        
      end subroutine
             
  end subroutine
 
  
  
   ! Solves the monoenergetic DKE using an algorithm for inverting
   ! block tridiagonal matrices. N_xi must be greater or equal to 2.
   ! 
   ! INPUTS:
   !     N_theta: Number of points in which the poloidal angle is discretized
   !      N_zeta: Number of points in which the toroidal angle is discretized
   !        N_xi: Number of Legendre modes
   !          nu: Collision frequency in s^{-1} divided by the speed in m/s.
   !         E_r: Radial electric field in kV/m divided by the speed in m/s.
   ! 
   ! OUTPUTS:
   !       Gamma: 3x3 matrix containing the monoenergetic coefficients
   !       Gamma_33_Spitzer: Monoenergetic coefficient of Spitzer conductivity.
   ! 
   subroutine Solve_BTD_DKE_Legendre_DF( N_theta, N_zeta, N_xi, nu, E_r, &
                                         Gamma, Gamma_33_Spitzer, M_xi, F1, F3 ) 
                                       
      integer, intent(in) :: N_theta, N_zeta, N_xi
      real, intent(in) :: nu, E_r
      real, intent(out) :: Gamma(3,3), Gamma_33_Spitzer 
      integer, intent(in) :: M_xi
      real, target, intent(out) :: F1(0:N_theta-1,0:N_zeta-1,0:M_xi), &  
                                   F3(0:N_theta-1,0:N_zeta-1,0:M_xi)  
      
      ! Intermediate block matrices for discretization and Forward elimination
      real :: L( N_theta*N_zeta, N_theta*N_zeta )   ! Lower  
      real :: D( N_theta*N_zeta, N_theta*N_zeta )   ! Diagonal  
      real :: U( N_theta*N_zeta, N_theta*N_zeta )   ! Upper
      
      ! Block matrices to be stored for Backward substitution
      real :: Lower( N_theta*N_zeta, N_theta*N_zeta, 1:M_xi )   ! Lower 
      real :: Delta( N_theta*N_zeta, N_theta*N_zeta, 0:M_xi )   ! Diagonal 
      real :: Upper( N_theta*N_zeta, N_theta*N_zeta, 0:M_xi )   ! Upper 
      
      ! Source terms for Backward substitution
      real, target :: Sigma_1(0:N_theta-1,0:N_zeta-1, 0:M_xi) ! Source for D11 and D31
      real, target :: Sigma_3(0:N_theta-1,0:N_zeta-1, 0:M_xi) ! Source for D13 and D33
      real, pointer :: pSigma_1(:,:), pSigma_3(:,:) 

      ! Distribution functions as vectors
      real, pointer :: pF1(:,:),  pF3(:,:) 

      ! Source terms for original BTD system
      real, target :: S1(0:N_theta-1,0:N_zeta-1, 0:N_xi), S3(0:N_theta-1,0:N_zeta-1, 0:N_xi) 
      real, pointer :: pS1(:,:), pS3(:,:)
      
      integer :: N_fs, ipiv(N_theta*N_zeta,0:M_xi)
      
      N_fs = N_theta * N_zeta  ! Number of points for discretizing the flux-surface

      ! *** Fixing coefficients and quantities required for the calculation
      call Set_DKE_Grids_and_coefficients( nu, E_r/psi_p, &
                                           N_theta, N_zeta, N_xi ) 

      ! *** Bounds remapping for distribution functions and source terms
      pF1(1:N_fs,0:M_xi) => F1     ; pF3(1:N_fs,0:M_xi) => F3
      pS1(1:N_fs,0:N_xi) => S1  ; pS3(1:N_fs,0:N_xi) => S3
      pSigma_1(1:N_fs,0:M_xi) => Sigma_1 ; pSigma_3(1:N_fs,0:M_xi) => Sigma_3
      
      Sigma_1 = 0 ; Sigma_3 = 0 
      S1 = -vm ; S1(0,0,0) = 0 ; Sigma_1(:,:,0:2) = S1(:,:,0:2)
      S3 = 0 ; S3(:,:,1) = B(0:N_theta-1,0:N_zeta-1) 
      Sigma_3(:,:,0:2) = S3(:,:,0:2)
 
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
        Gamma_11 = Integral( ["theta","zeta "], g * ( -2*vm(:,:,0) * U1(:,:,0) - 2*vm(:,:,2) * U1(:,:,2)/5 )  ) / V_prime
        Gamma_31 = Integral( ["theta","zeta "], g * 2* U1(:,:,1) * B/3 ) / V_prime 
        Gamma_13 = Integral( ["theta","zeta "], g * ( -2 * vm(:,:,0) * U3(:,:,0) - 2* vm(:,:,2) * U3(:,:,2) / 5 )  ) / V_prime
        Gamma_33 = Integral( ["theta","zeta "], g * 2 * B * U3(:,:,1) /3) / V_prime  
     
        ! ** Scaling to match DKES normalization
        Gamma_11 = Gamma_11 / psi_p**2 
        Gamma_31 = Gamma_31 / (psi_p * B00)  
        Gamma_13 = Gamma_13 / (psi_p * B00)
        Gamma_33 = Gamma_33 / B00**2           
        
        ! Onsager matrix including the relations D_{2j} = D_{1j}, D_{i2}= D_{i1}
        Gamma(1,:) = [ Gamma_11, Gamma_11, Gamma_13 ] 
        Gamma(2,:) = [ Gamma_11, Gamma_11, Gamma_13 ] 
        Gamma(3,:) = [ Gamma_31, Gamma_31, Gamma_33 ]        
        
        ! Spitzer Gamma_33 
        g_Spitzer = B / nu
        Gamma_33_Spitzer = Integral( ["theta","zeta "], g * 2 * B * g_Spitzer /3) / V_prime 
        Gamma_33_Spitzer = Gamma_33_Spitzer / B00**2          

      end subroutine

      ! Performs the forward elimination storing only the matrices
      ! associated at the modes k=0,1 and 2 to save memory. 
      subroutine Forward_elimination
        integer :: k
         
        integer :: c0, c1, c_rate, ipiv_temp(N_fs)
        real :: t0, t1, Delta_k(N_fs, N_fs), X(N_fs, N_fs) ! Auxiliary variables for calculations 
        real :: y1(N_fs,1), y3(N_fs,1)
        
        call system_clock(count_rate=c_rate) ! Get click rate  
        ! Initialize matrices as zero
        L = 0 ; D = 0 ; U = 0 ; Delta_k = 0
        
        call Block_L(N_xi, L)        ! L_{N_\xi}  
        call Block_D(N_xi, Delta_k)  ! Delta_{N_\xi}    
        call Solve_Matrix_System_LAPACK( Delta_k, L, .false., ipiv_temp, X ) ! X_{N_\xi}     
 
        ! Computation of Delta_k
        do k = N_xi-1, 0, -1
                                
           write(*,*) " Legendre Mode k = ", k  
           call system_clock(c0)
           call Block_L(k, L) ! L_k    
           call Block_D(k, D) ! D_k  
           call Block_U(k, U) ! U_k             
           call system_clock(c1)
           write(*,*) "      Time in Defining L, D, U = ", ( c1 - c0  ) / real(c_rate)               

           ! Schur complement Delta_k = D_k - U_k X_{k+1}
           call system_clock(c0)
           Delta_k = D       
           call dgemm_LAPACK( U, X, Delta_k, -1d0, 1d0 )      
           call system_clock(c1)
           write(*,*) "      Time in Defining Delta = ", ( c1 - c0  ) / real(c_rate)  
           
           ! Solve Delta_{k} * X_{k} = L_{k} for X_{k} for next iteration                         
           call system_clock(c0)        
           if(k>0) call Solve_Matrix_System_LAPACK( Delta_k, L, .false., ipiv_temp, X )
           call system_clock(c1)
           write(*,*) "      Time in Solving Delta * X = L ", ( c1 - c0  ) / real(c_rate) 
           
           ! Extract required matrices for sources and backward substitution
           if( k <= M_xi ) then
              if(k>0) Lower(:,:,k) = L
              Delta(:,:,k) = Delta_k
              Upper(:,:,k) = U
              ipiv(:,k) = ipiv_temp      ! LU factorization information
           end if 
           
        end do
        
        ! *** Source terms for equivalent lower triangular system are given by formula: 
        !   sigma^k = s^k - Uk y^{k+1},
        !   where Delta_{k+1} y^{k+1} = sigma^{k+1}
        pSigma_1(:,2) = pS1(:,2) ! S1 has only Legendre modes 0 and 2
        pSigma_3(:,1) = pS3(:,1) ! S3 only has Legendre mode 1        
        do k = 1, 0, -1           
           ! Solve Delta_{k+1} y^{k+1} = sigma^{k+1} for y^{k+1} reusing LU factorization
           call Solve_Matrix_System_LAPACK( Delta(:,:,k+1), pSigma_1(:,k+1:k+1), .true., ipiv(:,k+1), y1 )           
           if( k == 0 ) &       
           call Solve_Matrix_System_LAPACK( Delta(:,:,k+1), pSigma_3(:,k+1:k+1), .true., ipiv(:,k+1), y3 )              
           
           ! Compute sources for Backward substitution
           pSigma_1(:,k) = pS1(:,k) - matmul( Upper(:,:,k), y1(:,1) )            
           if( k == 0 ) &  
           pSigma_3(:,k) = - matmul( Upper(:,:,k), y3(:,1) )
        end do               
        
      end subroutine
      
      ! Solves the drift-kinetic equation once it is in lower block triangular
      ! form.
      subroutine Backward_substitution      
        integer :: k
        
        real :: y1(N_fs,1), y3(N_fs,1)
        
        ! Solve Delta_0 f^{0} = sigma^{0} reusing LU factorization when possible
        call Solve_Matrix_System_LAPACK( Delta(:,:,0), pSigma_1(:,0:0), .false., ipiv(:,0), pF1(:,0:0) )
        call Solve_Matrix_System_LAPACK( Delta(:,:,0), pSigma_3(:,0:0),  .true., ipiv(:,0), pF3(:,0:0) ) 
        do k = 1, M_xi
           y1(:,1) = pSigma_1(:,k)- matmul( Lower(:,:,k), pF1(:,k-1) )
           y3(:,1) = pSigma_3(:,k)- matmul( Lower(:,:,k), pF3(:,k-1) )
           call Solve_Matrix_System_LAPACK( Delta(:,:,k), y1, .true., ipiv(:,k), pF1(:,k:k) )
           call Solve_Matrix_System_LAPACK( Delta(:,:,k), y3, .true., ipiv(:,k), pF3(:,k:k) ) 
        end do
        
      end subroutine
      
      ! *** Computes the block matrix L_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function. 
      ! It requires to be given a zero matrix as it fills only the required
      ! positions
      subroutine Block_L(k, L)
        integer, intent(in) :: k
        real, intent(inout) :: L(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 
        
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
            
      end subroutine
       
      ! *** Computes the block matrix D_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function. 
      ! It requires to be given a zero matrix as it fills only the required
      ! positions 
      subroutine Block_D(k, D)
        integer, intent(in) :: k 
        real, intent(inout) :: D(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 

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
      end subroutine
       
      ! *** Computes the block matrix U_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function. 
      ! It requires to be given a zero matrix as it fills only the required
      ! positions
      subroutine Block_U(k, U)
        integer, intent(in) :: k
        real, intent(inout) :: U(N_fs, N_fs)
        
        integer :: i, j, i_row, i_col, ii, jj 
        
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
        
      end subroutine
             
  end subroutine
   
   
   
  subroutine Monoenergetic_lambda_function_NEW( N_lambda, N_theta, N_zeta, &
                                                N_xi, M_xi, F1, F3,        &
                                                lambda, lambda_c, Gamma_ij )                                          
      integer, intent(in) :: N_lambda, N_theta, N_zeta, N_xi, M_xi
      real, intent(in)  :: F1(0:N_theta-1,0:N_zeta-1,0:M_xi), &
                           F3(0:N_theta-1,0:N_zeta-1,0:M_xi) 
      real, intent(out) :: lambda(0:N_lambda), lambda_c, Gamma_ij(3,3,0:N_lambda)  
      
      real :: xi(0:N_theta,0:N_zeta,0:N_lambda) 
      real :: U1(0:N_theta,0:N_zeta,0:M_xi), U3(0:N_theta,0:N_zeta,0:M_xi) 
      real :: H1(0:N_theta,0:N_zeta,0:N_lambda,0:2) ! Weight function for integrals
      real :: H3(0:N_theta,0:N_zeta,0:N_lambda,0:2) ! Weight function for integrals
      real :: Bmin, Bmax 
      real :: Ik(0:N_theta, 0:N_zeta, 0:N_lambda, 0:M_xi ), &
              Jk(0:N_theta, 0:N_zeta, 0:N_lambda, 0:M_xi )
              
      integer :: i, j, k, kk, c0, c1, c_rate
      
      call system_clock(count_rate=c_rate) ! Get click rate
      call system_clock(c0)
      
      ! *** Legendre modes in periodic grid with repeated point.  
      do j = 0, N_zeta
         do i = 0, N_theta
            U1(i,j,:) = F1( modulo(i,N_theta), modulo(j,N_zeta), :)
            U3(i,j,:) = F3( modulo(i,N_theta), modulo(j,N_zeta), :)                
         end do
      end do                                
                           
      ! Define lambda grid
      Bmin =minval(B) ; Bmax = maxval(B) ; lambda_c = 1d0 /Bmax
      lambda = [( k*(1d0/Bmin)/N_lambda, k=0,N_lambda )]            
      
      ! *** Pitch-angle cosine as a function of lambda, B.
      xi = 0 
      do k = 0, N_lambda            
         where( 1 - lambda(k) * B >= 0 ) &
            xi(:,:,k) = sqrt( 1 - lambda(k) * B )
      end do
      
      ! *** Construction of H_j^{(0)}, H3^{(1)}, H1^{(2)}
      Ik = 0 ; Jk = 0 
      do k = 0, N_lambda
         do j = 0, N_zeta 
            do i = 0, N_theta
               if (  1 - lambda(k) * B(i,j) >= 0 ) &
                 call Compute_Legendre_Moments( xi(i,j,k), Ik(i,j,k,:), Jk(i,j,k,:) )
                 
            end do
         end do
      end do
           
      H1 = 0 ; H3 = 0
      ! Hj^{(0)} and Hj^{(2)}
      do kk = 0, M_xi, 2
         do k = 0, N_lambda
            H1(:,:,k,0) = H1(:,:,k,0) + U1(:,:,kk) * Ik(:,:,k,kk)
            H1(:,:,k,2) = H1(:,:,k,2) &
                        + U1(:,:,kk) * ( 3 * Jk(:,:,k,kk) - Ik(:,:,k,kk) )/2
                        
            H3(:,:,k,0) = H3(:,:,k,0) + U3(:,:,kk) * Ik(:,:,k,kk)
            H3(:,:,k,2) = H3(:,:,k,2) &
                        + U3(:,:,kk) * ( 3 * Jk(:,:,k,kk) - Ik(:,:,k,kk) )/2
         end do 
      end do 
      ! Hj^{(1)}  
      do kk = 1, M_xi, 2
         do k = 0, N_lambda
            H1(:,:,k,1) = H1(:,:,k,1) + U1(:,:,kk) * Jk(:,:,k,kk)
            H3(:,:,k,1) = H3(:,:,k,1) + U3(:,:,kk) * Jk(:,:,k,kk)
         end do 
      end do 
      
      Gamma_ij = 0
      do k = 0, N_lambda
         Gamma_ij(1,1,k) = Integral( ["theta","zeta "], g * ( -vm(:,:,0) * H1(:,:,k,0) - vm(:,:,2) * H1(:,:,k,2) )  ) / V_prime        
         Gamma_ij(1,3,k) = Integral( ["theta","zeta "], g * ( -vm(:,:,0) * H3(:,:,k,0) - vm(:,:,2) * H3(:,:,k,2) )  ) / V_prime
         
         Gamma_ij(3,1,k) = Integral( ["theta","zeta "], g * B * H1(:,:,k,1)  ) / V_prime
         Gamma_ij(3,3,k) = Integral( ["theta","zeta "], g * B * H3(:,:,k,1)  ) / V_prime         
      end do 
      Gamma_ij(1,1,:) = Gamma_ij(1,1,:) / (psi_p**2)
      Gamma_ij(1,3,:) = Gamma_ij(1,3,:) / (psi_p*B00)
      Gamma_ij(3,1,:) = Gamma_ij(3,1,:) / (psi_p*B00)
      Gamma_ij(3,3,:) = Gamma_ij(3,3,:) / (B00**2)

      call system_clock(c1)
      write(*,*) "      Time in Monoenergetic_lambda_function = ", ( c1 - c0  ) / real(c_rate)  
      
      contains
            
      ! *** Computes the integrals 2*\int_{x0}^{t} P_{k} \dd{\xi} for k<= M_xi
      ! where t >= 0. 
      ! For  odd k, x0 = -1   
      ! For even k, x0 = 0 (is equal to the integral starting at x0=-1) 
      function Integral_I(t) result(I_k)
        real, intent(in) :: t 
        real :: I_k(0:M_xi+3)
        
        real :: P(0:M_xi+4) ! Legendre polynomials
        integer :: k
        
        ! Legendre polynomials 
        P(0) = 1 ; P(1) = t
        do k = 2, M_xi+4
           P(k) = Legendre_3term( k, t, P(k-1), P(k-2))
        end do 
        
        ! Integral 2*\int_{0}^{t} P_{k} \dd{\xi} (up to a constant for odd k) 
        I_k(0) = 2*t  
        I_k(1:M_xi+3) = [( 2*(P(k+1) - P(k-1))/(2*k+1),  k=1,M_xi+3)]        
                
      end function 
      
      ! Computes the integrals:
      !  2*\int_{0}^{t} P_{k} \dd{\xi} = I(k)
      !  2*\int_{0}^{t} \xi P_{2k+1} \dd{\xi} = J(2k+1)
      !  2*\int_{0}^{t} \xi^2 P_{2k} \dd{\xi} = J(2k)
      ! for k=0,1,... M_xi.
      subroutine Compute_Legendre_Moments( t, I, J )
         real, intent(in) :: t
         real, intent(out) :: I(0:M_xi), J(0:M_xi)
         
         real :: II(0:M_xi+3), JJ(0:M_xi+2)
         integer :: k
         
         ! Integral 2*\int_{0}^{t} P_{k} \dd{\xi} (up to a constant for odd k) 
         II = Integral_I(t) 
                   
         ! Integral 2*\int_{0}^{t} \xi P_{2k+1} \dd{\xi} 
         do k = 1, M_xi+2, 2
            JJ(k) = ( (k+1)* II(k+1) + k * II(k-1) ) / (2*k + 1)          
         end do 
         
         ! Integral 2*\int_{0}^{t} \xi^2 P_{2k} \dd{\xi}
         JJ(0) = JJ(1)
         do k = 2, M_xi, 2
            JJ(k) = ( (k+1)*JJ(k+1) + k*JJ(k-1) )/(2*k+1)
         end do 
         
         ! Output vectors
         I(0:M_xi) = II(0:M_xi) ; J(0:M_xi) = JJ(0:M_xi) 
      
      end subroutine 
      
      ! Computes P_{k}(x) given k, x and P1 = P_{k-1}(x) and P2 = P_{k-2}(x)
      real function Legendre_3term( k, x, P1, P2) result(P) 
        integer, intent(in) :: k
        real, intent(in) :: x, P1, P2
        
        if( k==0 .or. x==1 ) then
          P = 1
        elseif (k==1) then
          P = x
        else
          P = ( (2*k-1)*x*P1 - (k-1)*P2 )/k
        end if
      end function                      
      
    end subroutine
    
    subroutine Monoenergetic_lambda_function( N_lambda, N_theta, N_zeta, &
                                             N_xi, M_xi, nu, E_r,       &
                                             lambda, lambda_c, &
                                             Gamma_ij, Gamma_33_Spitzer )                                          
      integer, intent(in) :: N_lambda, N_theta, N_zeta, N_xi, M_xi
      real, intent(in) :: nu, E_r
      real, intent(out) :: lambda(0:N_lambda), lambda_c  
      real, intent(out) :: Gamma_ij(3,3,0:N_lambda), Gamma_33_Spitzer  
      
      real :: xi(0:N_theta,0:N_zeta,0:N_lambda) 
      real :: U1(0:N_theta,0:N_zeta,0:M_xi), U3(0:N_theta,0:N_zeta,0:M_xi) 
      real :: F1(0:N_theta-1,0:N_zeta-1,0:M_xi), F3(0:N_theta-1,0:N_zeta-1,0:M_xi) 
      real :: H1(0:N_theta,0:N_zeta,0:N_lambda,0:2) ! Weight function for integrals
      real :: H3(0:N_theta,0:N_zeta,0:N_lambda,0:2) ! Weight function for integrals
      real :: Bmin, Bmax 
      real :: Ik(0:N_theta, 0:N_zeta, 0:N_lambda, 0:M_xi ), &
              Jk(0:N_theta, 0:N_zeta, 0:N_lambda, 0:M_xi )
              
      real :: M_Legendre(0:M_xi,0:M_xi) 
      integer :: i, j, k, kk, c0, c1, c_rate
      
      ! Solve the DKE for the distribution function and monoenergetic coefficients
      call Solve_BTD_DKE_Legendre_DF( N_theta, N_zeta, N_xi, nu, E_r, &
                                      Gamma_ij(:,:,0), Gamma_33_Spitzer, M_xi, F1, F3 ) 
      
      
      call system_clock(count_rate=c_rate) ! Get click rate
      call system_clock(c0) 
      do k = 0, M_xi
         write(*,*) sum( abs(F1(:,:,k)) )/size(F1(:,:,k)), &
                    sum( abs(F3(:,:,k)) )/size(F3(:,:,k))
      
      end do ; write(*,*)
      ! *** Legendre modes in periodic grid with repeated point.  
      do j = 0, N_zeta
         do i = 0, N_theta
            U1(i,j,:) = F1( modulo(i,N_theta), modulo(j,N_zeta), :)
            U3(i,j,:) = F3( modulo(i,N_theta), modulo(j,N_zeta), :)                
         end do
      end do           
                     
                           
      ! Define lambda grid
      Bmin =minval(B) ; Bmax = maxval(B) ; lambda_c = 1d0 /Bmax
      lambda = [( k*(1d0/Bmin)/N_lambda, k=0,N_lambda )]            
      
      ! *** Pitch-angle cosine as a function of lambda, B.
      xi = 0 
      do k = 0, N_lambda            
         where( 1 - lambda(k) * B >= 0 ) &
            xi(:,:,k) = sqrt( 1 - lambda(k) * B )
      end do
      
      ! *** Matrix for evaluating Legendre polynomials in a Chebyshev basis
      M_Legendre = Legendre_Chebyshev_Matrix(M_xi)  
      
      ! *** Construction of H_j^{(0)}, H3^{(1)}, H1^{(2)}
      Ik = 0 ; Jk = 0 
      do k = 0, N_lambda
         do j = 0, N_zeta 
            do i = 0, N_theta
               if (  1 - lambda(k) * B(i,j) >= 0 ) &
                 call Compute_Legendre_Moments( xi(i,j,k), Ik(i,j,k,:), Jk(i,j,k,:) )
                 
            end do
         end do
      end do
            
      
      H1 = 0 ; H3 = 0
      ! Hj^{(0)} and Hj^{(2)}
      do kk = 0, M_xi, 2
         do k = 0, N_lambda
            H1(:,:,k,0) = H1(:,:,k,0) + U1(:,:,kk) * Ik(:,:,k,kk)
            H1(:,:,k,2) = H1(:,:,k,2) &
                        + U1(:,:,kk) * ( 3 * Jk(:,:,k,kk) - Ik(:,:,k,kk) )/2
                        
            H3(:,:,k,0) = H3(:,:,k,0) + U3(:,:,kk) * Ik(:,:,k,kk)
            H3(:,:,k,2) = H3(:,:,k,2) &
                        + U3(:,:,kk) * ( 3 * Jk(:,:,k,kk) - Ik(:,:,k,kk) )/2
         end do 
      end do 
      ! Hj^{(1)}  
      do kk = 1, M_xi, 2
         do k = 0, N_lambda
            H1(:,:,k,1) = H1(:,:,k,1) + U1(:,:,kk) * Jk(:,:,k,kk)
            H3(:,:,k,1) = H3(:,:,k,1) + U3(:,:,kk) * Jk(:,:,k,kk)
         end do 
      end do 
      
      Gamma_ij = 0
      do k = 0, N_lambda
         Gamma_ij(1,1,k) = Integral( ["theta","zeta "], g * ( -vm(:,:,0) * H1(:,:,k,0) - vm(:,:,2) * H1(:,:,k,2) )  ) / V_prime        
         Gamma_ij(1,3,k) = Integral( ["theta","zeta "], g * ( -vm(:,:,0) * H3(:,:,k,0) - vm(:,:,2) * H3(:,:,k,2) )  ) / V_prime
         
         Gamma_ij(3,1,k) = Integral( ["theta","zeta "], g * B * H1(:,:,k,1)  ) / V_prime
         Gamma_ij(3,3,k) = Integral( ["theta","zeta "], g * B * H3(:,:,k,1)  ) / V_prime         
      end do 
      Gamma_ij(1,1,:) = Gamma_ij(1,1,:) / (psi_p**2)
      Gamma_ij(1,3,:) = Gamma_ij(1,3,:) / (psi_p*B00)
      Gamma_ij(3,1,:) = Gamma_ij(3,1,:) / (psi_p*B00)
      Gamma_ij(3,3,:) = Gamma_ij(3,3,:) / (B00**2)

      call system_clock(c1)
      write(*,*) "      Time in Monoenergetic_lambda_function = ", ( c1 - c0  ) / real(c_rate)  
      
      contains
      
      
      ! Computes the matrix  M that relates 
      ! g = M * f, where f and g are respectively the coefficients of
      ! a function in Legendre and Chebyshev (second kind) polynomials.
      ! f(x) = \sum_k f_k P_k(x) = \sum_k g_k T_k(x)
      ! Remark: T_k(x) = cos( k*acos(x) )
      function Legendre_Chebyshev_Matrix(M_xi) result(M)
         integer, intent(in) :: M_xi
         real :: M(0:M_xi,0:M_xi) 
         
         integer :: i, j
         real :: L 
         
         M = 0 
         ! First row 
         do j = 0, M_xi, 2 
            M(0,j) = Lambda_Legendre( j* 1d0/2)**2/pi            
         end do 
         
         ! Rest of the rows
         do j = 1, M_xi
            do i = 1, j
               if( mod(i+j,2) == 0 ) & 
                 M(i,j) = 2* Lambda_Legendre( (j-i)* 1d0/2 ) &
                           * Lambda_Legendre( (j+i)* 1d0/2 )/pi 
            end do
         end do 
         
      end function 
      
      real elemental function Lambda_Legendre(x) result(L)
         real, intent(in) :: x
         
         L = gamma(x+0.5)/gamma(x+1)
      
      end function
      
      ! *** Computes the integrals 2*\int_{x0}^{t} P_{k} \dd{\xi} for k<= M_xi
      ! where t >= 0. 
      ! For  odd k, x0 = -1   
      ! For even k, x0 = 0 (is equal to the integral starting at x0=-1) 
      function Integral_I(t) result(I_k)
        real, intent(in) :: t 
        real :: I_k(0:M_xi)
        
        real :: P(0:M_xi+1) ! Legendre polynomials
        integer :: k
        
        ! Legendre polynomials 
        P(0) = 1 ; P(1) = t
        do k = 2, M_xi+1
           P(k) = Legendre_3term( k, t, P(k-1), P(k-2))
        end do 
        !P(0) = 1 ; P(1) = t
        !P(2:M_xi) = [( Legendre_Chebyshev( k, t ), k=2,M_xi )]
        
        ! Integral 2*\int_{0}^{t} P_{k} \dd{\xi} (up to a constant for odd k) 
        I_k(0) = 2*t  
        I_k(1:M_xi) = [( 2*(P(k+1) - P(k-1))/(2*k+1),  k=1,M_xi)]        
                
      end function 
      
      ! Computes the integrals:
      !  2*\int_{0}^{t} P_{k} \dd{\xi} = I(k)
      !  2*\int_{0}^{t} \xi P_{2k+1} \dd{\xi} = J(2k+1)
      !  2*\int_{0}^{t} \xi^2 P_{2k} \dd{\xi} = J(2k)
      ! for k=0,1,... M_xi.
      subroutine Compute_Legendre_Moments( t, I, J )
         real, intent(in) :: t
         real, intent(out) :: I(0:M_xi), J(0:M_xi)
         
         real :: II(0:M_xi+1), JJ(0:M_xi+1)
         integer :: k
         
         ! Integral 2*\int_{0}^{t} P_{k} \dd{\xi} (up to a constant for odd k) 
         I = Integral_I(t)
         
         ! Auxiliar integrals for simpler recurrences 
         II(0:M_xi) = I ; II(M_xi+1) = 0 
         
         ! Integral 2*\int_{0}^{t} \xi P_{2k+1} \dd{\xi}
         do k = 1, M_xi, 2
            J(k) = ( (k+1)* II(k+1) + k * II(k-1) ) / (2*k + 1)          
         end do 
         
         ! Auxiliar integrals for simpler recurrences 
         JJ(0:M_xi) = J ; JJ(M_xi+1) = 0
         ! Integral 2*\int_{0}^{t} \xi^2 P_{2k} \dd{\xi}
         J(0) = J(1)
         do k = 2, M_xi, 2
            J(k) = ( (k+1)*JJ(k+1) + k*JJ(k-1) )/(2*k+1)
         end do 
      
      end subroutine 
      
      ! Computes P_{k}(x) given k, x and P1 = P_{k-1}(x) and P2 = P_{k-2}(x)
      real function Legendre_3term( k, x, P1, P2) result(P) 
        integer, intent(in) :: k
        real, intent(in) :: x, P1, P2
        
        if( k==0 .or. x==1 ) then
          P = 1
        elseif (k==1) then
          P = x
        else
          P = ( (2*k-1)*x*P1 - (k-1)*P2 )/k
        end if
      end function   
      
      real function Legendre_Chebyshev( k, x ) result( P )
        integer, intent(in) :: k
        real, intent(in) :: x
        
        real :: T(0:M_xi), beta(0:M_xi) 
        integer :: kk
        
        ! Chebyshev polynomials 
        T = [( cos( kk*acos(x) ) , kk = 0, M_xi)]
        
        ! Coefficients of P_k(x) in Chebyshev basis
        beta = M_Legendre(:,k) 
        
        ! Legendre polynomial in Chebyshev basis
        P = dot_product( beta, T )
         
      end function                     
      
   end subroutine
   
   ! *** Initializes the grids and coefficients/magnetic configuration
   ! related quantities required for solving the DKE and posterior computation
   ! of monoenergetic coefficients. 
   subroutine Set_DKE_Grids_and_coefficients( nu, E_psi, &
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
         call Grid_initialization( "zeta ", zeta, 1 )

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
               
        V_prime = Integral( ["theta","zeta "], g ) 
        B2_mean = Integral( ["theta","zeta "], g * B**2 ) / V_prime 
        
        ! *** Coefficients for DKE 
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
           a(:,:,k,0)  =  0.5*nu * k * (k+1)              
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
