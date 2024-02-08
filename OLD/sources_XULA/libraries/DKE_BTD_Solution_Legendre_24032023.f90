module DKE_BTD_Solution_Legendre 
   
   
   use Magnetic_configuration
   
   !use Lagrange_Interpolation
   use Finite_Differences
   use Barycentric_grids
   
   use Sparse_Matrices
   
   !use Linear_Algebra_PETSC
   use Linear_systems
   use plots   
   
   implicit none

   private

   public :: Compute_Monoenergetic_Database
   public :: DKE_zeta_Convergence
   public :: Solve_BTD_DKE_Legendre, Solve_BTD_DKE_Legendre_Reusing
   
   ! *** Collisionality and radial electric field
   real, save :: nu, E_psi!, eps
      
   ! Phase space grids
   integer, save :: N_theta, N_zeta, N_xi, N  
   logical, save :: Grids_and_Coefficients_DKE = .false. 
   integer, save :: q_theta = 5, q_zeta = 5 
   integer, save :: q1_theta , &!= ( q_theta - mod(q_theta,2) )/2, &
                    q1_zeta  , &!= ( q_zeta - mod(q_zeta,2) )/2, &
                    q2_theta , &!= ( q_theta + mod(q_theta,2) )/2, &
                    q2_zeta   != ( q_zeta + mod(q_zeta,2) )/2
                     
   real, save, allocatable :: theta(:), zeta(:)!, xi(:) 
   real, save, allocatable :: theta_P(:), zeta_P(:) ! auxiliar for FD with periodicity
   
   ! Discretization matrices for Fourier Pseudospectral
   real, save, allocatable :: D_theta(:,:), D_zeta(:,:) 
   
   ! Block matrices for Block TriDiagonal (BTD) solution
   real, save, allocatable :: A_Block(:,:,:), B_Block(:,:,:), C_Block(:,:,:)
   real, save, allocatable :: Beta_Block(:,:,:), Beta_Inv_Block(:,:,:)
   
   real, parameter :: pi = acos(-1d0)
      
   ! Source term and coefficients of the drift-kinetic equation
   ! DKE in Legendre basis: 
   !   a_theta(l-1) * u_theta(l-1) + a_zeta(l-1) * u_zeta(l-1) + a(l-1) + u(l-1)
   ! + a_theta(l) * u_theta(l) + a_zeta(l) * u_zeta(l) + a(l) + u(l)
   ! + a_theta(l+1) * u_theta(l+1) + a_zeta(l+1) * u_zeta(l+1)   + a(l+1) + u(l+1)
   ! = s(l).
   real, save, allocatable :: a_theta(:,:,:,:), a_zeta(:,:,:,:), a(:,:,:,:) 
   real, save, allocatable ::  vm(:,:,:)
   real, save, allocatable :: B(:,:), g(:,:), V_prime 

   
 contains
  
  
  subroutine DKE_zeta_Convergence( nu, E_rho, N_theta, N_zeta, N_xi, Gamma_11, Gamma_31, t_clock ) 
    real, intent(in) :: nu,  E_rho
    integer, intent(in) :: N_theta, N_zeta(:) , N_xi
    real, intent(out)  :: Gamma_11(size(N_zeta)), Gamma_31(size(N_zeta))
    real, intent(out)  :: t_clock(size(N_zeta))
    
    real :: t0, t1, rate, residual(0:N_xi)
    integer :: c0, c1, c_rate, M_zeta
    integer :: j
    
    ! Number of N_zeta's used
    M_zeta = size(N_zeta) 
    
    ! Initialize compiler internal clock 
    call system_clock(count_rate=c_rate) ; rate = real(c_rate)
    
    do j = 1, M_zeta       
       call system_clock(c0)
       call Solve_BTD_DKE_Legendre( N_theta, N_zeta(j), N_xi, &
                          nu, E_rho, Gamma_11(j), Gamma_31(j), residual )
       call system_clock(c1)
                  
       t0 = c0 / rate ; t1 = c1 / rate ! Instants of time in seconds
       t_clock(j) = t1 - t0            ! Elapsed time in seconds
       write(*,*) "Gamma_11 and Gamma_31 calculation elapsed time = ", t_clock(j)       
    end do
    
  end subroutine 
  

  
  subroutine Compute_Monoenergetic_Database( N_theta, N_zeta, N_xi, nu, E_rho, Gamma_11, Gamma_31 )
      integer, intent(in) :: N_theta, N_zeta, N_xi
      real, intent(in) :: nu(:), E_rho(:)
      real, intent(out) :: Gamma_11(size(nu),size(E_rho)), & 
                           Gamma_31(size(nu),size(E_rho))

      integer :: i, j, N_nu, N_E_rho
      real :: t_clock(size(nu),size(E_rho)), residual(size(nu),size(E_rho),0:N_xi)

      N_nu = size(nu) ; N_E_rho = size(E_rho) 
      do j = 1, N_E_rho
         do i = 1, N_nu
            
            call Solve_BTD_DKE_Legendre( N_theta, N_zeta, N_xi, nu(i), E_rho(j), &
                       Gamma_11(i,j), Gamma_31(i,j), residual(i,j,:) )            

         end do
      end do      
      
   end subroutine

   ! Solves the monoenergetic DKE using an algorithm for inverting
   ! block tridiagonal matrices. N_xi must be greater or equal to 2
   subroutine Solve_BTD_DKE_Legendre_Reusing( N_theta, N_zeta, N_xi, nu, E_rho, &
                                       Gamma_11, Gamma_31, Gamma_13, Gamma_33 )
      integer, intent(in) :: N_theta, N_zeta, N_xi
      real, intent(in) :: nu, E_rho
      real, intent(out) :: Gamma_11, Gamma_31, Gamma_13, Gamma_33
      
      ! Block matrices for discretization and forward elimination
      real :: L(N_theta*N_zeta,N_theta*N_zeta,1:2)   ! Lower  
      real :: D(N_theta*N_zeta,N_theta*N_zeta,0:2)   ! Diagonal  
      real :: U(N_theta*N_zeta,N_theta*N_zeta,0:2)   ! Upper
      
      real :: Lower(N_theta*N_zeta,N_theta*N_zeta,1:2)   ! Lower for Backward substitution
      real :: Delta(N_theta*N_zeta,N_theta*N_zeta,0:2)   ! Diagonal 
      
      ! Source term for backward substitution
      real, target :: Sigma_1(0:N_theta-1,0:N_zeta-1, 0:2)
      real, pointer :: pSigma_1(:,:) 
      real, target :: Sigma_3(0:N_theta-1,0:N_zeta-1, 0:2)
      real, pointer :: pSigma_3(:,:) 

      ! Distribution function 
      real, target :: F1(0:N_theta-1,0:N_zeta-1, 0:2), F3(0:N_theta-1,0:N_zeta-1, 0:2) 
      real, pointer :: pF1(:,:),  pF3(:,:) 

      ! Source term for original BTD system
      real, target :: S1(0:N_theta-1,0:N_zeta-1, 0:N_xi), S3(0:N_theta-1,0:N_zeta-1, 0:N_xi) 
      real, pointer :: pS1(:,:), pS3(:,:)
      
      integer :: N_fs 
      
      N_fs = N_theta * N_zeta
      N = N_fs * (N_xi+1) 

      ! *** Fixing coefficients and quantities required for the calculation
      call Set_DKE_Grids_and_coefficients( nu, E_rho/psi_p, &
                                         N_theta, N_zeta, N_xi, & 
                                         q_theta, q_zeta ) 

      ! *** Bounds remapping and source term
      pF1(1:N_fs,0:2) => F1
      pF3(1:N_fs,0:2) => F3
      pS1(1:N_fs,0:N_xi) => S1
      pS3(1:N_fs,0:N_xi) => S3
      pSigma_1(1:N_fs,0:2) => Sigma_1
      pSigma_3(1:N_fs,0:2) => Sigma_3
      
      S1 = -vm ; S1(0,0,0) = 0 ; Sigma_1 = S1(:,:,0:2)
      S3 = 0 ; S3(:,:,1) = - B(0:N_theta-1,0:N_zeta-1) ; Sigma_3 = S3(:,:,0:2)

      ! Solve the drift-kinetic equation with the BTD algorithm
      call Forward_elimination_reusing
      call Backward_substitution_reusing

      ! Once solved compute the monoenergetic coefficients      
      call Compute_monoenergetic_coefficients( Gamma_11, Gamma_31, Gamma_13, Gamma_33 )

    contains

      ! *** Subroutine to calculate the monoenergetic coefficients once
      ! the solution W(0:N_theta-1,0:N_zeta-1,0:N_xi) is known.
      subroutine Compute_monoenergetic_coefficients( Gamma_11, Gamma_31, Gamma_13, Gamma_33 )
        real, intent(out) :: Gamma_11, Gamma_31, Gamma_13, Gamma_33 
     
        real :: U1(0:N_theta,0:N_zeta,0:2), U3(0:N_theta,0:N_zeta,0:2)        
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
        Gamma_11 = Integral( ["theta_int","zeta_int"], g * ( -2 * vm(:,:,0) * U1(:,:,0) - 2* vm(:,:,2) * U1(:,:,2) / 5 )  ) / V_prime
        Gamma_31 = Integral( ["theta_int","zeta_int"], g * 2* U1(:,:,1) * B/3) / V_prime 
        Gamma_13 = Integral( ["theta_int","zeta_int"], g * ( -2 * vm(:,:,0) * U3(:,:,0) - 2* vm(:,:,2) * U3(:,:,2) / 5 )  ) / V_prime
        Gamma_33 = Integral( ["theta_int","zeta_int"], g * 2 * B * U3(:,:,1) /3) / V_prime  
     
        ! ** Scaling to match DKES normalization
        Gamma_11 =  Gamma_11 / psi_p**2 
        Gamma_31 =  Gamma_31 / (psi_p * B00) ! psi_p from vm^\rho = vm^\psi \psi'
        Gamma_13 =  Gamma_13 / (psi_p * B00)
        Gamma_33 = -Gamma_33 / B00          ! the source is -B * xi

   
      end subroutine

      ! Performs the forward elimination storing only the matrices
      ! associated at the modes k=0,1 and 2 to save memory. 
      subroutine Forward_elimination_reusing
        integer :: k
        
        real :: t0, t1 
        integer :: c0, c1, c_rate
        real :: Delta_k( N_theta*N_zeta, N_theta*N_zeta ) ! Auxiliary variable for calculations 
        
        call system_clock(count_rate=c_rate) ! Get click rate
        write(*,*) " Legendre Mode k = ", N_xi
        !L(:,:,2) = Block_L(N_xi)
        !Delta_k = Block_D(N_xi)  ! Delta_{N_\xi}
        
        L(:,:,2) = Block_L_NEW(N_xi)
        Delta_k  = Block_D_NEW(N_xi)  ! Delta_{N_\xi}
         
        D(:,:,2) = Invert_LU_LAPACK_NEW( Delta_k ) ! Delta_{N_\xi}^{-1} ! TO DO: PARALLELIZE THIS
        
        do k = N_xi-1, 0, -1
           
           write(*,*) " Legendre Mode k = ", k                        
           call system_clock(c0)
           ! Legendre mode of the discretization matrix
           !L(:,:,1) =  Block_L(k) ! L_k 
           !D(:,:,1) =  Block_D(k) ! D_k 
           !U(:,:,1) =  Block_U(k) ! U_k              
           
           L(:,:,1) =  Block_L_NEW(k) ! L_k 
           D(:,:,1) =  Block_D_NEW(k) ! D_k  
           U(:,:,1) =  Block_U_NEW(k) ! U_k            
           call system_clock(c1)
           write(*,*) "      Time in Defining L, D, U = ", ( c1 - c0  ) / real(c_rate)

           ! Schur complement Delta_k = D_k - U_k * Delta_{k+1}^{-1} L_{k+1}
           ! and its inverse
           call system_clock(c0)
           Delta_k = D(:,:,1) &
                   - matmul_LAPACK( U(:,:,1), matmul_LAPACK(D(:,:,2), L(:,:,2)) )   
           call system_clock(c1)
           write(*,*) "      Time in Defining Delta = ", ( c1 - c0  ) / real(c_rate)
           
           call system_clock(c0)
           D(:,:,1) = Invert_LU_LAPACK_NEW( Delta_k ) ! Delta_{k}^{-1} ! TO DO: PARALLELIZE THIS
           call system_clock(c1)
           write(*,*) "      Time in Inverting Delta = ", ( c1 - c0  ) / real(c_rate)

           ! Source terms for equivalent system 
           if( k <= 1 ) &  ! The source term S1 only has Legendre modes 0 and 2
               pSigma_1(:,k) = pS1(:,k) - matmul( U(:,:,1), matmul( Delta(:,:,k+1), pSigma_1(:,k+1) ) )
           
           if( k <=1 ) &  ! The source term S3 only has Legendre mode 1
               pSigma_3(:,k) = pS3(:,k) - matmul( U(:,:,1), matmul( Delta(:,:,k+1), pSigma_3(:,k+1) ) )

           ! Necessary matrices for next iterations
           call system_clock(c0)
           L(:,:,2) = L(:,:,1) 
           D(:,:,2) = D(:,:,1)
           call system_clock(c1)
           write(*,*) "      Time in Defining Next iteration = ", ( c1 - c0  ) / real(c_rate)
           
           ! *** Extracting matrices for modes 0, 1 and 2 for Backward substitution
           if( k <= 2 )             Delta(:,:,k) = D(:,:,1)
           if( 1 <= k .and. k <=2 ) Lower(:,:,k) = L(:,:,1)               

        end do

      end subroutine
      
      
      subroutine Backward_substitution_reusing       
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
         real :: L(N_theta*N_zeta, N_theta*N_zeta)
         
         integer :: i, j, i_row 
         type( Sparse_vector ) :: Lk_row 
         L = 0
         do j = 0, N_zeta-1
            do i = 0, N_theta-1
               
               ! Row as a sparse vector type
               Lk_row = DKE_row_mode_cols_Fourier(i,j,k, -1)
               ! Eliminating the k dependence in the index
               Lk_row % i =  Lk_row % i - (k-1) * N_theta * N_zeta
               
               i_row = 1 + i + j * N_theta
               L(i_row, Lk_row % i) = Lk_row % v
               
               deallocate( Lk_row % i, Lk_row % v )
               
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
        real :: D(N_theta*N_zeta, N_theta*N_zeta)
        
        integer :: i, j, i_row
        type( Sparse_vector ) :: Dk_row
        
        D = 0
        do j = 0, N_zeta-1
           do i = 0, N_theta-1
              
              ! Row as a sparse vector type
              Dk_row = DKE_row_mode_cols_Fourier(i,j,k, 0)
              
              ! Eliminating the k dependence in the index
              Dk_row % i =  Dk_row % i - k * N_theta * N_zeta
              i_row = 1 + i + j * N_theta
              D(i_row, Dk_row % i) = Dk_row % v
              
              deallocate( Dk_row % i, Dk_row % v )
              
           end do
        end do
        
        ! Elimination of nullspace
        if( k==0 ) D(1,1) = 1 ; if( k==0 ) D(1,2:N_theta*N_zeta) = 0
      end function
       
       
      ! *** Computes the block matrix L_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function.
      function Block_L_NEW(k) result(L)
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
      function Block_D_NEW(k) result(D)
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
        if( k==0 ) D(1,1) = 1 ; D(1,2:N_fs) = 0
      end function
       
       
       
       
      ! *** Computes the block matrix U_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1)
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function.
      function Block_U_NEW(k) result(U)
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
       
       
       
      ! *** Computes the k-th block matrix U_k(N_theta*N_zeta, N_theta*N_zeta)
      ! from the discretization
      ! L_k * f(:,:,k-1) + D_k * f(:,:,k) + U_k * f(:,:,k+1) 
      ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
      ! distribution function.
      function Block_U(k) result(U)
        integer, intent(in) :: k
        real :: U(N_theta*N_zeta, N_theta*N_zeta)
        
        integer :: i, j, i_row 
        type( Sparse_vector ) :: Uk_row 
        U = 0
        do j = 0, N_zeta-1
           do i = 0, N_theta-1
              
              ! Row as a sparse vector type
  	      Uk_row = DKE_row_mode_cols_Fourier(i,j,k, +1)     
  	      ! Eliminating the k dependence in the index
  	      Uk_row % i =  Uk_row % i - (k+1) * N_theta * N_zeta 
  	         
              i_row = 1 + i + j * N_theta
  	      U(i_row, Uk_row % i) = Uk_row % v
  	         
  	      deallocate( Uk_row % i, Uk_row % v )
  	      
           end do 
        end do 
        
        ! Elimination of nullspace 
        if( k == 0 ) U(1,:) = 0 
      end function 


      function DKE_row_mode_cols_Fourier(i,j,k, kk) result(F)
            integer, intent(in) :: i, j , k, kk
            type(Sparse_vector) :: F
  
            type(Sparse_vector) :: F_theta, F_zeta 
            real ::  F_ijk(1), F_theta_v(N_theta-1), F_zeta_v(N_zeta-1) 
            integer :: F_theta_i(N_theta-1), F_zeta_i(N_zeta-1) 
            integer :: ii, jj, N_row, i_row, s_theta(0:N_theta-1), s_zeta(0:N_zeta-1)
            
            ! *** Initialize sparse vectors for the summands of the DKE
            call Initialize_sparse_vector( F_theta, N_theta )
            call Initialize_sparse_vector(  F_zeta, N_zeta )
  
            ! *** Derivative along theta: a_theta * u_theta
            s_theta = [( ii, ii=0,N_theta-1 )] 
            F_theta % i = 1 + s_theta + j * N_theta + (k+kk)*N_theta*N_zeta 
            F_theta % v = a_theta(i,j,k,kk) * D_theta(i,:) 
  
            ! *** Derivative along zeta: a_zeta * u_zeta
            s_zeta = [( jj, jj=0,N_zeta-1 )] 
            F_zeta % i = 1 + i + s_zeta * N_theta + (k+kk) * N_theta*N_zeta 
            F_zeta % v = a_zeta(i,j,k,kk) * D_zeta(j,:) 
  
            ! *** Complete discretization row
            i_row = 1 + i + j * N_theta + (k+kk) * N_theta * N_zeta
            N_row = N_theta + N_zeta - 1           
            allocate( F % i( N_row ),  F % v( N_row ) )
            
            ! Diagonal term of the DKE
            F % i(1) = i_row 
            F_ijk = pack( F_theta % v, F_theta % i == i_row )   &
                  + pack(  F_zeta % v,  F_zeta % i == i_row ) 
            F % v(1) = F_ijk(1)  + a(i,j,k,kk)
            
            ! Non diagonal terms of the DKE
            F_theta_i = pack( F_theta % i, F_theta % i /= i_row ) 
            F_theta_v = pack( F_theta % v, F_theta % i /= i_row ) 
            
            F_zeta_i = pack( F_zeta % i, F_zeta % i /= i_row ) 
            F_zeta_v = pack( F_zeta % v, F_zeta % i /= i_row ) 
            
            F % i(2:N_row) = [ F_theta_i, F_zeta_i ] 
            F % v(2:N_row) = [ F_theta_v, F_zeta_v ] 
            
            deallocate( F_theta % i, F_theta % v, F_zeta % i,  F_zeta % v )
       end function

      
     end subroutine
     
    



   
   ! Solves the monoenergetic DKE using an algorithm for inverting
   ! block tridiagonal matrices 
   subroutine Solve_BTD_DKE_Legendre( N_theta, N_zeta, N_xi, nu, E_rho, &
                                       Gamma_11, Gamma_31, residual )
      integer, intent(in) :: N_theta, N_zeta, N_xi
      real, intent(in) :: nu, E_rho
      real, intent(out) :: Gamma_11, Gamma_31
      real, intent(out) :: residual(0:N_xi) 
      
      ! Block matrices
      real :: L(N_theta*N_zeta,N_theta*N_zeta,1:N_xi)   ! Lower  
      real :: D(N_theta*N_zeta,N_theta*N_zeta,0:N_xi)   ! Diagonal  
      real :: U(N_theta*N_zeta,N_theta*N_zeta,0:N_xi-1) ! Upper  

      ! Matrices from Forward elimination
      real :: Delta(N_theta*N_zeta,N_theta*N_zeta,0:N_xi)
      real :: Delta_Inv(N_theta*N_zeta,N_theta*N_zeta,0:N_xi)
      
      ! Source term for backward substitution
      real, target :: Sigma(0:N_theta-1,0:N_zeta-1, 0:N_xi)
      real, pointer :: pSigma(:,:) 

      ! Distribution function 
      real, target :: F(0:N_theta-1,0:N_zeta-1, 0:N_xi) 
      real, pointer :: pF(:,:) 

      ! Source term for original BTD system
      real, target :: S(0:N_theta-1,0:N_zeta-1, 0:N_xi) 
      real, pointer :: pS(:,:)

      real, target :: Res(N_theta*N_zeta, 0:N_xi)
      real, pointer :: pRes(:) 
      
      integer :: N

      N = N_theta * N_zeta * (N_xi+1) 
      
      ! Construction of BTD discretization matrix for equation
      ! L(k) F(k-1) + D(k) F(k) + U(k) F(k+1) = S(k) 
      call Set_Discretization_Matrix( N_theta, N_zeta, N_xi, nu, E_rho/psi_p, &
                                      L, D, U )

      ! Bounds remapping and source term
      pF(1:N_theta*N_zeta,0:N_xi) => F
      pS(1:N_theta*N_zeta,0:N_xi) => S
      pSigma(1:N_theta*N_zeta,0:N_xi) => Sigma
      S = -vm ; S(0,0,0) = 0 ; pSigma = pS

      call Forward_elimination
      call Backward_substitution

      ! Radial transport and bootstrap current monoenergetic coefficients
      call Compute_monoenergetic_coefficients( Gamma_11, Gamma_31 )      
      
      ! Residual computation 
      call Compute_residual 
    contains
      subroutine Compute_residual
         integer :: k
         
         residual(0) = &
            norm2(  matmul(D(:,:,0), pF(:,0))    &
                 +  matmul(U(:,:,0), pF(:,1))  - pS(:,0) )
         do k = 1, N_xi-1
            residual(k) = &
            norm2(  matmul(L(:,:,k), pF(:,k-1))  &
                 +  matmul(D(:,:,k), pF(:,k))    &
                 +  matmul(U(:,:,k), pF(:,k+1))  - pS(:,k) )         
         
         end do  
         residual(N_xi) = &
         norm2(  matmul(L(:,:,N_xi), pF(:,N_xi-1))  &
              +  matmul(D(:,:,N_xi), pF(:,N_xi))- pS(:,N_xi) )  
      
      
      end subroutine 

      subroutine Forward_elimination        
        integer :: k
        
        real :: t0, t1 
        integer :: c0, c1, c_rate
    
        call system_clock(count_rate=c_rate) ! Get click rate
        
        write(*,*) " Legendre Mode k = ", N_xi
        Delta(:,:,N_xi) = D(:,:,N_xi) 
        Delta_inv(:,:,N_xi) =  Invert_LU_LAPACK_NEW( Delta(:,:,N_xi) )
        !Delta_inv(:,:,N_xi) =  Invert_LU_LAPACK( Delta(:,:,N_xi) )
        !TODO: Delta_inv(:,:,N_xi) =  Invert_LU_PETSC( Delta(:,:,N_xi) )
        do k = N_xi-1, 0, -1
           write(*,*) " Legendre Mode k = ", k
           call system_clock(c0)
           Delta(:,:,k) = D(:,:,k) - matmul( U(:,:,k), matmul(Delta_inv(:,:,k+1), L(:,:,k+1)) )
           call system_clock(c1)
           write(*,*) "      Time in Defining Delta = ", ( c1 - c0  ) / real(c_rate)
           
           call system_clock(c0)
           Delta_inv(:,:,k) =  Invert_LU_LAPACK_NEW( Delta(:,:,k) )
           call system_clock(c1) 
           write(*,*) "      Time in Inverting Delta = ", ( c1 - c0  ) / real(c_rate)
           !Delta_inv(:,:,k) =  Invert_LU_LAPACK( Delta(:,:,k) )
          !TODO: Delta_inv(:,:,k) =  Invert_LU_PETSC( Delta(:,:,k) )
   
           if( k <= 1 ) & ! Avoid multiplications with zeros
           pSigma(:,k) = pS(:,k) - matmul( U(:,:,k), matmul(Delta_inv(:,:,k+1), pSigma(:,k+1)) )
        end do
        
      end subroutine

      subroutine Backward_substitution        
        integer :: k
        
        pF(:,0) = matmul( Delta_inv(:,:,0), pSigma(:,0) )
        do k = 1, N_xi
           pF(:,k) = matmul( Delta_inv(:,:,k), pSigma(:,k)- matmul( L(:,:,k), pF(:,k-1) ) ) 
        end do
        
      end subroutine

      ! *** Subroutine to calculate the monoenergetic coefficients once
      ! the solution W(0:N_theta-1,0:N_zeta-1,0:N_xi) is known.
      subroutine Compute_monoenergetic_coefficients( Gamma_11, Gamma_31 )
        real, intent(out) :: Gamma_11, Gamma_31 
     
        real :: U(0:N_theta,0:N_zeta,0:2)        
        real :: I_11(0:N_theta,0:N_zeta), I_31(0:N_theta,0:N_zeta)
        integer :: i, j 
     
        ! *** Legendre modes 1 and 2 for radial transport and bootstrap
        ! in periodic grid with repeated point.    
        do j = 0, N_zeta
           do i = 0, N_theta
              U(i,j,:) = F( modulo(i,N_theta), modulo(j,N_zeta), 0:2)
           end do
        end do
     
        ! *** Flux surface average        
        Gamma_11 = Integral( ["theta_int","zeta_int"], g * ( -2 * vm(:,:,0) * U(:,:,0) - 2* vm(:,:,2) * U(:,:,2) / 5 )  ) / V_prime
        Gamma_31 = Integral( ["theta_int","zeta_int"], g * 2* U(:,:,1) * B/3) / V_prime 
     
        ! ** Scaling to match DKES normalization
        Gamma_11 = Gamma_11 / psi_p**2 
        Gamma_31 = Gamma_31 / (psi_p * B00)

        !call Beidler_NF_2011_Normalization
   
      end subroutine
      
   end subroutine
   
   
   
  ! *** Initializes the discretization matrix for solving the DKE and also
  ! the array that relates each ( theta(i), zeta(j), xi(k) ) with its position
  ! in the state vector row_index(i,j,k) = 1 + i + j * N_theta + k * N_theta * N_zeta
   subroutine Set_Discretization_Matrix( N_theta, N_zeta, N_xi, nu, E_psi, &
                                         L, D, U )
       integer, intent(in) :: N_theta, N_zeta, N_xi
       real, intent(in) :: nu, E_psi 
       real, intent(out) :: L(N_theta*N_zeta, N_theta*N_zeta, 1:N_xi)
       real, intent(out) :: D(N_theta*N_zeta, N_theta*N_zeta, 0:N_xi)
       real, intent(out) :: U(N_theta*N_zeta, N_theta*N_zeta, 0:N_xi-1)
       
       real :: t0, t1

       call Set_DKE_Grids_and_coefficients( nu, E_psi, &
                                         N_theta, N_zeta, N_xi, & 
                                         q_theta, q_zeta ) 
       
       !call CPU_time(t0)
       call Compute_Block_L(L)
       call Compute_Block_D(D)
       call Compute_Block_U(U)     
       
       ! Elimination of nullspace      
       D(1,1,0) = 1 ; D(1,2:N_theta*N_zeta,0) = 0
       U(1,:,0) = 0       
       
       !call System_clock(t1)       
       !write(*,*) " Elapsed time in constructing BTD discretization matrix ", t1 - t0       
     
       contains
       
       ! *** Computes the block matrix L(N_theta*N_zeta, N_theta*N_zeta, 1:N_xi)
       ! from the discretization
       ! L(:,:,k) * f(:,:,k-1) + D(:,:,k) * f(:,:,k) + U(:,:,k) * f(:,:,k+1) 
       ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
       ! distribution function.
       subroutine Compute_Block_L(L)
         real, intent(out) :: L(N_theta*N_zeta, N_theta*N_zeta, 1:N_xi)
         
         integer :: i, j, k, i_row 
         type( Sparse_vector ) :: Lk_row 
         L = 0
         do k = 1, N_xi
            do j = 0, N_zeta-1
               do i = 0, N_theta-1
                  
                  ! Row as a sparse vector type
  		  Lk_row = DKE_row_mode_cols_Fourier(i,j,k, -1)     
  		  ! Eliminating the k dependence in the index
  		  Lk_row % i =  Lk_row % i - (k-1) * N_theta * N_zeta 
  		        
                  i_row = 1 + i + j * N_theta
  		  L(i_row, Lk_row % i, k) = Lk_row % v
  		        
  		  deallocate( Lk_row % i, Lk_row % v )
  		     
               end do 
            end do 
         end do 
         
       end subroutine 
       
       ! *** Computes the block matrix D(N_theta*N_zeta, N_theta*N_zeta, 0:N_xi)
       ! from the discretization
       ! L(:,:,k) * f(:,:,k-1) + D(:,:,k) * f(:,:,k) + U(:,:,k) * f(:,:,k+1) 
       ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
       ! distribution function.
       subroutine Compute_Block_D(D)
         real, intent(out) :: D(N_theta*N_zeta, N_theta*N_zeta, 0:N_xi)
         
         integer :: i, j, k, i_row 
         type( Sparse_vector ) :: Dk_row 
         D = 0
         do k = 0, N_xi
            do j = 0, N_zeta-1
               do i = 0, N_theta-1
                  
                  ! Row as a sparse vector type
  		  Dk_row = DKE_row_mode_cols_Fourier(i,j,k, 0)     
  		  ! Eliminating the k dependence in the index
  		  Dk_row % i =  Dk_row % i - k * N_theta * N_zeta 
  		          
                  i_row = 1 + i + j * N_theta
  		  D(i_row, Dk_row % i, k) = Dk_row % v
  		        
  		  deallocate( Dk_row % i, Dk_row % v )
  		     
               end do 
            end do 
         end do 
       end subroutine 
       
       ! *** Computes the block matrix U(N_theta*N_zeta, N_theta*N_zeta, 0:N_xi-1)
       ! from the discretization
       ! L(:,:,k) * f(:,:,k-1) + D(:,:,k) * f(:,:,k) + U(:,:,k) * f(:,:,k+1) 
       ! Where f(0:N_theta-1,0:N_zeta-1,k) is the k-th Legendre mode of the
       ! distribution function.
       subroutine Compute_Block_U(U)
         real, intent(out) :: U(N_theta*N_zeta, N_theta*N_zeta, 0:N_xi-1)
         
         integer :: i, j, k, i_row 
         type( Sparse_vector ) :: Uk_row 
         U = 0
         do k = 0, N_xi-1
            do j = 0, N_zeta-1
               do i = 0, N_theta-1
               
                  ! Row as a sparse vector type
  		  Uk_row = DKE_row_mode_cols_Fourier(i,j,k, +1)     
  		  ! Eliminating the k dependence in the index
  		  Uk_row % i =  Uk_row % i - (k+1) * N_theta * N_zeta 
  		        
                  i_row = 1 + i + j * N_theta
  		  U(i_row, Uk_row % i, k) = Uk_row % v
  		        
  		  deallocate( Uk_row % i, Uk_row % v )
  		     
               end do 
            end do 
         end do 
         
       end subroutine 
   
   
       function DKE_row_mode_cols_FD(i,j,k, kk) result(F)
            integer, intent(in) :: i, j , k, kk
            type(Sparse_vector) :: F
  
            type(Sparse_vector) :: F_theta, F_zeta 
            real ::  F_ijk(1), F_theta_v(q_theta), F_zeta_v(q_zeta) 
            integer :: F_theta_i(q_theta), F_zeta_i(q_zeta) 
            integer :: ii, jj, N_row, i_row, s_theta(0:q_theta), s_zeta(0:q_theta)
            
            ! *** Initialize sparse vectors for the summands of the DKE
            call Initialize_sparse_vector( F_theta, q_theta+1 )
            call Initialize_sparse_vector(  F_zeta,  q_zeta+1 )
  
            ! *** Derivative along theta: a_theta * u_theta
            ! Stencil in [-q_theta, N_theta + q_theta] 
            s_theta = [( ii + i - q1_theta, ii=0,q_theta )] 
            F_theta % i = 1 + modulo(s_theta, N_theta) + j * N_theta + (k+kk)*N_theta*N_zeta 
            F_theta % v = a_theta(i,j,k,kk) &
                        * [( Derivative( "theta_P", i+q_theta, s_theta(ii)+q_theta, 1 ), ii = 0, q_theta )]
  
            ! *** Derivative along zeta: a_zeta * u_zeta
            ! Stencil in [-q_zeta, N_zeta + q_zeta]
            s_zeta = [( jj + j - q1_zeta, jj=0,q_zeta )] 
            F_zeta % i = 1 + i + modulo(s_zeta, N_zeta) * N_theta + (k+kk) * N_theta*N_zeta 
            F_zeta % v = a_zeta(i,j,k,kk) &
                       * [( Derivative( "zeta_P", j+q_zeta, s_zeta(jj) + q_zeta, 1 ), jj=0,q_zeta )]
  
  
            ! *** Complete discretization row
            i_row = 1 + i + j * N_theta + (k+kk) * N_theta * N_zeta
            N_row = q_theta + q_zeta + 1         
            allocate( F % i( N_row ),  F % v( N_row ) )
            
            ! Diagonal term of the DKE
            F % i(1) = i_row 
            F_ijk = pack( F_theta % v, F_theta % i == i_row )   &
                  + pack(  F_zeta % v,  F_zeta % i == i_row ) 
            F % v(1) = F_ijk(1)  + a(i,j,k,kk)
            
            ! Non diagonal terms of the DKE
            F_theta_i = pack( F_theta % i, F_theta % i /= i_row ) 
            F_theta_v = pack( F_theta % v, F_theta % i /= i_row ) 
            
            F_zeta_i = pack( F_zeta % i, F_zeta % i /= i_row ) 
            F_zeta_v = pack( F_zeta % v, F_zeta % i /= i_row ) 
            
            F % i(2:N_row) = [ F_theta_i, F_zeta_i ] 
            F % v(2:N_row) = [ F_theta_v, F_zeta_v ] 
            
            deallocate( F_theta % i, F_theta % v, F_zeta % i,  F_zeta % v )
       end function
       
       function DKE_row_mode_cols_Fourier(i,j,k, kk) result(F)
            integer, intent(in) :: i, j , k, kk
            type(Sparse_vector) :: F
  
            type(Sparse_vector) :: F_theta, F_zeta 
            real ::  F_ijk(1), F_theta_v(N_theta-1), F_zeta_v(N_zeta-1) 
            integer :: F_theta_i(N_theta-1), F_zeta_i(N_zeta-1) 
            integer :: ii, jj, N_row, i_row, s_theta(0:N_theta-1), s_zeta(0:N_zeta-1)
            
            ! *** Initialize sparse vectors for the summands of the DKE
            call Initialize_sparse_vector( F_theta, N_theta )
            call Initialize_sparse_vector(  F_zeta, N_zeta )
  
            ! *** Derivative along theta: a_theta * u_theta
            ! Stencil in [-q_theta, N_theta + q_theta] 
            s_theta = [( ii, ii=0,N_theta-1 )] 
            F_theta % i = 1 + s_theta + j * N_theta + (k+kk)*N_theta*N_zeta 
            F_theta % v = a_theta(i,j,k,kk) * D_theta(i,:) 
  
            ! *** Derivative along zeta: a_zeta * u_zeta
            ! Stencil in [-q_zeta, N_zeta + q_zeta]
            s_zeta = [( jj, jj=0,N_zeta-1 )] 
            F_zeta % i = 1 + i + s_zeta * N_theta + (k+kk) * N_theta*N_zeta 
            F_zeta % v = a_zeta(i,j,k,kk) * D_zeta(j,:) 
  
            ! *** Complete discretization row
            i_row = 1 + i + j * N_theta + (k+kk) * N_theta * N_zeta
            N_row = N_theta + N_zeta - 1           
            allocate( F % i( N_row ),  F % v( N_row ) )
            
            ! Diagonal term of the DKE
            F % i(1) = i_row 
            F_ijk = pack( F_theta % v, F_theta % i == i_row )   &
                  + pack(  F_zeta % v,  F_zeta % i == i_row ) 
            F % v(1) = F_ijk(1)  + a(i,j,k,kk)
            
            ! Non diagonal terms of the DKE
            F_theta_i = pack( F_theta % i, F_theta % i /= i_row ) 
            F_theta_v = pack( F_theta % v, F_theta % i /= i_row ) 
            
            F_zeta_i = pack( F_zeta % i, F_zeta % i /= i_row ) 
            F_zeta_v = pack( F_zeta % v, F_zeta % i /= i_row ) 
            
            F % i(2:N_row) = [ F_theta_i, F_zeta_i ] 
            F % v(2:N_row) = [ F_theta_v, F_zeta_v ] 
            
            deallocate( F_theta % i, F_theta % v, F_zeta % i,  F_zeta % v )
       end function
       
       ! *** Gives the sparse row F (type(Sparse_vector)) associated
       ! to discretize Hirshman (DKES) drift-kinetic equation at the grid 
       ! point ( theta(i), zeta(j), xi(k) ) using finite differences of
       ! orders q_theta, q_zeta and q_xi respectively. It imposes periodicity
       ! through the stencils along theta and zeta. 
       function DKE_row_Fourier(i,j,k) result(F)
            integer, intent(in) :: i, j , k
            type(Sparse_vector) :: F
            
            type(Sparse_vector) :: F_m1, F_0, F_p1 
            
            if ( k == 0 ) then    
              ! *** Columns acting on Legendre mode 0     
              F_0 = DKE_row_mode_cols_Fourier(i,j,k, 0)           
              ! *** Columns acting on Legendre mode 1     
              F_p1 = DKE_row_mode_cols_Fourier(i,j,k, +1)        
              
              ! *** Complete row 
              allocate( F % i( size(F_0%i) + size(F_p1%i) ) )
              allocate( F % v( size(F_0%v) + size(F_p1%v) ) )
              
              F % i = [ F_0 % i, F_p1 % i ]
              F % v = [ F_0 % v, F_p1 % v ]
                          
              
              deallocate( F_0 % i , F_p1 % i)
              deallocate( F_0 % v , F_p1 % v)
              
            elseif ( k == N_xi ) then 
            
              ! *** Columns acting on Legendre mode k-1      
              F_m1 = DKE_row_mode_cols_Fourier(i,j,k, -1)     
              ! *** Columns acting on Legendre mode k     
              F_0 = DKE_row_mode_cols_Fourier(i,j,k, 0)         
              
              ! *** Complete row 
              allocate( F % i( size(F_m1%i) + size(F_0%i) ) )
              allocate( F % v( size(F_m1%v) + size(F_0%v) ) )
              F % i = [ F_m1 % i , F_0 % i ]
              F % v = [ F_m1 % v , F_0 % v ]
              
              deallocate(F_m1 % i , F_0 % i )
              deallocate(F_m1 % v , F_0 % v )      
              
              
            else 
            
              ! *** Columns acting on Legendre mode k-1      
              F_m1 = DKE_row_mode_cols_Fourier(i,j,k, -1)     
              ! *** Columns acting on Legendre mode k     
              F_0 = DKE_row_mode_cols_Fourier(i,j,k, 0)           
              ! *** Columns acting on Legendre mode k*1     
              F_p1 = DKE_row_mode_cols_Fourier(i,j,k, +1)        
              
              ! *** Complete row 
              allocate( F % i( size(F_m1%i) + size(F_0%i) + size(F_p1%i) ) )
              F % i = [ F_m1 % i , F_0 % i , F_p1 % i ]
              allocate( F % v( size(F_m1%v) + size(F_0%v) + size(F_p1%v) ) )
              F % v = [ F_m1 % v , F_0 % v , F_p1 % v ]
                            
              deallocate(F_m1 % i , F_0 % i , F_p1 % i)
              deallocate(F_m1 % v , F_0 % v , F_p1 % v)
            end if
  		   
       end function
       
   end subroutine 
   
   
  ! *** Initializes the grids and coefficients/magnetic configuration
  ! related quantities required for solving the DKE and posterior computation
  ! of monoenergetic coefficients. 
   subroutine Set_DKE_Grids_and_coefficients( nu_IN, E_psi_IN, &
                                         N_theta_IN, N_zeta_IN, N_xi_IN, & 
                                         q_theta_IN, q_zeta_IN ) 
      real, intent(in) :: nu_IN, E_psi_IN 
      integer, intent(in) :: N_theta_IN, N_zeta_IN, N_xi_IN
      integer,  intent(in) :: q_theta_IN, q_zeta_IN
      
      logical :: Grids_DKE = .false., Coefficients_DKE = .false. 
      
      ! *** If the variables have been allocated, deallocate them.
      if( Grids_and_Coefficients_DKE ) then
       deallocate( theta, zeta, theta_P, zeta_P ) 
       deallocate( a_theta, a_zeta, a, B, g, vm ) 
      end if 
      
      ! *** Setting the grid sizes 
      N_theta = N_theta_IN ; N_zeta = N_zeta_IN ; N_xi = N_xi_IN      
      
      ! *** Size of the system of equations
      N = N_theta * N_zeta * (N_xi+1)
      
      ! *** Setting the collisionality and radial electric field
      nu = nu_IN ; E_psi = E_psi_IN 
      
      ! *** Setting the orders from input
      q_theta = q_theta_IN ; q_zeta  = q_zeta_IN 
      q1_theta = ( q_theta - mod(q_theta,2) )/2 ; q2_theta = ( q_theta + mod(q_theta,2) )/2 
      q1_zeta  = ( q_zeta - mod(q_zeta,2) )/2   ;  q2_zeta = ( q_zeta + mod(q_zeta,2) )/2 
      
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
         
         ! *** Poloidal auxiliary extended angular grid for periodicity
         allocate( theta_P(-q_theta:N_theta+q_theta) )
         theta_P(0:N_theta) = theta
         theta_P(-q_theta:-1) = theta(0:q_theta-1) - ( theta(N_theta)  - theta(0) )
         theta_P(N_theta+1:N_theta+q_theta) = theta(1:q_theta) + ( theta(N_theta)  - theta(0) ) 
         
         ! *** Toroidal auxiliary extended angular grid for periodicity
         allocate( zeta_P(-q_zeta:N_zeta+q_zeta) )
         zeta_P(0:N_zeta) = zeta
         zeta_P(-q_zeta:-1) = zeta(0:q_zeta-1) - ( zeta(N_zeta)  - zeta(0) )
         zeta_P(N_zeta+1:N_zeta+q_zeta) = zeta(1:q_zeta) + ( zeta(N_zeta) - zeta(0) ) 
         
         ! *** Finite difference grids for discretization 
         call Grid_initialization( "theta", theta, q_theta )
         call Grid_initialization( "zeta", zeta, q_zeta )
         
         call Grid_initialization( "theta_P", theta_P, q_theta )
         call Grid_initialization( "zeta_P", zeta_P, q_zeta )
         
         ! *** For integrals using trapezoidal rule
         call Grid_initialization( "theta_int", theta, 1 )
         call Grid_initialization( "zeta_int", zeta, 1 )

         ! *** Discretization matrices for Pseudospectral Fourier
         if( allocated(D_theta) )  deallocate( D_theta ) 
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
         
         
         integer :: M1, M2, N1, N2 
         !integer, parameter :: M1 = N_theta - mod(N_theta,2) , M2 = N_theta + mod(N_theta,2) 
         !integer, parameter :: N1 = N_zeta - mod(N_zeta,2) , N2 = N_zeta + mod(N_zeta,2) 
         
         complex, allocatable ::  vds_l_modes(:,:), &
                                  g_inv_modes(:,:), &
                                 gB_inv_modes(:,:), &
                                      B_modes(:,:), &
                                     B2_modes(:,:), &
                                     dBdz_l_modes(:,:)
         
         
         
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
            / ( dldz_l * B**2 )                              ! magnetic radial drift
         dBdz_l = dBdzeta + iota * dBdtheta  ! Angular derivative of B along field lines
                
         V_prime  = Integral( ["theta_int","zeta_int"], g ) 
         B2_mean  = Integral( ["theta_int","zeta_int"], g * B**2) / V_prime      
             
         open(21, file="results/Magnetic_configuration.plt")
         write(21,*) " psi_p, iota, B_theta, B_zeta "
         write(21,*)   psi_p, iota, B_theta, B_zeta  
         write(21,*)  " theta(i), zeta(j), B(i,j), dBdtheta(i,j), dBdzeta(i,j) "
         do j = 0, N_zeta
            do i = 0, N_theta
               write(21,'(9999e)') theta(i), zeta(j), B(i,j), dBdtheta(i,j), dBdzeta(i,j)   
            end do
         end do
         close(21)    !; stop
         
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
         
           a(:,:,k,-1) = (k*(k-1)) * ( dBdzeta + iota * dBdtheta )       &
                       / ( 2*(2*k-1) * g * B**2 ) 
           a(:,:,k,+1) = -((k+1)*(k+2)) * ( dBdzeta + iota * dBdtheta )       &
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
         
         ! *** Fourier Spectra of the coefficients and source term
         M1 = N_theta - mod(N_theta,2) ; M2 = N_theta + mod(N_theta,2) 
         N1 =  N_zeta - mod(N_zeta,2)  ; N2 = N_zeta + mod(N_zeta,2) 
         allocate(  vds_l_modes(-M1/2:M2/2-1,-N1/2:N2/2-1), &
                    g_inv_modes(-M1/2:M2/2-1,-N1/2:N2/2-1), &
                   gB_inv_modes(-M1/2:M2/2-1,-N1/2:N2/2-1), &
                        B_modes(-M1/2:M2/2-1,-N1/2:N2/2-1), &
                       B2_modes(-M1/2:M2/2-1,-N1/2:N2/2-1), &
                   dBdz_l_modes(-M1/2:M2/2-1,-N1/2:N2/2-1) )
                       
         
         vds_l_modes = DFT_2D( theta, zeta, vds_l(0:N_theta-1,0:N_zeta-1) ) 
                       
         g_inv_modes = DFT_2D( theta, zeta, 1d0/g(0:N_theta-1,0:N_zeta-1) )             
         B2_modes = DFT_2D( theta, zeta, B(0:N_theta-1,0:N_zeta-1)**2 ) 
                         
         gB_inv_modes = DFT_2D( theta, zeta, 1d0/(g(0:N_theta-1,0:N_zeta-1)*B(0:N_theta-1,0:N_zeta-1)) )              
         B_modes = DFT_2D( theta, zeta, B(0:N_theta-1,0:N_zeta-1) )            
         dBdz_l_modes = DFT_2D( theta, zeta, dBdz_l(0:N_theta-1,0:N_zeta-1) )            
         
         
         call plot_contour( x=[(1d0*i,i=-M1/2,M2/2-1)], y=[(1d0*j,j=-N1/2,N2/2-1)], z=abs(  vds_l_modes ), x_label="$m$", y_label="$n$", &
                            path = "results/DKE_coefficients_Fourier_Spectra/vds_l_modes" , graph_type = "color") 
         
         call plot_contour( x=[(1d0*i,i=-M1/2,M2/2-1)], y=[(1d0*j,j=-N1/2,N2/2-1)], z=abs(  g_inv_modes ), x_label="$m$", y_label="$n$", &
                            path = "results/DKE_coefficients_Fourier_Spectra/g_inv_modes", graph_type = "color" )         
         call plot_contour( x=[(1d0*i,i=-M1/2,M2/2-1)], y=[(1d0*j,j=-N1/2,N2/2-1)], z=abs(  B2_modes ), x_label="$m$", y_label="$n$", &
                            path = "results/DKE_coefficients_Fourier_Spectra/B2_modes", graph_type = "color" ) 
         
         call plot_contour( x=[(1d0*i,i=-M1/2,M2/2-1)], y=[(1d0*j,j=-N1/2,N2/2-1)], z=abs(  gB_inv_modes ), x_label="$m$", y_label="$n$", &
                            path = "results/DKE_coefficients_Fourier_Spectra/gB_inv_modes", graph_type = "color" ) 
                            
                                    
         call plot_contour( x=[(1d0*i,i=-M1/2,M2/2-1)], y=[(1d0*j,j=-N1/2,N2/2-1)], z=abs(  B_modes ), x_label="$m$", y_label="$n$", &
                            path = "results/DKE_coefficients_Fourier_Spectra/B_modes", graph_type = "color" ) 
                                    
         call plot_contour( x=[(1d0*i,i=-M1/2,M2/2-1)], y=[(1d0*j,j=-N1/2,N2/2-1)], z=abs(  dBdz_l_modes ), x_label="$m$", y_label="$n$", &
                            path = "results/DKE_coefficients_Fourier_Spectra/dBdz_l_modes", graph_type = "color" ) 
         
         
         
       end subroutine
       
   end subroutine 
  


end module
