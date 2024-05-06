  module Barycentric_grids
  
  use Linear_systems
  
  implicit none
  
  
  
  
  contains
  
  integer function Factorial( n ) result (F)
      integer, intent(in) :: n
      
      integer :: i 
      if ( n <= 1 ) F = 1 
      if ( n > 1 ) F = product( [(i,i=1,n)] ) 
  
  end function
  
  function Barycentric_interpolant ( x, x_nodes, y_nodes, weights ) result(F)
     real, intent(in) :: x(0:), x_nodes(0:), y_nodes(0:)
     real, optional, intent(in) :: weights(0:size(x)-1)
     real :: F(0:size(x)-1)
  
     real :: w(0:size(x_nodes)-1), dx(0:size(x_nodes)-1)  
     integer :: i, j, N, M
     
     if(       present(weights) ) w = weights
     if( .not. present(weights) ) w = Barycentric_weights(x_nodes)
     
     M = size(x) - 1 ; N = size(x_nodes) - 1
     
     do i = 0, M
        dx = [( x(i)-x_nodes(j), j =0,N)]
        
        if( minval( abs(dx) ) > 0 ) then        
           dx = w / dx           
           F(i) = dot_product( dx, y_nodes ) / sum( dx )        
        else           
           F(i) = y_nodes(i)        
        end if      
     end do
     
  end function
  
  
  function Barycentric_derivatives( x, p, weights ) result(D)
     real, intent(in) :: x(0:)
     integer, intent(in) :: p
     real, optional, intent(in) :: weights(0:size(x)-1)
     real :: D(p,0:size(x)-1, 0:size(x)-1)
     
     
     real :: w(0:size(x)-1), A_m    
     integer :: i, j, r, k, m, N 
     
     if(       present(weights) ) w = weights
     if( .not. present(weights) ) w = Barycentric_weights(x)
          
     N = size(x) -1 
     do j = 0, N 
        do i = 0, N 
           do r = 1, p    
              
              if ( i /= j ) then
                D(r,i,j) = Factorial(r) *(-1)**(r+1) * w(j) / w(i)   &
                         / ( x(i) - x(j) )**r                         
                do m = 1, r-1
                
                   A_m = 0
                   do k =0, N
                      if(k/=i) A_m = A_m + (w(k)/w(i)) / ( x(i) - x(k) )**(r-m)                  
                   end do
                   A_m = A_m * (-1)**(r-m) * Factorial(r) / Factorial(m)
                   
                   D(r,i,j) = D(r,i,j) + A_m * D(m,i,j)
                end do                    
              
              else
                D(r,i,j) = 0     
              end if 
           
           end do
        end do     
     end do 
     
     do i = 0, N
        D(:,i,i) = [( -sum(D(r,i,:)), r=1,p )]        
     end do
  
  end function
  
  
  
  function Barycentric_weights(x) result( w )
      real, intent(in) :: x(0:)
      real :: w(0:size(x)-1)
      
      real :: dx(0:size(x)-2) 
      integer :: j, k, N
      
      N = size(x) - 1
      
      do j = 0, N     
          dx = [( x(j)-x(k), k=0, j-1 ), ( x(j)-x(k), k=j+1, N )]  
          w(j) = product( dx ) 
      end do
      
      w = 1 / w 
  
  end function
  
  function Barycentric_weights_OLD(x) result( w )
      real, intent(in) :: x(0:)
      real :: w(0:size(x)-1)
      
      integer :: j, k, N
      
      N = size(x) - 1
      
      w(0) = 1 
      do j = 1, N
         do k = 0, j - 1
            w(k) = w(k) * ( x(k) - x(j) )
         end do
         
          w(j) = product( [( x(j)-x(k), k=0,j-1 )] )            
      end do
      
      w = 1 / w 
  
  end function
  
  ! *** Gives the row of the differentiation matrix of order p that acts
  ! on periodic functions in the domain x \in [x_nodes(0),x_nodes(N)] 
  !
  !   -- INPUT:  x: Point in which the derivative is evaluated.
  !              x_nodes(0:N): Equispaced grid points
  !              p: Order of derivation.
  !   -- INPUT:  x(0:N): Equispaced grid points. 
  pure function Fourier_Derivative_Row( x, x_nodes, p ) result(D)
     real, intent(in) :: x, x_nodes(0:)
     integer, intent(in) :: p 
     real :: D(0:size(x_nodes)-2)
     
     real, parameter :: pi = acos(-1d0)
     complex :: F(0:size(x_nodes)-2), a 
     integer :: j, k, N, N1, N2 
     
     N = size(x_nodes) - 1 ; N1 = N - mod(N,2) ; N2 = N + mod(N,2) 
     a =2*(0,1)*pi/(x_nodes(N)-x_nodes(0))
       
     F = 0
     do j = 0, N-1 
        do k = -N1/2, N2/2-1 
           F(j) = F(j) + (a*k)**p * exp( a*k*(x-x_nodes(j)) )
        end do 
     end do  
     
     D = real(F) / N
  
  end function 
  
  ! *** Gives the row of the Fourier interpolation matrix that acts
  ! on periodic functions in the domain x \in [x_nodes(0),x_nodes(N)] 
  !
  !   -- INPUT:  x: Point in which the derivative is evaluated.
  !              x_nodes(0:N): Equispaced grid points
  
  !   -- INPUT:  x(0:N): Equispaced grid points. 
  pure function Fourier_Interpolation_Row( x, x_nodes ) result(D)
     real, intent(in) :: x, x_nodes(0:) 
     real :: D(0:size(x_nodes)-2)
     
     real, parameter :: pi = acos(-1d0)
     complex :: F(0:size(x_nodes)-2), a 
     integer :: j, k, N, N1, N2 
     
     N = size(x_nodes) - 1 ; N1 = N - mod(N,2) ; N2 = N + mod(N,2) 
     a =2*(0,1)*pi/(x_nodes(N)-x_nodes(0))
       
     F = 0
     do j = 0, N-1 
        do k = -N1/2, N2/2-1 
           F(j) = F(j) + exp( a*k*(x-x_nodes(j)) )
        end do 
     end do  
     
     D = real(F) / N
  
  end function 
   
  
  ! *** Gives the row of the Fourier interpolation matrix that acts
  ! on periodic functions in the domain x \in [x_nodes(0),x_nodes(N)] 
  !
  !   -- INPUT:  x(0:N): Equispaced grid points. 
  !              f(0:N-1): Image of the function at x(0:N-1).
  !   -- OUTPUT: DFT(-N1/2:N2/2-1): Discrete Fourier Modes 
  pure function DFT_1D( x, f ) result(DFT)
     real, intent(in) :: x(0:), f(0:size(x)-1-1) 
     complex :: DFT(-(size(x)-1-mod(size(x)-1,2))/2:(size(x)-1+mod(size(x)-1,2))/2-1)
     
     real, parameter :: pi = acos(-1d0)
     complex, parameter :: ii = (0,1)
     complex :: a 
     integer :: j, k, N, N1, N2 
     
     N = size(x) - 1 ; N1 = N - mod(N,2) ; N2 = N + mod(N,2) 
     a =2*ii*pi/(x(N)-x(0))
     
     do k = -N1/2, N2/2-1
        DFT(k) = sum( f * [( exp(-k*a*x(j)),j=0,N-1)] )/N     
     end do
  
  end function 
  pure function DFT_1Dc( x, f ) result(DFT)
     real, intent(in) :: x(0:) 
     complex, intent(in) :: f(0:size(x)-1-1) 
     complex :: DFT(-(size(x)-1-mod(size(x)-1,2))/2:(size(x)-1+mod(size(x)-1,2))/2-1)
     
     real, parameter :: pi = acos(-1d0)
     complex, parameter :: ii = (0,1)
     complex :: a 
     integer :: j, k, N, N1, N2 
     
     N = size(x) - 1 ; N1 = N - mod(N,2) ; N2 = N + mod(N,2) 
     a =2*ii*pi/(x(N)-x(0))
     
     do k = -N1/2, N2/2-1
        DFT(k) = sum( f *[( exp(-k*a*x(j)),j=0,N-1)] )/N     
     end do
  
  end function 
  
  pure function DFT_2D( x, y, f ) result(DFT)
     real, intent(in) :: x(0:), y(0:), f(0:size(x)-1-1, 0:size(y)-1-1) 
     complex :: DFT( -(size(x)-1-mod(size(x)-1,2))/2:(size(x)-1+mod(size(x)-1,2))/2-1, & 
                     -(size(y)-1-mod(size(y)-1,2))/2:(size(y)-1+mod(size(y)-1,2))/2-1 )
     
     integer :: i, j, M, N, N1, N2, M1, M2      
     complex :: DFT_y(0:size(x)-1-1, & 
                    -(size(y)-1-mod(size(y)-1,2))/2:(size(y)-1+mod(size(y)-1,2))/2-1  )
     
     M = size(x) - 1 ; M1 = M-mod(M,2) ; M2 = M+mod(M,2)
     N = size(y) - 1 ; N1 = N-mod(N,2) ; N2 = N+mod(N,2)
     
     do i = 0, M-1
        DFT_y(i,:) = DFT_1D( y, f(i,:) )     
     end do 
     
     do j = -N1/2, N2/2 -1     
        DFT(:,j) =  DFT_1Dc( x, DFT_y(:,j) )  
     end do 
     
  end function
   
  
  ! *** Gives the row of the Fourier interpolation matrix that acts
  ! on periodic functions in the domain x \in [x_nodes(0),x_nodes(N)] 
  !
  !   -- INPUT:  x(0:N): Equispaced grid points. 
  !              F_modes(-N1/2:N2/2-1): Discrete Fourier Modes. (N1+N2)/2 = N
  !   -- OUTPUT: F(0:N-1): Image of the function at x(0:N-1). 
  pure function IFT_1D( x, F_modes ) result(F)
     real, intent(in) :: x(0:) 
     complex, intent(in) :: F_modes(-(size(x)-1-mod(size(x)-1,2))/2:(size(x)-1+mod(size(x)-1,2))/2-1)
     real :: F(0:size(x)-1-1)
     
     real, parameter :: pi = acos(-1d0)
     complex, parameter :: ii = (0,1)
     complex :: a 
     integer :: j, k, N, N1, N2 
     
     N = size(x) - 1 ; N1 = N - mod(N,2) ; N2 = N + mod(N,2) 
     a =2*ii*pi/(x(N)-x(0))
     
     do j = 0, N-1
        F(j) = real( sum( F_modes * [( exp(k*a*x(j)),k=-N1/2,N2/2-1)]) ) 
     end do
  
  end function 
  
  pure function IFT_1Dc( x, F_modes ) result(F)
     real, intent(in) :: x(0:) 
     complex, intent(in) :: F_modes(-(size(x)-1-mod(size(x)-1,2))/2:(size(x)-1+mod(size(x)-1,2))/2-1)
     complex :: F(0:size(x)-1-1)
     
     real, parameter :: pi = acos(-1d0)
     complex, parameter :: ii = (0,1)
     complex :: a 
     integer :: j, k, N, N1, N2 
     
     N = size(x) - 1 ; N1 = N - mod(N,2) ; N2 = N + mod(N,2) 
     a =2*ii*pi/(x(N)-x(0))
     
     do j = 0, N-1
        F(j) =  sum( F_modes * [( exp(k*a*x(j)),k=-N1/2,N2/2-1)]) 
     end do
  
  end function 
  
  
  pure function IFT_2D( x, y, F_modes ) result(F)
     real, intent(in) :: x(0:), y(0:)
     complex, intent(in)  :: F_modes(-(size(x)-1-mod(size(x)-1,2))/2:(size(x)-1+mod(size(x)-1,2))/2-1, &
                        -(size(y)-1-mod(size(y)-1,2))/2:(size(y)-1+mod(size(y)-1,2))/2-1) 
     real :: F(0:size(x)-1-1,0:size(y)-1-1)
     
     complex :: IFT_x(0:size(x)-1-1, &
                      -(size(y)-1-mod(size(y)-1,2))/2:(size(y)-1+mod(size(y)-1,2))/2-1) 
     integer :: i, j, M, N, N1, N2
     
     M = size(x) - 1 ; N = size(y) - 1 ; N1 = N-mod(N,2) ; N2 = N+mod(N,2)
     
     do j = -N1/2, N2/2 - 1     
        IFT_x(:,j) = IFT_1Dc( x, F_modes(:,j) )
        
     end do
     
     do i = 0, M-1
        F(i,:) = IFT_1D( y, IFT_x(i,:) )
     end do
     
  end function 
  
  
  
  ! Gives the n+1 zeros of an orthogonal polynomial p_{n+1} which satisfies
  ! the three-term recurrence relation 
  ! u(n) p_{n+1}(x) = ( x - d(n) ) p_{n}(x) - l(n) p_{n-1}(x)
  ! where l(n) * u(n) > 0.
  !    -- Inputs: l, d, u : vectors of size n+1 containing the coefficients of the recurrence
  !                     n : number of zeros calculated minus one. 
  !    -- Output:       x : vector of size n+1 containing the zeros of p_{n+1}(x).
  subroutine Zeros_orthogonal_polynomial( l, d, u, n, x, w ) 
     real, intent(in) :: l(0:n), d(0:n), u(0:n)
     integer, intent(in) :: n
     real, intent(out) :: x(0:n), w(0:n) 
     
     real :: Jacobi(0:n,0:n), V(0:n,0:n) 
     integer :: k
     
     ! Jacobi matrix  
     Jacobi = 0 
     do k = 0, n
        Jacobi(k,k) = d(k)  
        if( k < n ) Jacobi(k,k+1) = sqrt( l(k+1) * u(k) )
        if( k < n ) Jacobi(k+1,k) = Jacobi(k,k+1)
     end do
     
     ! Zeros as eigenvalues of the Jacobi matrix  
     call Eigenvalues_Jacobi_LAPACK( Jacobi, x, V ) 
     w = V(0,:)**2 
     
  end subroutine 
  
  
  ! Gives the n+1 zeros of an orthonormal polynomial p_{n+1} which satisfies
  ! the three-term recurrence relation 
  ! l(n+1) p_{n+1}(x) = ( x - d(n) ) p_{n}(x) - l(n) p_{n-1}(x)
  ! where l(n) > 0.
  !    -- Inputs: l, d : vectors of size n+1 containing the coefficients of the recurrence
  !                  n : number of zeros calculated minus one. 
  !    -- Output:    x : vector of size n+1 containing the zeros of p_{n+1}(x).
  subroutine Zeros_orthonormal_polynomial( l, d, n, x, w ) 
     real, intent(in) :: l(0:n), d(0:n) 
     integer, intent(in) :: n
     real, intent(out) :: x(0:n), w(0:n) 
     
     real :: Jacobi(0:n,0:n), V(0:n,0:n) 
     integer :: k
     
     ! Jacobi matrix  
     Jacobi = 0 
     do k = 0, n
        Jacobi(k,k) = d(k)  
        if( k < n ) Jacobi(k,k+1) = l(k+1)  
        if( k < n ) Jacobi(k+1,k) = Jacobi(k,k+1)
     end do
     
     ! Zeros as eigenvalues of the Jacobi matrix  
     call Eigenvalues_Jacobi_LAPACK( Jacobi, x, V ) 
     w = V(0,:)**2 
     
  end subroutine 
  
  
  
  ! Given an integer N gives the ORTHONORMAL Maxwell polynomials up to degree N.
  ! Maxwell polynomials are orthogonal w.r.t to the weight function 
  ! exp(-x^2)  in the domain [0, infinity). 
  ! It has been tested that the tuning parameters "cut_off", "limit", "key"
  ! "abs_tol", "rel_tol" work up to 
  !  "N=~200" (for Visual Studio project using ifort in Javi's laptop and ifort version in XULA). 
  !  "N=~95" (for gfortran Ubuntu version in Javi's laptop). 
  ! 
  ! For bigger values of N, the algorithm fails due to ill-conditioning of the moment approach.
  ! (Perhaps in the future implement with a Chebyshev modified moment algorithm).
  subroutine Maxwell_polynomials( N, l, d, x, w )
     integer, intent(in) :: N
     real, intent(out) :: l(0:N+1), d(0:N)
     real, optional, intent(out) :: x(0:N), w(0:N)
     
     integer, parameter :: inf = 1, key = 6 
     real, dimension(18*N) :: alist, blist, rlist, elist !workspacesize=18*N
     integer, dimension(18*N) :: iord
     real, parameter :: pi = acos(-1d0), abs_tol = 0, rel_tol =1e-13, cut_off = 100
     real :: mu_0 , abs_err, ll
     integer :: i, k, N_evals, ierr, last, limit
     
     limit = 18*N
     mu_0 = sqrt(pi/4)     
     
     l = 1 ;  d = 0
     l(0) = 0 ; d(0) = 0.5/mu_0 ! Analytical
     
     do i = 0, N ; k = i
        ! *** Compute d(k) 
        if( i > 0 )  &
          ! Integral between [0, cut_off]  as approximation to [0, infinity) 
          call dqage( Integrand_d_k, 0d0, cut_off, abs_tol, rel_tol, key, limit,d(i), abs_err, &
                      N_evals, ierr, alist, blist, rlist, elist, iord, last)           
                     
        ! *** Compute l(k+1)**2  
        ! Integral between [0, cut_off]  as approximation to [0, infinity) 
        call dqage( Integrand_l_k1, 0d0, cut_off, abs_tol, rel_tol, key, limit, ll, abs_err , &
                    N_evals, ierr, alist, blist, rlist, elist, iord, last) 
        
        ! Compute l(k+1)
        l(i+1) = sqrt( ll )  
                 
     end do 
     
     ! *** Zeros and weights of Maxwell polynomials   
     if( present(x) .and. present(w) ) then
       call Zeros_orthonormal_polynomial( l, d, N, x, w )    
       w = w * mu_0 
     end if 
     
     contains

  function Orthon_polynom(x,n) result(P)
     real, intent(in) :: x
     integer, intent(in) :: n
     real :: P
     
     real :: P0, P1 
     integer :: k
     
     ! P_{-1} = 0, P_0 = 1/sqrt(mu_0)
     P0 = 0 ; P1 = 1/sqrt(mu_0) ; P = P1
     do k = 0, n-1
        P = ( x - d(k) ) * P1  - l(k) * P0
        P = P / l(k+1)
        
        P0 = P1
        P1 = P 
     end do
     
  
  end function 
  
  function Weight_function(x) result(w)
     real, intent(in) :: x
     real :: w
    
     w = exp(-x**2) 
  
  end function
  
  function Integrand_d_k(x) result(F)
     real, intent(in) :: x
     real :: F
     
     real :: w, p_k
     
     p_k = Orthon_polynom(x,k)
     w = Weight_function(x) 
     F = x * w * p_k **2
  
  end function 
  
  function Integrand_l_k1(x) result(F)
     real, intent(in) :: x
     real :: F
     
     real :: w, p_k1
     
     p_k1 = Orthon_polynom(x,k+1)
     w = Weight_function(x) 
     F = w * p_k1**2
  
  end function 
  
  end subroutine 
  
  
  
  
  end module
