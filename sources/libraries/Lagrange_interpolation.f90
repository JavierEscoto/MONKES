module Lagrange_interpolation

 implicit none 
 private
 
 interface Stencil
    module procedure Centered_Stencil, Upwind_Stencil
 end interface
 
 
 
 public :: Lagrange_polynomial, Lagrange_polynomials_FJEL ! *** FJEL Lagrange polynomials
 public :: Lagrange_derivative, Lagrange_derivatives      ! *** FJEL Lagrange derivatives
 public :: Lagrange_integrals                             ! *** FJEL Lagrange integrals
 public :: Interpolant, Interpolant_derivatives           ! *** FJEL Interpolant and derivatives
 public :: Interpolant_mask                               ! *** FJEL Interpolant with mask
 
 public :: Interpolated_value, Interpolated_derivatives  ! *** FJEL Interpolated value and derivatives
 
 public :: Lagrange_polynomials                           ! *** JAHR Lagrange polynomials, derivatives and integrals
 public :: Interpolant_JAHR   ! *** JAHR Piecewise Interpolation      
 public :: Interpolated_value_JAHR ! *** JAHR Interpolated value
 public :: Stencil
 
contains 
   
    
     
!*******************************************************************************************************************************************
!                                           Lagrange polynomials
!
!  Recurrence relation  to obtain derivatives and integral values of the Lagrange polynomials at xp from a set of nodes x(0:N)
!  
!
!      lagrange_j{N+1} (x) = (x-x0)(x-x1)(x-2)........(x-xn) / (xj- x0)(xj-x1)(xj-x2).......(xj-xn), j=0...N
!      lagrange_j{N+1} (x) = (x -xn)/(xj-xn) lagrange_j{N} (x) 
!
!      d^k lagrange_j{r+1} (x) /dx^k  = ( d^k lagrange_j{r} (x) /dx^k (x-xr) + k d^(k-1) lagrange_j{r} (x) /dx^(k-1) ) / (xj -xr ), r=0...N  
!
!      integral( lagrange_j{N+1}(x), from x0 to xp ) = integral ( l(xp) + dl/dx(xp) ( x -xp) +....... d^N l/dx^N(xp) * (x- xp)**N / N! ) 
!                                                    = sum_k (  -  dl^k/dx^k (xp) ( xj - xp )**(k+1) /(k+1)! )  
! 
!      k = -1 means integral value of lagrange_j{N+1} (x) from xj to xp (xj  nearest node from xp)
!      k =  0 means value of lagrange_j{N+1} (x) at xp 
!      k > 1    means derivative value of lagrange_j{N+1} (x) at xp 
!
! It returns Lagrange_polynomials(k,j) 
!
!     k: meaning given above 
!     j: corresponding wth L_j(x) (this polynomial equals to one in x_j and vanishes at the rest of nodes) 
!
! Author : Juan A Hernandez (juanantonio.hernandez@upm.es) 
!*****************************************************************************************************************************************
 pure function Lagrange_polynomials( x, xp ) !
   real, intent(in) :: x(0:), xp
   real ::  Lagrange_polynomials(-1:size(x)-1,0:size(x)-1) 


   integer :: j   ! node 
   integer :: r   ! recursive index 
   integer :: k   ! derivative 
   integer :: Nk  ! maximum order of the derivative 
   integer :: N, j1(1), jp 

   real :: d(-1:size(x)-1) 
   real :: f 

   Nk = size(x) - 1 
   N  = size(x) - 1 

   do j = 0, N 
      d(-1:Nk) = 0 
      d(0) = 1

 ! ** k derivative of lagrange(x) at xp 
      do  r = 0, N 

         if (r/=j) then 
          do k = Nk, 0, -1
            d(k) = ( d(k) *( xp - x(r) ) + k * d(k-1) ) /( x(j) - x(r) )  
          end do 
        endif 

      enddo 

!  ** integral of lagrange(x) from x(jp) to xp 
      f = 1 
      j1 = minloc( abs(x - xp) ) - 2
      
      jp = max(0, j1(1)) 
      do k=0, Nk 
         f = f * ( k + 1 ) 
         d(-1) = d(-1) -  d(k) * ( x(jp) - xp )**(k+1)
      enddo 

      Lagrange_polynomials(-1:Nk, j ) = d(-1:Nk) 

  end do 

end function 


!*******************************************************************************************************************************************
!                                           Lagrange polynomial
! Computes the grade q Lagrange polynomial l0 at a point x in [x0,xq] 
! given a set of nodes x0, x1,..., xq. 
!  
! Author : F. Javier Escoto (javier.escoto.lopez@gmail.com) 
!*****************************************************************************************************************************************
real pure function Lagrange_polynomial(x, x_nodes) result (L)
    real, intent(in) :: x, x_nodes(0:)
    integer ::  j, q
    real :: x0, xj
    
    q = size(x_nodes) - 1
    
    L = 1
    x0 = x_nodes(0)   
    do j = 1, q
        xj = x_nodes(j)
        L = L * ( x - xj ) / (x0 - xj)
    end do
    
end function 

!*******************************************************************************************************************************************
!                                           Lagrange polynomialS
! Computes the value at x of ALL the Lagrange polynomials l_j(x),  j=0,1...,q
!     given a set of nodes x_0, x_1, ..., x_q 
!
!  Author: Francisco Javier Escoto Lopez 
!          (javier.escoto.lopez@gmail.com)
! *******************************************************************************************************************
pure function Lagrange_polynomials_FJEL(x, x_nodes) result(L)
      real, intent(in) :: x, x_nodes(0:)
      real :: L(0:size(x_nodes)-1)
      real :: s(0:size(x_nodes)-1)
      integer :: j, q
      
      q = size(x_nodes)-1
      
   ! *** Lagrange_polynomials using permutations in x_nodes 
     do j = 0, q
        s = [ x_nodes(j:q), x_nodes(0:j-1) ]
         
        L(j) = Lagrange_polynomial( x, s )  
     end do
      
end function



!*******************************************************************************************************************************************
!                                           Lagrange derivative
! Computes the q derivatives and the value of the Lagrange polynomial l_0(x) at a point x in [x0,xq] 
! given a set of nodes x0, x1,..., xq. 
!  
!  Author: Francisco Javier Escoto Lopez 
!          (javier.escoto.lopez@gmail.com) 
!*****************************************************************************************************************************************
pure function Lagrange_derivative(x, x_nodes) result (dL)
    real, intent(in) :: x, x_nodes(0:)
    real ::  dL(0:size(x_nodes)-1)
    real :: d(-1:size(x_nodes)-1)
    real :: x0, xj
    integer :: j, k, q
    
    q = size(x_nodes)-1
    
    x0 = x_nodes(0) 
    d = 0  ; d(0) = 1
    ! *** Recursion to obtain the derivatives in descendent order
    do j = 1, q       
        xj = x_nodes(j)
        do k = q, 0, -1 
        d(k) = ( d(k) *( x - xj ) + k * d(k-1) ) /( x0 - xj )  
        end do
    end do
    dL = d(0:q)

end function


!*********************************************************************************************************
!                                           Lagrange derivativeS
! Computes the Lagrange polynomial and first q derivatives and the values for ALL the Lagrange polynomials l_j(x),
! j=0,1...,q at a point x in [x_0,x_q]  given a set of nodes x_0, x_1, ..., x_q 
!
!  Author: Francisco Javier Escoto Lopez 
!          (javier.escoto.lopez@gmail.com)
! ********************************************************************************************************
pure function Lagrange_derivatives(x, x_nodes) result(dL)
      real, intent(in) :: x, x_nodes(0:)
      real ::  dL(0:size(x_nodes)-1,0:size(x_nodes)-1)
      real :: s(0:size(x_nodes)-1)
      integer :: j, q
      
      q = size(x_nodes)-1
      
      ! *** Derivatives calculation by permutations on x_nodes
      do j = 0, q
         s = [ x_nodes(j:q), x_nodes(0:j-1) ]
         dL(:,j) = Lagrange_derivative(x, s)
      end do
      

end function


!*********************************************************************************************************
!                                           Lagrange derivativeS
! Computes the first q derivatives and the values for ALL the Lagrange polynomials l_j(x),  j=0,1...,q
! at a point x in [x_0,x_q]  given a set of nodes x_0, x_1, ..., x_q 
!
!  Author: Francisco Javier Escoto Lopez 
!          (javier.escoto.lopez@gmail.com)
! ********************************************************************************************************
pure function Lebesgue_functions(x, x_nodes) result(L)
      real, intent(in) :: x, x_nodes(0:)
      real :: L(0:size(x_nodes)-1)
      
      real :: dL(0:size(x_nodes)-1,0:size(x_nodes)-1) 
      integer :: j, q
      
      q = size(x_nodes)-1

      dL = Lagrange_derivatives(x, x_nodes) ! Lagrange derivatives
      ! *** Lebesgue functions as sums of Lagrange derivatives
      do j = 0, q
          L(j) = sum( abs( dL(:,j) ) )
      end do
      

end function

!*********************************************************************************************************
!                                           Lagrange Integrals
! Computes the integral for ALL the Lagrange polynomials l_j(x),  j=0,1...,q in the interval [x_0, x] (x<=x_q)
!   given a set of nodes { x_nodes(0),..., x_nodes(q) }
!
!  Author: Francisco Javier Escoto Lopez 
!          (javier.escoto.lopez@gmail.com)
! ********************************************************************************************************
pure function Lagrange_integrals(x, x0, x_nodes, derivatives) result(Integral)
      real, intent(in) :: x, x0, x_nodes(0:)
      real, optional, intent(in) :: derivatives(0:size(x_nodes)-1,0:size(x_nodes)-1)
      real :: Integral(0:size(x_nodes)-1)
      
      real ::  dL(0:size(x_nodes)-1,0:size(x_nodes)-1)
      real ::  Monomial(0:size(x_nodes)-1) 
      integer :: j, q, f, jq
      
      q = size(x_nodes)-1
            
      ! *** Derivatives (centered) of the Lagrange polynomials 
      if ( present(derivatives) ) then
          jq = ( q - mod(q,2) )/2 
          dL = derivatives 
      else
          ! *** Centered point of the stencil
          jq = ( q - mod(q,2) )/2 ! q/2 for even q and (q-1)/2 for odd q
          dL = Lagrange_derivatives( x_nodes(jq), x_nodes )
      end if
      
      ! *** Monomials of the Taylor expansion
      f = 1
      do j = 0, q
          f = f * (j+1)
          Monomial(j) = (  x - x_nodes(jq) )**(j+1) / f  &
                      - ( x0 - x_nodes(jq) )**(j+1) / f
      end do
      
      ! *** Integral of the Taylor expansion
      Integral = matmul( Monomial, dL)

end function


! ************************************************************
!                 STENCIL
!
!   Given a set of nodes x(0:), an integer i and an order
!   of interpolation q, computes the centered stencil S(0:q) that uses
!   the Lagrange polynomial  l_i(x)  
!
!  Author: Francisco Javier Escoto Lopez 
!          (javier.escoto.lopez@gmail.com)
! ***********************************************************

pure function Centered_Stencil(x,i,q) result(S)
     real, intent(in)   :: x(0:)
     integer,intent(in) :: i, q
                integer :: S(0:q)
                   
     integer :: N, first, j
     
     N = size(x) - 1
          
     ! *** Starting point of the stencil
     first = i - (q-mod(q,2))/2
     
     ! *** Stencil definition
     if( first < 0 ) then
         s = [(j,j=0,q)]       
     elseif ( first > N - q ) then
         s = [( N-q   + j, j=0,q )]   
     else
         s = [( first + j, j=0,q )]
     end if
                
end function

pure function Upwind_Stencil( x, i, q, sign, i_0 ) result(s)
     real, intent(in) :: x(0:)
     integer, intent(in) :: i, q, sign
     integer, optional, intent(in) :: i_0 ! first index value for the grid
     integer :: s(0:q)
                   
     integer :: j, N
     
     N = size(x) - 1
          
     ! *** Stencil definition
     if ( sign > 0 ) then
        if( i < q )  s = [( j, j = 0,q )]
        if( i >= q ) s = [( i+j-q, j=0,q )] 
     elseif ( sign < 0 ) then
        if( i <= N - q ) s = [( i+j  , j=0,q )]
        if( i  > N - q ) s = [( N-q+j, j = 0,q )] 
     else
        s = Stencil(x, i, q)
     end if
     
     ! *** For grids that do not start at i=0
     if ( present(i_0) ) s = s + i_0
end function


pure function Stencil_S(x,i,q,sign) result(S)
     real, intent(in)   :: x(0:)
     integer,intent(in) :: i, q, sign
                integer :: S(0:q)
                        
     ! *** Stencil centered or upwind
     if( sign == 0 ) then
         s = Stencil(x,i,q)    
     else
         s = Upwind_Stencil(x,i,q,sign) 
     end if
                
end function

! ************************************************************
!                 Stencil_1
!
!   Given a set of nodes x(1:), an integer i >= 1 and an order
!   of interpolation q, computes the stencil S= [s1, s2] that uses
!   the Lagrange polynomial  l_i(x)  
!
!  Author: Francisco Javier Escoto Lopez 
!          (javier.escoto.lopez@gmail.com)
! ***********************************************************
 pure function Stencil_1(x,i,q) result(S)   !pure
     real, intent(in)   :: x(:)
     integer,intent(in) :: i, q
     integer :: S(2)
                   
     integer :: S0(0:q)
     real ::  x0(0:size(x)-1)
     
     x0 = x 
     S0 = Stencil( x0, i-1, q )
     S = [S0(0), S0(q)] + 1
                
end function


!*************************************************************************************************************
!* Computes the piecewise polynomial interpolation of I(xp) from the data x(:) and y(:). 
!  The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!*************************************************************************************************************
function Interpolant_JAHR(x, y, degree, xp )
    real, intent(in) :: x(0:), y(0:), xp(0:) 
    integer, intent(in) :: degree
    real :: Interpolant_JAHR(0:degree, 0:size(xp)-1)
    
    integer :: N, M, s, i, j, k
    
!   maximum order of derivative and width of the stencil 
    integer :: Nk 
    
!   Lagrange coefficients and their derivatives at xp 
    real, allocatable :: Weights(:,:) 
    
    N = size(x) - 1
    M = size(xp) - 1 
    Nk = degree
    allocate( Weights(-1:Nk, 0:Nk))
    
do i=0, M 
        
    j = max(0, maxloc(x, 1, x < xp(i) ) - 1) 

    if( (j+1) <= N ) then ! the (j+1) cannot be accessed
        if( xp(i) > (x(j) + x(j + 1))/2 ) then   
            j = j + 1
        end if
    end if
   

    if (mod(Nk,2)==0) then 
        s = max( 0, min(j-Nk/2, N-Nk) )   ! For Nk=2 
    else 
        s = max( 0, min(j-(Nk-1)/2, N-Nk) )
    endif 

  
    Weights(-1:Nk, 0:Nk)=Lagrange_polynomials( x = x(s:s+Nk), xp = xp(i) )
     
    do k=0, Nk 
        Interpolant_JAHR(k, i) = sum ( Weights(k, 0:Nk) * y(s:s+Nk) ) 
    end do    
     
end do  
 
    deallocate(Weights)

end function


!*************************************************************************************************************
!* Computes the piecewise polynomial interpolation of I(xp) from the data x(:) and y(:). 
!  The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!*************************************************************************************************************
pure real function Interpolated_value_JAHR(x, y, xp, degree)!
    real, intent(in) :: x(0:), y(0:), xp
    integer, optional, intent(in) :: degree
    integer :: N, s, j
    
!   maximum order of derivative and width of the stencil 
    integer :: Nk !
    
!   Lagrange coefficients and their derivatives at xp 
    real, allocatable :: Weights(:,:) 

    N = size(x) - 1

    
    if(present(degree))then
        Nk = degree
    else
        Nk = 2
    end if

    allocate( Weights(-1:Nk, 0:Nk))
  
    j = max(0, maxloc(x, 1, x < xp ) - 1) 

    if( (j+1) <= N ) then ! if not the (j+1) cannot be accessed
        if( xp > (x(j) + x(j + 1))/2 ) then   
            j = j + 1
        end if
    end if
   
    ! *** First index of the stencil
    if (mod(Nk,2)==0) then 
        s = max( 0, min(j-Nk/2, N-Nk) )   ! For even Nk
    else 
        s = max( 0, min(j-(Nk-1)/2, N-Nk) )
    endif 

    Weights(-1:Nk, 0:Nk)  = Lagrange_polynomials( x = x(s:s+Nk), xp = xp )
    Interpolated_value_JAHR = sum ( Weights(0, 0:Nk) * y(s:s+Nk) ) 
    
    deallocate(Weights)

end function

! *****************************************************************
!       Interpolant(x, x_nodes, y_nodes, q)
! Given a set of nodes {x_nodes(0),...,x_nodes(N)} 
! and the data {y_nodes(0),...,y_nodes(N)}
! calculates the interpolant F(x) of order q
!  evaluated at {x(0),..., x(M)}  
! *****************************************************************
pure function Interpolant(x, x_nodes, y_nodes, q) result(F)
         real, intent(in) :: x(0:), x_nodes(0:), y_nodes(0:)
      integer, intent(in) :: q
      real :: F(0:size(x)-1)
      
      integer :: i, M     
      
      ! *** Number of evaluation points (M+1)
      M = size(x) - 1
    
      ! *** Interpolant at x(0:M)  
      F = [( Interpolated_value(x(i), x_nodes, y_nodes, q), i=0,M )]
      
end function




 function Interpolant_mask(x, x_nodes, y_nodes, q, mask) result(F)
         real, intent(in) :: x(0:), x_nodes(0:), y_nodes(0:)
      integer, intent(in) :: q
      logical, intent(in) :: mask(0:)
      real :: F(0:size(x)-1)
      
      integer :: i, M     
      real :: xx(0:count(mask)-1), yy(0:count(mask)-1), FF(0:count(mask)-1)
      
      F = 0
      ! *** Number of evaluation points (M+1)
      M = size(x)-1
      xx = pack(x_nodes, mask) 
      yy = pack(y_nodes, mask) 
      ! *** Interpolant at x(0:M)  
      FF = [( Interpolated_value(x(i), xx, yy, q), i=0,M )]
      
      where( mask ) F = FF
end function


! *****************************************************************
!       Interpolant_derivatives(x, x_nodes, y_nodes, q)
! Given a set of nodes {x_nodes(0),...,x_nodes(N)} 
! and the data {y_nodes(0),...,y_nodes(N)}
! calculates the interpolant F(x) of order q and its derivatives
!     evaluated at {x(0),..., x(M)}  
! *****************************************************************
pure function Interpolant_derivatives(x, x_nodes, y_nodes, q) result(F)
         real, intent(in) :: x(0:), x_nodes(0:), y_nodes(0:)
      integer, intent(in) :: q
      real :: F(0:q,0:size(x)-1)
       
      integer :: i, M     
      
      ! *** Number of evaluation points (M+1)
      M = size(x) - 1
    
      ! *** Interpolant at x(0:M)  
      do i = 0, M
      F(:,i) = Interpolated_derivatives(x(i), x_nodes, y_nodes, q)
      end do 
end function

! *****************************************************************
!       Interpolated_value
! Given a set of nodes { x_nodes(0),...,x_nodes(N) } 
! and the data { y_nodes(0),...,y_nodes(N) }
! calculates the interpolant evaluated at x
! *****************************************************************

real pure function Interpolated_value(x, x_nodes, y_nodes, q) result(F)
         real, intent(in) :: x, x_nodes(0:), y_nodes(0:)
      integer, intent(in) :: q
      
      real :: L(0:q) 
      integer :: js, q_eff        
      integer, allocatable :: s(:)    
      
      ! *** Stencil centered at the nearest point x(js)
      js = minloc( abs(x_nodes - x) , 1) - 1
      q_eff = minval( [q,size(x_nodes)-1] ) ; allocate( s(0:q_eff) )
      s = Stencil(x_nodes,js,q_eff) 
      
      ! *** Lagrange polynomials evaluated at x from stencil nodes
      L = Lagrange_polynomials_FJEL( x, x_nodes(s) ) 
      
      ! *** Interpolated value at x
      F = dot_product( y_nodes(s), L  )
end function


! *****************************************************************
!       Interpolated_derivatives
! Given a set of nodes { x_nodes(0),...,x_nodes(N) } 
! and the data { y_nodes(0),...,y_nodes(N) }
! calculates the interpolant and its derivatives evaluated at x
! *****************************************************************

pure function Interpolated_derivatives(x, x_nodes, y_nodes, q) result(F)
         real, intent(in) :: x, x_nodes(0:), y_nodes(0:)
      integer, intent(in) :: q
      real :: F(0:q)
                     
      real :: L(0:q,0:q) 
      integer :: s(0:q), js       
      
      ! *** Stencil centered at the nearest point x(js)
      js = minloc( abs(x_nodes - x) , 1) - 1
      s = Stencil(x_nodes,js,q) 
      
      ! *** Lagrange polynomials evaluated at x from stencil nodes
      L = Lagrange_derivatives( x, x_nodes(s) ) 
      
      ! *** Interpolated value at x
      F = matmul( L, y_nodes(s)  )
end function




end module
