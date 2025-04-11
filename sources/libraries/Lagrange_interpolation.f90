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
 public :: Linear_Interpolation                                
 
 public :: Interpolated_value, Interpolated_derivatives  ! *** FJEL Interpolated value and derivatives  
 public :: Stencil
 
contains 
   
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

 ! Linear interpolation of f(x) using y_nodes = f(x_nodes)   
 real pure function Linear_Interpolation(x, x_nodes, y_nodes) result(F)   !pure
     real, intent(in)   :: x, x_nodes(0:), y_nodes(0:) 
     real :: L(0:1), dx 
     integer :: S(0:1), j, N
     
     N = size(x_nodes) - 1
     
     j = minloc( abs(x - x_nodes), 1 ) -1    
     if( x_nodes(j) ==  x ) then
       F = y_nodes(j)
     else
     
       if( x_nodes(j) <  x ) then
               
         if ( j == N ) then 
           S = [ j-1, j ] 
         else
           S = [ j, j+1 ] 
         end if 
         
       elseif( x < x_nodes(j) ) then
         
         if ( j == 0 ) then  
           S = [ j, j+1 ] 
         else
           S = [ j-1, j ] 
         end if 
         
       end if
       
       L(0) = ( x-x_nodes(s(1)) )/( x_nodes(s(0))-x_nodes(s(1)) )
       L(1) = ( x-x_nodes(s(0)) )/( x_nodes(s(1))-x_nodes(s(0)) )
       
       F = dot_product( L, y_nodes(s) )
       
     end if
      
                
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
      integer :: s(0:q), js 
      
      ! *** Stencil centered at the nearest point x(js)
      js = minloc( abs(x_nodes - x) , 1) - 1 
      s = Stencil(x_nodes,js,q) 
      
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
