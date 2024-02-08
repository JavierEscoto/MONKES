module API_Example_Barycentric_grids

  use Barycentric_grids
  use Magnetic_configuration


contains

  ! Example showing how the p first derivatives
  ! of a function are approximated using the DFT


  subroutine Example_Convergence_DFT_1D


  end subroutine
  
  subroutine Example_DFT_1D
    integer, parameter :: N=1000, p=4, M=111, N1 = N-mod(N,2), N2 = N + mod(N,2) 
    real, parameter :: pi = acos(-1d0)
    real :: x_nodes(0:N), u_nodes(0:N-1), x(0:M), u(0:M,0:p), u_Exact(0:M)
    complex :: u_Modes(-N1/2:N2/2-1)
    integer :: i, k

    x_nodes = [(i*2*pi/N,i=0,N)]
    x = [(i*2*pi/M,i=0,M)]

    u_nodes = [( Exact(x_nodes(i)),i=0,N-1)]
    u_Exact = [( Exact( x(i) ),i=0,M )]
    u_Modes = DFT_1D( x_nodes, u_nodes ) 
    
    do k = 0, p
       do i = 0, M
          if( k==0 ) &
            u(i,k) = dot_product( Fourier_Interpolation_Row(x(i),x_nodes), u_nodes )
          if( k>0 ) &
            u(i,k) = dot_product( Fourier_Derivative_Row(x(i),x_nodes,k), u_nodes )
       end do
    end do
    
    open(21,file="results/Example_DFT_1D.plt")
    do i = 0, M
       write(21,*) " x "," u ", " u_Exact ", " Error ", " Derivatives "
       write(21,*)  x(i)  , u(i,0) ,  u_Exact(i) , u(i,0) - u_Exact(i) , u(i,1:p)
       
       end do 
    close(21)

    
    open(21,file="results/Example_DFT_1D_Modes.plt")
    do k = -N1/2, N2/2-1
       write(21,*) " Mode_number ", " Real part ", "Complex part ", " Modulus "
       write(21,*)  real(k), u_Modes(k) % Re, u_Modes(k) % Im, norm2( [ u_Modes(k) % Re, u_Modes(k) % Im ] )
       
    end do
    
    close(21)  

    contains

      real function Exact(x) result(f)
        real, intent(in) :: x

        integer, parameter :: N = 2000
        integer :: k
        real :: a(0:N)

        a = [ 1d0 - 2 * pi**2/6 , (-4d0*(-1)**k * exp(-0.01*k*1d0), k=1,N )]

        f = dot_product( a, [(cos(k*x),k=0,N )] ) 
        
        !f = exp(-sin(x**2) )
        !f = sqrt( sin(x)**2 )
        !if ( x <= pi ) f = x
        !if ( x > pi ) f = 2*pi-x

      end function
      
    
  end subroutine

  
end module 
