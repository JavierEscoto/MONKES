module Finite_differences

    use Lagrange_interpolation
    use Sparse_Matrices
    
    implicit none
   
    private 
	
    public :: Grid_Initialization, Initialize_grid
    public :: Derivative_weights, Upwind_Derivative_weights
    public :: Integral_weights, Antiderivative_weights
    public :: Grid_Stencil, Upwind_Grid_Stencil
    
    public :: Chebyshev_points 
    
    public :: Derivative_Delta
    public :: Derivative 
    
    public :: Derivative_ij, Derivative_F_ij, Derivative_B_ij
    
    public :: Backward_Derivative, Forward_Derivative
    
    public :: Integral_Delta
    public :: Integral, Antiderivative
    public :: Characterize_nodes, Index_loc
    
    public :: Stencil_Differential_operator 
    public :: Stencil_Differential_operator_mask
    
    public :: Inverse_Stencil, Upwind_Inverse_Stencil, Inverse_Stencil_mask
    
    public :: Periodic_Grid_Stencil, Periodic_Derivative_weights
    public :: Grid_bounds_1D, Periodic_Grid_bounds_1D
    public :: Interval_bounds
    public :: Extended_periodic_grid
    
    public :: Grids, Grid
    public :: Grid_Sparse_derivative_matrix, Grid_Sparse_derivative_2D_matrix
    public :: Sparse_derivative_1D_row, Sparse_upwind_derivative_1D_row
    public :: Sparse_derivative_2D_row, Sparse_derivative_2D_column
    
    interface Derivative

        module procedure Derivative_1D, Derivative_2D, Derivative_3D, &
                         Derivative_1D_mask, Derivative_2D_mask, &
                         Derivative_3D_mask
                         
        module procedure Derivative_FD_1D_ij, Derivative_FD_1D_ij_indices
        module procedure Derivative_FD_1D_upwind_ij, Derivative_FD_1D_upwind_ij_indices
        
                         
        module procedure Derivative_1D_point, Derivative_1D_mask_point
        module procedure Derivative_1D_indices_point
    end interface
    
    interface Backward_Derivative

        module procedure Backward_Derivative_1D_point
        module procedure Backward_Derivative_1D_mask_point
        
        module procedure Backward_Derivative_1D_indices_point
                                                  
    end interface
    
    interface Forward_Derivative

        module procedure Forward_Derivative_1D_point
        module procedure Forward_Derivative_1D_mask_point
                                                  
    end interface
    
    interface Derivative_Delta
        module procedure Derivative_1D_point_Delta_j
    end interface
    
    interface Integral_Delta
        module procedure Integral_1D_Delta_j, Integral_2D_Delta_ij
    end interface

    interface Integral

        module procedure Integral_1D, Integral_2D, &
                         Integral_1D_mask, Integral_1D_mask_2, &
                         Integral_2D_mask, &
                         Definite_Integral,  Definite_Integral_2
    end interface

    interface Antiderivative

        module procedure Antiderivative_1D, &
                         Antiderivative_2D, &
                         Antiderivative_3D, &
                         Antiderivative_1D_mask
    end interface

    interface Characterize_nodes

        module procedure Characterize_nodes_1D, &
                         Characterize_nodes_2D, &
                         Characterize_nodes_3D
    end interface


    
    interface Inverse_Stencil
    
        module procedure Centered_Inverse_Stencil
        module procedure Upwind_Inverse_Stencil
        module procedure Inverse_Stencil_mask
    
    end interface
    
    interface Interval_bounds
    
        module procedure Interval_bounds_1D, Interval_bounds_1D_mask
    
    end interface

  type Grid

    character(len=30) :: name  = "Empty"
    integer :: N, Order, q1, q2
    real, allocatable :: nodes(:)
    integer, allocatable :: Stencil(:,:)
    integer, allocatable :: Stencil_B(:,:)
    integer, allocatable :: Stencil_F(:,:)
    real, allocatable :: Derivatives(:, :, :)
    real, allocatable :: Derivatives_B(:, :, :)
    real, allocatable :: Derivatives_F(:, :, :)
    real, allocatable :: Integral_weights(:,:)
    real, allocatable :: Antiderivative_weights(:,:)

  end type



  integer, parameter :: Nmax = 20
  type (Grid), save :: Grids(1:Nmax)
  integer, save :: ind = 0 ! index to count the number of grids

  
    contains
    
    
    
  subroutine Initialize_grid( Mesh, x, q )
        class(Grid), intent(out) :: Mesh
        real, intent(in) :: x(0:)
        integer, intent(in) :: q

        integer :: N 
        
        N = size(x)-1 ! *** Grid size   

        
        if ( allocated( Mesh % nodes ) ) then   
            deallocate( Mesh % nodes, &
                        Mesh % Stencil, &
                        Mesh % Stencil_B, &
                        Mesh % Stencil_F, &
                        Mesh % Derivatives, &
                        Mesh % Derivatives_B, &
                        Mesh % Derivatives_F, &
                        Mesh % Integral_weights, &
                        Mesh % Antiderivative_weights  )
        end if

        ! *** Grid label, order, grid size
        Mesh % Order = q
        Mesh % N = N
        Mesh % q1 = ( q - mod(q,2) )/2 
        Mesh % q2 = ( q + mod(q,2) )/2

        ! *** Grid nodes
        allocate ( Mesh % nodes(0:N) )
        Mesh % nodes = x

        ! *** Grid Stencil
        allocate ( Mesh % Stencil(0:q,0:N) )
        Mesh % Stencil = Grid_Stencil( x, q )

        ! *** Backward Grid Stencil
        allocate ( Mesh % Stencil_B(0:q,0:N) )
        Mesh % Stencil_B = Upwind_Grid_Stencil( x, q, 1 )

        ! *** Forward Grid Stencil
        allocate ( Mesh % Stencil_F(0:q,0:N) )
        Mesh % Stencil_F = Upwind_Grid_Stencil( x, q, -1 )

        ! *** Centered derivative weights
        allocate( Mesh % Derivatives(0:q,0:q,0:N) )
        Mesh % Derivatives = Derivative_weights(x, q)

        ! *** Backward derivative weights
        allocate( Mesh % Derivatives_B(0:q,0:q,0:N) )
        Mesh % Derivatives_B = Upwind_Derivative_weights(x, q, 1)
        ! *** Forward derivative weights
        allocate( Mesh % Derivatives_F(0:q,0:q,0:N) )
        Mesh % Derivatives_F = Upwind_Derivative_weights(x, q, -1)

        ! *** Weights of integrals
        allocate( Mesh % Integral_weights(0:q,0:N) )
        Mesh % Integral_weights =  &
        Integral_weights(x, q, Mesh % Derivatives)

        ! *** Weights of Antiderivatives
        allocate( Mesh % Antiderivative_weights(0:q,0:N) )
        Mesh % Antiderivative_weights = &
        Antiderivative_weights(x, q, Mesh % Derivatives)
           
   end subroutine


   function Grid_Sparse_derivative_matrix(Mesh, p) result( M )
       class(Grid), intent(in) :: Mesh
       integer, intent(in) :: p
       type(Sparse_Matrix) :: M
       
       integer :: i, i_row, ii, m_i, N, s0, sq, q 
       
       N = Mesh % N
       q = Mesh % Order
       ! *** Initialization of sparse matrix of N+1 by N+1 rows and columns
       call Initialize_sparse_matrix( M, N+1, N+1 )
       ! *** Fill the discretization matrix by rows
       do i = 0, N 
        
          i_row = i + 1
          s0 = Mesh % Stencil(0,i)
          sq = Mesh % Stencil(q,i)
          
          ! *** Number of non zero elements
          m_i = sq - s0 + 1
          
          allocate( M % rows(i_row) % i(m_i), M % rows(i_row) % v(m_i) )
          M % rows(i_row) % i = 1 + [( ii, ii=s0,sq )]
          M % rows(i_row) % v = [( Mesh % Derivatives(p,ii,i), ii=0,sq-s0 )]
       
       end do
   
   
   end function
   

   function Grid_Sparse_derivative_2D_matrix_OLD( Mesh_x, Mesh_y, direction, p ) result( M )
       class(Grid), intent(in) :: Mesh_x,  Mesh_y 
       integer, intent(in) :: direction  ! Integer to select derivative along x or y ( if 1 -> along x )
       integer, intent(in) :: p
       type(Sparse_Matrix) :: M
       
       integer :: i, j, ii, jj, i_row, m_i, Nx, Ny, s(2), N, q 
       
       Nx = Mesh_x % N
       Ny = Mesh_y % N
       
       N = ( Nx + 1 ) * ( Ny + 1 )
       ! *** Initialization of sparse matrix of N by N rows and columns
       call Initialize_sparse_matrix( M, N, N )
       ! *** Fill the discretization matrix by rows
       do j = 0, Ny
          do i = 0, Nx 
             i_row = 1 + i + ( Nx + 1 ) * j
          
             if ( direction == 1 ) then
               q =  Mesh_x % Order
               s = [ Mesh_x % Stencil(0,i), Mesh_x % Stencil(q,i)] ! Stencil            
               m_i = s(2) - s(1) + 1 !  Number of non zero elements
                
               allocate( M % rows(i_row) % i(m_i), M % rows(i_row) % v(m_i) )
               M % rows(i_row) % i = [( 1 + ii + ( Nx + 1 ) * j , ii=s(1),s(2) )] 
               M % rows(i_row) % v = [( Mesh_x % Derivatives(p,ii,i), ii=0, s(2) - s(1) )]
             elseif ( direction == 2 ) then
               q =  Mesh_y % Order
               s = [ Mesh_y % Stencil(0,j), Mesh_y % Stencil(q,j)] ! Stencil
               m_i = s(2) - s(1) + 1 !  Number of non zero elements
             
               allocate( M % rows(i_row) % i(m_i), M % rows(i_row) % v(m_i) )
               M % rows(i_row) % i = [(  1 + i + ( Nx + 1 ) * jj, jj=s(1),s(2) )]
               M % rows(i_row) % v = [( Mesh_y % Derivatives(p,jj,j), jj=0, s(2) - s(1) )]
             end if            
       
          end do
       end do
   
   
   end function

   function Grid_Sparse_derivative_2D_matrix( Mesh_x, Mesh_y, direction, p ) result( M )
       class(Grid), intent(in) :: Mesh_x,  Mesh_y 
       integer, intent(in) :: direction  ! Integer to select derivative along x or y ( if 1 -> along x )
       integer, intent(in) :: p
       type(Sparse_Matrix) :: M
       
       integer :: i, j, i_row, Nx, Ny, N 
       
       Nx = Mesh_x % N
       Ny = Mesh_y % N
       
       N = ( Nx + 1 ) * ( Ny + 1 )
       ! *** Initialization of sparse matrix of N by N rows and columns
       call Initialize_sparse_matrix( M, N, N )
       ! *** Fill the discretization matrix by rows
       do j = 0, Ny
          do i = 0, Nx 
             i_row = 1 + i + ( Nx + 1 ) * j                         
             M % rows(i_row) = Sparse_derivative_2D_row( Mesh_x, Mesh_y, i, j, direction, p )
       
          end do
       end do
   
    
   end function

   function Sparse_derivative_1D_row( Mesh, i, p ) result( D )
       class(Grid), intent(in) :: Mesh 
       integer, intent(in) :: i, p
       type(Sparse_vector) :: D
       
       integer :: N, q
       
       N = Mesh % N
       q = Mesh % Order
       
       allocate( D % i(q+1), D % v(q+1) ) 
       D % i = Mesh % Stencil(:,i) + 1
       D % v = Mesh % Derivatives(p,:,i)      
   
   end function

   function Sparse_upwind_derivative_1D_row( Mesh, i, p, sign ) result( D )
       class(Grid), intent(in) :: Mesh 
       integer, intent(in) :: i, p, sign
       type(Sparse_vector) :: D
       
       integer :: N, q
       
       N = Mesh % N
       q =  Mesh % Order
       
       allocate( D % i(q+1), D % v(q+1) ) 
       if ( sign > 0 ) then
         D % i =  Mesh % Stencil_B(:,i) + 1
         D % v = Mesh % Derivatives_B(p,:,i)   
       elseif ( sign < 0 ) then
         D % i =  Mesh % Stencil_F(:,i) + 1
         D % v = Mesh % Derivatives_F(p,:,i)        
       else   
         D = Sparse_derivative_1D_row( Mesh, i, p )
       end if
   end function

   function Sparse_derivative_2D_row( Mesh_x, Mesh_y, i, j, direction, p ) result( D )
       class(Grid), intent(in) :: Mesh_x,  Mesh_y 
       integer, intent(in) :: i, j, direction  ! Integer to select derivative along x or y ( if 1 -> along x )
       integer, intent(in) :: p
       type(Sparse_vector) :: D
       
       integer :: ii, jj, m_i, Nx, Ny, s(2), q 
       
       Nx = Mesh_x % N
       Ny = Mesh_y % N       
     
       if ( direction == 1 ) then
         q =  Mesh_x % Order
         s = [ Mesh_x % Stencil(0,i), Mesh_x % Stencil(q,i)] ! Stencil            
         m_i = s(2) - s(1) + 1 !  Number of non zero elements
          
         allocate( D % i(m_i), D % v(m_i) )
         D % i = [( 1 + ii + ( Nx + 1 ) * j , ii=s(1),s(2) )] 
         D % v = [( Mesh_x % Derivatives(p,ii,i), ii=0, s(2) - s(1) )]
       elseif ( direction == 2 ) then
         q =  Mesh_y % Order
         s = [ Mesh_y % Stencil(0,j), Mesh_y % Stencil(q,j)] ! Stencil
         m_i = s(2) - s(1) + 1 !  Number of non zero elements
       
         allocate( D % i(m_i), D % v(m_i) )
         D % i = [(  1 + i + ( Nx + 1 ) * jj, jj=s(1),s(2) )]
         D % v = [( Mesh_y % Derivatives(p,jj,j), jj=0, s(2) - s(1) )]
       end if            
       
   end function

   function Sparse_derivative_2D_column( Mesh_x, Mesh_y, i, j, direction, p ) result( D )
       class(Grid), intent(in) :: Mesh_x,  Mesh_y 
       integer, intent(in) :: i, j, direction  ! Integer to select derivative along x or y ( if 1 -> along x )
       integer, intent(in) :: p
       type(Sparse_vector) :: D
       
       integer :: ii, jj, k, m_i, Nx, Ny, s(2), q 
       
       Nx = Mesh_x % N
       Ny = Mesh_y % N
       
       if ( direction == 1 ) then
         q =  Mesh_x % Order
         m_i = count( Mesh_x % Stencil(0,:) <= i .and. i <= Mesh_x % Stencil(q,:) )
         allocate( D % i(m_i), D % v(m_i) )        
         k = 1 
         do ii = 0, Nx
            s = [ Mesh_x % Stencil(0,ii), Mesh_x % Stencil(q,ii) ]
            if ( s(1) <= i .and. i <= s(2) ) then
                D % i(k) = 1 + ii + j*(Nx+1)
                D % v(k) = Mesh_x % Derivatives(p,i-s(1),ii)
                k = k +1 
            end if
         
         end do 
       elseif ( direction == 2 ) then
         q =  Mesh_y % Order
         
         m_i = count( Mesh_y % Stencil(0,:) <= j .and. j <= Mesh_y % Stencil(q,:) )
         allocate( D % i(m_i), D % v(m_i) )        
         k = 1 
         do jj = 0, Ny
            s = [ Mesh_y % Stencil(0,jj), Mesh_y % Stencil(q,jj) ]
            if ( s(1) <= j .and. j <= s(2) ) then
                D % i(k) = 1 + i + jj*(Nx+1)
                D % v(k) = Mesh_y % Derivatives(p,j-s(1),jj)
                k = k +1 
            end if
         
         end do 
         
       end if            
     
   
   
   end function

  ! *** Function to obtain the value of the centered derivative at the i-th grid point
  ! when the input function is a Delta Kronecker centered at the j-th grid point. 
  real function Derivative_ij( Mesh, i, j, p ) result(D)
    class(Grid), intent(in) :: Mesh
    integer, intent(in) :: i, j, p 
    
    integer :: q 
    
    q = Mesh % Order
    if( Mesh % Stencil(0,i) <= j .and. j <= Mesh % Stencil(q,i) ) &
        D = Mesh % Derivatives(p, j-Mesh % Stencil(0,i), i )    
    
  end function
 
 
  ! *** Function to obtain the value of the backward derivative at the i-th grid point
  ! when the input function is a Delta Kronecker centered at the j-th grid point. 
  real function Derivative_B_ij( Mesh, i, j, p ) result(D)
    class(Grid), intent(in) :: Mesh
    integer, intent(in) :: i, j, p 
    
    integer :: q 
    
    q = Mesh % Order
    if( Mesh % Stencil_B(0,i) <= j .and. j <= Mesh % Stencil_B(q,i) ) &
        D = Mesh % Derivatives_B(p, j-Mesh % Stencil_B(0,i), i )    
    
  end function
 
 
  ! *** Function to obtain the value of the backward derivative at the i-th grid point
  ! when the input function is a Delta Kronecker centered at the j-th grid point. 
  real function Derivative_F_ij( Mesh, i, j, p ) result(D)
    class(Grid), intent(in) :: Mesh
    integer, intent(in) :: i, j, p 
    
    integer :: q 
    
    q = Mesh % Order
    if( Mesh % Stencil_F(0,i) <= j .and. j <= Mesh % Stencil_F(q,i) ) &
        D = Mesh % Derivatives_F(p, j-Mesh % Stencil_F(0,i), i )    
    
  end function
  

    ! **********************************************************
    !      Chebyshev points: Gives the Chebyshev nodal points
    !      in an interval [a,b] in C(0:N)
    !
    ! ***********************************************************

    function Chebyshev_points( a, b, N ) result(C)
         real, intent(in) :: a, b
         integer, intent(in) :: N
         real :: C(0:N)

         real :: pi = acos(-1d0)
         integer :: i

         C = [( a + (b-a)*( 1-cos(pi*i/N) )/2 ,i=0,N)]

    end function

    ! *** Function to locate "string" in a vector of characters
    !     names(:) if it finds it gives back its position. Else
    !     it gives 0 as output
    integer function Index_loc( names, string )
        character(len=*), intent(in) :: names(:), string
        integer :: i

        Index_loc = 0
        do i = 1, size(names)
            if (trim(names(i)) == trim(string))then
                Index_loc = i
                exit
            end if
        end do

    end function


   subroutine Grid_initialization( name, x, q )
        character(len=*), intent(in) :: name
        real, intent(in) :: x(0:)
        integer, intent(in) :: q

        integer :: N, g
        
        N = size(x)-1 ! *** Grid size   

        ! *** Check if the grid is already initialized
        g = 0
        g = Index_loc( Grids(:) % name, name )
        
        if ( g == 0 ) then
            ind = ind + 1 ! *** Grid counter
            g = ind 
        else    
            deallocate( Grids(g) % nodes, &
                        Grids(g) % Stencil, &
                        Grids(g) % Stencil_B, &
                        Grids(g) % Stencil_F, &
                        Grids(g) % Derivatives, &
                        Grids(g) % Derivatives_B, &
                        Grids(g) % Derivatives_F, &
                        Grids(g) % Integral_weights, &
                        Grids(g) % Antiderivative_weights  )
        end if

        ! *** Grid label, order, grid size
        Grids(g) % name = name
        Grids(g) % Order = q
        Grids(g) % N = N
        Grids(g) % q1 = ( q - mod(q,2) )/2 
        Grids(g) % q2 = ( q + mod(q,2) )/2

        ! *** Grid nodes
        allocate ( Grids(g) % nodes(0:N) )
        Grids(g) % nodes = x

        ! *** Grid Stencil
        allocate ( Grids(g) % Stencil(0:q,0:N) )
        Grids(g) % Stencil = Grid_Stencil( x, q )

        ! *** Backward Grid Stencil
        allocate ( Grids(g) % Stencil_B(0:q,0:N) )
        Grids(g) % Stencil_B = Upwind_Grid_Stencil( x, q, 1 )

        ! *** Forward Grid Stencil
        allocate ( Grids(g) % Stencil_F(0:q,0:N) )
        Grids(g) % Stencil_F = Upwind_Grid_Stencil( x, q, -1 )

        ! *** Centered derivative weights
        allocate( Grids(g) % Derivatives(0:q,0:q,0:N) )
        Grids(g) % Derivatives = Derivative_weights(x, q)

        ! *** Backward derivative weights
        allocate( Grids(g) % Derivatives_B(0:q,0:q,0:N) )
        Grids(g) % Derivatives_B = Upwind_Derivative_weights(x, q, 1)
        ! *** Forward derivative weights
        allocate( Grids(g) % Derivatives_F(0:q,0:q,0:N) )
        Grids(g) % Derivatives_F = Upwind_Derivative_weights(x, q, -1)

        ! *** Weights of integrals
        allocate( Grids(g) % Integral_weights(0:q,0:N) )
        Grids(g) %               &
        Integral_weights = Integral_weights(x, q, Grids(g) % Derivatives)

        ! *** Weights of Antiderivatives
        allocate( Grids(g) % Antiderivative_weights(0:q,0:N) )
        Grids(g) % &
        Antiderivative_weights = &
        Antiderivative_weights(x, q, Grids(g) % Derivatives)
           
   end subroutine



! ***************************************************
!         Grid_Stencil( x, q ) 
! Computes the computational molecule or stencil of order q
! for the grid x 
! *****************************************************
pure function Grid_Stencil( x, q ) result(s)
     real, intent(in) :: x(0:)
     integer, intent(in) :: q
     integer :: s(0:q,0:size(x)-1)
     
     integer :: i, N 
     
     N = size(x) - 1
     do i = 0, N

        s(:,i) = Stencil(x,i,q)

     end do
end function

! ***************************************************
!         Upwind_Grid_Stencil( x, q ) 
! Computes the computational molecule or stencil of order q
! for the grid x 
! *****************************************************
pure function Upwind_Grid_Stencil( x, q, sign ) result(s)
     real, intent(in) :: x(0:)
     integer, intent(in) :: q, sign
     integer :: s(0:q,0:size(x)-1)
     
     integer :: i, N 
     
     N = size(x) - 1
     do i = 0, N

        s(:,i) = Stencil(x,i,q, sign)

     end do
end function

!***************************************************
!           Derivative_weights(x, q):
!Given a set of nodes x(i) for i = 0,1,..,N
!and an order of interpolation q, gives the
!weights for the computation of the first q derivatives
!
!******************************************


 pure function Derivative_weights(x, q) result(L)
       real, intent(in) :: x(0:)
       integer, intent(in) :: q
       real :: L(0:q,0:q,0:size(x)-1)

       integer :: i, N, s(0:q)

       N = size(x) - 1

       do i = 0, N

          s = Stencil(x,i,q)
          L(:,:,i)= Lagrange_derivatives( x(i), x(s) )

       end do

 end function

 pure function Upwind_Derivative_weights(x, q, sign) result(L)
       real, intent(in) :: x(0:)
       integer, intent(in) :: q, sign
       real :: L(0:q,0:q,0:size(x)-1)

       integer :: i, N, s(0:q)

       N = size(x) - 1

       do i = 0, N

          s = Stencil(x,i,q, sign)
          L(:,:,i)= Lagrange_derivatives( x(i), x(s) )

       end do

 end function

 !***************************************************
!           Integral_weights(x, q):
!Given a set of nodes x(i) for i = 0,1,..,N
!and an order of interpolation q, gives the
!weights for the computation of the integral
!
!******************************************
 pure function Integral_weights(x, q, Derivatives) result(Iw)
       real, intent(in) :: x(0:)
       integer, intent(in) :: q
       real, optional, intent(in) :: Derivatives(0:q,0:q,0:size(x)-1)
       real :: Iw(0:q,0:size(x)-1)


       real :: dL(0:q,0:q,0:size(x)-1)
       integer :: i, N, s(0:q), Ns, q1, q2

       N = size(x) - 1

       if (present(Derivatives)) then
          dL = Derivatives
       else
          dL = Derivative_weights(x, q)
       end if
       
       ! *** Weights of integrals

       Ns = (N-q)/q ! Number of stencils of order q
	   
       Iw = 0
       ! Weights between [x(0), x(Ns*q)]
       q1 = ( q - mod(q,2) )/2 ; q2 = ( q + mod(q,2) )/2
       do i = q1, Ns*q-q2, q
	   
          s = Stencil(x,i,q)
          Iw(:,i) = Lagrange_integrals( x(s(q)), x(s(0)), x(s), &
                                        dL(:,:, s(q1) ) )
       end do
	   
       ! Weights between [x(Ns*q),x(N-q)]
       s = Stencil(x, N-q-q2, q)
       Iw(:,N-q) = Lagrange_integrals( x(N-q), x(Ns*q), x(s), &
                                       dL(:,:, s(q1)) )
       ! Weights between [x(N-q),x(N)]
       s = Stencil(x, N, q)
       Iw(:,N) = Lagrange_integrals( x(N), x(N-q), x(s), &
                                      dL(:,:, s(q1)) )
 end function
 !***************************************************
!           Antiderivative_weights(x, q):
!Given a set of nodes x(i) for i = 0,1,..,N
!and an order of interpolation q, gives the
!weights for the computation of the antiderivative
!
!******************************************
 pure function Antiderivative_weights(x, q, Derivatives) result(Iw)
       real, intent(in) :: x(0:)
       integer, intent(in) :: q
       real, optional, intent(in) :: Derivatives(0:q,0:q,0:size(x)-1)
       real :: Iw(0:q,0:size(x)-1)

       real :: dL(0:q,0:q,0:size(x)-1)
       integer :: i, k, N, s(0:q), Ns, q1, q2

       if (present(Derivatives)) then
          dL = Derivatives
       else
          dL = Derivative_weights(x, q)
       end if
       
       N = size(x) - 1
       Ns = (N-q)/q ! Number of stencils of order q from x(0) to x(N-q)
       q1 = ( q - mod(q,2) )/2
       q2 = ( q + mod(q,2) )/2

       ! *** Weights of integrals
       ! Weights  between [x(0),x(Ns*q)]
       Iw(:,0) = 0
       do k = q1, Ns*q-q2, q
         s = Stencil(x,k,q)
         do i = s(1), s(q)
            Iw(:,i) = Lagrange_integrals( x(i), x(s(0)), x(s), dL(:,:, s(q1)) )
         end do

       end do

       ! Weights between [x(Ns*q),x(N-q)]
       s = Stencil(x, N-q-q2, q)
       do i = Ns*q+1, N-q
          Iw(:,i) = Lagrange_integrals( x(i), x(Ns*q), x(s), dL(:,:, s(q1)) )
       end do

       ! Weights between [x(N-q+1),x(N)]
       s = Stencil(x,N,q)
       do i = s(1), s(q)
          Iw(:,i) =  &
               Lagrange_integrals( x(i), x(s(0)), x(s), dL(:,:, s(q1)) )
       end do

 end function

    ! ************************************************
    !         Derivative_1D_point(name, x, u, p, i)
    !  Given the label stored in name of an initialized
    !  grid x(0:) and its nodal images u(0:) computes
    !  the derivative of order p at the point x(i).
    ! ************************************************

    real function Derivative_1D_point(name, u, p, i) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p , i       ! Derivative order

        integer :: g, q ! Grid number and Order
        integer, allocatable :: s(:) ! Grid number and Stencil
        
        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        if( maxval( abs( u( Grids(g) % Stencil(:, i) ) ) ) > 0 ) then
           
           ! *** Stencil at x(i).
           q = Grids(g) % Order
           allocate( s(0:q) )
           s =  Grids(g) % Stencil(:, i)
           
           ! *** Derivative as linear combination of Derivatives weights
           D = dot_product( Grids(g) % Derivatives(p,:,i), u(s) )
           deallocate(s)
        else
           D = 0
        end if 
        
    end function

    real function Backward_Derivative_1D_point(name, u, p, i) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p , i       ! Derivative order

        integer :: g, q ! Grid number and Order
        integer, allocatable :: s(:) ! Grid number and Stencil
        
        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        if( maxval( abs( u( Grids(g) % Stencil_B(:, i) ) ) ) > 0 ) then
           
           ! *** Stencil at x(i).
           q = Grids(g) % Order
           allocate( s(0:q) )
           s =  Grids(g) % Stencil_B(:, i)
           
           ! *** Derivative as linear combination of Derivatives weights
           D = dot_product( Grids(g) % Derivatives_B(p,:,i), u(s) )
           deallocate(s)
        else
           D = 0
        end if 
        
    end function

    real function Forward_Derivative_1D_point(name, u, p, i) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p , i       ! Derivative order

        integer :: g, q ! Grid number and Order
        integer, allocatable :: s(:) ! Grid number and Stencil
        
        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        if( maxval( abs( u( Grids(g) % Stencil_F(:, i) ) ) ) > 0 ) then
           
           ! *** Stencil at x(i).
           q = Grids(g) % Order
           allocate( s(0:q) )
           s =  Grids(g) % Stencil_F(:, i)
           
           ! *** Derivative as linear combination of Derivatives weights
           D = dot_product( Grids(g) % Derivatives_F(p,:,i), u(s) )
           deallocate(s)
        else
           D = 0
        end if 
        
    end function

    ! ************************************************
    !         Derivative_1D_point_Delta_j(name, j, p, i)
    !  Given the label stored in name of an initialized
    !  grid x(0:), computes the derivative of order p at the point x(i)
    !  when u = Delta_{jk} e_k is a Delta Kronecker on position j.
    ! ************************************************
    real function Derivative_1D_point_Delta_j(name, j, p, i) result(D)
        character(len=*), intent(in) :: name
        integer, intent(in) :: j, p, i       ! Derivative order

        integer :: g, q ! Grid number and Order
        integer, allocatable :: s(:) ! Grid number and Stencil
        
        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        ! *** Stencil at x(i).
        q = Grids(g) % Order
        allocate( s(0:q) )
        
        s =  Grids(g) % Stencil(:, i)
        if( s(0) <= j .and. j <= s(q) ) then
          ! If j belongs to the stencil of i then D /= 0
          D = Grids(g) % Derivatives(p,j-s(0),i)
        else          
          ! If j does not belong to the stencil of i then D = 0
          D = 0
        end if
        
        deallocate(s)
    end function
    
    
    ! ******************************************************************
    !  
    function Derivative_1D_Delta_j( name, u, p, j ) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p, j ! Derivative order and Delta index
        real :: D(0:size(u)-1)
        
        integer :: g, q, i! Grid number and Order
        integer :: s(2) ! Grid number and Stencil
        
        ! *** Grid selection and grid size
        g = Index_loc( Grids(:) % name, trim(name) )
        q = Grids(g) % Order
        
        ! *** Derivative image when u is a Delta Kronecker at j
        D = 0
        s = Stencil_Differential_operator( Grids(g) % nodes, j, q ) - 1 
        D(s(1):s(2)) = [( Derivative_Delta( name, j, p, i), i=s(1),s(2) )]
            
    end function

    ! ************************************************
    !         Derivative_1D(name, x, u, p)
    !  Given the label stored in name of an initialized
    !  grid x(0:) and its nodal images u(0:) computes
    !  the derivative of order p.
    ! ************************************************

    function Derivative_1D(name, u, p, periodicity) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p        ! Derivative order
        character(len=*), optional, intent(in) :: periodicity
        real :: D(0:size(u)-1)

        integer :: i, j, q, N, g, q1, q2
        integer, allocatable :: s(:) ! Stencil
        real, allocatable :: z(:), dL(:,:)  ! nodes for periodic region

        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        if ( g == 0 ) then
            stop
        end if

        ! *** Order and size of the grid.
        q = Grids(g) % Order
        N = Grids(g) % N

        if( q < p ) then
            stop
        end if

        allocate ( s(0:q) )

        ! *** Derivatives as a linear combination of
        !     Lagrange derivatives.
        do i = 0, N
           s = Stencil( Grids(g) % nodes, i, q )
           ! *** Avoid operations with zeros to speed up
           if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
              D(i) = 0
           else
              D(i) = dot_product( Grids(g) % Derivatives(p,:,i), u(s) )
           end if
        end do

        ! *** Derivatives at the periodic region
        if( present( periodicity ).and. &
                   ( periodicity == "periodic") ) then
            q1 = ( q - mod(q,2) )/2
            q2 = ( q + mod(q,2) )/2
        allocate( z(0:q), dL(0:q,0:q)  )
            do i = 1, q-1

                ! *** Periodic stencil for
                !     N-q2+i, ... ,N-1,0,...,i
                s = mod( [(j,j=N-q+i,N+i)] , N)

                ! *** Periodic grid for derivatives
                z = Periodic_grid( Grids(g) % nodes, s, i)

                ! *** Derivative at x(j)
                j = mod(N-q2+i,N)
                if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
                    D(j) = 0
                else
                    dL = Lagrange_derivatives( z(q1), z )
                    D(j) = dot_product( dL(p,:), u(s) )
                end if

                ! *** Derivative at x(N) == x(0)
                if ( j == 0 ) then
                   s(q-i) = N
                   if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
                       D(N) = 0
                   else
                       D(N) = dot_product( dL(p,:), u(s) )
                   end if

                end if

            end do
        deallocate(z, dL)

        end if
        deallocate(s)
        
    end function

    ! ***************************************
    !        Periodic_grid( x,s, i)
    !  Gives the nodal points z(0:q) associated
    !  to the periodic stencil at positions 0 and N
    ! ************************************************
    function Periodic_grid( x,s, i) result(z)
         real, intent(in) :: x(0:)
         integer, intent(in) :: s(0:), i
         real :: z(0:size(s)-1)

         integer :: j, q, N

         q = size(s) - 1
         N = size(x) - 1
         z(0)=0
         do j = 1, q-i-1
             z(j) = x(s(j)) - x(s(0))
         end do
         z(q-i) = z(q-i-1) + x(N) - x(N-1)
         do j = q-i+1, q
             z(j) = z(q-i) + x(s(j)) - x(0)
         end do

    end function



    ! ************************************************
    !         Derivative_1D_mask(name, x, u, mask)
    !  Given the label stored in name of an initialized
    !  grid x(0:) and its nodal images u(0:) computes
    !  the derivative of order p on the regions of x(0:)
    !  in which the logical mask(0:) is .true. .
    ! ************************************************

    function Derivative_1D_mask(name, u, p, mask) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p        ! Derivative order
        logical, intent(in) :: mask(0:size(u)-1)
        real :: D(0:size(u)-1)

        real :: x(0:size(u)-1)
        integer :: i, j, q, q1, q2, N, g ! g = index to select the grid
        integer :: N_lower, N_upper, N_j  ! number of bounds and size of regions
        logical :: upper_bound(0:size(u)-1), lower_bound(0:size(u)-1)
        integer, allocatable :: s(:), i_lower(:), i_upper(:) ! Stencil
        real, allocatable :: dL(:,:)  ! Derivatives at each point

        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        if ( g == 0 )  stop

        ! *** Order, size and nodes of the grid.
        q = Grids(g) % Order  ;  N = Grids(g) % N ; x = Grids(g) % nodes

        if( q < p ) stop

        ! *** Characterization of grid nodes in lower and upper bounds

        call Characterize_nodes(x, mask, lower_bound, upper_bound)

        ! *** Count lower and upper bounds and determine their indices
        N_lower = count(lower_bound) ; N_upper = count(upper_bound)

        allocate( i_lower(N_lower), i_upper(N_upper) )
        
        i_lower = pack( [(i, i =0,N)] , lower_bound )
        i_upper = pack( [(i, i =0,N)] , upper_bound )

        ! *** Derivatives as a linear combination of
        !     Lagrange derivatives in the allowed regions.
        allocate ( s(0:q),  dL(0:q,0:q) )
        D = 0
        do j = 1, N_lower

           N_j = i_upper(j) - i_lower(j)

           if( (N_j < q) .and. (0 < N_j) ) then ! For regions with less than q+1 points

              s(0:N_j) = [( i_lower(j) + i, i=0, N_j )]
              do i = s(0), s(N_j)
                 dL = 0
                 dL(0:N_j,0:N_j) = Lagrange_derivatives( x(i), x( s(0:N_j) ) )
                 D(i) = dot_product(dL(p,0:N_j), u(s(0:N_j)))
              end do

           elseif( (N_j >= q) .and. (0 < N_j) ) then

              q1 = ( q - mod(q,2))/2 ; q2 = ( q + mod(q,2))/2

              s = [( i_lower(j) + i , i=0,q )]
              do i = i_lower(j), i_lower(j) + q1
                 dL = Lagrange_derivatives( x(i), x(s) )
                 D(i) = dot_product(dL(p,:), u(s))
              end do

              do i = i_lower(j) + q1 , i_upper(j) - q2
                 s = Stencil(x,i,q)
                 dL = Grids(g) % Derivatives(:,:,i)
                 D(i) = dot_product(dL(p,:), u(s))
              end do

              s = [( i_upper(j) - q + i , i=0,q )]
              do i = i_upper(j) - q2 + 1, i_upper(j)
                 dL = Lagrange_derivatives( x(i), x(s) )
                 D(i) = dot_product(dL(p,:), u(s))
              end do

           end if

        end do

        deallocate(s, dL, i_lower, i_upper )

    end function

    real function Derivative_1D_mask_point(name, u, p, mask, i) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p        ! Derivative order
        logical, intent(in) :: mask(0:size(u)-1)
        integer, intent(in) :: i         ! Point of evaluation of the derivative
        
        real :: x(0:size(u)-1)
        integer :: j, k, q, q1, q2, N, g ! g = index to select the grid
        integer :: N_lower, N_upper, N_j  ! number of bounds and size of regions
        logical :: upper_bound(0:size(u)-1), lower_bound(0:size(u)-1)
        integer, allocatable :: s(:), i_lower(:), i_upper(:) ! Stencil
        real, allocatable :: dL(:,:)  ! Derivatives at each point
        
        if ( mask(i) ) then
        
           ! *** Grid selection
           g = Index_loc( Grids(:) % name, trim(name) )
		   
           ! *** Order, size and nodes of the grid.
           q = Grids(g) % Order  ;  N = Grids(g) % N ; x = Grids(g) % nodes
		   
           ! *** Characterization of grid nodes in lower and upper bounds
           call Characterize_nodes(x, mask, lower_bound, upper_bound)
		   
           ! *** Count lower and upper bounds and determine their indices
           N_lower = count(lower_bound) ; N_upper = count(upper_bound)
		   
           allocate( i_lower(N_lower), i_upper(N_upper) )
           i_lower = pack( [(k, k=0,N)] , lower_bound )
           i_upper = pack( [(k, k=0,N)] , upper_bound )
		   
           ! *** Derivatives as a linear combination of
           !     Lagrange derivatives in the allowed regions.
           allocate ( s(0:q),  dL(0:q,0:q) )
           do j = 1, N_lower
              if ( i_lower(j) <= i .and. i <= i_upper(j) ) then
                  N_j = i_upper(j) - i_lower(j)   ! Length of the region
		   
                  if( N_j < q ) then 
		             ! *** For regions with less than q+1 points
                     s(0:N_j) = [( i_lower(j) + k, k=0, N_j )]
                     dL = 0
                     dL(0:N_j,0:N_j) = Lagrange_derivatives( x(i), x( s(0:N_j) ) )
                     D = dot_product( dL(p,0:N_j), u(s(0:N_j)) )
		          
                  else
		             ! *** For regions with more than or with q+1 points
                     q1 = ( q - mod(q,2))/2 ; q2 = ( q + mod(q,2))/2
                     
                     if ( i <= i_lower(j) + q1 ) then
                     
                        s = [( i_lower(j) + k , k=0,q )]
                        dL = Lagrange_derivatives( x(i), x(s) )
                        D = dot_product( dL(p,:), u(s) )
                                                
                     elseif ( i_upper(j) - q2 + 1 <= i ) then
                     
                        s = [( i_upper(j) - q + k , k=0,q )]
                        dL = Lagrange_derivatives( x(i), x(s) )
                        D = dot_product( dL(p,:), u(s) )
                        
                     else 
                     
                         s = Grids(g) % Stencil(:,i)
                        dL = Grids(g) % Derivatives(:,:,i)
                         D = dot_product(dL(p,:), u(s))
                     end if
		          
                  end if
                  
              end if 
              
           end do
		   
           deallocate(s, dL, i_lower, i_upper )
        
        else
        
          D = 0
          
        end if

    end function

    real function Derivative_1D_indices_point(name, u, p, indices, i) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p, indices(2), i         ! Point of evaluation of the derivative
        
        real :: x(0:size(u)-1)
        integer :: i0, i1
        integer :: q, q1, q2, N, N_j, g
        integer, allocatable :: s(:) 
        real, allocatable :: dL(:,:)  ! Derivatives at each point
        
        i0 = indices(1) ; i1 = indices(2) 
        if ( i0 <= i .and. i <= i1 .and.  p <= i1 - i0 ) then
           ! *** Grid selection
           g = Index_loc( Grids(:) % name, trim(name) )
		   
           ! *** Order, size and nodes of the grid.
           q = Grids(g) % Order  ;  N = Grids(g) % N ; x = Grids(g) % nodes
		   q1 = ( q - mod(q,2) ) / 2
		   q2 = ( q + mod(q,2) ) / 2		     
		   
		   
           N_j = min( i1 - i0, q ) ! size of the interval grid points 
           allocate ( s(0:N_j),  dL(0:N_j,0:N_j) )
                      
		   s(0:N_j) = i0 + Stencil( x(i0:i1), i - i0, N_j )  ! stencil
		   
           ! *** Reusing Lagrange derivatives' weights if possible
		   if ( i0 <= i - q1 .and. i + q2 <= i1 ) then
		     dL = Grids(g) % Derivatives(:,:,i)	  		   
		   else
		     dL = Lagrange_derivatives( x(i), x(s) )
		   end if 
		   
		   ! *** Derivatives as a linear combination of
           !     Lagrange derivatives' weights
		   D = dot_product( dL(p,:), u(s(:)) ) 		   		   
           
           deallocate(s, dL )
        
        else        
          D = 0  ! if not sufficient points for p-th derivative or x(i) not in
                 ! the interval [x(i0), x(i1)]        
        end if

    end function

    real function Backward_Derivative_1D_indices_point(name, u, p, indices, i) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p, indices(2), i         ! Point of evaluation of the derivative
        
        real :: x(0:size(u)-1)
        integer :: i0, i1
        integer :: q, q1, q2, N, N_j, g
        integer, allocatable :: s(:) 
        real, allocatable :: dL(:,:)  ! Derivatives at each point
        
        i0 = indices(1) ; i1 = indices(2) 
        if ( i0 <= i .and. i <= i1 .and.  p <= i1 - i0 ) then 
           ! *** Grid selection
           g = Index_loc( Grids(:) % name, trim(name) )
		   
           ! *** Order, size and nodes of the grid.
           q = Grids(g) % Order  ;  N = Grids(g) % N ; x = Grids(g) % nodes
		   q1 = q  ;  q2 = 0  	   
		   
           N_j = min( i1 - i0, q ) ! size of the interval grid points 
           allocate ( s(0:N_j),  dL(0:N_j,0:N_j) )
           
           ! Backward stencil starting at i0
		   s = Stencil( x(i0:i1), i - i0, N_j, 1, i0 )  
		   
           ! *** Reusing Lagrange derivatives' weights if possible
		   if ( i0 <= i - q1 .and. i + q2 <= i1 ) then
		     dL = Grids(g) % Derivatives_B(:,:,i)	  		   
		   else
		     dL = Lagrange_derivatives( x(i), x(s) )
		   end if 		   
		   
		   ! *** Derivatives as a linear combination of
           !     Lagrange derivatives' weights
		   D = dot_product( dL(p,:), u(s) ) 		   		   
           
           deallocate(s, dL )
        
        else        
          D = 0  ! if not sufficient points for p-th derivative or x(i) not in
                 ! the interval [x(i0), x(i1)]        
        end if

    end function

    real function Backward_Derivative_1D_mask_point(name, u, p, mask, i) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p        ! Derivative order
        logical, intent(in) :: mask(0:size(u)-1)
        integer, intent(in) :: i         ! Point of evaluation of the derivative
        
        real :: x(0:size(u)-1)
        integer :: j, k, q, q1, q2, N, g ! g = index to select the grid
        integer :: N_lower, N_upper, N_j  ! number of bounds and size of regions
        logical :: upper_bound(0:size(u)-1), lower_bound(0:size(u)-1)
        integer, allocatable :: s(:), i_lower(:), i_upper(:) ! Stencil
        real, allocatable :: dL(:,:)  ! Derivatives at each point
        
        if ( mask(i) ) then
           D = 0
           ! *** Grid selection
           g = Index_loc( Grids(:) % name, trim(name) )
		   
           ! *** Order, size and nodes of the grid.
           q = Grids(g) % Order  ;  N = Grids(g) % N ; x = Grids(g) % nodes
		   
           ! *** Characterization of grid nodes in lower and upper bounds
           call Characterize_nodes(x, mask, lower_bound, upper_bound)
		   
           ! *** Count lower and upper bounds and determine their indices
           N_lower = count(lower_bound) ; N_upper = count(upper_bound)
		   
           allocate( i_lower(N_lower), i_upper(N_upper) )
           i_lower = pack( [(k, k=0,N)] , lower_bound )
           i_upper = pack( [(k, k=0,N)] , upper_bound )
		   
           ! *** Derivatives as a linear combination of
           !     Lagrange derivatives in the allowed regions.
           allocate ( s(0:q),  dL(0:q,0:q) )
           do j = 1, N_lower
              if ( i_lower(j) <= i .and. i <= i_upper(j) ) then
                  N_j = i_upper(j) - i_lower(j)   ! Length of the region
		   
                  if( N_j < q ) then 
		             ! *** For regions with less than q+1 points
                     s(0:N_j) = [( i_lower(j) + k, k=0, N_j )]
                     dL = 0
                     dL(0:N_j,0:N_j) = Lagrange_derivatives( x(i), x( s(0:N_j) ) )
                     D = dot_product( dL(p,0:N_j), u(s(0:N_j)) )
		          
                  else
		             ! *** For regions with more than or with q+1 points
                     q1 = q ; q2 = 0
                     
                     if ( i <= i_lower(j) + q1 ) then
                     
                        s = [( i_lower(j) + k , k=0,q )]
                        dL = Lagrange_derivatives( x(i), x(s) )
                        D = dot_product( dL(p,:), u(s) )
                                                
                     elseif ( i_upper(j) - q2 + 1 <= i ) then
                     
                        s = [( i_upper(j) - q + k , k=0,q )]
                        dL = Lagrange_derivatives( x(i), x(s) )
                        D = dot_product( dL(p,:), u(s) )
                        
                     else 
                     
                         s = Grids(g) % Stencil_B(:,i)
                        dL = Grids(g) % Derivatives_B(:,:,i)
                         D = dot_product(dL(p,:), u(s))
                     end if
		          
                  end if
                  
              end if 
              
           end do
		   
           deallocate(s, dL, i_lower, i_upper )
        
        else
        
          D = 0
          
        end if

    end function

    real function Forward_Derivative_1D_mask_point(name, u, p, mask, i) result(D)
        character(len=*), intent(in) :: name
        real, intent(in) :: u(0:)
        integer, intent(in) :: p        ! Derivative order
        logical, intent(in) :: mask(0:size(u)-1)
        integer, intent(in) :: i         ! Point of evaluation of the derivative
        
        real :: x(0:size(u)-1)
        integer :: j, k, q, q1, q2, N, g ! g = index to select the grid
        integer :: N_lower, N_upper, N_j  ! number of bounds and size of regions
        logical :: upper_bound(0:size(u)-1), lower_bound(0:size(u)-1)
        integer, allocatable :: s(:), i_lower(:), i_upper(:) ! Stencil
        real, allocatable :: dL(:,:)  ! Derivatives at each point
        
        if ( mask(i) ) then
           ! *** Grid selection
           g = Index_loc( Grids(:) % name, trim(name) )
		   
           ! *** Order, size and nodes of the grid.
           q = Grids(g) % Order  ;  N = Grids(g) % N ; x = Grids(g) % nodes
		   
           ! *** Characterization of grid nodes in lower and upper bounds
           call Characterize_nodes(x, mask, lower_bound, upper_bound)
		   
           ! *** Count lower and upper bounds and determine their indices
           N_lower = count(lower_bound) ; N_upper = count(upper_bound)
		   
           allocate( i_lower(N_lower), i_upper(N_upper) )
           i_lower = pack( [(k, k=0,N)] , lower_bound )
           i_upper = pack( [(k, k=0,N)] , upper_bound )
		   
           ! *** Derivatives as a linear combination of
           !     Lagrange derivatives in the allowed regions.
           allocate ( s(0:q),  dL(0:q,0:q) )
           do j = 1, N_lower
              if ( i_lower(j) <= i .and. i <= i_upper(j) ) then
                  N_j = i_upper(j) - i_lower(j)   ! Length of the region
		   
                  if( N_j < q ) then 
		             ! *** For regions with less than q+1 points
                     s(0:N_j) = [( i_lower(j) + k, k=0, N_j )]
                     dL = 0
                     dL(0:N_j,0:N_j) = Lagrange_derivatives( x(i), x( s(0:N_j) ) )
                     D = dot_product( dL(p,0:N_j), u(s(0:N_j)) )
		          
                  else
		             ! *** For regions with more than or with q+1 points
                     q1 = 0 ; q2 = q
                     
                     if ( i <= i_lower(j) + q1 ) then
                     
                        s = [( i_lower(j) + k , k=0,q )]
                        dL = Lagrange_derivatives( x(i), x(s) )
                        D = dot_product( dL(p,:), u(s) )
                                                
                     elseif ( i_upper(j) - q2 + 1 <= i ) then
                     
                        s = [( i_upper(j) - q + k , k=0,q )]
                        dL = Lagrange_derivatives( x(i), x(s) )
                        D = dot_product( dL(p,:), u(s) )
                        
                     else 
                     
                         s = Grids(g) % Stencil_F(:,i)
                        dL = Grids(g) % Derivatives_F(:,:,i)
                         D = dot_product(dL(p,:), u(s))
                     end if
		          
                  end if
                  
              end if 
              
           end do
		   
           deallocate(s, dL, i_lower, i_upper )
        
        else
        
          D = 0
          
        end if

    end function

    ! ***************************************************************
    !       Characterize_nodes_1D: Given a set of nodes x(0:N) and a logical
    !      mask(0:N), gives as output the logicals lower_bound(0:N) and
    !      upper_bound(0:N) which are true for the points which are lower
    !      and upper limits of the allowed subregions of x in which mask
    !      is true. 
    ! ******************************************************************
    subroutine Characterize_nodes_1D(x, mask, lower_bound, upper_bound)
        real, intent(in) :: x(0:)
        logical, intent(in) :: mask(0:size(x)-1)
        logical, intent(out):: lower_bound(0:size(x)-1), &
                               upper_bound(0:size(x)-1)
        integer :: i, N

        N = size(x)-1

        ! *** Characterization of grid nodes in lower and upper bounds
        lower_bound = .false. ; upper_bound = .false.
        
        ! *** Characterization for x(0) and x(N) edges
        lower_bound(0) = mask(0) .and. mask(1)
        upper_bound(N) = mask(N) .and. mask(N-1)

        ! *** Characterization for inner grid points 0 < i < N
        do i = 1, N-1
           lower_bound(i) = mask(i) .and. (.not. mask(i-1)) !.and. mask(i+1)
           upper_bound(i) = mask(i) .and. (.not. mask(i+1)) !.and. mask(i-1)
        end do

    end subroutine


    subroutine Characterize_nodes_2D( x, y, mask, lower_bound, upper_bound )
          real, intent(in) :: x(0:), y(0:)
          logical, intent(in) :: mask(0:size(x)-1, 0:size(y)-1)
          logical, intent(out):: lower_bound(0:size(x)-1,0:size(y)-1,2), &
                                 upper_bound(0:size(x)-1,0:size(y)-1,2)

          integer :: i, j, Nx, Ny

          Nx = size(x)-1   ; Ny = size(y)-1
          ! *** Characterize nodes along x
          do j = 0, Ny
             call Characterize_nodes( x, mask(:,j), &
                                      lower_bound(:,j,1), &
                                      upper_bound(:,j,1) )
          end do

          ! *** Characterize nodes along y
          do i = 0, Nx
             call Characterize_nodes( y, mask(i,:), &
                                      lower_bound(i,:,2), &
                                      upper_bound(i,:,2) )
          end do

    end subroutine


 subroutine Characterize_nodes_3D( x, y, z, mask, lower_bound, upper_bound )
     real, intent(in) :: x(0:), y(0:), z(0:)
     logical, intent(in) :: mask(0:size(x)-1,0:size(y)-1,0:size(z)-1)
     logical, intent(out):: lower_bound(0:size(x)-1,0:size(y)-1,0:size(z)-1,3), &
                            upper_bound(0:size(x)-1,0:size(y)-1,0:size(z)-1,3)

     integer :: i, j, k,  Nx, Ny, Nz

     Nx = size(x)-1   ; Ny = size(y)-1 ; Nz = size(z) - 1

     ! *** Characterize nodes along x and y
     do k = 0, Nz
        call Characterize_nodes( x, y, mask(:,:,k), &
                                 lower_bound(:,:,k,1:2), &
                                 upper_bound(:,:,k,1:2) )
     end do

     ! *** Characterize nodes along z

     do j = 0, Ny
        do i = 0, Nx

           call Characterize_nodes( z, mask(i,j,:), &
                                    lower_bound(i,j,:,3), &
                                    upper_bound(i,j,:,3) )
        end do
     end do

 end subroutine     

    ! *****************************************************
    !     Derivative_2D_mask(coords, name, u, p,mask)
    !   Computes the derivative of order p in the nodal points where
    !   mask = .true.
    !   along the direction given by "name" for a function
    !   dependent on the coordinates "coords" and stored on u(:,:)
    ! ***************************************************
    function Derivative_2D_mask(coords, name, u, p, mask) result(D)
        character(len=*), intent(in) :: coords(2), name
        real, intent(in) :: u(0:,0:)
        integer, intent(in) :: p        ! Derivative order
        logical, intent(in) :: mask(0:,0:)
        real :: D(0:size(u,dim=1)-1,0:size(u,dim=2)-1)

        integer :: i, j, Nx, Ny
        integer :: g ! index for the grid

        Nx = size(u,dim=1)-1
        Ny = size(u,dim=2)-1

        ! *** Grid selection
        g = Index_loc( coords, trim(name) )
        if ( g == 0 ) then
            stop
        end if

        ! *** Derivative along x
        if ( g == 1 ) then
          do j = 0, Ny
             D(:,j) = Derivative_1D_mask( name, u(:,j), p, mask(:,j) )
          end do
        end if

        ! *** Derivative along y
        if ( g == 2 ) then
          do i = 0, Nx
             D(i,:) = Derivative_1D_mask( name, u(i,:), p, mask(i,:) )
          end do
        end if

    end function


    ! *****************************************************
    !     Derivative_2D(coords, name, u, p)
    !   Computes the derivative of order p in the nodal points
    !   along the direction given by "name" for a function
    !   dependent on the coordinates "coords" and stored on u(:,:)
    ! ***************************************************
    function Derivative_2D(coords, name, u, p, periodicity) result(D)
        character(len=*), intent(in) :: coords(2), name
        real, intent(in) :: u(0:,0:)
        integer, intent(in) :: p        ! Derivative order
        character(len=*), optional, intent(in) :: periodicity
        real :: D(0:size(u,dim=1)-1,0:size(u,dim=2)-1)

        integer :: i, j, Nx, Ny
        integer :: g ! index for the grid

        Nx = size(u,dim=1)-1
        Ny = size(u,dim=2)-1  

        ! *** Grid selection
        g = Index_loc( coords, trim(name) )
        if ( g == 0 ) then
            stop
        end if

        ! *** Derivative along x
        if ( g == 1 ) then
          do j = 0, Ny
             if(present(periodicity)) then
                 D(:,j) = Derivative_1D( name, u(:,j), p, periodicity )
             else
                 D(:,j) = Derivative_1D( name, u(:,j), p )
             end if

          end do
        end if

        ! *** Derivative along y
        if ( g == 2 ) then
          do i = 0, Nx
            if(present(periodicity)) then
             D(i,:) = Derivative_1D( name, u(i,:), p, periodicity )
            else
             D(i,:) = Derivative_1D( name, u(i,:), p )
            end if
          end do
        end if

    end function



    ! *****************************************************
    !     Derivative_3D(coords, name, u, p)
    !   Computes the derivative of order p in the nodal points
    !   along the direction given by "name" for a function
    !   dependent on the coordinates "coords" and stored on u(:,:,:)
    ! ***************************************************
    function Derivative_3D(coords, name, u, p, periodicity) result(D)
        character(len=*), intent(in) :: coords(3), name
        real, intent(in) :: u(0:,0:,0:)
        integer, intent(in) :: p       ! Derivative order
        character(len=*), optional, intent(in) :: periodicity
        real :: D(0:size(u,dim=1)-1,0:size(u,dim=2)-1,0:size(u,dim=3)-1)

        integer :: i, k, Nx, Ny, Nz
        integer :: g ! index for the grid

        Nx = size(u,dim=1)-1
        Ny = size(u,dim=2)-1
        Nz = size(u,dim=3)-1

        ! *** Grid selection
        g = Index_loc( coords, trim(name) )
        if ( g == 0 ) then
            stop
        end if

        ! *** Derivative along x or y
        if ( (g==1).or.(g==2) ) then
          do k = 0, Nz
             if(present(periodicity)) then
                 D(:,:,k) = Derivative_2D( coords(1:2), name, u(:,:,k), p, &
                                           periodicity )
             else
                 D(:,:,k) = Derivative_2D( coords(1:2), name, u(:,:,k), p )
             end if

          end do
        ! *** Derivative along z
        elseif (g==3) then
          do i = 0, Nx

            if(present(periodicity)) then
                D(i,:,:) = Derivative_2D( coords(2:3), name, u(i,:,:), p, &
                                          periodicity )
            else
                D(i,:,:) = Derivative_2D( coords(2:3), name, u(i,:,:), p )
            end if

          end do
        end if

    end function


    ! *****************************************************
    !     Derivative_3D_mask(coords, name, u, p)
    !   Computes the derivative of order p in the nodal points in which
    !   mask(i,j,k) = .true. along the direction given by
    !   "name" for a function dependent on the coordinates "coords"
    !   and stored on u(:,:,:). At the nodes in which mask(i,j,k) = .false.
    !   D(i,j,k) = 0.
    ! ***************************************************
    function Derivative_3D_mask(coords, name, u, p, mask) result(D)
        character(len=*), intent(in) :: coords(3), name
        real, intent(in) :: u(0:,0:,0:)
        integer, intent(in) :: p       ! Derivative order
        logical, intent(in) :: mask(0:,0:,0:)
        real :: D(0:size(u,dim=1)-1,0:size(u,dim=2)-1,0:size(u,dim=3)-1)

        integer :: i, k, Nx, Ny, Nz
        integer :: g ! index for the grid

        Nx = size(u,dim=1)-1
        Ny = size(u,dim=2)-1
        Nz = size(u,dim=3)-1

        ! *** Grid selection
        g = Index_loc( coords, trim(name) )
        if ( g == 0 ) then
            stop
        end if

        ! *** Derivative along x or y
        if ( (g==1).or.(g==2) ) then
          do k = 0, Nz

             D(:,:,k) = Derivative_2D_mask( coords(1:2), name, u(:,:,k), p, &
                                            mask(:,:,k) )
          end do
        ! *** Derivative along z
        elseif (g==3) then
          do i = 0, Nx

            D(i,:,:) = Derivative_2D_mask( coords(2:3), name, u(i,:,:), p, &
                                           mask(i,:,:) )
          end do
        end if

    end function

real function Definite_integral ( name, x, u, bounds, formula ) result(F)
        character(len=*), intent(in) :: name
        real, intent(in) :: x(0:), u(0:)
        real, intent(in) :: bounds(2)
        character(len=*), optional, intent(in) :: formula
        
        logical :: mask(0:size(x)-1)
        real :: FF(0:size(x)-1)
        real ::  dz_1, dz_N, a, b
        real, allocatable :: Fz(:), z(:), w(:)
        integer :: Nz

        ! *** Logical condition for Integral
        if( present(formula) .and. trim(formula) == "open" ) then
            mask = ( bounds(1) < x ) .and. ( x < bounds(2) )
        else
            mask = ( bounds(1) <= x ) .and. ( x <= bounds(2) )
        end if
        Nz = count(mask) ! number of elements of the grid
                
        allocate( Fz(Nz), z(Nz), w(Nz) )
        ! *** Integral in inner region ( bounds(1), bounds(2) )
        FF = Integral_1D_mask( name, x, u, mask )
        Fz = pack( FF, mask )
         F = Fz(1)
        
        ! *** Integral at the edges for open formulas
        a = bounds(1)  ; b = bounds(2) 
          z = pack( x, mask)   ;  w = pack( u, mask) 
          dz_1 =  z(2)  - z(1) ;  dz_N =  z(Nz) - z(Nz-1)
        if ( present(formula) .and. trim(formula) == "open" &
             .and. Nz > 1 ) then
        
          F = F &
          + ( -w(1)*dz_1**2 + w(1)*(a-z(2))**2  - w(2)*(a-z(1))**2 ) &
          / ( 2*dz_1 ) &
          - ( w(Nz-1)*(b-z(Nz))**2 - w(Nz)*(b-z(Nz-1))**2 + w(Nz)*dz_N**2  )  &
          / (2*dz_N) 
          
        elseif ( present(formula) .and. trim(formula) == "open" &
                 .and.  Nz == 1) then
          F = F + w(1) * ( b - a  )
        end if
        deallocate( Fz, z, w )

end function





! *** Integrates u in the grid x using Newton Cotes quadratures of order q
 real function Definite_Integral_2( x, u, q ) result(Integral_1D)
        real, intent(in) :: x(0:), u(0:)
        integer, intent(in) :: q

        integer :: i, N, s(0:q)
        integer :: Ns, q1, q2

        real :: Iw(0:q,0:size(x)-1)

        ! *** Size of the grid.
        N = size(x) - 1
        Ns = (N-q)/q

        if ( N < q ) then
            stop
        end if

        q1 = ( q - mod(q,2) ) /2 ; q2 = ( q + mod(q,2) ) /2

        Iw = Integral_weights(x,q)
        ! *** Integral between [x(0), x(Ns*q)]
        Integral_1D =  0
        do i = q1, Ns*q-q2, q
            s = Stencil(x,i,q)
           ! *** Avoid operations with zeros to speed up
           if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
               Integral_1D = Integral_1D
           else
               Integral_1D = Integral_1D &
                 + dot_product( Iw(:,i), u(s) )
           end if
        end do

        ! *** Integral between [x(Ns*q),x(N-q)]
        s = Stencil(x, N-q-q2, q)
        ! *** Avoid operations with zeros to speed up
        if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
            Integral_1D = Integral_1D
        else
            Integral_1D = Integral_1D &
              + dot_product( Iw(:,N-q), u(s) )
        end if

        ! *** Integral between [x(N-q),x(N)]
        s = Stencil(x, N, q)
        ! *** Avoid operations with zeros to speed up
        if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
            Integral_1D = Integral_1D
        else
            Integral_1D = Integral_1D &
              + dot_product( Iw(:,N),  u(s) )
        end if

 end function
! ********************************************************
!            Integral_1D: integrates u in x, where x is
!      labeled by the character name given by a previous initialization
!
! ********************************************************
    real function Integral_1D( name, x, u )
        character(len=*), intent(in) :: name
        real, intent(in) :: x(0:), u(0:)

        integer, allocatable :: s(:)
        integer :: i, q, N, Ns, q1, q2, g

        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        
        ! *** Order and size of the grid.
        q = Grids(g) % Order
        N = size(x) - 1
        Ns = (N-q)/q


        q1 = ( q - mod(q,2) ) /2 ; q2 = ( q + mod(q,2) ) /2

        allocate( s(0:q) )

        ! *** Integral between [x(0), x(Ns*q)]
        Integral_1D =  0
        do i = q1, Ns*q-q2, q
            s = Stencil(x,i,q)
           ! *** Avoid operations with zeros to speed up
           if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
               Integral_1D = Integral_1D
           else
               Integral_1D = Integral_1D &
                 + dot_product( Grids(g) % Integral_weights(:,i), u(s) )
           end if
        end do

        ! *** Integral between [x(Ns*q),x(N-q)]
        s = Stencil(x, N-q-q2, q)
        ! *** Avoid operations with zeros to speed up
        if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
            Integral_1D = Integral_1D
        else
            Integral_1D = Integral_1D &
              + dot_product( Grids(g) % Integral_weights(:,N-q), u(s) )
        end if


        ! *** Integral between [x(N-q),x(N)]
        s = Stencil(x, N, q)
        ! *** Avoid operations with zeros to speed up
        if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
            Integral_1D = Integral_1D
        else
            Integral_1D = Integral_1D &
              + dot_product( Grids(g) % Integral_weights(:,N),  u(s) )
        end if


        deallocate( s )

    end function
 
 
! ********************************************************
!          Integral_1D_Delta_j: integrates u = delta_{ij} e_i in x, where x is
!      labelled by the character name by a previous initialization
!
! ********************************************************
    real function Integral_1D_Delta_j( name, x, j ) result(Integral)
        character(len=*), intent(in) :: name
        real, intent(in) :: x(0:)
        integer, intent(in) :: j

        integer, allocatable :: s(:)
        integer :: i, q, N, q1, q2, g, Ns
        
        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )

        ! *** Order and size of the grid and other integer quantities.
        q = Grids(g) % Order  ;  N = size(x) - 1
        q1 = Grids(g) % q1    ; q2 = Grids(g) % q2
        Ns = ( N - q ) / q
        
        allocate(s(0:q))
        
        Integral =  0
        ! *** Integral between [x(0), x(Ns*q)]
        do i = q1, Ns*q-q2, q
           s = Grids(g) % Stencil(:,i)
           if ( s(0) <= j .and. j <= s(q) ) &
           Integral = Integral + Grids(g) % Integral_weights(j-s(0),i)
        end do
        
        ! *** Integral between [x(Ns*q),x(N-q)]
        s = Grids(g) % Stencil(:, N-q-q2)
        if ( s(0) <= j .and. j <= s(q) ) &
        Integral = Integral + Grids(g) % Integral_weights( j-s(0), N-q ) 

        ! *** Integral between [x(N-q),x(N)]
        s = Grids(g) % Stencil(:, N)
        if ( s(0) <= j .and. j <= s(q) ) &
        Integral = Integral + Grids(g) % Integral_weights( j-s(0), N ) 

        deallocate( s )

    end function
 
 
    real function Integral_1D_Delta_j_OLD( name, x, j ) result(Integral)
        character(len=*), intent(in) :: name
        real, intent(in) :: x(0:)
        integer, intent(in) :: j

        integer, allocatable :: s(:)
        integer :: i, k, q, N, q1, g, ss(2), Ns
        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )

        ! *** Order of the grid.
        q = Grids(g) % Order
        N = Grids(g) % N
        q1 = ( q - mod(q,2) )/2
        Ns = N - q / q
        ! *** Inverse stencil for integrals (interval of indices affected by j)
        ss = Inverse_Stencil_Integral( x, j, q )
        
        allocate(s(0:q))
        
        Integral =  0
        k=0
        ! *** Integral as integral weights (Lagrange integrals) 
        do i = ss(1) + q1, ss(2), q
           s = Grids(g) % Stencil(:, i )
           
           ! *** Index k to choose the non zero integral weight
           if ( j < N -2*q ) k = i
           if ( N -2*q <= j .and. j <= N-q ) then
               if( i == ss(1) + q1 ) k = i
               if( i  > ss(1) + q1 ) k = N-q
           end if
           if ( j > N-q ) k = N
           Integral = Integral + Grids(g) % Integral_weights(j - s(0),k)
        end do
        
        
        deallocate( s )

    end function
 
 
 
! ***Given the values of a function u(0:N) over nodal points x(0:N) and a
!    logical condition mask(0:N). Computes the integral over the subregions
!    where mask = .true. stored over a vector of size N+1.
!    In forbidden regions (mask = .false.) it gives back the 0 value.

    function Integral_1D_mask( name, x, u, mask ) result(F)
        character(len=*), intent(in) :: name
        real, intent(in) :: x(0:), u(0:)
        logical, intent(in) :: mask(0:size(x)-1)
        real :: F(0:size(u)-1)

        logical :: lower_bound(0:size(x)-1), upper_bound(0:size(x)-1)
        real, allocatable :: dL(:)
        integer, allocatable :: i_lower(:), i_upper(:), s(:)
        integer :: i, j, k, q, N, Ns, q1, q2, g, N_lower, N_upper, N_j

        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        if ( g == 0 ) then
          stop
        end if

        ! *** Order and size of the grid.
        q = Grids(g) % Order
        N = size(x) - 1
        Ns = (N-q)/q

        if ( N < q ) then
            stop
        end if

        ! *** Characterization of grid nodes in lower and upper bounds
        call Characterize_nodes(x, mask, lower_bound, upper_bound)

        ! *** Count lower and upper bounds and determine their indices
        N_lower = count(lower_bound) ; N_upper = count(upper_bound)

        allocate( i_lower(N_lower), i_upper(N_upper), dL(0:q), s(0:q) )

        ! *** Indices for lower and upper bounds
        !j = 1  ; k = 1
        !do i = 0, N
        !   if( lower_bound(i) ) then
        !      i_lower(j) = i ; j = j + 1
        !   elseif( upper_bound(i) ) then
        !      i_upper(k) = i ; k = k + 1
        !   end if
        !end do
        
        i_lower = pack( [(i, i =0,N)] , lower_bound )
        i_upper = pack( [(i, i =0,N)] , upper_bound )

        q1 = ( q - mod(q,2) ) /2 ; q2 = ( q + mod(q,2) ) /2

        F = 0
        ! *** Loop on lower bounds
        do j = 1, N_lower

           ! *** Number of steps in the subregion
           N_j =  i_upper(j)-i_lower(j)

           ! If the number of points is lesser than the order
           if( (N_j < q) .and. (0 < N_j) ) then
             s(0:N_j) = [( i_lower(j) + k, k=0,N_j)]

             if( maxval(abs(u(s(0:N_j)))) > epsilon(0d0) ) then
               dL(0:N_j) = Lagrange_integrals( x(s(N_j)), x(s(0)), x(s(0:N_j)) )
               F(i_lower(j)) = dot_product( u(s(0:N_j)) , dL(0:N_j) )
             end if
             

           elseif( N_j >= q ) then
              Ns = (N_j - q)/q

              F(i_lower(j)) = 0
              ! *** Integral between x(i_lower(j)), x(i_lower(j) + Ns*q)
              do i = i_lower(j) + q1, i_lower(j) + Ns*q - q2, q
                 s = Stencil(x,i,q)
                 
                 if( maxval(abs(u(s))) > epsilon(0d0) ) then
                    dL = Lagrange_integrals( x(s(q)), x(s(0)), x(s), &
                                            Grids(g) % Derivatives(:,:,i) )
				    
				    F(i_lower(j)) = F(i_lower(j)) + dot_product( dL, u(s) )
                 end if
              end do

              ! *** Integral between [x(i_lower(j) + Ns*q), x(i_upper(j)-q)]
              s = Stencil(x, i_lower(j) + Ns*q + q1, q)
              if( maxval(abs(u(s))) > epsilon(0d0) ) then
                dL = Lagrange_integrals( x(i_upper(j)-q), x(i_lower(j) + Ns*q), x(s), &
                                    Grids(g) % Derivatives(:,:, s(q1)) )
			    
                F(i_lower(j)) = F(i_lower(j)) + dot_product( dL, u(s) )
              end if
              ! *** Integral between [x(i_upper(j)-q),x(i_upper(j))]
              s = Stencil(x, i_upper(j)- q2, q)
              if( maxval(abs(u(s))) > epsilon(0d0) ) then
                dL = Lagrange_integrals( x(i_upper(j)), x(i_upper(j)-q), x(s), &
                                  Grids(g) % Derivatives(:,:, s(q1)) )

                F(i_lower(j)) = F(i_lower(j)) + dot_product( dL, u(s) )
              end if
           end if

           F(i_lower(j)+1:i_upper(j)) = F(i_lower(j))

        end do

        deallocate(dL, s)


    end function




   function Integral_1D_mask_2( coords, name, u, mask ) result(F)
        character(len=*), intent(in) :: coords(2), name
        real, intent(in) :: u(0:,0:)
        logical, intent(in) :: mask(0:size(u,dim=1)-1,0:size(u,dim=2)-1)
        real :: F(0:size(u,dim=1)-1,0:size(u,dim=2)-1)

        integer :: i, j, g, Nx, Ny

        Nx = size(u,dim=1)-1 ; Ny = size(u,dim=2)-1

        ! *** Grid selection
        g = Index_loc( coords, trim(name) )
        if ( g == 0 ) then
          stop
        end if

        ! *** Integral along x
        if ( g == 1 ) then
            ! *** Grid selection
            g = Index_loc( Grids(:) % name, trim(name) )
            do j = 0, Ny
                F(:,j) = Integral_1D_mask( name, Grids(g) % nodes, &
                                           u(:,j), mask(:,j) )
            end do
         elseif ( g == 2 ) then
            ! *** Grid selection
            g = Index_loc( Grids(:) % name, trim(name) )
            do i = 0, Nx
                F(i,:) = Integral_1D_mask( name, Grids(g) % nodes, &
                                           u(i,:), mask(i,:) )
            end do
         end if

    end function



    function Integral_2D_mask( coords, u, mask ) result(F)
        character(len=*), intent(in) :: coords(2)
        real, intent(in) :: u(0:,0:)
        logical, intent(in) :: mask(0:size(u,dim=1)-1,0:size(u,dim=2)-1)
        real :: F(0:size(u,dim=1)-1,0:size(u,dim=2)-1)

        real :: Ix(0:size(u,dim=1)-1,0:size(u,dim=2)-1)
        integer :: gx, gy

        ! *** Grid selection
        gx = Index_loc( Grids(:) % name, trim(coords(1)) )
        gy = Index_loc( Grids(:) % name, trim(coords(2)) )

        if ( (gx == 0).or.(gy == 0) ) then
            stop
        end if

        ! *** Integral along x
        Ix = Integral_1D_mask_2(coords, coords(1), u, mask)

        ! *** Integral along y
        F = Integral_1D_mask_2(coords, coords(2), Ix, mask)

    end function

 ! ***************************************************
  real function Integral_2D( coords,  u )
        character(len=*), intent(in) :: coords(2)
        real, intent(in) :: u(0:,0:)

        integer :: gx, gy, Nx, Ny
        integer :: i
        real :: Iy(0:size(u,dim=1)-1)

        ! *** Grid selection
        gx = Index_loc( Grids(:) % name, trim(coords(1)) )
        gy = Index_loc( Grids(:) % name, trim(coords(2)) )

        if ( (gx == 0).or.(gy == 0) ) then
            stop
        end if

        Nx = Grids(gx) % N ; Ny = Grids(gy) % N

        ! *** Integral along y
        do i = 0, Nx
            Iy(i) = Integral_1D( coords(2), Grids(gy) % nodes, u(i,:)  )
        end do

        ! *** Integral along x
        Integral_2D = Integral_1D( coords(1), Grids(gx) % nodes, Iy  )

  end function


 ! ***************************************************
  real function Integral_2D_Delta_ij( coords,  ij ) result(Integral)
        character(len=*), intent(in) :: coords(2)
        integer, intent(in) :: ij(2)

        integer :: gx, gy
        real :: Iy
               
        ! *** Grid selection
        gx = Index_loc( Grids(:) % name, trim(coords(1)) )
        gy = Index_loc( Grids(:) % name, trim(coords(2)) )

        ! *** Integral along y for Delta Kronecker at ij(2)
        Iy = Integral_1D_Delta_j( coords(2), Grids(gy) % nodes, ij(2)  )
        
        ! *** Integral along x for Delta Kronecker at ij(1) 
        Integral = Iy * Integral_1D_Delta_j( coords(1), Grids(gx) % nodes, ij(1)  )
  end function


  ! *******************************************************************
  !        Antiderivative_1D_mask
  ! *******************************************************************
   function Antiderivative_1D_mask( name, x, u, mask ) result( F )
        character(len=*), intent(in) :: name
        real, intent(in) :: x(0:), u(0:)
        logical, intent(in) :: mask(0:)
        real :: F(0:size(u)-1)

        real, allocatable :: dL(:)
        logical :: lower_bound(0:size(x)-1), upper_bound(0:size(x)-1)
        integer, allocatable :: s(:), i_lower(:),i_upper(:)
        integer :: i, j, k, q, N, Ns, q1, q2, g, N_j, N_lower, N_upper

        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        if ( g == 0 ) then
            stop
        end if

        ! *** Order and size of the grid.
        q = Grids(g) % Order  ;  N = size(x) - 1

        if ( N < q ) then
           stop
        end if

        ! *** Characterization of grid nodes in lower and upper bounds
        call Characterize_nodes(x, mask, lower_bound, upper_bound)

        ! *** Count lower and upper bounds and determine their indices
        N_lower = count(lower_bound) ; N_upper = count(upper_bound)

        allocate( i_lower(N_lower), i_upper(N_upper), dL(0:q), s(0:q) )

        ! *** Loop to obtain indices of lower and upper bounds
        
        i_lower = pack( [(i, i =0,N)] , lower_bound )
        i_upper = pack( [(i, i =0,N)] , upper_bound )

        q1 = ( q - mod(q,2) ) /2    ;    q2 = ( q + mod(q,2) ) /2

        F = 0
        ! *** Loop on lower bounds
        do j = 1, N_lower

           ! *** Number of steps in the subregion
           N_j =  i_upper(j)-i_lower(j)

           ! If the number of points is lesser than the order
           if( (N_j < q) .and. (0 < N_j) ) then

             s(0:N_j) = [( i_lower(j) + k, k=0,N_j)]

             if( maxval( abs(u(s(0:N_j))) ) > epsilon(0d0) ) then
                do i = 1, N_j
                  dL(0:N_j) = Lagrange_integrals( x(s(i)), x(s(0)), x(s(0:N_j)) )
			    
                  F(s(i)) = dot_product( u(s(0:N_j)) , dL(0:N_j) )
                end do
             end if
           elseif( (N_j >= q) .and. (0 < N_j) ) then
              ! *** Number of stencils from i_lower(j) to i_upper(j) - q
              Ns = (N_j - q)/q

              ! *** Antiderivative between x(i_lower(j)), x(i_lower(j) + Ns*q)
              do k = i_lower(j) + q1, i_lower(j) + Ns*q - q2, q
                 s = Stencil(x,k,q)   ; F(s) = F(s(0))
                 if( maxval( abs(u(s)) ) > epsilon(0d0) ) then
                   do i = 1, q
                      dL = Lagrange_integrals( x(s(i)), x(s(0)), x(s), &
                                           Grids(g) % Derivatives(:,:,s(q1)) )
				   
				      F(s(i)) = F(s(i)) + dot_product( dL, u(s) )
                   end do
                 end if
              end do

              ! *** Antiderivative between [x(i_lower(j) + Ns*q), x(i_upper(j)-q)]

              s = Stencil(x, i_lower(j) + Ns*q + q1, q)
              F(s) = F(s(0))
              if( maxval( abs(u(s)) ) > epsilon(0d0) ) then
                do i = 1, q
                   dL = Lagrange_integrals( x(s(i)), x(s(0)), x(s),         &
                                     Grids(g) % Derivatives(:,:, s(q1)) )
			    
                   F(s(i)) = F(s(i)) + dot_product( dL, u(s) )
                end do
              end if
              
              ! *** Antiderivative between [x(i_upper(j)-q),x(i_upper(j))]
              s = Stencil(x, i_upper(j)- q2, q)
              F(s) = F(s(0))
              if( maxval( abs(u(s)) ) > epsilon(0d0) ) then
                do i = 1, q
                   dL = Lagrange_integrals( x(s(i)), x(s(0)), x(s), &
                                     Grids(g) % Derivatives(:,:, s(q1)) )
                   F(s(i)) = F(s(i)) + dot_product( dL, u(s) )
                end do
              end if
           end if

        end do

        deallocate( dL, s, i_lower, i_upper )

   end function

   ! *****************************************************
   !     Antiderivative_1D
   ! *****************************************************

   function Antiderivative_1D( name, x, u ) result( F )
        character(len=*), intent(in) :: name
        real, intent(in) :: x(0:), u(0:)
        real :: F(0:size(u)-1)

        real, allocatable :: dL(:)
        integer, allocatable :: s(:)
        integer :: i, q, N, Ns
        integer :: q1, q2, g, k

        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )
        if ( g == 0 ) then
            stop
        end if

        ! *** Order and size of the grid.
        q = Grids(g) % Order
        N = size(x) - 1

        if ( N < q ) then
            stop
        end if

        q1 = ( q - mod(q,2) ) /2
        q2 = ( q + mod(q,2) ) /2

        ! *** Number of stencils between [x(0),x(N-q)]
        Ns = ( N - q ) / q

        allocate( dL(0:q), s(0:q) )

        ! *** Antiderivative between [x(0),x(Ns*q)]
        F(0) = 0
        do k = q1, Ns*q-q2, q
          s = Stencil(x,k,q)
          ! *** Avoid operations with zeros to speed up
          if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then

              F(s(0):s(q)) = F(s(0))

          else

            do i = s(1), s(q)
               dL = Grids(g) % Antiderivative_weights(:,i)
               F(i) = F(s(0)) + dot_product( dL,  u(s))
            end do

          end if
        end do

       ! *** Antiderivative between [x(Ns*q),x(N-q)]
       s = Stencil(x, N-q-q2, q)
       ! *** Avoid operations with zeros to speed up
       if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
            F(Ns*q+1:N-q) = F(Ns*q)

       else
          do i = Ns*q+1, N-q
             dL = Grids(g) % Antiderivative_weights(:,i)
             F(i) = F(Ns*q) + dot_product( dL ,  u(s)  )
          end do
       end if

       ! *** Antiderivative between [x(N-q+1),x(N)]
       s = Stencil(x,N,q)
       ! *** Avoid operations with zeros to speed up
       if ( maxval( abs(u(s)) ) < epsilon(0d0) ) then
          do i = s(1), s(q)
              F(i) = F(N-q)
          end do
       else
          do i = s(1), s(q)
             dL = Grids(g) % Antiderivative_weights(:,i)
             F(i) = F(N-q) + dot_product( dL,  u(s))
          end do

       end if


        deallocate( dL, s )

   end function

   ! *****************************************************
   !     Antiderivative_2D
   ! *****************************************************

   function Antiderivative_2D ( coords, name, u) result(F)
        character(len=*), intent(in) :: coords(2), name
        real, intent(in) :: u(0:,0:)
        real :: F(0:size(u,dim=1)-1,0:size(u,dim=2)-1)

        integer :: i, j, g
        integer :: Nx, Ny

        Nx = size(u,dim=1)-1
        Ny = size(u,dim=2)-1

        ! *** Grid selection
        g = Index_loc( coords, trim(name) )
        if ( g == 0 ) then
            stop
        end if

        ! *** Antiderivative in x
        if ( g == 1 ) then
            ! *** Grid selection
            g = Index_loc( Grids(:) % name, trim(name) )
            do j = 0, Ny
                F(:,j) = Antiderivative_1D( name, Grids(g) % nodes, u(:,j) )
            end do

        ! *** Antiderivative in y
        elseif ( g == 2 ) then
            ! *** Grid selection
            g = Index_loc( Grids(:) % name, trim(name) )
            do i = 0, Nx
                F(i,:) = Antiderivative_1D( name, Grids(g) % nodes, u(i,:) )
            end do
        end if


   end function


   ! *****************************************************
   !     Antiderivative_3D
   ! *****************************************************

   function Antiderivative_3D ( coords, name, u ) result(F)
        character(len=*), intent(in) :: coords(3), name
        real, intent(in) :: u(0:,0:,0:)
        real :: F(0:size(u,dim=1)-1,0:size(u,dim=2)-1,0:size(u,dim=3)-1)

        integer :: i, k, g
        integer :: Nx, Ny, Nz

        Nx = size(u,dim=1)-1
        Ny = size(u,dim=2)-1
        Nz = size(u,dim=3)-1

        ! *** Grid selection
        g = Index_loc( coords, name )
        if ( g == 0 ) then
            stop
        end if

        ! *** Antiderivative in x or y
        if ( (g == 1).or.(g == 2) ) then

            do k = 0, Nz
                F(:,:,k) = Antiderivative_2D( coords(1:2), name, &
                                              u(:,:,k) )
            end do

        ! *** Antiderivative in z
        elseif ( g == 3 ) then

            do i = 0, Nx
                F(i,:,:) = Antiderivative_2D( coords(2:3) , name, &
                                              u(i,:,:) )
            end do

        end if

   end function







! ************************************************************
!                 Stencil_Differential_operator
!
!   Given a set of grid points x(:) and the index i, gives back the set
!  of indices s(0:q) which use the value of a function evaluated at x(i)
!  to compute derivatives in the x(s(1)), x(s(1)+1), ..., x(s(2)) 
!
!  Author: Francisco Javier Escoto Lopez 
!          (javier.escoto.lopez@gmail.com)
! ***********************************************************
 pure function Stencil_Differential_operator_OLD_SLOW(x,i,q) result(s)    
     real, intent(in)   :: x(:)
     integer,intent(in) :: i, q
     integer :: s(2)
                   
     integer :: q1, q2, first, last, i_0, i_N
     
     q1 = ( q - mod(q,2) )/2 ; q2 = ( q + mod(q,2) )/2
     i_0 = lbound(x,1)       ; i_N = ubound(x,1)
     ! First and last points that use "i" in their stencil (far from bounds)
     first = i - q2 ! first point for which i belongs to the stencil
     last  = i + q1 ! last point for which i belongs to the stencil
      
     if( first >= i_0 + q .and. last <= i_N - q ) then     
        s = [ first, last ]     
     else
        s = [ max(i_0,first-q), min(i_N, last+q)]![ i_0, i_N ]!
     end if
     
     !if( i >  i_0 + q ) s(1) = first
     !if( i <= i_0 + q ) s(1) = i_0
     !
     !if( i <  i_N - q ) s(2) = last
     !if( i >= i_N - q ) s(2) = i_N
             
end function

 pure function Stencil_Differential_operator(x,i,q) result(s)    
     real, intent(in)   :: x(:)
     integer,intent(in) :: i, q
     integer :: s(2)
                   
     integer :: q1, q2, first, last, i_0, i_N
     
     q1 = ( q - mod(q,2) )/2 ;  q2 = ( q + mod(q,2) )/2
     i_0 = lbound(x,1)       ; i_N = ubound(x,1)
     ! First and last points that use "i" in their stencil (far from bounds)
     first = i - q2 ! first point for which i belongs to the stencil
     last  = i + q1 ! last point for which i belongs to the stencil
           
     if( i >  i_0 + q ) s(1) = first
     if( i <= i_0 + q ) s(1) = i_0
     
     if( i <  i_N - q ) s(2) = last
     if( i >= i_N - q ) s(2) = i_N
             
end function

! j \in [0, N]
pure function Inverse_Stencil_Integral( x, j, q ) result(s)
     real, intent(in) :: x(0:)
     integer,intent(in) :: j, q  
     integer :: s(2)
     
     integer :: Ns, r, N
     
     N = size(x) - 1  
     Ns = ( N - q )/q  ! number of stencils in [0, N -q]
     r = j/q   ! Stencil number for j where r \in [0, Ns]\cup[ Ns, floor(N,q)]
     if ( j <= Ns * q ) then  
       
       if ( mod(j,q) == 0 ) then
         ! *** When j = r*q two stencils are affected
         s(1) = max( r * q - q, 0 )
         if( j /= Ns*q ) s(2) =  r * q + q
         if( j == Ns*q ) s(2) =  N - q 
       
       else
         ! *** When j = r*q + p ( 0 < p < q ) one stencil is affected
         s(1) =  r * q  
         if ( j <  N - 2*q ) s(2) =  r * q + q
         if ( j >= N - 2*q ) s(2) =  N - q 
         
       end if
     
     elseif ( Ns*q < j .and. j <= N - q ) then
        s(1) = N - 2*q
        if( j == N-q ) s(2) = N
        if( j <  N-q ) s(2) = N-q
        
     else
         
       ! *** When j > N-q one stencil is affected
       s(1) = N - q
       s(2) = N
     
     end if

end function


! j \in [0, N]
pure function Inverse_Stencil_Integral_OLD( x, j, q ) result(s)
     real, intent(in) :: x(0:)
     integer,intent(in) :: j, q  
     integer :: s(2)
     
     integer :: Ns, r, N
     
     N = size(x) - 1  
     Ns = ( N - q )/q  ! number of stencils in [0, N -q]
     r = j/q   ! Stencil number for j where r \in [0, Ns]\cup[ Ns, floor(N,q)]
     if ( r <= Ns ) then ! For j \in [0, N-q + c]
       
       if ( mod(j,q) == 0 ) then
         ! *** When j = r*q two stencils are affected
         s(1) = max( r * q - q, 0 )
         s(2) =  r * q + q
       
       else
         ! *** When j = r*q + p ( 0 < p < q ) one stencil is affected
         s(1) =  r * q  
         if ( j <  N - 2*q ) s(2) =  r * q + q
         if ( j >= N - 2*q ) s(2) =  N - q 
         
       end if
     
     else
         
       ! *** When j>N-q one stencil is affected
       s(1) = N - q
       s(2) = N
     
     end if

end function

 function Stencil_Differential_operator_mask(x,i,q,mask) result(s)    
     real, intent(in)   :: x(:)
     integer,intent(in) :: i, q  
     logical, intent(in)   :: mask(:)
     integer :: s(2)
                   
     integer :: N_lower, N_upper, j, N
     integer, allocatable :: j_lower(:), j_upper(:)
     logical :: lower_bound(size(x)), upper_bound(size(x))
     
     ! *** Characterization of grid nodes in lower and upper bounds
     call Characterize_nodes(x, mask, lower_bound, upper_bound)

     ! *** Count lower and upper bounds and determine their indices
     N_lower = count(lower_bound) ; N_upper = count(upper_bound)
     
     allocate( j_lower(N_lower), j_upper(N_upper) )
     
     N = size(x) 
     j_lower = pack( [(i, i=1,N)] , lower_bound )
     j_upper = pack( [(i, i=1,N)] , upper_bound )
     
     do j = 1, N_lower
        
        if ( j_lower(j) <= i .and. i <= j_upper(j) ) then
        
          s = Stencil_Differential_operator( x(j_lower(j):j_upper(j)), &
                                             i - j_lower(j) + 1, &
                                             min(q, j_upper(j) - j_lower(j)) )
          s = s + j_lower(j) - 1
        end if
     end do
     
     deallocate( j_lower, j_upper )
             
end function









 pure function Centered_Inverse_Stencil(x,i,q) result(s)    
     real, intent(in)   :: x(0:)
     integer,intent(in) :: i, q
     integer :: s(2)
                   
     integer :: q1, q2, first, last , N
     
     q1 = ( q - mod(q,2) )/2 ;  q2 = ( q + mod(q,2) )/2
     N = size(x) - 1
     ! First and last points that use "i" in their stencil (far from edges)
     first = i - q2 ! first point for which i belongs to the stencil
     last  = i + q1 !  last point for which i belongs to the stencil
           
     if( i >   q ) s(1) = first
     if( i <=  q ) s(1) = 0
     
     if( i <  N - q ) s(2) = last
     if( i >= N - q ) s(2) = N
             
 end function
 
 



 
 pure function Upwind_Inverse_Stencil(x,i,q,sign) result(s)    
     real, intent(in)   :: x(0:)
     integer,intent(in) :: i, q, sign
     integer :: s(2)
                   
     integer :: N 
     
     N = size(x) - 1
     ! *** Inverse Stencil definition
     if ( sign > 0 ) then
       s = [ i, min(N, i + q) ]
       if ( i <= q ) s(1) = 0
     elseif ( sign < 0 ) then
       s = [ max(0, i - q), i ]
       if( i >= N-q ) s(2) = N         
     elseif ( sign == 0 ) then
       s =  Centered_Inverse_Stencil(x, i, q)
     end if             
     
end function




 function Inverse_Stencil_mask( x, i, q, mask, sign ) result(s)    
     real, intent(in)   :: x(0:)
     integer,intent(in) :: i, q  
     logical, intent(in) :: mask(0:size(x)-1)
     integer, optional, intent(in) :: sign  
     integer :: s(2)
                   
     integer :: N_lower, N_upper, j, N, N_j
     integer, allocatable :: j_lower(:), j_upper(:)
     logical :: lower_bound(0:size(x)-1), upper_bound(0:size(x)-1)
     
     ! *** Characterization of grid nodes in lower and upper bounds
     call Characterize_nodes(x, mask, lower_bound, upper_bound)

     ! *** Count lower and upper bounds and determine their indices
     N_lower = count(lower_bound) ; N_upper = count(upper_bound)
     
     allocate( j_lower(N_lower), j_upper(N_upper) )
     
     N = size(x) - 1
     j_lower = pack( [(i, i=0,N)] , lower_bound )
     j_upper = pack( [(i, i=0,N)] , upper_bound )
     
     do j = 1, N_lower
        
        if ( j_lower(j) <= i .and. i <= j_upper(j) ) then
          N_j = j_upper(j) - j_lower(j)
          
          if ( present(sign) ) then
             s = Upwind_Inverse_Stencil( x(j_lower(j):j_upper(j)), &
                                         i - j_lower(j), min(q,N_j), sign )
                                            
          else
             s = Centered_Inverse_Stencil( x(j_lower(j):j_upper(j)), &
                                            i - j_lower(j), min(q,N_j) )
          
          end if
          s = s + j_lower(j) 
        end if
     end do
     
     deallocate( j_lower, j_upper )
             
end function


 

 
  ! *** Grid_bounds_1D: Given a grid x(1:N) and a logical mask(1:N)
  ! that delimites the allowed subregions of x, gives back the indices values
  ! i_lower(1:N), i_upper(1:N) which delimite each region.
  subroutine Grid_bounds_1D( x, mask, i_lower, i_upper ) 
     real, intent(in) :: x(:) 
     logical, intent(in) :: mask(size(x))
     integer, intent(out) :: i_lower( size(x) ), i_upper( size(x) )
     
     logical :: lower_bound( size(x) ), upper_bound( size(x) )
     integer :: i, j, N, N_j
     integer, allocatable :: j_lower(:), j_upper(:)
     
     N = size(x) 
     ! *** Characterization of grid nodes in lower and upper bounds 
     ! using logicals
     lower_bound = .false. ; upper_bound = .false.
     lower_bound(1) = mask(1) .and. mask(2)  
     upper_bound(N) = mask(N) .and. mask(N-1)    
     do i = 2, N-1 
        lower_bound(i) = mask(i) .and. (.not. mask(i-1)) 
        upper_bound(i) = mask(i) .and. (.not. mask(i+1)) 
     end do
     
     ! *** Indices for lower and upper bounds of size N_j
     N_j = count(lower_bound)
     
     allocate( j_lower(N_j), j_upper(N_j) ) 
     j_lower = pack( [(j, j =1,N)] , lower_bound )
     j_upper = pack( [(j, j =1,N)] , upper_bound )  
     
     ! *** Indices for lower and upper bounds of size N
     i_lower = - 1 ; i_upper = -1   ! -1 value in not allowed regions
     
     do j = 1, N_j
        do i = j_lower(j), j_upper(j)
           i_lower(i) = j_lower(j)
           i_upper(i) = j_upper(j)
        end do
     end do
        
     deallocate( j_lower, j_upper ) 
  end subroutine


 
  ! *** Interval_bounds_1D: Given a grid x(0:N) and the bounds of an interval
  ! [a, b] = [bounds(1), bounds(2)]. Gives back the minimum and maximum indices values
  ! i_lower, i_upper such that a <= x(i_lower) .and. x(i_upper) <= b.
  subroutine Interval_bounds_1D( x, bounds, i_lower, i_upper, i0 ) 
     real, intent(in) :: x(0:), bounds(2) 
     integer, intent(out) :: i_lower, i_upper
     integer, optional, intent(in) :: i0
     
     logical :: mask(0:size(x)-1)
     integer :: j, N, N_j
     integer, allocatable :: allowed(:)
     
     N = size(x) - 1
     
     ! *** Mask to define points belonging in the interval
     mask = bounds(1) <= x .and. x <= bounds(2)
     
     N_j = count(mask)  ! number of points in the interval
     
     if ( N_j > 0 ) then
       allocate( allowed(N_j) ) 
       allowed = pack( [(j, j =0,N)] , mask )
       i_lower = allowed(1)
       i_upper = allowed(N_j)     
       deallocate( allowed )
     else
       i_lower = -1
       i_upper = -2
     end if
     
     if(present(i0)) i_lower = i_lower + i0
     if(present(i0)) i_upper = i_upper + i0
      
  end subroutine


  ! *** Interval_bounds_1D: Given a grid x(0:N) and the bounds of an interval
  ! [a, b] = [bounds(1), bounds(2)]. Gives back the minimum and maximum indices values
  ! i_lower, i_upper such that a <= x(i_lower) .and. x(i_upper) <= b.
  subroutine Interval_bounds_1D_mask( x, bounds, mask, i_lower, i_upper, i0 ) 
     real, intent(in) :: x(0:), bounds(2) 
     logical, intent(in) :: mask(0:size(x)-1)
     integer, intent(out) :: i_lower, i_upper
     integer, optional, intent(in) :: i0
     
     logical :: mask_interval(0:size(x)-1)
     integer :: j, N, N_j
     integer, allocatable :: allowed(:)
     
     N = size(x) - 1
     
     ! *** Mask to define points belonging in the interval
     mask_interval = bounds(1) <= x .and. x <= bounds(2)
     mask_interval = mask_interval .and. mask
     
     N_j = count(mask_interval)  ! number of points in the interval
     
     if ( N_j > 0 ) then
       allocate( allowed(N_j) ) 
       allowed = pack( [(j, j =0,N)] , mask_interval )
       i_lower = allowed(1)
       i_upper = allowed(N_j)     
       deallocate( allowed )
     else
       i_lower = -1
       i_upper = -2
     end if
     
     if(present(i0)) i_lower = i_lower + i0
     if(present(i0)) i_upper = i_upper + i0
      
  end subroutine



  ! *** Periodic_Grid_bounds_1D: Given a periodic grid x(1:N) and a logical mask(1:N)
  ! that delimites the allowed subregions of x, gives back the indices values
  ! i_lower(1:N), i_upper(1:N) which delimite each region.
  subroutine Periodic_Grid_bounds_1D( x, mask, i_lower, i_upper ) 
     real, intent(in) :: x(:) 
     logical, intent(in) :: mask(size(x))
     integer, intent(out) :: i_lower( size(x) ), i_upper( size(x) )
      
     integer :: N
     
     N = size(x)
     
     call Grid_bounds_1D( x, mask, i_lower, i_upper ) 
     
     if( mask(N).and.mask(1) ) i_lower(1:i_upper(1)) = i_lower(N)
     if( mask(N).and.mask(1) ) i_upper(i_lower(N):N) = i_upper(1)
     
  end subroutine

  
!***************************************************
!           Periodic_Derivative_weights(x, q, sign):
! Given a set of nodes x(i) for i = 0,1,..,N
! and an order of interpolation q, gives the
! weights for the computation of the first q derivatives in periodic domains.
! The periodicity must be such that x(0) = x(N-1) + ( x(1) - x(0) ). ( no repeated points )
!  That is, if we extend x(N) := x(N-1) + ( x(1) - x(0) ), the periodicity
! is such that x(0) = x(N)
! In case that sign is present it uses upwind stencils
! for the Lagrange polynomials involved on the 
! weights computation.
!  These weights must be used with the appropriate stencil.
!******************************************

pure function Periodic_Derivative_weights(x, q, sign) result(L)
       real, intent(in) :: x(0:)
       integer, intent(in) :: q
       integer, optional, intent(in) :: sign
       real :: L(0:q,0:q,0:size(x)-2)

       integer :: i, N, s(0:q), q1, q2
       real, allocatable :: z(:)
       real :: dx(0:size(x)-2)

       N = size(x) -1    
       if( present(sign).and. sign > 0 )      q1 = q                   ! backward 
       if( present(sign).and. sign < 0 )      q1 = 0                   ! forward                                    
       if( .not.present(sign).or. sign == 0 ) q1 = ( q - mod(q,2) )/2  ! centered  
       q2 = q - q1
              
       ! *** Extended grid for periodicity
       allocate( z(-q1:N+q2) )
       dx = x(1:N-1) - x(0:N-2) 
       
       z(0:N-1) = x
       do i = -1, - q1, - 1
          z(i) = z(i+1) - dx( modulo(i,N) ) 
       end do
       do i = N-1, N + q2-1
          z(i+1) = z(i) + dx( modulo(i,N) ) 
       end do
       
       do i = 0, N -1
          ! *** Centered or upwind stencil selection
          if( present(sign) )       s = Stencil(z,i+q1,q,sign) -q1 ! upwind
          if( .not. present(sign) ) s = Stencil(z,i+q1,q)      -q1 ! centered
          
          ! *** Lagrange polynomials derivatives from stencil
          L(:,:,i)= Lagrange_derivatives( z(i), z(s) )

       end do
       
       deallocate(z) 
 end function

! *** Given a grid x(0:N-1) gives back the extended periodic grid
! z(-q1+N+q2) for which z(i+1) = z(i) + dx(modulo(i,N)).
pure function Extended_periodic_grid(x,q1,q2) result(z)
     real, intent(in) :: x(0:)
     integer, intent(in) :: q1, q2
     real :: z(-q1:size(x)+q2)
     
     real :: dx(0:size(x)-2)
     integer :: N, i
     
     N = size(x) 
     
     dx = x(1:N-1) - x(0:N-2) 
       
     z(0:N-1) = x
     do i = -1, - q1, - 1
        z(i) = z(i+1) - dx( modulo(i,N) ) 
     end do
     do i = N, N + q2-1
        z(i+1) = z(i) + dx( modulo(i,N) ) 
     end do

end function


pure function Periodic_Grid_Stencil( x, q, sign ) result(s)
     real, intent(in) :: x(0:)
     integer, intent(in) :: q
     integer, optional, intent(in) :: sign
     integer :: s(0:q,0:size(x)-1)
     
     integer :: i, N, q1, q2
     real, allocatable :: z(:)
     
     N = size(x)    
     if( present(sign).and. sign > 0 )      q1 = q                   ! backward 
     if( present(sign).and. sign < 0 )      q1 = 0                   ! forward                                    
     if( .not.present(sign).or. sign == 0 ) q1 = ( q - mod(q,2) )/2  ! centered  
     q2 = q - q1
     
     ! *** Extended grid for periodicity
     allocate( z(-q1:N+q2) )
     z=0
     do i = 0, N-1
        if( present(sign) )       s(:,i) = Stencil(z,i+q1,q,sign)   -q1
        if( .not. present(sign) ) s(:,i) = Stencil(z,i+q1,q)  -q1 
        s(:,i) = modulo(  s(:,i), N )
     end do
     
     deallocate(z)
end function 





    
    ! ************************************************
    !     Derivative_FD_1D_ij(name, i, j, p) 
    !  Gives back the value of the derivative of order p at x_i when 
    ! u = delta_kronecker(x_j) centered at x_j. The derivative uses
    ! a centered stencil
    ! This is equivalent to the
    ! element (i,j) of the discretization matrix of the derivative of order
    ! p.
    ! ************************************************
    real function Derivative_FD_1D_ij(name, i, j, p) result(D)
        character(len=*), intent(in) :: name
        integer, intent(in) :: i, j, p 

        integer :: q, g
        integer, allocatable :: s(:) ! Stencil

        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )

        ! *** Order and size of the grid.
        q = Grids(g) % Order

        allocate ( s(0:q) )
        ! *** Image of the derivative
        s = Grids(g) % Stencil(:,i)
        if ( s(0) <= j .and. j <= s(q) ) then ! if j belongs to the stencil
           D = Grids(g) % Derivatives(p, j-s(0), i)
        else ! if j does not belong to the stencil
           D = 0
        end if
         
    end function
    
    ! ************************************************
    !     Derivative_FD_1D_ij_indices(name, i, j, p, indices) 
    !  Gives back the value of the derivative of order p at x_i when 
    ! u = delta_kronecker(x_j) centered at x_j and the grid nodes are
    ! delimited by indices(1) and indices(2). That is, we only use the 
    ! nodes Grids(g) % nodes(indices(1):indices(2)). 
    ! The derivative uses a centered stencil of the greatest possible 
    ! order up to Grids(g) % Order. 
    ! This is equivalent to the
    ! element (i,j) of the discretization matrix of the derivative of order
    ! p.
    ! ************************************************
    real function Derivative_FD_1D_ij_indices(name, i, j, p, indices) result(D)
        character(len=*), intent(in) :: name
        integer, intent(in) :: i, j, p, indices(2)

        integer :: q, g, q1, q2 
        integer, allocatable :: s(:) ! Stencil
        real, allocatable :: dL(:,:) ! Lagrange derivatives
        integer :: j0, j1, N_j
        logical :: allowed, in_stencil
        
        D = 0
        j0 = indices(1) ; j1 = indices(2)         
        allowed = j0 <= min(i,j) .and. max(i,j) <= j1 .and. j1-j0 >= p
        if ( allowed ) then
           ! *** Grid selection
           g = Index_loc( Grids(:) % name, trim(name) )
		   
           ! *** Order and size of the grid.
           q = Grids(g) % Order ; N_j = min( q, j1 - j0 )
           q1 = ( q - mod(q,2) )/2 ; q2 = ( q + mod(q,2) )/2		    
		           
           ! *** Centered stencil at i - i0 in the allowed region
           allocate ( s(0:N_j), dL(0:N_j,0:N_j) ) 
		   s = j0 + Stencil( Grids(g) % nodes(j0:j1), i - j0, N_j ) 
           
           in_stencil = s(0) <= j .and. j <= s(N_j)
           if ( in_stencil ) then
             ! *** Reusing Lagrange derivatives' weights if possible
		     if ( j0 <= i - q1 .and. i + q2 <= j1 ) then
		       dL = Grids(g) % Derivatives(:,:,i)	  		   
		     else
		       dL = Lagrange_derivatives( Grids(g) % nodes(i), Grids(g) % nodes(s) )
		     end if 			
		     ! *** If both in allowed region and stencil     
		     D = dL(p,j-s(0))	
           end if 
                     
        end if    
           
    end function
    
 
 
   
    ! ************************************************
    !      Derivative_FD_1D_upwind_ij(name, i, j, p, sign) 
    !  Gives back the value of the derivative of order p at x_i when 
    ! u = delta_kronecker(x_j) centered at x_j. The derivative uses
    ! an upwind stencil
    ! This is equivalent to the
    ! element (i,j) of the discretization matrix of the derivative of order
    ! p.
    ! ************************************************
    real function Derivative_FD_1D_upwind_ij(name, i, j, p, sign) result(D)
        character(len=*), intent(in) :: name
        integer, intent(in) :: i, j, p, sign 

        integer :: q, g
        integer, allocatable :: s(:) ! Stencil
        real, allocatable :: dL(:) ! Derivative weights

        ! *** Grid selection
        g = Index_loc( Grids(:) % name, trim(name) )

        ! *** Order and size of the grid.
        q = Grids(g) % Order

        ! Stencil for upwind scheme
        allocate ( s(0:q) )
        if( sign > 0 )  s = Grids(g) % Stencil_B(:,i)
        if( sign < 0 )  s = Grids(g) % Stencil_F(:,i)
        
        ! Weights for upwind scheme
        allocate ( dL(0:q) )
        if( sign > 0 )  dL = Grids(g) % Derivatives_B(p,:,i)
        if( sign < 0 )  dL = Grids(g) % Derivatives_F(p,:,i)
        
        ! *** Image of the derivative
        if ( s(0) <= j .and. j <= s(q) ) then ! if j belongs to the stencil
           D = dL(j-s(0))
        else ! if j does not belong to the stencil
           D = 0
        end if
        
    end function
    
    
    
    ! ************************************************
    !     Derivative_FD_1D_upwind_ij_indices(name, i, j, p, indices) 
    !  Gives back the value of the derivative of order p at x_i when 
    ! u = delta_kronecker(x_j) centered at x_j and the grid nodes are
    ! delimited by indices(1) and indices(2). That is, we only use the 
    ! nodes Grids(g) % nodes(indices(1):indices(2)). 
    ! The derivative uses a centered stencil of the greatest possible 
    ! order up to Grids(g) % Order. 
    ! This is equivalent to the
    ! element (i,j) of the discretization matrix of the derivative of order
    ! p.
    ! ************************************************
    real function Derivative_FD_1D_upwind_ij_indices(name, i, j, p, sign, indices) result(D)
        character(len=*), intent(in) :: name
        integer, intent(in) :: i, j, p, sign, indices(2)

        integer :: q, g, q1, q2 
        integer, allocatable :: s(:) ! Stencil
        real, allocatable :: dL(:,:) ! Lagrange derivatives
        integer :: j0, j1, N_j
        logical :: allowed, in_stencil
        
        D = 0
        j0 = indices(1) ; j1 = indices(2)         
        allowed = j0 <= min(i,j) .and. max(i,j) <= j1 .and. j1-j0 >= p
        if ( allowed ) then
           ! *** Grid selection
           g = Index_loc( Grids(:) % name, trim(name) )
		   
           ! *** Order and size of the grid.
           q = Grids(g) % Order ; N_j = min( q, j1 - j0 )
           
           ! *** Choosing between backward and forward
           if( sign > 0 ) q1 = q
           if( sign < 0 ) q1 = 0           
           q2 = ( q - q1 )		    		          
		                    
           ! *** Centered stencil at i - i0 in the allowed region
           allocate ( s(0:N_j), dL(0:N_j,0:N_j) ) 
		   s = Stencil( Grids(g) % nodes(j0:j1), i - j0, N_j, sign, j0 ) 
            
           in_stencil = s(0) <= j .and. j <= s(N_j)
           if ( in_stencil ) then 
             ! *** Reusing Lagrange derivatives' weights if possible
		     if ( j0 <= i - q1 .and. i + q2 <= j1 ) then
		       if( sign > 0 ) dL = Grids(g) % Derivatives_B(:,:,i)	  		   
		       if( sign < 0 ) dL = Grids(g) % Derivatives_F(:,:,i)	  		   
		     else
		       dL = Lagrange_derivatives( Grids(g) % nodes(i), Grids(g) % nodes(s) )
		     end if 			
		     ! *** If both in allowed region and stencil 
		     D = dL(p,j-s(0))	
		     
           end if 
                     
        end if    
           
    end function
    
    
    


end module

