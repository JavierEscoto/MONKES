module Sparse_Matrices
   
   implicit none
   
   
   type Sparse_vector
       real, allocatable :: v(:)     ! Non zero components of the vector
       integer, allocatable :: i(:)  ! Indices for Non zero components of the vector           
   end type 
   
   type Sparse_Matrix
      integer :: m, n        
      type(Sparse_vector), allocatable :: rows(:), columns(:)  
       contains
       
       procedure :: Initialize_sparse_matrix
   end type
   
 contains
 
 function Sparse_Identity(N) result(A)
    integer, intent(in) :: N
    type(Sparse_Matrix) :: A
    
    integer :: i
    
    call Initialize_sparse_matrix(A, N, N)
    
    do i = 1, N
       allocate( A % rows(i) % i(1), A % rows(i) % v(1) )
       A % rows(i) % i = i
       A % rows(i) % v = 1    
    end do
 
 end function
 
 function Sparse_Matmul( A, x ) result(y)
     class(Sparse_Matrix), intent(in) :: A 
     real, intent(in) :: x( A % n ) 
     real :: y( A % m )
     
     integer :: i
     
     do i = 1, A % m 
        y(i) = dot_product( A % rows(i) % v, x(A % rows(i) % i) )        
     end do
     
 
 end function 
 
 function Sparse_Diagonal(N, D) result(A)
    integer, intent(in) :: N
    real, intent(in) :: D(N)
    type(Sparse_Matrix) :: A
    
    integer :: i
    
    call Initialize_sparse_matrix(A, N, N)
    
    do i = 1, N
       allocate( A % rows(i) % i(1), A % rows(i) % v(1) )
       A % rows(i) % i = i
       A % rows(i) % v = D(i)   
    end do
 
 end function
 
 
 ! *** Extracts the submatrix constituded by the positions
 ! (i,j)=( row_index(ii), column_index(jj) ) (ii=1,..M, jj=1,..N) and stores
 ! them in the MxN matrix B.
 subroutine Extract_submatrix( A, row_index, column_index, B )
    type(Sparse_Matrix), intent(in) :: A
    integer, intent(in) :: row_index(:), column_index(:) 
    real :: B( size(row_index), size(column_index) )
    
    integer :: i, j, k, ii, jj, M, N
    
    M = size(row_index) ; N = size(column_index)
    
    B = 0
    do ii = 1, M ; i = row_index(ii)        
       do jj = 1, N ; j = column_index(jj)  
          
          do k = 1, size( A % rows(i) % i )          
             if ( A % rows(i) % i(k) == j ) B(ii,jj) = A % rows(i) % v(k)     
          end do
       
       end do
    end do 
 
 end subroutine 
 
 subroutine Fill_Sparse_Matrix_columns( A )
     class(Sparse_Matrix), intent(inout) :: A
     
     integer :: i, ii, j, k, N_nonzero  
     
     allocate ( A % columns(A % n) )
     do j = 1, A % n
        
        ! *** Count non zero rows for column j
        N_nonzero = 0
        do i = 1, A % m
           do k = 1, size( A % rows(i) % i )                
              if( A % rows(i) % i(k) == j ) N_nonzero = N_nonzero + 1
           end do           
        end do
        
	    ! *** Extract non zero rows for column j
        if ( N_nonzero > 0 ) then
        
          allocate( A % columns(j) % i(N_nonzero), &
                    A % columns(j) % v(N_nonzero) )                    
          ii = 0 
          do i = 1, A % m             
             do k = 1, size( A % rows(i) % i )
                
                if( A % rows(i) % i(k) == j ) then
                   ii = ii + 1 
                   A % columns(j) % i(ii) = i 
                   A % columns(j) % v(ii) = A % rows(i) % v(k)                 
                end if
             
             end do             
          end do
          
        else
          ! *** If its a zero column vector
          allocate( A % columns(j) % i(1), A % columns(j) % v(1) )
          A % columns(j) % i(1) = j ; A % columns(j) % v(1) = 0                 
        end if
     
     end do
 
 end subroutine 
 
 
 function Allocate_sparse_vector(N) result(Sp_vector)
    integer, intent(in) :: N
    type(Sparse_vector) :: Sp_vector
    
    allocate( Sp_vector % i(N) )
    allocate( Sp_vector % v(N) )
     
    Sp_vector % i = 0
    Sp_vector % v = 0
 end function
 
 real function Sp_dot_product(x, y) result(C)
    class(Sparse_vector), intent(in) :: x, y
    
    integer :: i, j
    
    C = 0
    do j = 1, size( y % i ) 
       do i = 1, size( x % i ) 
         
         if( x % i(i) == y % i(j) ) C = C + x % v(i) * y % v(j)
       
       end do
    end do
    
 end function
 
 function Sum_sparse_vector(x, y) result(z)
    class(Sparse_vector), intent(in) :: x, y
    type(Sparse_vector) :: z
    
    integer :: m_x, m_y, m_z, i, j, k
    integer, allocatable :: i_temp(:) 
    logical, allocatable :: repeated(:) 
    
    m_x = size( x % i ) 
    m_y = size( y % i ) 
    
    ! *** Temporal array of indices
    allocate( i_temp(m_x + m_y), repeated(m_y) ) 
    i_temp = [ x % i, -1 ]
    repeated = .false.
    do i = 1, m_y       
       ! *** Assign non repeated indices
       if ( count( x % i == y % i(i) ) == 0 ) i_temp(m_x+i)  = y % i(i)
       if ( count( x % i == y % i(i) ) /= 0 ) repeated(i) = .true.    
    end do
    
    ! *** Dimension of the sparse sum vector z = x + y
    m_z = m_x + count( .not. repeated )
    
    allocate( z % i(m_z),  z % v(m_z) )     
    
    z % v(1:m_x) = x % v
    z % i(1:m_x) = x % i
    
    j = 1
    do i = 1, m_y
       if( .not. repeated(i) ) then
          z % v(m_x+j) = y % v(i) 
          z % i(m_x+j) = y % i(i) 
          j = j + 1
       else
          k = minloc( abs( x % i - y % i(i) ), dim=1 )
          z % v(k) = z % v(k) + y % v(i) 
       end if 
    end do 
    
    deallocate( i_temp, repeated )
     
 end function
  
  
 subroutine Sum_sparse_vector_sub(x, y, z)
    class(Sparse_vector), intent(in) :: x, y
    class(Sparse_vector), intent(out) :: z
    
    integer :: m_x, m_y, m_z, i, j, k
    integer, allocatable :: i_temp(:) 
    logical, allocatable :: repeated(:) 
    
    m_x = size( x % i ) 
    m_y = size( y % i ) 
    
    ! *** Temporal array of indices
    allocate( i_temp(m_x + m_y), repeated(m_y) ) 
    i_temp = [ x % i, -1 ]
    repeated = .false.
    do i = 1, m_y       
       ! *** Assign non repeated indices
       if ( count( x % i == y % i(i) ) == 0 ) i_temp(m_x+i)  = y % i(i)
       if ( count( x % i == y % i(i) ) /= 0 ) repeated(i) = .true.    
    end do
    
    ! *** Dimension of the sparse sum vector z = x + y
    m_z = m_x + count( .not. repeated )
    
    allocate( z % i(m_z),  z % v(m_z) )     
    
    z % v(1:m_x) = x % v
    z % i(1:m_x) = x % i
    
    j = 1
    do i = 1, m_y
       if( .not. repeated(i) ) then
          z % v(m_x+j) = y % v(i) 
          z % i(m_x+j) = y % i(i) 
          j = j + 1
       else
          k = minloc( abs( x % i - y % i(i) ), dim=1 )
          z % v(k) = z % v(k) + y % v(i) 
       end if 
    end do 
    
    deallocate( i_temp, repeated )
     
 end subroutine
  
  
 ! *** Initialize a sparse matrix of m rows and n columns
 subroutine Initialize_sparse_matrix( A, m, n)
     class(Sparse_Matrix) :: A
     integer, intent(in) :: m, n

     A % m = m ; A % n = n  
     if( allocated( A % rows ) ) deallocate( A % rows )
     if( allocated( A % columns ) ) deallocate( A % columns )
     ! *** Matrix size   
     allocate( A % rows(m) )
     allocate( A % columns(n) )
     
 end subroutine
 
 subroutine Initialize_sparse_vector( x, n )
     class(Sparse_vector), intent(out) :: x
     integer, intent(in) :: n     
     
     if( allocated( x % i ) ) deallocate( x % i )
     if( allocated( x % v ) ) deallocate( x % v )     
        
     allocate( x % i(n), x % v(n) )
 
 end subroutine 
 
 ! *** Given a vector v(1:n) fills the Sparse column type Sparse_v with the non
 ! zero elements of v in Sparse_v % v
 subroutine Fill_Sparse_vector( v, Sparse_v )
     real, intent(in) :: v(:)
     class(Sparse_vector), intent(out) :: Sparse_v
 
     integer :: m_j, i
     
     ! *** Number of non zero rows of v
     m_j = count(v/=0)        
     allocate( Sparse_v % v(m_j), Sparse_v % i(m_j) )   
       
     ! *** Index position of non zero rows of v
     Sparse_v % i = pack( [( i, i=1,size(v) )], v/=0 ) 
     
     ! *** Values of non zero rows of v
     Sparse_v % v = pack( v, v/=0 ) 
 
 end subroutine


 subroutine Fill_Sparse_Matrix_row( i, v, A )
     integer, intent(in) :: i ! Row to be filled
     real, intent(in) :: v(:)
     class(Sparse_Matrix) :: A
 
     call Fill_Sparse_vector( v, A % rows(i) )
 
 end subroutine
 
 subroutine Test_Fill_Sparse_columns
    integer, parameter :: M=150000, N = M
    
    integer :: i, j, N_nonzero
    real :: A(M,N), t0, t1
    type(Sparse_Matrix) :: A_sp
    
    call Initialize_sparse_matrix( A_sp, M, N )
    
    call random_number(A)
    
   
    A = 0
    
    ! *** Fill sparse_matrix by rows
    do i = 1, M    
       N_nonzero = count( A(i,:) /= 0 )
       write(*,*) " N_nonzero ", N_nonzero       
       allocate( A_sp % rows(i) % i(N_nonzero), &
                 A_sp % rows(i) % v(N_nonzero) )
       
       A_sp % rows(i) % i = pack( [(j,j=1,N)], A(i,:) /= 0 )
       A_sp % rows(i) % v = pack( A(i,:), A(i,:) /= 0 )    
    end do
    
    call cpu_time(t0) 
    call Fill_Sparse_Matrix_columns(A_sp) 
    call cpu_time(t1) 
    write(*,*) " cpu_time Fill_Sparse_Matrix_columns", t1- t0
    
    do j = 1, N
      !write(*,*) " Column ", j
      !write(*,*) " A(:,j) " , A(:,j)
      !write(*,*) " A_sp % columns(j) % i " , A_sp % columns(j) % i
      !write(*,*) " A_sp % columns(j) % v " , A_sp % columns(j) % v
      !write(*,*) " maxval( abs( pack(A(:,j),A(:,j)/=0) - A_sp % columns(j) % v ) ) " 
      write(*,*)   maxval( abs( pack(A(:,j),A(:,j)/=0) - A_sp % columns(j) % v ) )  == 0
      !write(*,*) " maxval( abs( pack( [(i,i=1,M)], A(:,j)/=0) - A_sp % columns(j) % i ) ) " 
      write(*,*)   maxval( abs( pack( [(i,i=1,M)], A(:,j)/=0) - A_sp % columns(j) % i ) )  == 0
      !write(*,*)  
    
    end do
 
 
 end subroutine
 
 subroutine Test_Fill_Sparse_vector
    integer, parameter :: N = 20000000
    real :: v(N)
    type(Sparse_vector) :: Sp_v
    
    integer :: k, m_j
    real :: t0, t1
    
    call cpu_time(t0)
    ! *** Toy vector
    v=0 ; v(14:16) = [1,2,3]
    
    ! *** Sparse vector    
    call Fill_Sparse_vector( v, Sp_v )
    
    call cpu_time(t1)
    write(*,*)  " Time in Fill_Sparse_vector of dimension", &
       N, " = ", t1-t0 
    m_j = count(v/=0)
    write(*,*) " m_j ", m_j
    write(*,*)  " k, Sp_v % i(k), Sp_v % v(k), v(Sp_v % i(k))  " 
    do k = 1, m_j
       write(*,*) k, Sp_v % i(k), Sp_v % v(k), v(Sp_v % i(k))
    end do
    
 end subroutine









end module
