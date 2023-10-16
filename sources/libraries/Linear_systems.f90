module Linear_systems 
  
 use Sparse_Matrices
 implicit none 


contains 
  
 function matmul_LAPACK( A, B ) result(C)
    real, intent(in) :: A(:,:), B(:,:)
    real :: C( size(A,dim=1), size(B,dim=2) )
    integer :: M, N, K
    
    M = size(A,dim=1) ; K = size(A,dim=1) ; N = size(B,dim=2)    
    call dgemm('N','N', M, N, K, 1d0, A, M, B, K, 0d0, C, M ) 
 
 end function
 
 function Invert_LU_LAPACK( A ) result(A_inv) 
    real, intent(in) :: A(:,:) 
    real :: A_inv(size(A,dim=1),size(A,dim=2)) 
        
    real :: A_LU(size(A,dim=1),size(A,dim=2))
    logical, parameter :: timing = .true.  
    integer :: N, info, ipiv(size(A,dim=1)) 
    integer :: c0, c1, c_rate ! for clock time      
    
    ! Initialize compiler internal clock 
    if(timing) call system_clock(count_rate=c_rate)  
    
    N = size(A,dim=1)  
    
    A_LU = A 
    ! *** LU Factorization
    if(timing) write(*,*)  "         *** Clock time for Invert_LU_LAPACK "
    if(timing) call system_clock(c0)
    call DGETRF( N, N, A_LU, N, ipiv, info )    
    if(timing) call system_clock(c1)
    
    ! *** Instants of time in seconds
    if(timing)  & 
    write(*,*) "           Time in LU factorization = ", ( c1 - c0  ) / real(c_rate)
    
    ! *** Inverse of the factorized matrix
    if(timing) call system_clock(c0)
    A_inv = Identity(N)     
    call DGETRS( 'N', N, N, A_LU, N, ipiv, A_inv, N, info) ! Computes solution for given source     
    if(timing) call system_clock(c1)
    
    ! Instants of time in seconds
    if(timing)  & 
    write(*,*) "           Time in inverse calculation from LU = ", ( c1 - c0  ) / real(c_rate)
 end function
 

 function Solve_System_LU_LAPACK( A, B ) result(X) 
    real, intent(in) :: A(:,:), B(:,:)
    real :: X(size(A,dim=1),size(A,dim=2)) 
        
    real :: A_LU(size(A,dim=1),size(A,dim=2))
    logical, parameter :: timing = .true.  
    integer :: N, info, ipiv(size(A,dim=1)) 
    integer :: c0, c1, c_rate ! for clock time      
    
    ! Initialize compiler internal clock 
    if(timing) call system_clock(count_rate=c_rate)  
    
    N = size(A,dim=1)  
    
    A_LU = A 
    ! *** LU Factorization
    if(timing) write(*,*)  "         *** Clock time for Solve_System_LU_LAPACK "
    if(timing) call system_clock(c0)
    call DGETRF( N, N, A_LU, N, ipiv, info )    
    if(timing) call system_clock(c1)
    
    ! *** Instants of time in seconds
    if(timing)  & 
    write(*,*) "           Time in LU factorization = ", ( c1 - c0  ) / real(c_rate)
    
    ! *** Solution of A*X=B of the factorized matrix
    if(timing) call system_clock(c0)
    X = B     
    call DGETRS( 'N', N, N, A_LU, N, ipiv, X, N, info) ! Computes solution for given source     
    if(timing) call system_clock(c1)
    
    ! Instants of time in seconds
    if(timing)  & 
    write(*,*) "           Time in solving system from LU = ", ( c1 - c0  ) / real(c_rate)
 end function

 subroutine Solve_Matrix_System_LAPACK( A, B, factorized, ipiv, X )  
    real, intent(inout) :: A(:,:) 
    real, intent(in) ::    B(:,:)
    logical, intent(in) :: factorized
    integer, intent(inout) :: ipiv(size(A,dim=1)) 
    real, intent(out) :: X( size(B,dim=1), size(B,dim=2) ) 
    
         
    logical, parameter :: timing = .true.  
    integer :: N, N_rhs, info
    integer :: c0, c1, c_rate ! for clock time      
    
    ! Initialize compiler internal clock 
    if(timing) call system_clock(count_rate=c_rate)  
    
    N = size(A,dim=1) ; N_rhs = size(B, dim=2)
     
    if( .not. factorized ) then
      ! *** LU Factorization
      if(timing) write(*,*)  "         *** Clock time for Solve_System_LU_LAPACK "
      if(timing) call system_clock(c0)
      call DGETRF( N, N, A, N, ipiv, info )    
      if(timing) call system_clock(c1)
      
      ! *** Instants of time in seconds
      if(timing)  & 
      write(*,*) "           Time in LU factorization = ", ( c1 - c0  ) / real(c_rate)
    end if
    
    ! *** Solution of A*X=B of the factorized matrix
    if(timing) call system_clock(c0)
    X = B     
    call DGETRS( 'N', N, N_rhs, A, N, ipiv, X, N, info) ! Computes solution for given source     
    if(timing) call system_clock(c1)
    
    ! Instants of time in seconds
    if(timing)  & 
    write(*,*) "           Time in solving system from LU = ", ( c1 - c0  ) / real(c_rate)
 end subroutine
  
pure function Identity(N) result(Id)
   integer, intent(in) :: N
   real :: Id(N,N)
   
   integer :: i 
   
   Id = 0   
   do i = 1, N           
      Id(i,i) = 1           
   end do   
   
end function

 
end module 
