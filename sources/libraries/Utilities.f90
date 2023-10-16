

module Utilities



 implicit none

    type Pointer_vector
       real, pointer :: p
    end type

    contains

 function Maxloc_2D( U ) result(indices)
     real, intent(in) :: U(:,:)
     integer :: indices(2) 
     
     integer :: i, j, l_max, Nx, Ny
     
     real, target :: UU(size(U,dim=1), size(U,dim=2))
     real, pointer :: pUU(:)
     
     Nx = size(U,dim=1) ; Ny = size(U,dim=2) 
     ! *** Bounds remapping of input
     UU = U
     pUU(1:Nx*Ny) => UU 
     
     l_max = maxloc(pUU,dim=1) 
     
     do j = 1, Ny
        do i = 1, Nx
           if ( U(i,j) == pUU(l_max)  ) indices = [i,j]
        end do
     end do
     
 end function

!*********************************
! mantissa, exponent
!*********************************
 subroutine mantissa_exponent(x, m, e)
      real, intent(in) :: x
      real, intent(out) :: m
      integer, intent(out) ::  e

      real :: xx

      if (x < 0.) then

         xx = -x
         e = int(log10(xx))
         if (e > 0) e = e + 1
         m = - xx * 1.e1 ** (-e)

      else if (x > 0.) then

         e = int(log10(x))
         if (e > 0) e = e + 1
         m = x * 1.e1 ** (-e)

      else

         e = 0.
         m = 0.

      endif

end subroutine




!*****************************************************************
!*
!*****************************************************************
subroutine Data_pointer( N1, N2, U, pU )
     integer, intent(in) :: N1, N2
     real, target, intent(in) :: U(N1, N2)
     real, pointer, intent(out) :: pU(:,:)

     pU => U

end subroutine

!*****************************************************************
!*      Points a vector of explicit shape U(N1*N2) with the pointer
!      pU(:).  If this subroutine is called with an explicit expression,
!      then TKR rule can be broken to allow U to be, for example an
!      array of rank 3 U(1:M,1:N,1:K) with M*N*K = N1*N2.
!*****************************************************************
subroutine Vector_pointer( N1, N2, U, pU )
     integer, intent(in) :: N1, N2
     real, target, intent(in) :: U(N1*N2)
     real, pointer, intent(out) :: pU(:)

     pU => U

end subroutine






! *****************************************************************
!
!     Function that gives the associated index l(Indices) between the indices
!    of the target U(:,:,...,:) and the pointer pU(:) given the values of the
!    indices stored on Indices(:) and their maximum values Bounds(:)
!
!
pure integer function Vector_pointer_index(Indices,Bounds) result(L)
           integer, intent(in) :: Indices(:), Bounds(:)
           integer :: i, N, Factor(size(Bounds))

           N = size(Bounds)

           Factor(1) = 1
           do i = 2, N
              Factor(i) = Bounds(i-1) * Factor(i-1)
           end do

           L = 1 + dot_product( Indices-1, Factor)

end function


function Reshape_pointer(N, A) result(v)
     integer, intent(in) :: N
     real, target, intent(in) :: A(N)
     real, pointer :: v(:)

     v => A

end function




!*****************************************************************************
!* It searches the element x in an ordered list V ( x = V(ix) )
!* If the element x is not in the list it returns -1
!*****************************************************************************
subroutine Binary_search(x, V, ix)
    real, intent (in) :: x, V(0:)
    integer, intent (out) :: ix

  integer:: m, upper_limit, lower_limit, N
  N = size(V) - 1

  lower_limit = 0; upper_limit = N;

  do while( upper_limit > lower_limit )

     m = ( upper_limit + lower_limit ) / 2

     if ( x > V(m) ) then
                          lower_limit = m + 1
     else
                          upper_limit = m
     endif

  end do

  !if (x == V(m) ) then
  !      ix = m
  !else
  !      ix = -1
  !endif
  ix = lower_limit


end subroutine




!*****************************************************************************
!* It calculates trace of matrix A
!*****************************************************************************
real function trace(A)
  real, intent(in) :: A(:,:)

  integer :: i
  real :: S

  S = 0
  do i=1, size(A, dim=1)

      S = S + A(i,i)

  end do

  trace = S


end function

!*****************************************************************************
!* It allocates a Vandermonde matrix of dimension M
!*****************************************************************************
function Vandermonde(M) result(A)
  integer, intent(in) :: M
  real, allocatable :: A(:,:)

  integer :: i, j

  allocate ( A(M, M) )

  do i=1, M; do j=1, M
      A(i,j) = ( i / real(M) ) **(j-1)
  end do; end do


end function



!****************************************************************************************************
!* It looks for the given name in an database (array of names) and if returns its position
!****************************************************************************************************
logical function  equal(a, b)
   class(*), intent(in) :: a, b


   select type (a)

       type is (integer)
            select type (b)
            type is (integer)
            equal = (a == b)
            end select


       type is (character(len=*))
            select type (b)
            type is (character(len=*))
            equal = (trim(a) == trim(b))
            end select

       type is (real)
            select type (b)
            type is (real)
            equal = (a == b)
            end select

       class default
             !write(*,*) "ERROR "
             stop

   end select




end function

!****************************************************************************************************
!* It looks for the given element x in a set of values. If it exists, returns its position
!****************************************************************************************************
  function  look_forg(x, set)
   class(*), intent(in) :: x, set(:)
   integer :: look_forg


   integer :: i, N

   N = size(set)

   look_forg = 0

   do i=1, N
      if ( equal( x, set(i) ) ) then
                                look_forg = i
                                exit
      endif

   end do


end function


!*****************************************************************************
!* It splits a string into N different substrings based on the delimiter
!*****************************************************************************
subroutine split( string, delimiter, substrings, N  )
   character(len=*), intent(in) :: string
   character(len=1), intent(in) :: delimiter
   character(len=*), intent(out) :: substrings(:)
   integer, intent(out) :: N


      integer :: i, index, k

    i = 1
    k = 1
    substrings = "a"
    do while (i<=len(string) )

        index = scan( string(i:), delimiter)

        if (index>0) then
                          substrings(k) = string(i:i+index-2)
                          i = i + index
                          k = k + 1
        else
             substrings(k) = string(i:len(string))
             N = k
             exit
        endif

    end do

    !write(*,'(A,A)') " strings = ", string
    !
    !do k=1, N
    !   write(*,*) " substrings = ", substrings(k)
    !end do
    !read(*,*)

end subroutine

!***********************************************************************
!   Given a real vector v(0:N) gives back a vector x(0:N) with the values
!   of v sorted in crescent order from i=0,1,...,N
!***********************************************************************

function Sort_vector( v ) result( x )
      real, intent(in) :: v(0:)
      real :: x(0:size(v)-1)

      integer :: i, N, min_index

      ! *** Size of v
      N = size(v) - 1

      ! *** Sorting of the vector
      x(0) = minval(v)
      do i = 1, N

         min_index = minloc( v, dim = 1, mask = v > x(i-1) ) - 1
         x(i) = v(min_index)

      end do

end function


! given a vector u(0:N) it reorders it and gives back a vector sorted in
! crescent order v(0:N) and optionally it gives an integer array
! indices(0:N,2) which relates both as
!    v(i) = u(indices(i,1))    ; u(i) = v(indices(i,2))

subroutine Sort_vector_indices( u , v , indices, i0 ) 
      real, intent(in) :: u(:)          ! unsorted input
      real, intent(out) :: v(size(u))   !   sorted output
      integer, optional, intent(out) :: indices(size(u),2)
      integer, optional, intent(in) :: i0 ! if we want the indices to
                                          ! start at i0

      integer :: i, N, min_index

      ! *** Size of v
      N = size(v)
      ! *** Sorting of the vector
      v(1) = minval(u)
      min_index = minloc(u, dim = 1)
      if( present(indices) ) indices(1,1) = min_index
      if( present(indices) ) indices(min_index, 2) = 1
      
      do i = 2, N
         min_index = minloc( u, dim = 1, mask = u > v(i-1) ) 
         v(i) = u(min_index)
         if( present(indices) ) indices(i,1) = min_index
         if( present(indices) ) indices(min_index,2) = i
         
      end do
      
      ! *** If is required that the indices start at i0 instead of 1
      if( present(indices) .and. present(i0) ) indices = indices - 1 + i0
end subroutine

! given a vector u(0:N) it reorders it and gives back a vector sorted in
! crescent order v(0:N) and optionally it gives an integer array
! indices(0:N,2) which relates both as
!    v(i) = u(indices(i,1))    ; u(i) = v(indices(i,2))

subroutine Sort_vector_indices_0( u , v , indices ) 
      real, intent(in) :: u(0:)             ! unsorted input
      real, intent(out) :: v(0:size(u)-1)   ! sorted output
      integer, optional, intent(out) :: indices(0:size(u)-1,2)

      integer :: i, N, min_index

      ! *** Size of v
      N = size(v)-1
      ! *** Sorting of the vector
      v(0) = minval(u)
      min_index = minloc(u, dim = 1) - 1
      if( present(indices) ) indices(1,1) = min_index
      if( present(indices) ) indices(min_index, 2) = 0
      
      do i = 1, N
         min_index = minloc( u, dim = 1, mask = u > v(i-1) ) - 1 
         v(i) = u(min_index)
         if( present(indices) ) indices(i,1) = min_index
         if( present(indices) ) indices(min_index,2) = i
         
      end do
      
end subroutine

subroutine Sort_vector_indices_OLD( u , v , indices ) 
      real, intent(in) :: u(:)             ! unsorted input
      real, intent(out) :: v(size(u))   ! sorted output
      integer, optional, intent(out) :: indices(size(u),2)

      integer :: i, N, min_index

      ! *** Size of v
      N = size(v)
      ! *** Sorting of the vector
      v(1) = minval(u)
      min_index = minloc(u, dim = 1)
      if( present(indices))          indices(1,1) = min_index
      if( present(indices)) indices(min_index, 2) = 1
      
      do i = 2, N
         min_index = minloc( u, dim = 1, mask = u > v(i-1) ) 
         v(i) = u(min_index)
         if( present(indices))         indices(i,1) = min_index
         if( present(indices)) indices(min_index,2) = i
         
      end do
      
end subroutine


!*******************************************************************************
!    Counts the number of lines of a file given by the string input "path".
!*******************************************************************************

integer function Count_lines(path) result(lines)
        character(len=*), intent(in) :: path

        integer :: ios
        
        open(10,file=path)
        lines = 0

        do
            read(10,*,iostat=ios)
            if (ios/=0) exit
            lines = lines + 1
        end do
        close(10)

end function


!*************************************************************************
!    Read_lines. Reads the N first data lines of a file
!*************************************************************************

 subroutine Read_lines( path, N, U , N_header)
        character(len=*), intent(in) :: path
        integer, intent(in) ::  N
        integer,optional, intent(in) :: N_header
        real, intent(out) :: U(:, :) 

        integer :: i
        open(21,file=path)

        if ( present(N_header) ) then
             ! *** Skip header lines
             do i = 1, N_header
                 read(20,*)
             end do
         end if

        do i = 1, N
            read(21,*) U(i,:)
        end do
        close(21)
 end subroutine




end module





