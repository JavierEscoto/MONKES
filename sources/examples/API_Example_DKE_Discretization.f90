module API_Example_DKE_Discretization

  use Barycentric_grids
  use Magnetic_configuration

  implicit none

  private

  public :: Test_Fourier_Discretization_BTD
contains

  subroutine Test_Fourier_Discretization_BTD
    integer, parameter :: N_theta = 23, N_zeta = 55, N_xi = 400, N_fs = N_theta * N_zeta
    integer, parameter :: N_m = 4, N_n = 4
    real, parameter :: pi = acos(-1d0)
    integer :: i, j, k
    real :: m, n
    real :: B_max 
    real :: theta(0:N_theta), zeta(0:N_zeta), D_theta(0:N_theta-1,0:N_theta-1), D_zeta(0:N_zeta-1,0:N_zeta-1)
    real :: e_mn(0:N_theta-1,0:N_zeta-1), B(0:N_theta-1,0:N_zeta-1), dlnBdl(0:N_theta-1,0:N_zeta-1)
    real :: b_nabla_zeta(0:N_theta-1,0:N_zeta-1), dBdtheta(0:N_theta-1,0:N_zeta-1), dBdzeta(0:N_theta-1,0:N_zeta-1)
    real :: Le_mn(0:N_theta-1,0:N_zeta-1), De_mn(0:N_theta-1,0:N_zeta-1), Du_mn(0:N_theta-1,0:N_zeta-1)
    real :: Le_mn_exact(0:N_theta-1,0:N_zeta-1), De_mn_exact(0:N_theta-1,0:N_zeta-1), Ue_mn_exact(0:N_theta-1,0:N_zeta-1)
    real :: LB(0:N_theta-1,0:N_zeta-1), DB(0:N_theta-1,0:N_zeta-1), UB(0:N_theta-1,0:N_zeta-1)
    real :: LB_exact(0:N_theta-1,0:N_zeta-1), DB_exact(0:N_theta-1,0:N_zeta-1), UB_exact(0:N_theta-1,0:N_zeta-1)
    theta = [(i*2*pi/N_theta, i=0,N_theta)]
    zeta = [(j*2*pi/(Np*N_zeta), j=0,N_zeta)]

    m = 17.15!( N_theta - 1 )/2+1 
    n = ( N_zeta - 1 )/2 
    do j = 0, N_zeta -1
       D_zeta(j,:) = Fourier_Derivative_Row( zeta(j), zeta, 1 )
    end do
    do i = 0, N_theta -1
       D_theta(i,:) = Fourier_Derivative_Row( theta(i), theta, 1 )
    end do
    
    do j = 0, N_zeta-1
       do i = 0, N_theta-1
          B(i,j) = Magnetic_Field_Boozer( theta(i), zeta(j) )  
          b_nabla_zeta(i,j) = B(i,j) / abs( B_zeta + iota * B_theta )
          dBdtheta(i,j) = Magnetic_Field_Boozer_Derivative_theta( theta(i), zeta(j) ) 
          dBdzeta(i,j) = Magnetic_Field_Boozer_Derivative_zeta( theta(i), zeta(j) ) 
          dlnBdl(i,j) = ( 1d0 / abs( B_zeta + iota * B_theta ) ) & 
                      * ( iota *  dBdtheta(i,j) +  dBdzeta(i,j) )

        end do
     end do
     B_max = 1!maxval(B)
     B = B / B_max
     dBdtheta = dBdtheta/B_max
     dBdzeta = dBdzeta/B_max
     
     open(31,file="Test_DKE_DFT_Discretization.dat")
     write(31,'(9999A25)') "k ", "max Error L", "min Error L", " mean L " , " mean Error " , " mean Error /  mean L" 
     do k = 0, N_xi
        
        LB = Operator_L(k, (B)**m )
        !LB_exact = b_nabla_zeta *  ( iota *dBdtheta + dBdzeta ) &
            ! + 0.5*(k-1) * dlnBdl * B
        !LB_exact = 2* k * LB_exact /(2*k-1)
        LB_exact = ( iota *dBdtheta + dBdzeta ) * m * (B)**(m-1) 


             
        write(31,'(9999e)') ; write(31,'(9999e)') 
        write(31,'(9999e)')  real(k), maxval(abs(LB - LB_exact)), &
                      minval(abs(LB - LB_exact)),    &
                      sum( abs(LB_exact) ) /(N_theta*N_zeta), &
                      sum(abs(LB - LB_exact) ) / (N_theta*N_zeta), &
                      sum(abs(LB - LB_exact) ) / sum( abs(LB_exact) )

        do j =0, N_zeta-1 ; do i = 0, N_theta -1
            write(313,'(9999e)') theta(i), zeta(j), LB(i,j), LB_exact(i,j), LB(i,j) - LB_exact(i,j)
        end do ; end do 
     end do
     close(31)
   contains

     function Operator_L(k,F) result(L)
       integer, intent(in) :: k
       real, intent(in) :: F(0:N_theta-1,0:N_zeta-1)
       real :: L(0:N_theta-1,0:N_zeta-1)    

       integer :: i, j

       do j = 0, N_zeta -1 ; do i = 0, N_theta -1 
       L(i,j) = b_nabla_zeta(i,j) &
                * ( iota * Derivative_theta(i,j, F)  &
                    + Derivative_zeta(i,j, F) ) &
                + 0.5*(k-1) * dlnBdl(i,j) * F(i,j)
       L(i,j) =  k * L(i,j) / (2*k-1)
       L(i,j) =  ( iota * Derivative_theta(i,j, F)  &
                    + Derivative_zeta(i,j, F) ) 
       end do ; end do 
     end function


     real function Derivative_theta(i,j, U) result(D)
       integer, intent(in) :: i, j
       real, intent(in) :: U(0:N_theta-1,0:N_zeta-1)
       
       D = dot_product( D_theta(i,:), U(:,j)  )
       
     end function

     
     real function Derivative_zeta(i,j, U) result(D)
       integer, intent(in) :: i, j
       real, intent(in) :: U(0:N_theta-1,0:N_zeta-1)
       
       D = dot_product( D_zeta(j,:), U(i,:)  )
       
     end function
           
    
    


       
    
  end subroutine

  

end module
