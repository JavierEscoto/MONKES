program main_monkes
  
  use Magnetic_configuration
  !use API_Example_Barycentric_grids
  use API_Example_DKE_BTD_Solution_Legendre 
  
implicit none 
  
  character(len=20) :: input_case="none"
  
  write(*,*) " ****************************************************** "
  write(*,*) " This is MONKES: "
  write(*,*) " MONoenergetic Kinetic Equation solver for Stellarators "
  write(*,*)  
  write(*,*) " Author: F. Javier Escoto Lopez (javier.escoto.lopez@gmail.com) "
  write(*,*) " ****************************************************** "
  write(*,*)  
    
  
  ! Select input type: BOOZER_XFORM output or DKES input
  call Select_Input(input_case) 
  call Select_surface  ! Select flux surface (needed for BOOZER_XFORM output)
  
  select case (input_case)
  
     case ("boozmn.nc")         
     call read_boozer_xform_output(s)
     
     case ("ddkes2.data")
     call read_ddkes2data("ddkes2.data")
     
     case default
     write(*,*) " There is no BOOZER_XFORM output 'boozmn.nc' &
                  nor ddkes2.data in this folder. "
     
   end select

  !call Test_Fourier_Discretization_BTD
  !call Monoenergetic_convergence_xi
  !call Adjoint_Monoenergetic_Database_Input
  !call Test_Monoenergetic_lambda_function

  call Monoenergetic_Database_Input_NEW
  !call Monoenergetic_Database_Input
  
  
  !call Monoenergetic_Database_Input_Orthonormal
  !call Monoenergetic_Database_Example_1_Orthonormal
  !call Monoenergetic_Database_Example_1 
  !call DKE_zeta_Convergence_Example_1
  !call Example_DFT_1D
  contains
  
  subroutine Select_Input(input_case) 
    character(len=20), intent(out) :: input_case  
    
    integer :: ierr     
    
    open(1, file= "boozmn.nc", status="old", iostat=ierr) ; close(1) 
    
    if( ierr == 0 ) then
      input_case = "boozmn.nc"
    else
      input_case = "ddkes2.data" 
    end if     
    
  end subroutine
  
  
  
end program
