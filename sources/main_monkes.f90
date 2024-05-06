program main_monkes
  
  use Magnetic_configuration
  use API_Example_DKE_BTD_Solution_Legendre 
  
implicit none 
  
  character(len=30) :: input_case="none"
  
  write(*,*) " ****************************************************** "
  write(*,*) " This is MONKES: "
  write(*,*) " MONoenergetic Kinetic Equation Solver "
  write(*,*)  
  write(*,*) " Author: F. Javier Escoto Lopez (javier.escoto.lopez@gmail.com) "
  write(*,*) " ****************************************************** "
  write(*,*)  
      
  ! *** Select input type: BOOZER_XFORM output or DKES input
  call Select_Input(input_case) 
  call Select_surface  ! Select flux surface (needed for BOOZER_XFORM output)
  
  !input_case = "VMEC.nc"  ! ad-hoc for testing (TO BE ERASED)
  select case (input_case)
  
     case ("monkes_input")         
     call read_MONKES_input_configuration("monkes_input.configuration")
     
     case ("boozmn.nc")         
     call read_boozer_xform_output(s)
     
     case ("ddkes2.data")
     call read_ddkes2data("ddkes2.data")
     
     case ("VMEC.nc")
     call read_VMEC_output(s) 
     
     case default
     write(*,*) " ERROR: There is no BOOZER_XFORM output 'boozmn.nc' &
                  nor 'ddkes2.data' in this folder. "
     
  end select
  
  ! *** Select the case: Either a monoenergetic database for specified collisionalities 
  ! or neoclassical flows (WORK IN PROGRESS).
  call Monoenergetic_Database_Input
  
  
  ! *** Test of speed grid using Maxwell polynomials (WORK IN PROGRESS)
  call Monoenergetic_Database_Maxwell_points ; stop
  
  
  contains
  
  subroutine Select_Input(input_case) 
    character(len=20), intent(out) :: input_case  
    
    integer :: ierr     
    
    open(1, file= "monkes_input.configuration", status="old", iostat=ierr) ; close(1)  
    if( ierr == 0 ) then
        input_case = "monkes_input"
    else
    
      open(1, file= "boozmn.nc", status="old", iostat=ierr) ; close(1)       
      if( ierr == 0 ) then
        input_case = "boozmn.nc"
      else    
        
        open(1, file= "VMEC.nc", status="old", iostat=ierr) ; close(1)       
        if( ierr == 0 ) then
          input_case = "VMEC.nc"
        else
          input_case = "ddkes2.data" 
        end if
      end if     
    end if      
    
  end subroutine
  
  
  
end program
