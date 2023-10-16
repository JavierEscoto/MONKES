module API_Example_Coarray
 


implicit none


contains 


  subroutine Test_Hello_world_Coarray
    integer :: i  
    write(*,*) " *** Test_Hello_world_Coarray "
    
    ! Distribute information to other images
    do i = 1, num_images()
      
       write (*,*) "Hello from image", this_image(), "of", num_images()
    end do
    
    sync all ! Barrier to make sure the data has arrived
    
    
    stop
  end subroutine




end module
