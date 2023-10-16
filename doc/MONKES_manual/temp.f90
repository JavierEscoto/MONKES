  subroutine Monoenergetic_Database_Scan( nu, E_r, N_theta, N_zeta, N_xi, DD, DD_33_Sp  )
     real, intent(in) :: nu(:), E_r(:) 
     integer, intent(in) :: N_theta(:), N_zeta(:), N_xi(:) 
     real, optional, intent(out) :: DD( 3, 3, size(nu), size(E_r), size(N_theta), size(N_zeta), size(N_xi) )
     real, optional, intent(out) :: DD_33_Sp( size(nu), size(E_r), size(N_theta), size(N_zeta), size(N_xi) )
     
     real :: D( 3, 3, size(nu), size(E_r), size(N_theta), size(N_zeta), size(N_xi) )
     real :: D_33_Sp( size(nu), size(E_r), size(N_theta), size(N_zeta), size(N_xi) )
  
     integer :: i, j, ii, jj, kk, N_nu, N_E_r, M_theta, M_zeta, M_xi    
     character(len=500) :: file_path
     real :: t_clock, rate, t0, t1
     integer :: c0, c1, c_rate 
     
     ! Initialize compiler internal clock for computing wall-clock time
     call system_clock(count_rate=c_rate) ; rate = real(c_rate)
     
     ! Vectors sizes
     N_nu = size(nu) ; N_E_r = size(E_r)
     M_theta = size(N_theta) ; M_zeta = size(N_zeta) ; M_xi = size(N_xi)
     
     ! Location for output
     file_path = "monkes_Monoenergetic_Database.dat"     
     ! Open output file and write header
     open(21, file=trim(file_path))
     write(21,'(9999A25)') " nu/v [m^-1]", " E_r/v [kV s /m^2]", &
                             " N_theta ", " N_zeta ", " N_xi ", &
                             " D_11 ", " D_31 ", &
                             " D_13 ", " D_33 ", &
                             " D_33_Spitzer ",   &
                             " Wall-clock time [s] ",   &
                             " CPU time [s] "
     do j = 1, N_E_r ! Loop electric field value
        do i = 1, N_nu ! Loop collisionality value   
           do kk = 1, M_xi ! Loop number of Legendre modes
              do ii = 1, M_theta ! Loop number of theta points
                 do jj = 1, M_zeta ! Loop number of zeta points                  
                 
                    call system_clock(c0) ; call cpu_time(t0) 
                    call Solve_BTD_DKE_Legendre( N_theta(ii),          &
                                                 N_zeta(jj),           &
                                                 N_xi(kk),             &
                                                 nu(i), E_r(j),        &
                                                 D(:,:,i,j,ii,jj,kk),  &
                                                 D_33_Sp(i,j,ii,jj,kk) )                 
                    call system_clock(c1) ; call cpu_time(t1)
                    
                    ! Wall-clock time in seconds
                    t_clock = ( c1 - c0 ) / rate   
                                  
                    ! Writing results on "monkes_Monoenergetic_Database.plt"
                    write(21,'(9999e)') nu(i), E_r(j), &
                                        real(N_theta(ii)), &
                                        real(N_zeta(jj)), &
                                        real(N_xi(kk)), &
                                        D(1,1,i,j,ii,jj,kk), &
                                        D(3,1,i,j,ii,jj,kk), &
                                        D(1,3,i,j,ii,jj,kk), &
                                        D(3,3,i,j,ii,jj,kk), &
                                        D_33_Sp(i,j,ii,jj,kk), &
                                        t_clock, t1-t0 
                 end do  
              end do  
           end do  
        end do  
     end do
     close(21)
     
     ! Monoenergetic coefficients as output
     if ( present(DD) ) DD = D
     if ( present(DD_33_Sp) ) DD_33_Sp = D_33_Sp
  end subroutine 
