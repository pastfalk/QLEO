!> Adjust wavenumber interval to unstable wavenumber range
!! \param krange_full full wavenumber grid
!! \param Bksq_full magnetic energy spectrum over the full wavenumber grid
!! \param omega_full frequency spectrum over the full wavenumber grid
!! \param gamma_full growth rate spectrum over the full wavenumber grid
!! \param krange_adapt adjusted wavenumber grid which covers unstable wavenumber range only
!! \param Bksq_adapt magnetic energy spectrum over adjusted wavenumber grid
!! \param increment_adapt frequency increment to determine the initial guesses for subsequent wavenumbers of adjusted wavenumber grid
!! \param omega_start_adapt initial frequency guess for the first wavenumber of the adjusted grid
subroutine adapt_krange(krange_full,Bksq_full,omega_full,gamma_full,&
     & krange_adapt,Bksq_adapt,increment_adapt,omega_start_adapt)
  use param_mod
  implicit none
  integer :: ik, ikold
  complex :: omega_start_adapt
  complex :: increment_adapt
  real, dimension(nk) :: krange_adapt, Bksq_adapt
  real, dimension(nk) :: krange_full, Bksq_full, omega_full, gamma_full
  real :: kstart, kend, dk
  real :: om1, om2, ga1, ga2
  real, dimension(nk-1,0:3) :: c_Bk, c_om, c_ga


  !locate unstable modes with lowest and highest wavenumber in considered wavenumber interval
  
  do ik=1,nk
     if(gamma_full(ik).gt.10.0**(-5)) then
        kstart=krange_full(ik)
        omega_start_adapt=omega_full(ik)+i*gamma_full(ik)
        exit
     endif
  enddo

  do ik=nk,1,-1
     if(gamma_full(ik).gt.10.0**(-5)) then
        kend=krange_full(ik)
        exit
     endif
  enddo


  !define adjusted wavenumber grid

  dk=(kend-kstart)/(nk-1.0)
  do ik=1,nk
     krange_adapt(ik)=kstart+(ik-1)*dk
  enddo


  !spline-interpolate magnetic energy, frequency, and growth rate spectrum

  call splcoeff_spec(krange_full,Bksq_full,c_Bk,nk)
  call splcoeff_spec(krange_full,omega_full,c_om,nk)
  call splcoeff_spec(krange_full,gamma_full,c_ga,nk)


  !sample the frequency, growth rate, and magnetic energy spectrum on the adusted wavenumber grid using the corresponding spline interpolations

  do ik=1,nk
     
     ikold=0

     do ikold=1,nk-1
        if((   (krange_full(1).lt.krange_full(2)).and.(krange_full(ikold).gt.krange_adapt(ik))).or.&
             &((krange_full(1).gt.krange_full(2)).and.(krange_full(ikold).lt.krange_adapt(ik)))) then
           exit
        endif
     enddo
     
     Bksq_adapt(ik)=c_Bk(ikold-1,3)*krange_adapt(ik)**3 + c_Bk(ikold-1,2)*krange_adapt(ik)**2 +&
          &   c_Bk(ikold-1,1)*krange_adapt(ik)    + c_Bk(ikold-1,0)

     if(ik.eq.1) then

        om1=c_om(ikold-1,3)*krange_adapt(ik)**3 + c_om(ikold-1,2)*krange_adapt(ik)**2 +&
          &   c_om(ikold-1,1)*krange_adapt(ik)    + c_om(ikold-1,0)

        ga1=c_ga(ikold-1,3)*krange_adapt(ik)**3 + c_ga(ikold-1,2)*krange_adapt(ik)**2 +&
          &   c_ga(ikold-1,1)*krange_adapt(ik)    + c_ga(ikold-1,0)

        om2=c_om(ikold-1,3)*krange_adapt(ik+1)**3 + c_om(ikold-1,2)*krange_adapt(ik+1)**2 +&
          &   c_om(ikold-1,1)*krange_adapt(ik+1)    + c_om(ikold-1,0)

        ga2=c_ga(ikold-1,3)*krange_adapt(ik+1)**3 + c_ga(ikold-1,2)*krange_adapt(ik+1)**2 +&
          &   c_ga(ikold-1,1)*krange_adapt(ik+1)    + c_ga(ikold-1,0)

        increment_adapt=(om2-om1)+i*(ga2-ga1)

    endif

  enddo


end subroutine adapt_krange
