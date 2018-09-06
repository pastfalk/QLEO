!> Scans through the wavenumber interval and computes corresponding frequencies (based on dispersion relation solver LEOPARD, see Astfalk & Jenko, JGR 2017)
!! \param krange wavenumber grid
!! \param omega_start initial frequency guess for the first wavenumber
!! \param increment frequency increment to determine the initial guesses for subsequent wavenumbers, before subroutine polyfit() takes over
!! \param om real part of the approximate root of dispersion relation obtained from Muller iteration
!! \param ga imaginary part of the approximate root of dispersion relation obtained from Muller iteration
subroutine disp_rel(krange,omega_start,increment,om,ga)
  use param_mod
  implicit none
  complex :: omega_start, increment
  integer :: iarb, ik
  complex, dimension (nk) :: solution
  real, dimension (nk) :: krange
  real :: gamma_max
  real, dimension(nk) :: om, ga
  real, allocatable, dimension(:,:,:,:,:) :: splcoeff1, splcoeff2

  !spline-interpolate the velocity distributions

  allocate(splcoeff1(npara(1)-1,nperp(1)-1,4,3,narb))
  allocate(splcoeff2(npara(1)-1,nperp(1)-1,4,3,narb))

  do iarb=1,narb
     call splcoeff_dist(iarb,splcoeff1(:,:,:,:,iarb),splcoeff2(:,:,:,:,iarb))
  enddo

  !scan through wavenumber interval
  do ik=1,nk

     if(ik.gt.1) then

        if(aimag(solution(ik-1)).eq.0.0) then
           solution(ik)=0.0
           cycle
        endif
     endif

     !use Muller method to iterate root of dispersion relation

     call muller(omega_start,krange(ik),solution(ik),splcoeff1, splcoeff2)

     if ((ik .ge. 3).and.(ik .lt. nk))  then

        !if three subsequent solutions omega(k) are found, use quadratic polynomial fit to guess next starting frequency for Muller iteration
        call polyfit(krange(ik-2:ik+1),solution(ik-2:ik),omega_start)

     else

        !for the first two solution omega(k) guess next starting frequency for Muller iteration by raising the computed omega by an increment
        omega_start=solution(ik)+increment

     end if

     write(57,'(F12.8,E20.10,E20.10)') krange(ik), real(solution(ik)), aimag(solution(ik))

  enddo

  gamma_max=maxval(aimag(solution))

  write(*,*) '************'
  write(*,*) 'k_1 /omega:', krange(1), solution(1)
  write(*,*) 'k_nk/omega:', krange(nk), solution(nk)
  write(*,*) 'gamma_max:', gamma_max
  write(*,*) '************'

  om=real(solution)
  ga=aimag(solution)

  omega_start=solution(1)
  increment=solution(2)-solution(1)

  deallocate(splcoeff1,splcoeff2)


end subroutine disp_rel
