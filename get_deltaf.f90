!> Carry out the wavenumber integration in the quasilinear equation for the velocity distribution's time evolution 
!! \param krange wavenumber grid
!! \param omega frequency spectrum
!! \param gamma growth rate spectrum
!! \param c_om spline coefficients of interpolated frequency spectrum
!! \param c_ga spline coefficients of interpolated growth rate spectrum
!! \param c_Bk spline coefficients of interpolated magnetic energy spectrum
!! \param vpa parallel velocity
!! \param sgn signum in the integrand's denominator which determines whether right-hand or left-hand modes are considered
!! \param contr1 solutions of the wavenumber integration
!! \param contr2 solutions of the wavenumber integration
subroutine get_deltaf(krange,c_om,c_ga,c_Bk,vpa,sgn,contr1,contr2)
  use param_mod
  implicit none
  integer :: n
  real, dimension (nk) :: krange
  real :: k1, k2
  complex, dimension(3) :: k0
  complex, dimension(0:2) :: cubic_coeff 
  integer :: sgn
  real, dimension (nk-1,4) :: t1_sum
  real, dimension (nk-1,2) :: t2_sum
  real :: t1a,t2a,t3a,t4a,t5a,t6a
  real :: t1b,t2b,t3b,t4b,t5b,t6b
  real, dimension(4) :: contr1
  real, dimension(2) :: contr2
  real, dimension (nk-1,0:3) :: c_om, c_ga
  real, dimension (nk-1,0:3) :: c_Bk
  real :: vpa

  !scan through the wavenumber grid and integrate piecewise-analytically the different terms in the wavenumber integral

  do n=1,nk-1

     k1=krange(n)
     k2=krange(n+1)

     if( (c_ga(n,0).eq.0.0).and.(c_ga(n,1).eq.0.0).and.(c_ga(n,2).eq.0.0).and.(c_ga(n,3).eq.0.0)) then

        t1_sum(n,:)=0.0
        t2_sum(n,:)=0.0

     else

        cubic_coeff(0)=(c_om(n,0)-sgn*mu(1)*q(1)+i*c_ga(n,0))/&
             & (c_om(n,3)+i*c_ga(n,3))

        cubic_coeff(1)=(c_om(n,1)-vpa+i*c_ga(n,1))/&
             & (c_om(n,3)+i*c_ga(n,3))

        cubic_coeff(2)=(c_om(n,2)+i*c_ga(n,2))/&
             & (c_om(n,3)+i*c_ga(n,3))

        call cubic_sol(cubic_coeff,k0)

        call Term1(k1,k2,k0, c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa, sgn,t1a)
        call Term4(k1,k2,k0, c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa, sgn,t4a)
        call Term5(k1,k2,k0, c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa, sgn,t5a)


        if(vpa.eq.0.0) then

           t2a=0.0
           t3a=0.0
           t6a=0.0

        else

           call Term2(k1,k2,k0, c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa, sgn,t2a)
           call Term3(k1,k2,k0, c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa, sgn,t3a)
           call Term6(k1,k2,k0, c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa,t6a)

        endif


        cubic_coeff(0)=(-c_om(n,0)+sgn*mu(1)*q(1)+i*c_ga(n,0))/&
             & (-c_om(n,3)+i*c_ga(n,3))

        cubic_coeff(1)=(-c_om(n,1)-vpa+i*c_ga(n,1))/&
             & (-c_om(n,3)+i*c_ga(n,3))

        cubic_coeff(2)=(-c_om(n,2)+i*c_ga(n,2))/&
             & (-c_om(n,3)+i*c_ga(n,3))

        call cubic_sol(cubic_coeff,k0)

        call Term1(k1,k2,k0,-c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa,-sgn,t1b)
        call Term4(k1,k2,k0,-c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa,-sgn,t4b)
        call Term5(k1,k2,k0,-c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa,-sgn,t5b)


        if(vpa.eq.0.0) then

           t2b=0.0
           t3b=0.0
           t6b=0.0

        else

           call Term2(k1,k2,k0,-c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa,-sgn,t2b)
           call Term3(k1,k2,k0,-c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa,-sgn,t3b)
           call Term6(k1,k2,k0,-c_om(n,:), c_ga(n,:),c_Bk(n,:),vpa,t6b)

        endif

        t1_sum(n,1)=t1a+t1b
        t1_sum(n,2)=t2a+t2b
        t1_sum(n,3)=t3a+t3b
        t1_sum(n,4)=t4a+t4b

        t2_sum(n,1)=t5a+t5b
        t2_sum(n,2)=t6a+t6b

     endif

  enddo

  !add up the analytic solutions of the integrals over the whole wavenumber grid

  contr1=0.0
  contr2=0.0

  do n=1,nk-1

     contr1(1)=contr1(1)+t1_sum(n,1)
     contr1(2)=contr1(2)+t1_sum(n,2)
     contr1(3)=contr1(3)+t1_sum(n,3)
     contr1(4)=contr1(4)+t1_sum(n,4)

     contr2(1)=contr2(1)+t2_sum(n,1)
     contr2(2)=contr2(2)+t2_sum(n,2)

  enddo

  if(krange(1).gt.krange(2)) then
     contr1=-contr1
     contr2=-contr2
  endif

  contr1=-0.25 * contr1 * q(1)**2 * mu(1)**2
  contr2=-0.25 * contr2 * q(1)**2 * mu(1)**2


end subroutine get_deltaf
