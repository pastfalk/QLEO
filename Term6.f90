!> Analytic solution of the wavenumber integration of term T22 = dB_k^2 * k /(omega-k*vpa-sgn*mu*q)^2 over interval [k1,k2]
!! \param k1 lower limit of the considered wavenumber interval
!! \param k2 upper limit of the considered wavenumber interval
!! \param k0 roots of the integrand's denominator
!! \param c_om spline coefficients of interpolated frequency spectrum
!! \param c_ga spline coefficients of interpolated growth rate spectrum
!! \param c_Bk spline coefficients of interpolated magnetic energy spectrum
!! \param vpa parallel velocity
!! \param t6 result of the integration
subroutine Term6(k1,k2,k0,c_om,c_ga,c_Bk,vpa,t6)
  use param_mod
  implicit none
  real :: k1, k2
  complex, dimension(3) :: k0
  real, dimension(0:3) :: c_om, c_ga, c_Bk
  complex :: tk1, tk2
  real :: t6
  real, dimension(0:3) :: B
  integer :: sigma
  complex :: t6pole, contr_pole
  real :: vpa

  if(vpa.eq.0.0) then
     t6=0.0
  else

     B(0)= c_Bk(0)

     B(1)= c_Bk(1)

     B(2)= c_Bk(2)

     B(3)= c_Bk(3)

     tk1= (B(3)*k0(1)**4 + B(2)*k0(1)**3 + B(1)*k0(1)**2 + B(0)*k0(1))/&
          & ( (k0(1)-k1)*(k0(1)-k0(2))**2 *(k0(1)-k0(3))**2 )

     tk1=tk1+&
          & (B(3)*k0(2)**4 + B(2)*k0(2)**3 + B(1)*k0(2)**2 + B(0)*k0(2))/&
          & ( (k0(2)-k1)*(k0(2)-k0(1))**2 *(k0(2)-k0(3))**2 )

     tk1=tk1+&
          & (B(3)*k0(3)**4 + B(2)*k0(3)**3 + B(1)*k0(3)**2 + B(0)*k0(3))/&
          & ( (k0(3)-k1)*(k0(3)-k0(1))**2 *(k0(3)-k0(2))**2 )

     tk1=tk1+&
          & (-B(3)*k0(1)**3 *(k0(1)*k0(2)+k0(1)*k0(3)-2*k0(2)*k0(3))-&
          &   0.5*B(2)*k0(1)**2 *(k0(1)**2 + k0(1)*k0(2)+ k0(1)*k0(3) -3*k0(2)*k0(3))-&
          &   B(1)*k0(1)*(k0(1)**2 -k0(2)*k0(3))-&
          &   0.5*B(0)*(3*k0(1)**2 -k0(1)*k0(2)-k0(1)*k0(3)-k0(2)*k0(3)))*&
          & (log( k0(1)*conjg(k0(1))-2*real(k0(1))*k1 + k1**2)-&
          & 2*i*atan(aimag(k0(1))/(k1-real(k0(1)))) )/&
          & ((k0(1)-k0(2))**3 *(k0(1)-k0(3))**3 )

     tk1=tk1+&
          & (-B(3)*k0(2)**3 *(k0(1)*k0(2)+k0(2)*k0(3)-2*k0(1)*k0(3))-&
          &   0.5*B(2)*k0(2)**2 *(k0(2)**2 + k0(1)*k0(2)+ k0(2)*k0(3) -3*k0(1)*k0(3))-&
          &   B(1)*k0(2)*(k0(2)**2 -k0(1)*k0(3))-&
          &   0.5*B(0)*(3*k0(2)**2 -k0(1)*k0(2)-k0(1)*k0(3)-k0(2)*k0(3)))*&
          & (log( k0(2)*conjg(k0(2))-2*real(k0(2))*k1 + k1**2)-&
          & 2*i*atan(aimag(k0(2))/(k1-real(k0(2)))) )/&
          & ((k0(2)-k0(1))**3 *(k0(2)-k0(3))**3 )

     tk1=tk1+&
          & (-B(3)*k0(3)**3 *(k0(1)*k0(3)+k0(2)*k0(3)-2*k0(1)*k0(2))-&
          &   0.5*B(2)*k0(3)**2 *(k0(3)**2 + k0(1)*k0(3)+ k0(2)*k0(3) -3*k0(1)*k0(2))-&
          &   B(1)*k0(3)*(k0(3)**2 -k0(1)*k0(2))-&
          &   0.5*B(0)*(3*k0(3)**2 -k0(1)*k0(2)-k0(1)*k0(3)-k0(2)*k0(3)))*&
          & (log( k0(3)*conjg(k0(3))-2*real(k0(3))*k1 + k1**2)-&
          & 2*i*atan(aimag(k0(3))/(k1-real(k0(3)))) )/&
          & ((k0(3)-k0(1))**3 *(k0(3)-k0(2))**3 )

     tk2= (B(3)*k0(1)**4 + B(2)*k0(1)**3 + B(1)*k0(1)**2 + B(0)*k0(1))/&
          & ( (k0(1)-k2)*(k0(1)-k0(2))**2 *(k0(1)-k0(3))**2 )

     tk2=tk2+&
          & (B(3)*k0(2)**4 + B(2)*k0(2)**3 + B(1)*k0(2)**2 + B(0)*k0(2))/&
          & ( (k0(2)-k2)*(k0(2)-k0(1))**2 *(k0(2)-k0(3))**2 )

     tk2=tk2+&
          & (B(3)*k0(3)**4 + B(2)*k0(3)**3 + B(1)*k0(3)**2 + B(0)*k0(3))/&
          & ( (k0(3)-k2)*(k0(3)-k0(1))**2 *(k0(3)-k0(2))**2 )

     tk2=tk2+&
          & (-B(3)*k0(1)**3 *(k0(1)*k0(2)+k0(1)*k0(3)-2*k0(2)*k0(3))-&
          &   0.5*B(2)*k0(1)**2 *(k0(1)**2 + k0(1)*k0(2)+ k0(1)*k0(3) -3*k0(2)*k0(3))-&
          &   B(1)*k0(1)*(k0(1)**2 -k0(2)*k0(3))-&
          &   0.5*B(0)*(3*k0(1)**2 -k0(1)*k0(2)-k0(1)*k0(3)-k0(2)*k0(3)))*&
          & (log( k0(1)*conjg(k0(1))-2*real(k0(1))*k2 + k2**2)-&
          & 2*i*atan(aimag(k0(1))/(k2-real(k0(1)))) )/&
          & ((k0(1)-k0(2))**3 *(k0(1)-k0(3))**3 )

     tk2=tk2+&
          & (-B(3)*k0(2)**3 *(k0(1)*k0(2)+k0(2)*k0(3)-2*k0(1)*k0(3))-&
          &   0.5*B(2)*k0(2)**2 *(k0(2)**2 + k0(1)*k0(2)+ k0(2)*k0(3) -3*k0(1)*k0(3))-&
          &   B(1)*k0(2)*(k0(2)**2 -k0(1)*k0(3))-&
          &   0.5*B(0)*(3*k0(2)**2 -k0(1)*k0(2)-k0(1)*k0(3)-k0(2)*k0(3)))*&
          & (log( k0(2)*conjg(k0(2))-2*real(k0(2))*k2 + k2**2)-&
          & 2*i*atan(aimag(k0(2))/(k2-real(k0(2)))) )/&
          & ((k0(2)-k0(1))**3 *(k0(2)-k0(3))**3 )

     tk2=tk2+&
          & (-B(3)*k0(3)**3 *(k0(1)*k0(3)+k0(2)*k0(3)-2*k0(1)*k0(2))-&
          &   0.5*B(2)*k0(3)**2 *(k0(3)**2 + k0(1)*k0(3)+ k0(2)*k0(3) -3*k0(1)*k0(2))-&
          &   B(1)*k0(3)*(k0(3)**2 -k0(1)*k0(2))-&
          &   0.5*B(0)*(3*k0(3)**2 -k0(1)*k0(2)-k0(1)*k0(3)-k0(2)*k0(3)))*&
          & (log( k0(3)*conjg(k0(3))-2*real(k0(3))*k2 + k2**2)-&
          & 2*i*atan(aimag(k0(3))/(k2-real(k0(3)))) )/&
          & ((k0(3)-k0(1))**3 *(k0(3)-k0(2))**3 )


     !if one of the denominator's roots lies within the considered wavenumber interval add contribution from the pole

     contr_pole=0.0

     if(k1.lt.k2) then

        if((real(k0(1)).ge.k1).and.(real(k0(1)).lt.k2)) then

           call Term6_pole(k0(1),k0(2),k0(3),B,t6pole)

           if(aimag(k0(1)).lt.0.0) then
              sigma=-1
           else if(aimag(k0(1)).eq.0.0) then
              sigma=0
           else
              sigma=1
           endif

           contr_pole=contr_pole+i*pi*sigma*t6pole

        endif

        if((real(k0(2)).ge.k1).and.(real(k0(2)).lt.k2)) then

           call Term6_pole(k0(2),k0(1),k0(3),B,t6pole)

           if(aimag(k0(2)).lt.0.0) then
              sigma=-1
           else if(aimag(k0(2)).eq.0.0) then
              sigma=0
           else
              sigma=1
           endif

           contr_pole=contr_pole+i*pi*sigma*t6pole

        endif

        if((real(k0(3)).ge.k1).and.(real(k0(3)).lt.k2)) then

           call Term6_pole(k0(3),k0(1),k0(2),B,t6pole)

           if(aimag(k0(3)).lt.0.0) then
              sigma=-1
           else if(aimag(k0(3)).eq.0.0) then
              sigma=0
           else
              sigma=1
           endif

           contr_pole=contr_pole+i*pi*sigma*t6pole

        endif

     else

        if((real(k0(1)).lt.k1).and.(real(k0(1)).ge.k2)) then

           call Term6_pole(k0(1),k0(2),k0(3),B,t6pole)

           if(aimag(k0(1)).lt.0.0) then
              sigma=1
           else if(aimag(k0(1)).eq.0.0) then
              sigma=0
           else
              sigma=-1
           endif

           contr_pole=contr_pole+i*pi*sigma*t6pole

        endif

        if((real(k0(2)).lt.k1).and.(real(k0(2)).ge.k2)) then

           call Term6_pole(k0(2),k0(1),k0(3),B,t6pole)

           if(aimag(k0(2)).lt.0.0) then
              sigma=1
           else if(aimag(k0(2)).eq.0.0) then
              sigma=0
           else
              sigma=-1
           endif

           contr_pole=contr_pole+i*pi*sigma*t6pole

        endif

        if((real(k0(3)).lt.k1).and.(real(k0(3)).ge.k2)) then

           call Term6_pole(k0(3),k0(1),k0(2),B,t6pole)

           if(aimag(k0(3)).lt.0.0) then
              sigma=1
           else if(aimag(k0(3)).eq.0.0) then
              sigma=0
           else
              sigma=-1
           endif

           contr_pole=contr_pole+i*pi*sigma*t6pole

        endif

     endif


     t6=aimag(((tk2-tk1)+contr_pole)/(c_om(3)+i*c_ga(3))**2)


  endif


end subroutine Term6
