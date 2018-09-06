!> Analytic solution of the wavenumber integration of term T12 = dB_k^2 * omega^* /k /(omega-k*vpa-sgn*mu*q) over interval [k1,k2]
!! \param k1 lower limit of the considered wavenumber interval
!! \param k2 upper limit of the considered wavenumber interval
!! \param k0 roots of the integrand's denominator
!! \param c_om spline coefficients of interpolated frequency spectrum
!! \param c_ga spline coefficients of interpolated growth rate spectrum
!! \param c_Bk spline coefficients of interpolated magnetic energy spectrum
!! \param vpa parallel velocity
!! \param sgn signum in the integrand's denominator which determines whether right-hand or left-hand modes are considered
!! \param t2 result of the integration
subroutine Term2(k1,k2,k0,c_om,c_ga,c_Bk,vpa,sgn,t2)
  use param_mod
  implicit none
  real :: k1, k2
  complex, dimension(3) :: k0
  real, dimension(0:3) :: c_om, c_ga, c_Bk
  complex :: tk1, tk2
  real :: t2
  real :: vpa
  integer :: sgn
  complex, dimension(0:2) :: D
  complex, dimension(0:6) :: B
  complex, dimension(3) :: lambda
  complex, dimension(0:3) :: alpha
  integer :: sigma
  complex :: t2pole, contr_pole

  if(vpa.eq.0.0) then
     t2=0.0
  else

     B(0)=  (c_om(0)-i*c_ga(0))*c_Bk(0)

     B(1)=  (c_om(1)-i*c_ga(1))*c_Bk(0)+(c_om(0)-i*c_ga(0))*c_Bk(1)

     B(2)=  (c_om(2)-i*c_ga(2))*c_Bk(0)+(c_om(1)-i*c_ga(1))*c_Bk(1)+&
          & (c_om(0)-i*c_ga(0))*c_Bk(2)

     B(3)=  (c_om(3)-i*c_ga(3))*c_Bk(0)+(c_om(2)-i*c_ga(2))*c_Bk(1)+&
          & (c_om(1)-i*c_ga(1))*c_Bk(2)+(c_om(0)-i*c_ga(0))*c_Bk(3)

     B(4)=  (c_om(3)-i*c_ga(3))*c_Bk(1)+(c_om(2)-i*c_ga(2))*c_Bk(2)+&
          & (c_om(1)-i*c_ga(1))*c_Bk(3)

     B(5)=  (c_om(3)-i*c_ga(3))*c_Bk(2)+(c_om(2)-i*c_ga(2))*c_Bk(3)

     B(6)=  (c_om(3)-i*c_ga(3))*c_Bk(3)

     D(0)=(c_om(0)-sgn*q(1)*mu(1)+i*c_ga(0))/&
          & (c_om(3)+i*c_ga(3))

     D(1)=(c_om(1)-vpa+i*c_ga(1))/&
          & (c_om(3)+i*c_ga(3))

     D(2)=(c_om(2)+i*c_ga(2))/(c_om(3)+i*c_ga(3))

     lambda(1)=B(4)-D(1)*B(6)-D(2)*(B(5)-D(2)*B(6))

     lambda(2)=B(3)-D(0)*B(6)-D(1)*(B(5)-D(2)*B(6))

     lambda(3)=B(2)-D(0)*(B(5)-D(2)*B(6))

     alpha(0)=B(0)

     alpha(1)=B(1)-D(0)*lambda(1)

     alpha(2)=lambda(3)-D(1)*lambda(1)

     alpha(3)=lambda(2)-D(2)*lambda(1)

     tk1=   B(6)/3.0 * k1**3 +&
          & 0.5*(B(5)-D(2)*B(6))*k1**2 +&
          & lambda(1)*k1+&
          & (-B(0)*log(k1))/(k0(1)*k0(2)*k0(3))

     tk1=tk1+0.5*&
          ( alpha(3)*k0(1)**2 + alpha(2)*k0(1) +&
          & alpha(1)+B(0)/k0(1))*&
          & (log( k0(1)*conjg(k0(1))-2*real(k0(1))*k1 + k1**2)-&
          & 2*i*atan(aimag(k0(1))/(k1-real(k0(1)))) )/&
          & ((k0(1)-k0(2))*(k0(1)-k0(3)))

     tk1=tk1+0.5*&
          ( alpha(3)*k0(2)**2 + alpha(2)*k0(2) +&
          & alpha(1)+B(0)/k0(2))*&
          & (log( k0(2)*conjg(k0(2))-2*real(k0(2))*k1 + k1**2)-&
          & 2*i*atan(aimag(k0(2))/(k1-real(k0(2)))) )/&
          & ((k0(2)-k0(1))*(k0(2)-k0(3)))

     tk1=tk1+0.5*&
          ( alpha(3)*k0(3)**2 + alpha(2)*k0(3) +&
          & alpha(1)+B(0)/k0(3))*&
          & (log( k0(3)*conjg(k0(3))-2*real(k0(3))*k1 + k1**2)-&
          & 2*i*atan(aimag(k0(3))/(k1-real(k0(3)))) )/&
          & ((k0(3)-k0(1))*(k0(3)-k0(2)))

     tk2=   B(6)/3.0 * k2**3 +&
          & 0.5*(B(5)-D(2)*B(6))*k2**2 +&
          & lambda(1)*k2+&
          & (-B(0)*log(k2))/(k0(1)*k0(2)*k0(3))

     tk2=tk2+0.5*&
          ( alpha(3)*k0(1)**2 + alpha(2)*k0(1) +&
          & alpha(1)+B(0)/k0(1))*&
          & (log( k0(1)*conjg(k0(1))-2*real(k0(1))*k2 + k2**2)-&
          & 2*i*atan(aimag(k0(1))/(k2-real(k0(1)))) )/&
          & ((k0(1)-k0(2))*(k0(1)-k0(3)))

     tk2=tk2+0.5*&
          ( alpha(3)*k0(2)**2 + alpha(2)*k0(2) +&
          & alpha(1)+B(0)/k0(2))*&
          & (log( k0(2)*conjg(k0(2))-2*real(k0(2))*k2 + k2**2)-&
          & 2*i*atan(aimag(k0(2))/(k2-real(k0(2)))) )/&
          & ((k0(2)-k0(1))*(k0(2)-k0(3)))

     tk2=tk2+0.5*&
          ( alpha(3)*k0(3)**2 + alpha(2)*k0(3) +&
          & alpha(1)+B(0)/k0(3))*&
          & (log( k0(3)*conjg(k0(3))-2*real(k0(3))*k2 + k2**2)-&
          & 2*i*atan(aimag(k0(3))/(k2-real(k0(3)))) )/&
          & ((k0(3)-k0(1))*(k0(3)-k0(2)))


     !if one of the denominator's roots lies within the considered wavenumber interval add contribution from the pole

     contr_pole=0.0

     if(k1.lt.k2) then

        if((real(k0(1)).ge.k1).and.(real(k0(1)).lt.k2)) then

           call Term2_pole(k0(1),k0(2),k0(3),alpha,t2pole)

           if(aimag(k0(1)).lt.0.0) then
              sigma=-1
           else if(aimag(k0(1)).eq.0.0) then
              sigma=0
           else
              sigma=1
           endif

           contr_pole=contr_pole+i*pi*sigma*t2pole

        endif

        if((real(k0(2)).ge.k1).and.(real(k0(2)).lt.k2)) then

           call Term2_pole(k0(2),k0(1),k0(3),alpha,t2pole)

           if(aimag(k0(2)).lt.0.0) then
              sigma=-1
           else if(aimag(k0(2)).eq.0.0) then
              sigma=0
           else
              sigma=1
           endif

           contr_pole=contr_pole+i*pi*sigma*t2pole

        endif

        if((real(k0(3)).ge.k1).and.(real(k0(3)).lt.k2)) then

           call Term2_pole(k0(3),k0(1),k0(2),alpha,t2pole)

           if(aimag(k0(3)).lt.0.0) then
              sigma=-1
           else if(aimag(k0(3)).eq.0.0) then
              sigma=0
           else
              sigma=1
           endif

           contr_pole=contr_pole+i*pi*sigma*t2pole

        endif

     else

        if((real(k0(1)).lt.k1).and.(real(k0(1)).ge.k2)) then

           call Term2_pole(k0(1),k0(2),k0(3),alpha,t2pole)

           if(aimag(k0(1)).lt.0.0) then
              sigma=1
           else if(aimag(k0(1)).eq.0.0) then
              sigma=0
           else
              sigma=-1
           endif

           contr_pole=contr_pole+i*pi*sigma*t2pole

        endif

        if((real(k0(2)).lt.k1).and.(real(k0(2)).ge.k2)) then

           call Term2_pole(k0(2),k0(1),k0(3),alpha,t2pole)

           if(aimag(k0(2)).lt.0.0) then
              sigma=1
           else if(aimag(k0(2)).eq.0.0) then
              sigma=0
           else
              sigma=-1
           endif

           contr_pole=contr_pole+i*pi*sigma*t2pole

        endif

        if((real(k0(3)).lt.k1).and.(real(k0(3)).ge.k2)) then

           call Term2_pole(k0(3),k0(1),k0(2),alpha,t2pole)

           if(aimag(k0(3)).lt.0.0) then
              sigma=1
           else if(aimag(k0(3)).eq.0.0) then
              sigma=0
           else
              sigma=-1
           endif

           contr_pole=contr_pole+i*pi*sigma*t2pole

        endif

     endif


     t2=aimag(((tk2-tk1)+contr_pole)/(c_om(3)+i*c_ga(3)))


  endif


end subroutine Term2
