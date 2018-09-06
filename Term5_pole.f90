!> Compute pole contribution to the integral over T21 in the case where one of the denominator's roots lies within the considered wavenumber interval
!! \param k01 root of the integrand's denominator
!! \param k02 root of the integrand's denominator
!! \param k03 root of the integrand's denominator
!! \param alpha coefficients required for computation of pole contribution
!! \param t5pole pole contribution
subroutine Term5_pole(k01,k02,k03,alpha,t5pole)
  use param_mod
  implicit none
  complex, dimension(0:5) :: alpha
  complex :: k01, k02, k03
  complex :: t5pole
  
  t5pole=(5*alpha(5)*k01**4 + 4*alpha(4)*k01**3 + 3*alpha(3)*k01**2+&
       &  2*alpha(2)*k01    +   alpha(1))/&
       &((k01-k02)*(k01-k03))**2-&
       & 2*(2*k01-k02-k03)*&
       & (alpha(5)*k01**5 + alpha(4)*k01**4 + alpha(3)*k01**3+&
       &  alpha(2)*k01**2 + alpha(1)*k01    + alpha(0))/&
       & ((k01-k02)*(k01-k03))**3


end subroutine Term5_pole
