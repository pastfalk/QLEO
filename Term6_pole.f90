!> Compute pole contribution to the integral over T22 in the case where one of the denominator's roots lies within the considered wavenumber interval
!! \param k01 root of the integrand's denominator
!! \param k02 root of the integrand's denominator
!! \param k03 root of the integrand's denominator
!! \param alpha coefficients required for computation of pole contribution
!! \param t6pole pole contribution
subroutine Term6_pole(k01,k02,k03,alpha,t6pole)
  use param_mod
  implicit none
  real, dimension(0:3) :: alpha
  complex :: k01, k02, k03
  complex :: t6pole
  
  t6pole=(4*alpha(3)*k01**3 + 3*alpha(2)*k01**2+&
       &  2*alpha(1)*k01    +   alpha(0))/&
       &((k01-k02)*(k01-k03))**2 -&
       & 2*(2*k01-k02-k03)*&
       & (alpha(3)*k01**4 + alpha(2)*k01**3+&
       &  alpha(1)*k01**2 + alpha(0)*k01)/&
       &((k01-k02)*(k01-k03))**3


end subroutine Term6_pole
