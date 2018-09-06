!> Compute pole contribution to the integral over T14 in the case where one of the denominator's roots lies within the considered wavenumber interval
!! \param k01 root of the integrand's denominator
!! \param k02 root of the integrand's denominator
!! \param k03 root of the integrand's denominator
!! \param alpha coefficients required for computation of pole contribution
!! \param t4pole pole contribution
subroutine Term4_pole(k01,k02,k03,alpha,t4pole)
  use param_mod
  implicit none
  complex, dimension(0:2) :: alpha
  complex :: k01, k02, k03
  complex :: t4pole

  t4pole=(alpha(2)*k01**2 + alpha(1)*k01 + alpha(0))/&
       &((k01-k02)*(k01-k03))


end subroutine Term4_pole
