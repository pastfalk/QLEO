!> Compute pole contribution to the integral over T13 in the case where one of the denominator's roots lies within the considered wavenumber interval
!! \param k01 root of the integrand's denominator
!! \param k02 root of the integrand's denominator
!! \param k03 root of the integrand's denominator
!! \param alpha coefficients required for computation of pole contribution
!! \param t3pole pole contribution
subroutine Term3_pole(k01,k02,k03,alpha,t3pole)
  use param_mod
  implicit none
  complex, dimension(0:3) :: alpha
  complex :: k01, k02, k03
  complex :: t3pole
    
  t3pole=(alpha(3)*k01**2 + alpha(2)*k01 + alpha(1) +alpha(0)/k01)/&
       &((k01-k02)*(k01-k03))


end subroutine Term3_pole
