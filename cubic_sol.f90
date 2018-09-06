!> Compute complex roots of cubic polynomial with coefficients a(0:2), i.e. 0= x^3 + a(2)*x^2 + a(1)*x + a(0)
!! \param a complex coefficients of cubic equation
!! \param sol complex roots of the cubic equation
subroutine cubic_sol(a,sol)
  use param_mod
  implicit none
  complex, dimension(0:2) :: a
  complex, dimension(3) :: sol
  complex :: c1,c2

  c1=1.0/(3.0*2.0**(1.0/3.0)) *&
       & ( -2*a(2)**3 + 9*a(2)*a(1) - 27*a(0) +&
       &   3*sqrt(3.0)*sqrt(-a(2)**2 *a(1)**2 + 4*a(1)**3 + 4*a(2)**3 *a(0)-&
       &   18*a(2)*a(1)*a(0)+27*a(0)**2))**(1.0/3.0)

  c2=(-a(2)**2 + 3*a(1))/(9.0*c1)


  sol(1)= -a(2)/3.0 + c1 - c2

  sol(2)= -a(2)/3.0 + 0.5*(-1.0-sqrt(3.0)*i )*c1 -&
       & 0.5*(-1.0+sqrt(3.0)*i )*c2

  sol(3)= -a(2)/3.0 + 0.5*(-1.0+sqrt(3.0)*i )*c1 -&
       & 0.5*(-1.0-sqrt(3.0)*i )*c2


end subroutine cubic_sol

