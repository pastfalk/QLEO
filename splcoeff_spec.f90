!> Computes the spline coefficients for the given array of input values
!! \param krange wavenumber grid which the array of input values is spline-interpolated over
!! \param spec array of values which are interpolated with cubic splines
!! \param splcoeff array of spline coefficients
!! \param nk number of grid points
subroutine splcoeff_spec(krange,spec,splcoeff,nk)
  implicit none
  real, dimension (nk) :: krange
  real, dimension (nk) :: spec
  real, dimension (nk-1,0:3) :: splcoeff
  real, dimension (nk,0:3) :: splcoeff_dummy
  integer :: ik, nk

  !s_i = a_0 + a_1*(k-k_i) + a_2*(k-k_i)^2 + a_3*(k-k_i)^3
  
  splcoeff_dummy(:,0)=spec

  call spline_interpol(krange,splcoeff_dummy(:,0),splcoeff_dummy(:,1),&
          & splcoeff_dummy(:,2),splcoeff_dummy(:,3),nk)

  !s_i = c_0 + c_1*k + c_2*k^2 + c_3*k^3

  do ik=1,nk-1
     splcoeff(ik,3)=splcoeff_dummy(ik,3)
     splcoeff(ik,2)=splcoeff_dummy(ik,2)-3*splcoeff_dummy(ik,3)*krange(ik)
     splcoeff(ik,1)=splcoeff_dummy(ik,1)-2*splcoeff_dummy(ik,2)*krange(ik)+&
          & 3*splcoeff_dummy(ik,3)*krange(ik)**2
     splcoeff(ik,0)=splcoeff_dummy(ik,0)-  splcoeff_dummy(ik,1)*krange(ik)+&
          & splcoeff_dummy(ik,2)*krange(ik)**2-  splcoeff_dummy(ik,3)*krange(ik)**3
  enddo


end subroutine splcoeff_spec
