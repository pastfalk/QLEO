!> Compute total magnetic energy from magnetic energy spectrum
!! \param krange_full full wavenumber grid
!! \param Bksq_full magnetic energy spectrum over the full wavenumber grid
!! \param EB total magnetic energy
subroutine compute_EB(krange_full,Bksq_full,EB)
  use param_mod
  implicit none
  integer :: ik
  real, dimension (nk) :: krange_full
  real, dimension (nk) :: Bksq_full
  real :: EB
  real, dimension (nk-1,4) :: c_Bk

  !interpolate magnetic energy spectrum with cubic splines
  call splcoeff_spec(krange_full,Bksq_full,c_Bk,nk)

  !compute total magnetic energy by (piecewise-)analytically evaluating the wavenumber integral over the magnetic energy spectrum

  EB=0.0

  do ik=1,nk-1
     
     EB=EB+ c_Bk(ik,4)/4.0 *(krange_full(ik+1)**4 - krange_full(ik)**4)+&
          & c_Bk(ik,3)/3.0 *(krange_full(ik+1)**3 - krange_full(ik)**3)+&
          & c_Bk(ik,2)/2.0 *(krange_full(ik+1)**2 - krange_full(ik)**2)+&
          & c_Bk(ik,1)/1.0 *(krange_full(ik+1)**1 - krange_full(ik)**1)
  
  enddo
  
  if(krange_full(1).gt.krange_full(2)) then
     EB=-EB
  endif


end subroutine compute_EB
