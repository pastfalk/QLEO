!> Print out dispersion relation and magnetic energy spectrum at selected time steps
!! \param tstep current time step
!! \param time simulation time
!! \param krange wavenumber grid
!! \param Bksq magnetic energy spectrum
!! \param omega frequency spectrum
!! \param gamma growth rate spectrum
subroutine print_omega(tstep,time,krange,Bksq,omega,gamma)
  use param_mod
  implicit none
  integer :: ik
  integer :: tstep
  real :: time
  real, dimension(nk) :: krange, Bksq, omega, gamma
  character(len=6) :: aux1
  character(len=1024) :: filename

  write (aux1,'(I6)') tstep
  aux1 = adjustl(aux1)

  filename='distribution/evolution/omega/'//trim(aux1)//'.dat'

  open(unit=72,status='unknown',file=filename)
  write(72,*) time

  do ik=1,nk
     write(72,'(F12.8,E20.10,E20.10,E20.10)') krange(ik), Bksq(ik), omega(ik), gamma(ik)
  enddo

  close(72)


end subroutine print_omega
