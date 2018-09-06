!> Initiate simulation by setting up the wavenumber grid and the initial magnetic field spectrum
!! \param restart index of the distribution files from which the simulation shall be restarted
!! \param Bstart amplitude of the initial magnetic energy perturbation
!! \param kstart start value of the requested wavenumber interval
!! \param kend final value of the requested wavenumber interval
!! \param krange wavenumber grid
!! \param Bksq initial magnetic energy spectrum
!! \param time simulation time
!! \param increment frequency increment suggested by the code user to determine the initial guesses for subsequent wavenumbers before subroutine polyfit() takes over
!! \param omega_start first initial frequency guess for the initial wavenumber provided by input.dat
subroutine initiate(restart,Bstart,kstart,kend,krange,Bksq,time,increment,omega_start)
  use param_mod
  implicit none
  integer :: ik
  real :: time
  complex :: omega_start
  complex :: increment
  real, dimension(nk) :: krange, Bksq, omega, gamma
  real :: dk, kend, kstart, Bstart
  integer :: restart
  character(len=5) :: aux1
  character(len=1024) :: filename
  integer :: iarb
  character(len=2) :: aux2

  if(restart.eq.0) then

     dk=(kend-kstart)/(nk-1.0)

     do ik=1,nk
        krange(ik)=kstart+(ik-1)*dk
        Bksq(ik)=Bstart
     enddo

     do iarb=1,narb

        write (aux2,'(I2)') iarb
        aux2 = adjustl(aux2)      

        open(unit=57,status='replace',file='omega.dat')
        close(57)
        open(unit=57,status='replace',file='anis_'//trim(aux2)//'.dat')
        close(57)

     enddo

     time=0.0
     restart=1

  else

     write (aux1,'(I5)') restart
     aux1 = adjustl(aux1)

     filename='distribution/evolution/omega/'//trim(aux1)//'.dat'
     open(unit=17,status='old',file=filename)

     read(17,*) time

     do ik=1,nk
        read(17,*) krange(ik), Bksq(ik), omega(ik), gamma(ik)
     enddo

     close(17)

     omega_start=omega(1)+i*gamma(1)
     increment=(omega(2)-omega(1))+i*(gamma(2)-gamma(1))

  endif


end subroutine initiate
