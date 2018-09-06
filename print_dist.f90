!> Print out velocity distribution and velocity derivatives of the distribution at selected time steps
!! \param tstep current time step
!! \param iarb index of species among the particle species with arbitrary velocity distribution
subroutine print_dist(tstep,iarb)
  use param_mod
  implicit none
  integer :: ipara, iperp
  integer :: tstep, iarb
  character(len=6) :: aux1
  character(len=2) :: aux2
  character(len=1024) :: filename

  write (aux1,'(I6)') tstep
  aux1 = adjustl(aux1)

  write (aux2,'(I2)') iarb
  aux2 = adjustl(aux2)

  filename='distribution/evolution/dist/'//trim(aux1)//'_'//trim(aux2)//'.dat'

  open(unit=72,status='unknown',file=filename)

  do ipara=1,npara(1)
     do iperp=1,nperp(1)     
        write(72,'(F12.8,F12.8,E20.10)') vpara(ipara,iarb), vperp(iperp,iarb), distribution(ipara,iperp,iarb) 
     enddo
  enddo

  close(72)

  filename='distribution/evolution/delfdelpa/'//trim(aux1)//'_'//trim(aux2)//'.dat'

  open(unit=72,status='unknown',file=filename)

  do ipara=nhalf(1),npara(1)
     do iperp=1,nperp(1)     
        write(72,'(F12.8,F12.8,E20.10,E20.10,E20.10)') vpara(ipara,iarb), vperp(iperp,iarb),&
             & delfdelpa(ipara-nhalf(1)+1,iperp,iarb), delfdelpapa(ipara-nhalf(1)+1,iperp,iarb),&
             & delfdelpape(ipara-nhalf(1)+1,iperp,iarb)  
     enddo
  enddo

  close(72)

  filename='distribution/evolution/delfdelpe/'//trim(aux1)//'_'//trim(aux2)//'.dat'

  open(unit=72,status='unknown',file=filename)

  do ipara=nhalf(1),npara(1)
     do iperp=1,nperp(1)     
        write(72,'(F12.8,F12.8,E20.10,E20.10)') vpara(ipara,iarb), vperp(iperp,iarb),&
             & delfdelpe(ipara-nhalf(1)+1,iperp,iarb), delfdelpepe(ipara-nhalf(1)+1,iperp,iarb)  
     enddo
  enddo

  close(72)


end subroutine print_dist
