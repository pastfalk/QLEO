!> Print out change in velocity distribution due to different terms in quasilinear equation at selected time steps
!! \param tstep current time step
!! \param df estimated total change in velocity distribution at current time step 
!! \param contr1 solutions of wavenumber integration
!! \param contr2 solutions of wavenumber integration
!! \param iarb index of species among the particle species with arbitrary velocity distribution
subroutine print_df(tstep,df,contr1,contr2,iarb)
  use param_mod
  implicit none
  integer :: ipara, iperp, iarb
  integer :: tstep
  character(len=6) :: aux1
  character(len=2) :: aux2
  character(len=1024) :: filename
  real, dimension(npara_max,nperp_max) :: df
  real, dimension(npara_max,4) :: contr1
  real, dimension(npara_max,2) :: contr2
  
  write (aux1,'(I6)') tstep
  aux1 = adjustl(aux1)

  write (aux2,'(I2)') iarb
  aux2 = adjustl(aux2)

  filename='distribution/evolution/df/'//trim(aux1)//'_'//trim(aux2)//'.dat'

  open(unit=72,status='unknown',file=filename)

  do ipara=nhalf(iarb),npara(iarb)
     do iperp=1,nperp(iarb)     
        write(72,'(F12.8,F12.8,E20.10)') vpara(ipara,iarb), vperp(iperp,iarb), df(ipara-nhalf(iarb)+1,iperp)
     enddo
  enddo

  close(72)

  filename='distribution/evolution/T1/'//trim(aux1)//'_'//trim(aux2)//'.dat'

  open(unit=72,status='unknown',file=filename)

  do ipara=nhalf(iarb),npara(iarb)
     write(72,'(F12.8,E20.10,E20.10,E20.10,E20.10)') vpara(ipara,iarb), contr1(ipara-nhalf(iarb)+1,1),&
          & contr1(ipara-nhalf(iarb)+1,2), contr1(ipara-nhalf(iarb)+1,3), contr1(ipara-nhalf(iarb)+1,4)
  enddo

  close(72)

  filename='distribution/evolution/T2/'//trim(aux1)//'_'//trim(aux2)//'.dat'

  open(unit=72,status='unknown',file=filename)
  
  do ipara=nhalf(iarb),npara(iarb)
     write(72,'(F12.8,E20.10,E20.10)') vpara(ipara,iarb), contr2(ipara-nhalf(iarb)+1,1),&
          & contr2(ipara-nhalf(iarb)+1,2)
  enddo

  close(72)


end subroutine print_df
