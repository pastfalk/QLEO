!> Reads the velocity distribution data from the provided distribution files
!! \param restart index of the distribution files from which the simulation shall be restarted. If restart=0, distributions are read from the 'distribution/' directory and not from the restart files.
subroutine read_dist(restart)
  use param_mod
  implicit none
  integer :: ios
  real :: start_pe
  real :: vpa, vpe
  real :: dist_value
  integer :: ipara, iperp, iarb
  integer :: n
  integer :: restart
  character(len=1024), allocatable, dimension (:) :: filename

  npara_max=0
  nperp_max=0

  if(narb.ne.0) then

     allocate(npara(narb),nperp(narb))
     allocate(filename(narb))

     !determine dimension of the velocity grid for each velocity distribution

     do iarb=1,narb

        if(restart.eq.0) then

           write(filename(iarb),'(A25,I1,A4)') 'distribution/distribution',iarb,'.dat'

        else

           if(restart.lt.10) then
              write(filename(iarb),'(A28,I1,A1,I1,A4)') 'distribution/evolution/dist/',restart,'_',iarb,'.dat'
           else if((restart.ge.10).and.(restart.lt.100)) then
              write(filename(iarb),'(A28,I2,A1,I1,A4)') 'distribution/evolution/dist/',restart,'_',iarb,'.dat'
           else if((restart.ge.100).and.(restart.lt.1000)) then
              write(filename(iarb),'(A28,I3,A1,I1,A4)') 'distribution/evolution/dist/',restart,'_',iarb,'.dat'
           else if((restart.ge.1000).and.(restart.lt.10000)) then
              write(filename(iarb),'(A28,I4,A1,I1,A4)') 'distribution/evolution/dist/',restart,'_',iarb,'.dat'
           else if((restart.ge.10000).and.(restart.lt.100000)) then
              write(filename(iarb),'(A28,I5,A1,I1,A4)') 'distribution/evolution/dist/',restart,'_',iarb,'.dat'
           else
              write(filename(iarb),'(A28,I6,A1,I1,A4)') 'distribution/evolution/dist/',restart,'_',iarb,'.dat'
           endif

        endif

        open(unit=17,status='old',file=filename(iarb))

        npara(iarb)=0
        nperp(iarb)=0

        read(17,*,iostat=ios) vpa, start_pe, dist_value
        rewind(17)

        do while(.true.)

           read(17,*,iostat=ios) vpa, vpe, dist_value

           if (ios.ne.0) exit

           if (vpe.eq.start_pe) then
              npara(iarb)=npara(iarb)+1
              nperp(iarb)=0
           endif

           nperp(iarb)=nperp(iarb)+1

        enddo


        if(npara(iarb).gt.npara_max) then
           npara_max=npara(iarb)
           nperp_max=nperp(iarb)
        endif

        close(17)

        write(*,*) iarb, npara(iarb),nperp(iarb)
        
     enddo

     !read distribution data from the files

     allocate(distribution(npara_max,nperp_max,narb))
     allocate(vpara(npara_max,narb))
     allocate(vperp(nperp_max,narb))
     
     do iarb=1,narb

        open(unit=17,status='old',file=filename(iarb))

        do ipara=1,npara(iarb)

           do iperp=1,nperp(iarb)

              read(17,*) vpa, vpe, dist_value
              distribution(ipara,iperp,iarb)=dist_value
              if(ipara.eq.1) vperp(iperp,iarb)=vpe

           enddo

           vpara(ipara,iarb)=vpa

        enddo

        close(17)

     enddo

  endif

  allocate(nhalf(narb))

  do iarb=1,narb
     if(sym(iarb)) then
        nhalf(iarb)=int(0.5*(npara(iarb)+1))
     else
        nhalf(iarb)=1
     endif
  enddo


end subroutine read_dist
