!> Initializes the setup and advances the velocity distribution function in time by coupling the dispersion solver LEOPARD with the quasilinear kinetic system of differential equations
program main
  use param_mod
  use omp_lib
  implicit none
  integer :: ipara,iperp,iarb
  integer :: ik
  real :: kstart, kend
  real, allocatable, dimension (:) :: krange_adapt, omega_adapt, gamma_adapt, Bksq_adapt
  real, allocatable, dimension (:) :: krange_full, Bksq_full, omega_full, gamma_full
  real :: delta_f
  real, allocatable, dimension (:,:,:) :: f1
  real, allocatable, dimension (:,:,:) :: old_dist
  real :: C01, C10, C02, C20, C11
  real :: time, dt, tmax
  integer :: tstep
  integer :: sgn
  complex :: omega_start_adapt, increment_adapt
  complex :: omega_start_full, increment_full
  real :: Bstart
  real, allocatable, dimension(:) :: bpara, bperp
  real :: start, finish
  real, allocatable, dimension(:,:) :: contr1
  real, allocatable, dimension(:,:) :: contr2
  real, dimension(4) :: dummy1
  real, dimension(2) :: dummy2
  real, allocatable, dimension(:,:) :: c_om, c_ga, c_Bk
  integer :: restart
  integer :: print_interval
  real :: EB, Etot
  real :: gamma_max
  character(len=2) :: aux

  !read input data and velocity distributions and set up all arrays

  call omp_set_num_threads(32)
  call read_data(omega_start_full,increment_full,kstart,kend,Bstart,tmax,dt,sgn,restart)
  call read_dist(restart)

  tstep=restart
  print_interval=1

  allocate(krange_adapt(nk),omega_adapt(nk),gamma_adapt(nk),Bksq_adapt(nk))
  allocate(krange_full(nk),omega_full(nk),gamma_full(nk),Bksq_full(nk))
  allocate(delfdelpa(npara_max,nperp_max,narb),delfdelpe(npara_max,nperp_max,narb))
  allocate(delfdelpape(npara_max,nperp_max,narb))
  allocate(delfdelpepe(npara_max,nperp_max,narb),delfdelpapa(npara_max,nperp_max,narb))
  allocate(f1(npara_max,nperp_max,narb))
  allocate(old_dist(npara_max,nperp_max,narb))
  allocate(contr1(npara_max,4),contr2(npara_max,2))
  allocate(c_om(nk-1,0:3), c_ga(nk-1,0:3),c_Bk(nk-1,0:3))
  allocate(bpara(narb),bperp(narb))

  call initiate(restart,Bstart,kstart,kend,krange_full,Bksq_full,time,increment_full,omega_start_full)

  !loop over time

  do while(.true.)

     start=omp_get_wtime()

     do iarb=1,narb

        old_dist=distribution

        !print out current velocity distribution (needed for restart of simulation)

        if((modulo(tstep,print_interval).eq.0.0).and.(tstep.ne.restart)) then
           call print_dist(tstep,iarb)
        endif

     enddo

     write(*,*) '******************************', time
     write(*,*) '******************************* disp full'

     open(unit=57,status='old', position="append",file='omega.dat')

     write(57,*) '******************************', time, 'full'

     !compute dispersion relation over the full wavenumber range

     call disp_rel(krange_full,omega_start_full,increment_full,omega_full,gamma_full)

     !print out dispersion relation and magnetic energy spectrum for the full wavenumber range (needed for restart of simulation)

     if(modulo(tstep,print_interval).eq.0.0) then
        call print_omega(tstep,time,krange_full,Bksq_full,omega_full,gamma_full)
     endif


     !adjust wavenumber range to unstable modes only

     call adapt_krange(krange_full,Bksq_full,omega_full,gamma_full,&
          & krange_adapt,Bksq_adapt,increment_adapt,omega_start_adapt)

     write(*,*) '******************************* disp adapt'

     write(57,*) '******************************', time, 'adapt'


     !compute dispersion relation for the adjusted wavenumber range

     call disp_rel(krange_adapt,omega_start_adapt,increment_adapt, omega_adapt,gamma_adapt)

     write(*,*) '*******************************'


     !compute derivatives of velocity distribution function

     call deriv_dist


     !interpolate frequency, growth rate, and magnetic energy spectrum with cubic splines

     call splcoeff_spec(krange_adapt,omega_adapt,c_om,nk)
     call splcoeff_spec(krange_adapt,gamma_adapt,c_ga,nk)
     call splcoeff_spec(krange_adapt,Bksq_adapt, c_Bk,nk)


     !compute the change of the velocity distribution, df, based on the kinetic quasilinear equation

     do iarb=1,narb

        !$omp parallel do private(ipara,dummy1,dummy2)

        do ipara=nhalf(iarb),npara(iarb)

           call get_deltaf(krange_adapt,c_om,c_ga,c_Bk,vpara(ipara,iarb),sgn,dummy1,dummy2)

           contr1(ipara-nhalf(iarb)+1,:)=dummy1*dt
           contr2(ipara-nhalf(iarb)+1,:)=dummy2*dt

        enddo

        !$omp end parallel do



        !iperp=1

        do ipara=nhalf(iarb),npara(iarb)

           C01=   contr1(ipara-nhalf(iarb)+1,1) -&
                & contr1(ipara-nhalf(iarb)+1,2)*vpara(ipara,iarb)-&
                & contr1(ipara-nhalf(iarb)+1,3)*vpara(ipara,iarb)+&
                & contr1(ipara-nhalf(iarb)+1,4)*(vpara(ipara,iarb)**2)

           C10=   2*contr1(ipara-nhalf(iarb)+1,2)-&
                & 2*contr1(ipara-nhalf(iarb)+1,4)*vpara(ipara,iarb)

           C02=   contr1(ipara-nhalf(iarb)+1,1) -&
                & contr1(ipara-nhalf(iarb)+1,2)*vpara(ipara,iarb) -&
                & contr1(ipara-nhalf(iarb)+1,3)*vpara(ipara,iarb) +&
                & contr1(ipara-nhalf(iarb)+1,4)*vpara(ipara,iarb)**2

           C20= 0.0

           C11= 0.0

           delta_f=C01*delfdelpe(ipara-nhalf(iarb)+1,1,iarb)+&
                &  C10*delfdelpa(ipara-nhalf(iarb)+1,1,iarb)+&
                &  C02*delfdelpepe(ipara-nhalf(iarb)+1,1,iarb)+&
                &  C20*delfdelpapa(ipara-nhalf(iarb)+1,1,iarb)+&
                &  C11*delfdelpape(ipara-nhalf(iarb)+1,1,iarb)

           f1(ipara-nhalf(iarb)+1,1,iarb)=delta_f

        enddo

        !iperp>1

        !$omp parallel do private(ipara,iperp,delta_f,C01,C10,C02,C20,C11)

        do ipara=nhalf(iarb),npara(iarb)

           do iperp=2,nperp(iarb)

              C01=(  contr1(ipara-nhalf(iarb)+1,1) -&
                   & contr1(ipara-nhalf(iarb)+1,2)*vpara(ipara,iarb)-&
                   & contr1(ipara-nhalf(iarb)+1,3)*vpara(ipara,iarb)+&
                   & contr1(ipara-nhalf(iarb)+1,4)*(vpara(ipara,iarb)**2 - vperp(iperp,iarb)**2))/&
                   & vperp(iperp,iarb)+&
                   & contr2(ipara-nhalf(iarb)+1,1)*vperp(iperp,iarb)-&
                   & contr2(ipara-nhalf(iarb)+1,2)*vperp(iperp,iarb)*vpara(ipara,iarb)

              C10=   2*contr1(ipara-nhalf(iarb)+1,2)-&
                   & 2*contr1(ipara-nhalf(iarb)+1,4)*vpara(ipara,iarb)+&
                   & contr2(ipara-nhalf(iarb)+1,2)*vperp(iperp,iarb)**2

              C02=   contr1(ipara-nhalf(iarb)+1,1) -&
                   & contr1(ipara-nhalf(iarb)+1,2)*vpara(ipara,iarb) -&
                   & contr1(ipara-nhalf(iarb)+1,3)*vpara(ipara,iarb) +&
                   & contr1(ipara-nhalf(iarb)+1,4)*vpara(ipara,iarb)**2

              C20=   contr1(ipara-nhalf(iarb)+1,4)*vperp(iperp,iarb)**2

              C11=   contr1(ipara-nhalf(iarb)+1,2)*vperp(iperp,iarb)+&
                   & contr1(ipara-nhalf(iarb)+1,3)*vperp(iperp,iarb)-&
                   & 2*contr1(ipara-nhalf(iarb)+1,4)*vperp(iperp,iarb)*vpara(ipara,iarb)

              delta_f=C01*delfdelpe(ipara-nhalf(iarb)+1,iperp,iarb)+&
                   &  C10*delfdelpa(ipara-nhalf(iarb)+1,iperp,iarb)+&
                   &  C02*delfdelpepe(ipara-nhalf(iarb)+1,iperp,iarb)+&
                   &  C20*delfdelpapa(ipara-nhalf(iarb)+1,iperp,iarb)+&
                   &  C11*delfdelpape(ipara-nhalf(iarb)+1,iperp,iarb)

              f1(ipara-nhalf(iarb)+1,iperp,iarb)=delta_f

           enddo
        enddo

        !$omp end parallel do




        if(modulo(tstep,print_interval).eq.0.0) then
           call print_df(tstep,f1(:,:,iarb),contr1,contr2,iarb)
        endif


        !update the magnetic energy spectrum and the velocity distribution function 

        do ik=1,nk
           Bksq_full(ik)=Bksq_full(ik)+2*gamma_full(ik)*Bksq_full(ik)*dt
        enddo


        if(sym(iarb)) then

           do ipara=nhalf(iarb),npara(iarb)
              do iperp=1,nperp(iarb)

                 distribution(ipara,iperp,iarb)=old_dist(ipara,iperp,iarb)+f1(ipara-nhalf(iarb)+1,iperp,iarb)
                 distribution(npara(iarb)+1-ipara,iperp,iarb)=old_dist(npara(iarb)+1-ipara,iperp,iarb)+&
                      & f1(ipara-nhalf(iarb)+1,iperp,iarb)

              enddo
           enddo

        else

           do ipara=1,npara(iarb)
              do iperp=1,nperp(iarb)

                 distribution(ipara,iperp,iarb)=old_dist(ipara,iperp,iarb)+f1(ipara,iperp,iarb)

              enddo
           enddo

        endif


     enddo


     !compute beta parameters

     call get_beta(bpara,bperp)

     if(narb.gt.1) then
        do iarb=1,narb
           write(*,*) iarb, '- beta para:', bpara(iarb)
           write(*,*) iarb, '- beta perp:', bperp(iarb)
        enddo
     else
        write(*,*) 'beta para:', bpara(1)
        write(*,*) 'beta perp:', bperp(1)
     endif

     write(*,*) '************'

     !compute total magnetic energy and check energy conservation 

     call compute_EB(krange_full,Bksq_full,EB)

     Etot=EB

     do iarb=1,narb
        Etot=Etot+bperp(iarb)+0.5*bpara(iarb)
     enddo

     write(*,*) 'Etot:', Etot
     write(*,*) 'Bmax:', maxval(Bksq_full), sum(Bksq_full)
     write(*,*) '************'

     gamma_max=maxval(gamma_adapt)

     do iarb=1,narb

        write (aux,'(I2)') iarb
        aux = adjustl(aux)

        open(unit=157,status='old', position="append",file='anis_'//trim(aux)//'.dat')
        write(157,'(F12.8,E20.10,F12.8,F12.8,E20.10)') time, EB, bpara(iarb), bperp(iarb), gamma_max
        close(157)

     enddo

     finish=omp_get_wtime()
     write(*,*) '*******************************'
     write(*,*) 'Time elapsed: ', finish-start
     write(*,*) ' '
     time=time+dt

     tstep=tstep+1

     if(time.gt.tmax) then
        exit
     endif


     !adjust time step to avoid numerical instability in later stage of simulation

     if((time.ge.45.0).and.(time.lt.50.55)) then
        dt=0.05
        print_interval=4
     else if((time.ge.50.55).and.(time.lt.59.19)) then
        dt=0.01
        print_interval=20
     else if((time.ge.59.19).and.(time.lt.63.59)) then
        dt=0.005
        print_interval=40
     else if((time.ge.63.59).and.(time.lt.70.95)) then
        dt=0.002
        print_interval=80
     else if((time.ge.70.95).and.(time.lt.86.796)) then
        dt=0.001
        print_interval=160
     else if((time.ge.86.796).and.(time.lt.96.956)) then
        dt=0.0005
        print_interval=320
     else if(time.ge.96.956) then
        dt=0.0002
        print_interval=320
     endif

     close(57)

  enddo

  deallocate(krange_adapt,omega_adapt,gamma_adapt,Bksq_adapt)
  deallocate(krange_full,omega_full,gamma_full,Bksq_full)
  deallocate(vpara,vperp,distribution,f1,old_dist)
  deallocate(delfdelpa,delfdelpe,delfdelpape)
  deallocate(delfdelpapa,delfdelpepe)
  deallocate(q,mu,dens, drift,beta_para,beta_perp,beta_ratio)
  deallocate(contr1,contr2)
  deallocate(c_om,c_ga,c_Bk)
  deallocate(bpara,bperp)


end program main
