!> Compute velocity derivatives of distribution function using local exponential fits
subroutine deriv_dist
  use param_mod
  implicit none
  integer :: ipara, iperp, iarb
  real :: b1, c1
  real :: b2, c2
  real :: vpa1, vpa2, vpa3
  real :: vpe1, vpe2, vpe3
  real :: dist1, dist2, dist3

  do iarb=1,narb

     !$omp parallel do private(ipara,iperp,b1,c1,vpa1,vpa2,vpa3,b2,c2,vpe1,vpe2,vpe3,dist1,dist2,dist3)

     do ipara=nhalf(iarb),npara(iarb)
        do iperp=1,nperp(iarb)

           dist2=distribution(ipara,iperp,iarb)
           vpa2=vpara(ipara,iarb)
           vpe2=vperp(iperp,iarb)

           if(ipara.eq.nhalf(iarb)) then

              vpa3=vpara(ipara+1,iarb)
              vpa1=vpara(ipara+2,iarb)

              dist3=distribution(ipara+1,iperp,iarb)
              dist1=distribution(ipara+2,iperp,iarb)

              c1= 0.5*((vpa2**2-vpa1**2)/  (log(dist1) - log(dist2) ) -&
                   &     (vpa3**2-vpa1**2)/  (log(dist1) - log(dist3) )  )/&
                   & (  (vpa1-vpa3)/  ( log(dist1) - log(dist3)) -&
                   &    (vpa1-vpa2)/  ( log(dist1) - log(dist2)) )

              b1= ((vpa2-c1)**2-(vpa1-c1)**2)/(log(dist1)-log(dist2))

           else if(ipara.eq.npara(iarb)) then

              vpa3=vpara(ipara-2,iarb)
              vpa1=vpara(ipara-1,iarb)

              dist3=distribution(ipara-2,iperp,iarb)
              dist1=distribution(ipara-1,iperp,iarb)

              c1= 0.5*((vpa2**2-vpa1**2)/  (log(dist1) - log(dist2) ) -&
                   &     (vpa3**2-vpa1**2)/  (log(dist1) - log(dist3) )  )/&
                   & (  (vpa1-vpa3)/  ( log(dist1) - log(dist3)) -&
                   &    (vpa1-vpa2)/  ( log(dist1) - log(dist2)) )

              b1= ((vpa2-c1)**2-(vpa1-c1)**2)/(log(dist1)-log(dist2))

           else

              vpa1=vpara(ipara-1,iarb)
              vpa3=vpara(ipara+1,iarb)

              dist1=distribution(ipara-1,iperp,iarb)
              dist3=distribution(ipara+1,iperp,iarb)

              if(dist1.ne.dist3) then

                 c1= 0.5*((vpa2**2-vpa1**2)/  (log(dist1) - log(dist2) ) -&
                      &     (vpa3**2-vpa1**2)/  (log(dist1) - log(dist3) )  )/&
                      & (  (vpa1-vpa3)/  ( log(dist1) - log(dist3)) -&
                      &    (vpa1-vpa2)/  ( log(dist1) - log(dist2)) )

                 b1= ((vpa2-c1)**2-(vpa1-c1)**2)/(log(dist1)-log(dist2))

              else

                 c1=vpa2
                 b1=(vpa3-vpa1)**2

              endif


           endif

           if(iperp.eq.1) then

              vpe1=vperp(iperp+1,iarb)
              dist1=distribution(ipara,iperp+1,iarb)

              b2= (vpe1**2-vpe2**2)/(log(dist2)-log(dist1))

           else if(iperp.eq.nperp(iarb)) then

              vpe3=vperp(iperp-2,iarb)
              vpe1=vperp(iperp-1,iarb)

              dist3=distribution(ipara,iperp-2,iarb)
              dist1=distribution(ipara,iperp-1,iarb)

              c2= 0.5*((vpe2**2-vpe1**2)/  (log(dist1) - log(dist2) ) -&
                   &     (vpe3**2-vpe1**2)/  (log(dist1) - log(dist3) )  )/&
                   & (  (vpe1-vpe3)/  ( log(dist1) - log(dist3)) -&
                   &    (vpe1-vpe2)/  ( log(dist1) - log(dist2)) )

              b2= ((vpe2-c2)**2-(vpe1-c2)**2)/(log(dist1)-log(dist2))

           else

              vpe1=vperp(iperp-1,iarb)
              vpe3=vperp(iperp+1,iarb)

              dist1=distribution(ipara,iperp-1,iarb)
              dist3=distribution(ipara,iperp+1,iarb)

              c2= 0.5*((vpe2**2-vpe1**2)/  (log(dist1) - log(dist2) ) -&
                   &     (vpe3**2-vpe1**2)/  (log(dist1) - log(dist3) )  )/&
                   & (  (vpe1-vpe3)/  ( log(dist1) - log(dist3)) -&
                   &    (vpe1-vpe2)/  ( log(dist1) - log(dist2)) )

              b2= ((vpe2-c2)**2-(vpe1-c2)**2)/(log(dist1)-log(dist2))

           endif


           if((ipara.eq.nhalf(iarb)).and.(ipara.ne.1)) then
              delfdelpa(ipara-nhalf(iarb)+1,iperp,iarb)=0.0
           else
              delfdelpa(ipara-nhalf(iarb)+1,iperp,iarb)=-2*(vpa2-c1)/b1 *dist2
           endif

           if(iperp.eq.1) then
              delfdelpe(ipara-nhalf(iarb)+1,iperp,iarb)=-2.0/b2 *dist2
              delfdelpepe(ipara-nhalf(iarb)+1,iperp,iarb)=-2.0/b2 *dist2
           else
              delfdelpe(ipara-nhalf(iarb)+1,iperp,iarb)=-2*(vpe2-c2)/b2 *dist2
              delfdelpepe(ipara-nhalf(iarb)+1,iperp,iarb)=2*(2*(vpe2-c2)**2/b2**2 -1.0/b2 )*dist2
           endif

           if((ipara.eq.nhalf(iarb)).or.(iperp.eq.1)) then
              delfdelpape(ipara-nhalf(iarb)+1,iperp,iarb)=0.0
           else
              delfdelpape(ipara-nhalf(iarb)+1,iperp,iarb)=4*(vpe2-c2)/b2 *(vpa2-c1)/b1  *dist2
           endif

           delfdelpapa(ipara-nhalf(iarb)+1,iperp,iarb)=2*(2*(vpa2-c1)**2/b1**2 -1.0/b1 )*dist2

        enddo

     enddo

     !$omp end parallel do

  enddo


end subroutine deriv_dist
