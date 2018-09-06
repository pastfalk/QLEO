!> Compute the parallel and perpendicular beta by taking second velocity moments of the distribution function using the distribution's spline interpolation
!! \param bpara parallel beta
!! \param bperp perpendicular beta
subroutine get_beta(bpara,bperp)
  use param_mod
  implicit none
  real, dimension(npara(1)-1,nperp(1)-1,4,4) :: splcoeff
  real, dimension(narb) :: bpara, bperp
  integer :: ipara,iperp,iarb
  real :: h1,h2,h3,h4


  do iarb=1,narb

     !interpolate velocity distribution with cubic splines
     call splcoeff_beta(splcoeff,iarb)

     !compute parallel beta by (piecewise-)analytically evaluating the corresponding velocity moment integral

     bpara(iarb)=0.0

     do iperp=1,nperp(1)-1

        h1=0.0
        h2=0.0
        h3=0.0
        h4=0.0

        do ipara=1,npara(1)-1

           h1=h1+splcoeff(ipara,iperp,1,1)*1.0/6.0*(vpara(ipara+1,iarb)**6 - vpara(ipara,iarb)**6)
           h1=h1+splcoeff(ipara,iperp,2,1)*1.0/5.0*(vpara(ipara+1,iarb)**5 - vpara(ipara,iarb)**5)
           h1=h1+splcoeff(ipara,iperp,3,1)*1.0/4.0*(vpara(ipara+1,iarb)**4 - vpara(ipara,iarb)**4)
           h1=h1+splcoeff(ipara,iperp,4,1)*1.0/3.0*(vpara(ipara+1,iarb)**3 - vpara(ipara,iarb)**3)

           h2=h2+splcoeff(ipara,iperp,1,2)*1.0/6.0*(vpara(ipara+1,iarb)**6 - vpara(ipara,iarb)**6)
           h2=h2+splcoeff(ipara,iperp,2,2)*1.0/5.0*(vpara(ipara+1,iarb)**5 - vpara(ipara,iarb)**5)
           h2=h2+splcoeff(ipara,iperp,3,2)*1.0/4.0*(vpara(ipara+1,iarb)**4 - vpara(ipara,iarb)**4)
           h2=h2+splcoeff(ipara,iperp,4,2)*1.0/3.0*(vpara(ipara+1,iarb)**3 - vpara(ipara,iarb)**3)

           h3=h3+splcoeff(ipara,iperp,1,3)*1.0/6.0*(vpara(ipara+1,iarb)**6 - vpara(ipara,iarb)**6)
           h3=h3+splcoeff(ipara,iperp,2,3)*1.0/5.0*(vpara(ipara+1,iarb)**5 - vpara(ipara,iarb)**5)
           h3=h3+splcoeff(ipara,iperp,3,3)*1.0/4.0*(vpara(ipara+1,iarb)**4 - vpara(ipara,iarb)**4)
           h3=h3+splcoeff(ipara,iperp,4,3)*1.0/3.0*(vpara(ipara+1,iarb)**3 - vpara(ipara,iarb)**3)

           h4=h4+splcoeff(ipara,iperp,1,4)*1.0/6.0*(vpara(ipara+1,iarb)**6 - vpara(ipara,iarb)**6)
           h4=h4+splcoeff(ipara,iperp,2,4)*1.0/5.0*(vpara(ipara+1,iarb)**5 - vpara(ipara,iarb)**5)
           h4=h4+splcoeff(ipara,iperp,3,4)*1.0/4.0*(vpara(ipara+1,iarb)**4 - vpara(ipara,iarb)**4)
           h4=h4+splcoeff(ipara,iperp,4,4)*1.0/3.0*(vpara(ipara+1,iarb)**3 - vpara(ipara,iarb)**3)

        enddo

        bpara(iarb)=bpara(iarb)+1.0/5.0 *(vperp(iperp+1,iarb)**5 - vperp(iperp,iarb)**5)*h1
        bpara(iarb)=bpara(iarb)+1.0/4.0 *(vperp(iperp+1,iarb)**4 - vperp(iperp,iarb)**4)*h2
        bpara(iarb)=bpara(iarb)+1.0/3.0 *(vperp(iperp+1,iarb)**3 - vperp(iperp,iarb)**3)*h3
        bpara(iarb)=bpara(iarb)+1.0/2.0 *(vperp(iperp+1,iarb)**2 - vperp(iperp,iarb)**2)*h4

     enddo

     bpara(iarb)=4.0*pi*bpara(iarb)/mu(iarb)*dens(iarb)


     !compute perpendicular beta by (piecewise-)analytically evaluating the corresponding velocity moment integral

     bperp(iarb)=0.0

     do iperp=1,nperp(1)-1

        h1=0.0
        h2=0.0
        h3=0.0
        h4=0.0

        do ipara=1,npara(1)-1

           h1=h1+splcoeff(ipara,iperp,1,1)*1.0/4.0*(vpara(ipara+1,iarb)**4 - vpara(ipara,iarb)**4)
           h1=h1+splcoeff(ipara,iperp,2,1)*1.0/3.0*(vpara(ipara+1,iarb)**3 - vpara(ipara,iarb)**3)
           h1=h1+splcoeff(ipara,iperp,3,1)*1.0/2.0*(vpara(ipara+1,iarb)**2 - vpara(ipara,iarb)**2)
           h1=h1+splcoeff(ipara,iperp,4,1)        *(vpara(ipara+1,iarb)    - vpara(ipara,iarb)   )

           h2=h2+splcoeff(ipara,iperp,1,2)*1.0/4.0*(vpara(ipara+1,iarb)**4 - vpara(ipara,iarb)**4)
           h2=h2+splcoeff(ipara,iperp,2,2)*1.0/3.0*(vpara(ipara+1,iarb)**3 - vpara(ipara,iarb)**3)
           h2=h2+splcoeff(ipara,iperp,3,2)*1.0/2.0*(vpara(ipara+1,iarb)**2 - vpara(ipara,iarb)**2)
           h2=h2+splcoeff(ipara,iperp,4,2)        *(vpara(ipara+1,iarb)    - vpara(ipara,iarb)   )

           h3=h3+splcoeff(ipara,iperp,1,3)*1.0/4.0*(vpara(ipara+1,iarb)**4 - vpara(ipara,iarb)**4)
           h3=h3+splcoeff(ipara,iperp,2,3)*1.0/3.0*(vpara(ipara+1,iarb)**3 - vpara(ipara,iarb)**3)
           h3=h3+splcoeff(ipara,iperp,3,3)*1.0/2.0*(vpara(ipara+1,iarb)**2 - vpara(ipara,iarb)**2)
           h3=h3+splcoeff(ipara,iperp,4,3)        *(vpara(ipara+1,iarb)    - vpara(ipara,iarb)   )

           h4=h4+splcoeff(ipara,iperp,1,4)*1.0/4.0*(vpara(ipara+1,iarb)**4 - vpara(ipara,iarb)**4)
           h4=h4+splcoeff(ipara,iperp,2,4)*1.0/3.0*(vpara(ipara+1,iarb)**3 - vpara(ipara,iarb)**3)
           h4=h4+splcoeff(ipara,iperp,3,4)*1.0/2.0*(vpara(ipara+1,iarb)**2 - vpara(ipara,iarb)**2)
           h4=h4+splcoeff(ipara,iperp,4,4)        *(vpara(ipara+1,iarb)    - vpara(ipara,iarb)   )

        enddo

        bperp(iarb)=bperp(iarb)+1.0/7.0 *(vperp(iperp+1,iarb)**7 - vperp(iperp,iarb)**7)*h1
        bperp(iarb)=bperp(iarb)+1.0/6.0 *(vperp(iperp+1,iarb)**6 - vperp(iperp,iarb)**6)*h2
        bperp(iarb)=bperp(iarb)+1.0/5.0 *(vperp(iperp+1,iarb)**5 - vperp(iperp,iarb)**5)*h3
        bperp(iarb)=bperp(iarb)+1.0/4.0 *(vperp(iperp+1,iarb)**4 - vperp(iperp,iarb)**4)*h4

     enddo


     bperp(iarb)=2.0*pi*bperp(iarb)/mu(iarb)*dens(iarb)


  enddo


end subroutine get_beta
