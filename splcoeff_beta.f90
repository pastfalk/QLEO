!> Computes the spline coefficients for both parallel and perpendicular interpolation of the velocity distribution data needed for the determination of the parallel and perpendicular beta
!! \param splcoeff1c array of spline coefficients
subroutine splcoeff_beta(splcoeff1c,iarb)
  use param_mod
  use omp_lib
  implicit none
  integer :: ipara,iperp,iarb
  real, dimension (npara(1),nperp(1),0:3) :: splcoeff1a
  real, dimension (npara(1),nperp(1),0:3,4) :: splcoeff1b
  real, dimension (npara(1)-1,nperp(1)-1,4,4) :: splcoeff1c

  !Interpolate distribution over perpendicular velocity while scanning through the parallel velocity grid

  !$omp parallel do private(ipara,iperp)

  do ipara=1,npara(1)

     do iperp=1,nperp(1)
        splcoeff1a(ipara,iperp,0)=distribution(ipara,iperp,iarb)
     enddo

     call spline_interpol(vperp(:,iarb),splcoeff1a(ipara,:,0),splcoeff1a(ipara,:,1),&
          & splcoeff1a(ipara,:,2),splcoeff1a(ipara,:,3),nperp(1))

  enddo

  !$omp end parallel do


  !Interpolate coefficients of perpendicular spline interpolation over parallel velocity

  !$omp parallel do private(ipara,iperp)

  do iperp=1,nperp(1)-1

     do ipara=1,npara(1)

        splcoeff1b(ipara,iperp,0,1)= splcoeff1a(ipara,iperp,3)
        splcoeff1b(ipara,iperp,0,2)= splcoeff1a(ipara,iperp,2)-3*splcoeff1a(ipara,iperp,3)*vperp(iperp,iarb)
        splcoeff1b(ipara,iperp,0,3)= splcoeff1a(ipara,iperp,1)-2*splcoeff1a(ipara,iperp,2)*vperp(iperp,iarb)+&
             & 3*splcoeff1a(ipara,iperp,3)*vperp(iperp,iarb)*vperp(iperp,iarb)
        splcoeff1b(ipara,iperp,0,4)= splcoeff1a(ipara,iperp,0)-splcoeff1a(ipara,iperp,1)*vperp(iperp,iarb)+&
             & splcoeff1a(ipara,iperp,2)*vperp(iperp,iarb)*vperp(iperp,iarb)-&
             & splcoeff1a(ipara,iperp,3)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)
     enddo

     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,0,1),splcoeff1b(:,iperp,1,1),&
          & splcoeff1b(:,iperp,2,1),splcoeff1b(:,iperp,3,1),npara(1))
     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,0,2),splcoeff1b(:,iperp,1,2),&
          & splcoeff1b(:,iperp,2,2),splcoeff1b(:,iperp,3,2),npara(1))
     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,0,3),splcoeff1b(:,iperp,1,3),&
          & splcoeff1b(:,iperp,2,3),splcoeff1b(:,iperp,3,3),npara(1))
     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,0,4),splcoeff1b(:,iperp,1,4),&
          & splcoeff1b(:,iperp,2,4),splcoeff1b(:,iperp,3,4),npara(1))

     do ipara=1,npara(1)-1

        splcoeff1c(ipara,iperp,1,1)=splcoeff1b(ipara,iperp,3,1)
        splcoeff1c(ipara,iperp,2,1)=splcoeff1b(ipara,iperp,2,1)-3*splcoeff1b(ipara,iperp,3,1)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,3,1)=splcoeff1b(ipara,iperp,1,1)-2*splcoeff1b(ipara,iperp,2,1)*vpara(ipara,iarb)+&
             &  3*splcoeff1b(ipara,iperp,3,1)*vpara(ipara,iarb)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,4,1)=splcoeff1b(ipara,iperp,0,1)-splcoeff1b(ipara,iperp,1,1)*vpara(ipara,iarb)+&
             & splcoeff1b(ipara,iperp,2,1)*vpara(ipara,iarb)*vpara(ipara,iarb)-&
             & splcoeff1b(ipara,iperp,3,1)*vpara(ipara,iarb)*vpara(ipara,iarb)*vpara(ipara,iarb)

        splcoeff1c(ipara,iperp,1,2)=splcoeff1b(ipara,iperp,3,2)
        splcoeff1c(ipara,iperp,2,2)=splcoeff1b(ipara,iperp,2,2)-3*splcoeff1b(ipara,iperp,3,2)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,3,2)=splcoeff1b(ipara,iperp,1,2)-2*splcoeff1b(ipara,iperp,2,2)*vpara(ipara,iarb)+&
             &  3*splcoeff1b(ipara,iperp,3,2)*vpara(ipara,iarb)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,4,2)=splcoeff1b(ipara,iperp,0,2)-splcoeff1b(ipara,iperp,1,2)*vpara(ipara,iarb)+&
             & splcoeff1b(ipara,iperp,2,2)*vpara(ipara,iarb)*vpara(ipara,iarb)-&
             & splcoeff1b(ipara,iperp,3,2)*vpara(ipara,iarb)*vpara(ipara,iarb)*vpara(ipara,iarb)

        splcoeff1c(ipara,iperp,1,3)=splcoeff1b(ipara,iperp,3,3)
        splcoeff1c(ipara,iperp,2,3)=splcoeff1b(ipara,iperp,2,3)-3*splcoeff1b(ipara,iperp,3,3)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,3,3)=splcoeff1b(ipara,iperp,1,3)-2*splcoeff1b(ipara,iperp,2,3)*vpara(ipara,iarb)+&
             &  3*splcoeff1b(ipara,iperp,3,3)*vpara(ipara,iarb)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,4,3)=splcoeff1b(ipara,iperp,0,3)-splcoeff1b(ipara,iperp,1,3)*vpara(ipara,iarb)+&
             & splcoeff1b(ipara,iperp,2,3)*vpara(ipara,iarb)*vpara(ipara,iarb)-&
             & splcoeff1b(ipara,iperp,3,3)*vpara(ipara,iarb)*vpara(ipara,iarb)*vpara(ipara,iarb)

        splcoeff1c(ipara,iperp,1,4)=splcoeff1b(ipara,iperp,3,4)
        splcoeff1c(ipara,iperp,2,4)=splcoeff1b(ipara,iperp,2,4)-3*splcoeff1b(ipara,iperp,3,4)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,3,4)=splcoeff1b(ipara,iperp,1,4)-2*splcoeff1b(ipara,iperp,2,4)*vpara(ipara,iarb)+&
             &  3*splcoeff1b(ipara,iperp,3,4)*vpara(ipara,iarb)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,4,4)=splcoeff1b(ipara,iperp,0,4)-splcoeff1b(ipara,iperp,1,4)*vpara(ipara,iarb)+&
             & splcoeff1b(ipara,iperp,2,4)*vpara(ipara,iarb)*vpara(ipara,iarb)-&
             & splcoeff1b(ipara,iperp,3,4)*vpara(ipara,iarb)*vpara(ipara,iarb)*vpara(ipara,iarb)

     enddo

  enddo

  !$omp end parallel do


end subroutine splcoeff_beta
