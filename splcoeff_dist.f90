!> Computes the spline coefficients for both parallel and perpendicular interpolation of the velocity distribution data
!! \param iarb index of species among the particle species with arbitrary velocity distribution which is to be interpolated
!! \param splcoeff1c array of spline coefficients for which distribution is first interpolated in parallel direction and then coefficients are interpolated again in perpendicular direction
!! \param splcoeff2c array of spline coefficients for which distribution is first interpolated in perpendicular direction and then coefficients are interpolated again in parallel direction
subroutine splcoeff_dist(iarb,splcoeff1c,splcoeff2c)
  use param_mod
  use omp_lib
  implicit none
  integer :: iarb
  integer :: ipara,iperp
  real, dimension (npara_max,nperp_max,0:3) :: splcoeff1a
  real, dimension (npara_max,nperp_max,0:3,3) :: splcoeff1b
  real, dimension (npara_max,nperp_max,0:3) :: splcoeff2a
  real, dimension (npara_max,nperp_max,0:3,3) :: splcoeff2b
  real, dimension (npara_max-1,nperp_max-1,4,3) :: splcoeff1c
  real, dimension (npara_max-1,nperp_max-1,4,3) :: splcoeff2c

  !Interpolate distribution over perpendicular velocity while scanning through the parallel velocity grid

  !$omp parallel do private(ipara,iperp)

  do ipara=1,npara(iarb)

     do iperp=1,nperp(iarb)
        splcoeff1a(ipara,iperp,0)=distribution(ipara,iperp,iarb)
     enddo

     call spline_interpol(vperp(:,iarb),splcoeff1a(ipara,:,0),splcoeff1a(ipara,:,1),&
          & splcoeff1a(ipara,:,2),splcoeff1a(ipara,:,3),nperp(iarb))

  enddo

  !$omp end parallel do


  !Interpolate coefficients of perpendicular spline interpolation over parallel velocity after taking first derivative w.r.t. to perpendicular velocity

  !$omp parallel do private(ipara,iperp)

  do iperp=1,nperp(iarb)-1

     do ipara=1,npara(iarb)

        splcoeff1b(ipara,iperp,0,1)= 3*splcoeff1a(ipara,iperp,3)
        splcoeff1b(ipara,iperp,0,2)= 2*splcoeff1a(ipara,iperp,2)-6*splcoeff1a(ipara,iperp,3)*vperp(iperp,iarb)
        splcoeff1b(ipara,iperp,0,3)= splcoeff1a(ipara,iperp,1)-2*splcoeff1a(ipara,iperp,2)*vperp(iperp,iarb)+&
             & 3*splcoeff1a(ipara,iperp,3)*vperp(iperp,iarb)*vperp(iperp,iarb)
     enddo

     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,0,1),splcoeff1b(:,iperp,1,1),&
          & splcoeff1b(:,iperp,2,1),splcoeff1b(:,iperp,3,1),npara(iarb))
     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,0,2),splcoeff1b(:,iperp,1,2),&
          & splcoeff1b(:,iperp,2,2),splcoeff1b(:,iperp,3,2),npara(iarb))
     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,0,3),splcoeff1b(:,iperp,1,3),&
          & splcoeff1b(:,iperp,2,3),splcoeff1b(:,iperp,3,3),npara(iarb))

     do ipara=1,npara(iarb)-1

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

     enddo

  enddo

  !$omp end parallel do


  !Now, interpolate distribution over parallel velocity while scanning through the perpendicular velocity grid

  !$omp parallel do private(ipara,iperp)

  do iperp=1,nperp(iarb)

     do ipara=1,npara(iarb)
        splcoeff2a(ipara,iperp,0)=distribution(ipara,iperp,iarb)
     enddo

     call spline_interpol(vpara(:,iarb),splcoeff2a(:,iperp,0),splcoeff2a(:,iperp,1),&
          & splcoeff2a(:,iperp,2),splcoeff2a(:,iperp,3),npara(iarb))

  enddo

  !$omp end parallel do


  !Interpolate coefficients of parallel spline interpolation over perpendicular velocity after taking first derivative w.r.t. to parallel velocity

  !$omp parallel do private(ipara,iperp)

  do ipara=1,npara(iarb)-1

     do iperp=1,nperp(iarb)
        splcoeff2b(ipara,iperp,0,1)=3*splcoeff2a(ipara,iperp,3)
        splcoeff2b(ipara,iperp,0,2)=2*splcoeff2a(ipara,iperp,2)-6*splcoeff2a(ipara,iperp,3)*vpara(ipara,iarb)
        splcoeff2b(ipara,iperp,0,3)=splcoeff2a(ipara,iperp,1)-2*splcoeff2a(ipara,iperp,2)*vpara(ipara,iarb)+&
             & 3*splcoeff2a(ipara,iperp,3)*vpara(ipara,iarb)*vpara(ipara,iarb)

     enddo

     call spline_interpol(vperp(:,iarb),splcoeff2b(ipara,:,0,1),splcoeff2b(ipara,:,1,1),&
          & splcoeff2b(ipara,:,2,1),splcoeff2b(ipara,:,3,1),nperp(iarb))
     call spline_interpol(vperp(:,iarb),splcoeff2b(ipara,:,0,2),splcoeff2b(ipara,:,1,2),&
          & splcoeff2b(ipara,:,2,2),splcoeff2b(ipara,:,3,2),nperp(iarb))
     call spline_interpol(vperp(:,iarb),splcoeff2b(ipara,:,0,3),splcoeff2b(ipara,:,1,3),&
          & splcoeff2b(ipara,:,2,3),splcoeff2b(ipara,:,3,3),nperp(iarb))

     do iperp=1,nperp(iarb)-1

        splcoeff2c(ipara,iperp,1,1)=splcoeff2b(ipara,iperp,3,1)
        splcoeff2c(ipara,iperp,2,1)=splcoeff2b(ipara,iperp,2,1)-3*splcoeff2b(ipara,iperp,3,1)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,3,1)=splcoeff2b(ipara,iperp,1,1)-2*splcoeff2b(ipara,iperp,2,1)*vperp(iperp,iarb)+&
             &  3*splcoeff2b(ipara,iperp,3,1)*vperp(iperp,iarb)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,4,1)=splcoeff2b(ipara,iperp,0,1)-splcoeff2b(ipara,iperp,1,1)*vperp(iperp,iarb)+&
             & splcoeff2b(ipara,iperp,2,1)*vperp(iperp,iarb)*vperp(iperp,iarb)-&
             & splcoeff2b(ipara,iperp,3,1)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)

        splcoeff2c(ipara,iperp,1,2)=splcoeff2b(ipara,iperp,3,2)
        splcoeff2c(ipara,iperp,2,2)=splcoeff2b(ipara,iperp,2,2)-3*splcoeff2b(ipara,iperp,3,2)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,3,2)=splcoeff2b(ipara,iperp,1,2)-2*splcoeff2b(ipara,iperp,2,2)*vperp(iperp,iarb)+&
             &  3*splcoeff2b(ipara,iperp,3,2)*vperp(iperp,iarb)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,4,2)=splcoeff2b(ipara,iperp,0,2)-splcoeff2b(ipara,iperp,1,2)*vperp(iperp,iarb)+&
             & splcoeff2b(ipara,iperp,2,2)*vperp(iperp,iarb)*vperp(iperp,iarb)-&
             & splcoeff2b(ipara,iperp,3,2)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)

        splcoeff2c(ipara,iperp,1,3)=splcoeff2b(ipara,iperp,3,3)
        splcoeff2c(ipara,iperp,2,3)=splcoeff2b(ipara,iperp,2,3)-3*splcoeff2b(ipara,iperp,3,3)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,3,3)=splcoeff2b(ipara,iperp,1,3)-2*splcoeff2b(ipara,iperp,2,3)*vperp(iperp,iarb)+&
             &  3*splcoeff2b(ipara,iperp,3,3)*vperp(iperp,iarb)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,4,3)=splcoeff2b(ipara,iperp,0,3)-splcoeff2b(ipara,iperp,1,3)*vperp(iperp,iarb)+&
             & splcoeff2b(ipara,iperp,2,3)*vperp(iperp,iarb)*vperp(iperp,iarb)-&
             & splcoeff2b(ipara,iperp,3,3)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)

     enddo

  enddo

  !$omp end parallel do


end subroutine splcoeff_dist
