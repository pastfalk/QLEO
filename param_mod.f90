!> Module containing all necessary parameters
module param_mod

implicit none

complex, parameter :: i=(0.0,1.0)
real, parameter :: pi= 4.0*atan(1.0)
integer :: nk
real, allocatable, dimension (:) :: q,mu
real, allocatable, dimension (:) :: dens
real, allocatable, dimension (:) :: drift
integer :: Nspecies
real :: theta
real, allocatable, dimension (:,:,:) :: distribution
real, allocatable, dimension (:,:) :: vpara,vperp
real, allocatable, dimension (:,:,:) :: delfdelpa
real, allocatable, dimension (:,:,:) :: delfdelpe
real, allocatable, dimension (:,:,:) :: delfdelpape

real, allocatable, dimension (:,:,:) :: delfdelpapa
real, allocatable, dimension (:,:,:) :: delfdelpepe
real :: delta 
real :: rf_error
real :: eps_error
real, allocatable, dimension (:) :: beta_para
real, allocatable, dimension (:) :: beta_perp
real, allocatable, dimension (:) :: beta_ratio

integer :: narb
integer, allocatable, dimension (:) :: mode
integer, allocatable, dimension (:) :: npara
integer, allocatable, dimension (:) :: nperp
integer, allocatable, dimension (:) :: nhalf


end module param_mod
