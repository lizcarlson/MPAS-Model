! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_AOT_DA_velocity
!
!> \brief MPAS ocean AOT data assimilation for horizontal velocity
!> \author  Todd Ringler
!> \date    5 January 2012
!> \details
!>  This module contains the routine for updating the horizontal velocity
!>  tendency with the AOT data assimilation algorithm (current version: AOT
!>  only on the normalVelocity, no interpolation)
!
!-----------------------------------------------------------------------
!/
!&AOT_data_assimilation
!    config_AOT_DA_on = .false.
!    config_AOT_DA_mu_vel = 0.0
!    config_AOT_DA_mu_TS = 0.0

module ocn_AOT_DA_velocity

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use ocn_constants

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_AOT_DA_velocity_tend, &
             ocn_AOT_DA_velocity_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!   real (kind=RKIND) :: velocityTend_mu 


!***********************************************************************

contains

!***********************************************************************
!
!  ocn_AOT_DA_velocity_tend
!
!> \brief   Computes AOT data assimilation for horizontal velocity
!> \author  Todd Ringler
!> \date    5 January 2012
!> \details
!>  This routine computes the horizontal velocity
!>  tendency with regards to the AOT data assimilation algorithm
!
!-----------------------------------------------------------------------

   subroutine ocn_AOT_DA_velocity_tend(meshPool, normalVelocity_interp, & 
                                        normalVelocity_interp_obs, &
                                        velocityTend_mu, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity_interp, &        !< Input: velocity
         normalVelocity_interp_obs    !< Input: observed velocity

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), intent(in) :: velocityTend_mu

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, k, nEdges
      integer, dimension(:), pointer :: nEdgesArray
      integer, dimension(:), pointer :: maxLevelEdgeTop
      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------
      err = 0

!      if ( .not. rayleighFrictionOn ) return ! EDIT THIS LINE FOR AOT DATA ASSIMILATION

      call mpas_timer_start('AOT DA horizontal velocity')

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      nEdges = nEdgesArray( 1 )

!      if ( rayleighFrictionOn ) then
         !$omp do schedule(runtime) private(k)
         do iEdge = 1, nEdges
           do k = 1, maxLevelEdgeTop(iEdge)

              tend(k,iEdge) = tend(k,iEdge) - velocityTend_mu * (normalVelocity_interp(k,iEdge) - normalVelocity_interp_obs(k,iEdge))

           enddo
         enddo
         !$omp end do
!      endif

!      if ( rayleighBottomFrictionOn ) then
!         !$omp do schedule(runtime) private(k)
!         do iEdge = 1, nEdges
!
!           k = maxLevelEdgeTop(iEdge)
!           tend(k,iEdge) = tend(k,iEdge) - rayleighDampingCoef * normalVelocity(k,iEdge)
!
!         enddo
!         !$omp end do
!      endif

      call mpas_timer_stop('AOT DA horizontal velocity')

   !--------------------------------------------------------------------

   end subroutine ocn_AOT_DA_velocity_tend!}}}

!***********************************************************************
!
!  ocn_vel_forcing_rayleigh_init
!
!> \brief   Initializes ocean Rayleigh friction
!> \author  Todd Ringler
!> \date    5 January 2012
!> \details
!>  This routine initializes quantities related to
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_AOT_DA_velocity_init(velocityTend_mu, err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag
      real (kind=RKIND), intent(out) :: velocityTend_mu
      logical, pointer :: config_AOT_DA_on
      real (kind=RKIND), pointer :: config_AOT_DA_mu_vel
      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_AOT_DA_mu_vel', config_AOT_DA_mu_vel)
      call mpas_pool_get_config(ocnConfigs, 'config_AOT_DA_on', config_AOT_DA_on)

      velocityTend_mu = 0.0

      if (config_AOT_DA_on) then
          velocityTend_mu = config_AOT_DA_mu_vel
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_AOT_DA_velocity_init!}}}

!***********************************************************************

end module ocn_AOT_DA_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
