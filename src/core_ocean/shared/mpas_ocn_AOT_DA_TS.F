! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_AOT_DA_TS
!
!> \brief MPAS AOT data assimilation for TS
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This module contains routines for assimilating temperature and salinity
!>  observations via the AOT algorithm.
!
!-----------------------------------------------------------------------

module ocn_AOT_DA_TS

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
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

   public :: ocn_AOT_DA_TS_compute, &
             ocn_AOT_DA_TS_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_AOT_DA_TS_compute
!
!> \brief   computes AOT DA term for TS 
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine computes the AOT DA term for temperature and salinity. 
!
!-----------------------------------------------------------------------

   subroutine ocn_AOT_DA_TS_compute(nTracers, nCellsSolve, maxLevelCell, &
                  tracer_interp, temperature_obs, salinity_obs,  &
                    tracer_mu, tracer_tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         maxLevelCell

!      ! two dimensional arrays
!      real (kind=RKIND), dimension(:,:), intent(in) :: &
!         layerThickness

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracer_interp
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         temperature_obs, salinity_obs

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve
      real (kind=RKIND), intent(in) :: tracer_mu

!      real (kind=RKIND), dimension(:,:,:), allocatable :: &
!         tracer_interp_obs

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
        tracer_tend

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, iLevel, iTracer
      integer, pointer :: indexTemperature, indexSalinity 
      err = 0

!      type (mpas_pool_type) :: statePool !< Input: State information
!      type (mpas_pool_type), pointer :: tracersPool            ! tracers and their tendency
!
!      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
!      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
!
!      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
!      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      !$omp do schedule(runtime) private(iLevel, iTracer)
      do iCell=1,nCellsSolve
        do iLevel=1,maxLevelCell(iCell)
             tracer_tend(1, iLevel, iCell) = tracer_tend(1, iLevel, iCell)  - tracer_mu * &
                                                   (tracer_interp(1, iLevel, iCell) -  &
                                                    temperature_obs(iLevel, iCell))
        enddo
      enddo

      do iCell=1,nCellsSolve
        do iLevel=1,maxLevelCell(iCell)
             tracer_tend(2, iLevel, iCell) = tracer_tend(2, iLevel, iCell)  - tracer_mu * &
                                                   (tracer_interp(2, iLevel, iCell) -  &
                                                    salinity_obs(iLevel, iCell))
        enddo
      enddo
      !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_AOT_DA_TS_compute!}}}

!***********************************************************************
!
!  routine ocn_AOT_DA_TS_init
!
!> \brief   Initializes ocean AOT DA for TS
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine initializes fields required for data assimilation of TS
!
!-----------------------------------------------------------------------

   subroutine ocn_AOT_DA_TS_init(tracer_mu, err)!{{{


      integer, intent(out) :: err !< Output: error flag
      real (kind=RKIND), intent(out) :: tracer_mu
      logical, pointer :: config_AOT_DA_on
      real (kind=RKIND), pointer :: config_AOT_DA_mu_TS
      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_AOT_DA_mu_TS', config_AOT_DA_mu_TS)
      call mpas_pool_get_config(ocnConfigs, 'config_AOT_DA_on', config_AOT_DA_on)

      tracer_mu = 0.0

      if (config_AOT_DA_on) then
          tracer_mu = config_AOT_DA_mu_TS
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_AOT_DA_TS_init!}}}

!***********************************************************************

end module ocn_AOT_DA_TS
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
