! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_debug
!
!> \brief MPAS ocean debug
!> \author Sara Calandrini
!> \date   11/01/2015
!> \details
!>  This module contains routines for initializing debug tracer
!
!-----------------------------------------------------------------------

module ocn_tracer_debug

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use mpas_log
   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

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

   public :: ocn_tracer_debug_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer, parameter ::  &
      debug_tracer_cnt = 2

contains

!***********************************************************************
!
!  routine ocn_tracer_debug_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Sara Calandrini
!> \date    10/28/2022
!> \details
!>  This routine initializes fields required for tracer surface flux restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_debug_init(domain,err)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: verticalMeshPool

      ! three dimensional pointers
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        debugTracers

      integer :: iCell, k, i

      !
      ! get tracers pools
      !

      err = 0

      if (config_use_debugTracers) then

      !
      ! Get tracer group so we can get the number of tracers in it
      !

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)

      do i=1,debug_tracer_cnt
         do iCell=1,nCellsOwned
            if (latCell(iCell)<0.0) then
               do k=1,maxLevelCell(iCell)
                  debugTracers(i,k,iCell) = 0.0_RKIND
               end do
            else
               do k=1,maxLevelCell(iCell)
                  debugTracers(i,k,iCell) = 1.0_RKIND
               end do
            end if
         end do
      end do

      end if  !  config_use_debugTracers)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_debug_init!}}}

end module ocn_tracer_debug

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

